#include "rng.h"

#include <climits>
#include <cmath>
#include <cstdint>
#include <limits>
#include <thread>
#include <vector>

#include "gtest/gtest.h"

namespace {

using rng_h::RNG;
using rng_h::uint;

// All tests use fixed seeds, so they are deterministic.  Statistical
// tolerances are still set at about five standard errors so the tests stay
// valid if the generator or the seeds are ever changed.

struct Moments {
  double mean, var;
};

template <typename F>
Moments SampleMoments(int n, F draw) {
  double sum = 0.0, sum_sq = 0.0;
  for (int i = 0; i < n; ++i) {
    const double x = draw();
    sum += x;
    sum_sq += x * x;
  }
  const double mean = sum / n;
  return {mean, sum_sq / n - mean * mean};
}

// Checks the sample mean and variance against the distribution's mean and
// variance.  'kurt' is the distribution's (non-excess) kurtosis bound used
// to size the tolerance on the variance.
void ExpectMoments(const Moments& m, int n, double mean, double var,
                   double kurt = 10.0) {
  EXPECT_NEAR(m.mean, mean, 5.0 * std::sqrt(var / n));
  EXPECT_NEAR(m.var, var, 5.0 * var * std::sqrt(kurt / n));
}

const int kN = 400000;

// ___________________________________________________________________________
// The core generator

// A straightforward 32-bit reference implementation of the generator.
struct ReferenceKiss {
  uint32_t z, w, jsr, jcong;
  uint32_t Next() {
    z = 36969 * (z & 65535) + (z >> 16);
    w = 18000 * (w & 65535) + (w >> 16);
    jcong = 69069 * jcong + 1234567;
    jsr ^= jsr << 13;
    jsr ^= jsr >> 17;
    jsr ^= jsr << 5;
    return (((z << 16) + w) ^ jcong) + jsr;
  }
};

TEST(Kiss, MatchesReferenceImplementation) {
  ReferenceKiss ref = {12345, 65435, 34221, 12345};
  RNG rng(12345, 65435, 34221, 12345);
  for (int i = 0; i < 100000; ++i) {
    ASSERT_EQ(ref.Next(), rng.rand_int32()) << "at draw " << i;
  }
}

// Marsaglia's 1999 posting includes a test program and its expected output.
// The MWC and CONG components must reproduce those values exactly.  (SHR3
// intentionally differs from the posting; see rng.h.)
TEST(Kiss, ComponentsMatchMarsagliaTestValues) {
  RNG rng(12345, 65435, 34221, 12345);
  for (int i = 0; i < 1000256; ++i) rng.rand_int32();
  uint32_t k = 0;
  for (int i = 0; i < 1000000; ++i) k = rng.CONG();
  EXPECT_EQ(1529210297u, k);
  for (int i = 0; i < 1000000; ++i) k = rng.MWC();
  EXPECT_EQ(904977562u, k);
}

TEST(Kiss, RandIntCombinesTwo32BitDraws) {
  RNG a(1, 2, 3, 4), b(1, 2, 3, 4);
  for (int i = 0; i < 1000; ++i) {
    const uint64_t lo = b.rand_int32();
    const uint64_t hi = b.rand_int32();
    ASSERT_EQ(lo | (hi << 32), a.rand_int());
  }
}

TEST(Seeding, SameSeedSameSequence) {
  RNG a(42), b(42), c(43);
  bool differs = false;
  for (int i = 0; i < 1000; ++i) {
    const uint32_t x = a.rand_int32();
    ASSERT_EQ(x, b.rand_int32());
    differs = differs || (x != c.rand_int32());
  }
  EXPECT_TRUE(differs);
}

TEST(Seeding, InitRestartsSequence) {
  RNG a(7);
  std::vector<uint32_t> first;
  for (int i = 0; i < 100; ++i) first.push_back(a.rand_int32());
  a.init(7);
  for (int i = 0; i < 100; ++i) ASSERT_EQ(first[i], a.rand_int32());
}

TEST(Seeding, DefaultConstructedGeneratorsDiffer) {
  RNG a, b;
  int same = 0;
  for (int i = 0; i < 100; ++i) same += (a.rand_int32() == b.rand_int32());
  EXPECT_LT(same, 5);
}

// With all-zero state the multiply-with-carry and shift-register generators
// would be stuck at zero forever, leaving only the weak congruential one.
TEST(Seeding, DegenerateStateIsRepaired) {
  RNG zero(0, 0, 0, 0);
  EXPECT_NE(0u, zero.MWC());
  EXPECT_NE(0u, zero.SHR3());
  RNG fixed_points(0x9068ffffu, 0x464fffffu, 0, 0);
  const uint32_t first = fixed_points.MWC();
  EXPECT_NE(first, fixed_points.MWC());
  // A seed of zero must still give a healthy stream.
  RNG seeded(0);
  ExpectMoments(SampleMoments(kN, [&] { return seeded.rand_halfclosed01(); }),
                kN, 0.5, 1.0 / 12.0);
}

// Each of the 32 output bits should be set about half the time.
TEST(Kiss, BitsAreBalanced) {
  RNG rng(2024);
  int counts[32] = {};
  for (int i = 0; i < kN; ++i) {
    const uint32_t x = rng.rand_int32();
    for (int b = 0; b < 32; ++b) counts[b] += (x >> b) & 1;
  }
  for (int b = 0; b < 32; ++b) {
    EXPECT_NEAR(counts[b], kN / 2.0, 5.0 * std::sqrt(kN / 4.0)) << "bit " << b;
  }
}

// ___________________________________________________________________________
// Uniform variates

TEST(Uniform, RangesAndMoments) {
  RNG rng(1);
  for (int i = 0; i < kN; ++i) {
    const double closed = rng.rand_closed01();
    ASSERT_TRUE(closed >= 0.0 && closed <= 1.0);
    const double open = rng.rand_open01();
    ASSERT_TRUE(open > 0.0 && open < 1.0);
    const double half_closed = rng.rand_halfclosed01();
    ASSERT_TRUE(half_closed >= 0.0 && half_closed < 1.0);
    const double half_open = rng.rand_halfopen01();
    ASSERT_TRUE(half_open > 0.0 && half_open <= 1.0);
    const int32_t i31 = rng.rand_int31();
    ASSERT_GE(i31, 0);
  }
  const double var = 1.0 / 12.0;
  ExpectMoments(SampleMoments(kN, [&] { return rng.rand_closed01(); }),
                kN, 0.5, var);
  ExpectMoments(SampleMoments(kN, [&] { return rng.rand_open01(); }),
                kN, 0.5, var);
  ExpectMoments(SampleMoments(kN, [&] { return rng.rand_halfclosed01(); }),
                kN, 0.5, var);
  ExpectMoments(SampleMoments(kN, [&] { return rng.rand_halfopen01(); }),
                kN, 0.5, var);
  ExpectMoments(SampleMoments(kN, [&] { return rng.uniform(-2.0, 6.0); }),
                kN, 2.0, 64.0 / 12.0);
}

// The interval endpoints are honoured even for the most extreme raw draws.
// (Dividing a 64-bit integer by 2^64 in floating point would round to 1.0.)
TEST(Uniform, EndpointsAreExact) {
  const double two53 = 9007199254740992.0;
  const uint64_t max = std::numeric_limits<uint64_t>::max();
  // These mirror the conversions in rng.h for the largest raw value.
  EXPECT_LT((static_cast<double>(max >> 11)) / two53, 1.0);
  EXPECT_LT((static_cast<double>(max >> 12) + 0.5) / (two53 / 2), 1.0);
  EXPECT_EQ((static_cast<double>(max >> 11) + 1.0) / two53, 1.0);
  EXPECT_GT((static_cast<double>(uint64_t{0} >> 12) + 0.5) / (two53 / 2), 0.0);
}

// ___________________________________________________________________________
// Continuous distributions

TEST(Normal, Moments) {
  RNG rng(2);
  ExpectMoments(SampleMoments(kN, [&] { return rng.normal(); }), kN, 0.0, 1.0);
  ExpectMoments(SampleMoments(kN, [&] { return rng.normal(3.0, 2.0); }),
                kN, 3.0, 4.0);
}

TEST(Normal, ShapeMatchesNormalCdf) {
  RNG rng(3);
  // Bin edges at -3, -2, ..., 3 and the tail beyond the ziggurat's base
  // strip (|x| > 3.4426), which exercises the rejection path.
  int below[7] = {};
  int tail = 0;
  for (int i = 0; i < kN; ++i) {
    const double x = rng.normal();
    for (int e = 0; e < 7; ++e) below[e] += (x < e - 3.0);
    tail += (std::abs(x) > 3.442619855899);
  }
  for (int e = 0; e < 7; ++e) {
    const double p = 0.5 * std::erfc(-(e - 3.0) / std::sqrt(2.0));
    EXPECT_NEAR(below[e], kN * p, 5.0 * std::sqrt(kN * p * (1 - p)))
        << "P(X < " << e - 3 << ")";
  }
  const double p_tail = std::erfc(3.442619855899 / std::sqrt(2.0));
  EXPECT_NEAR(tail, kN * p_tail, 5.0 * std::sqrt(kN * p_tail));
  EXPECT_GT(tail, 0);
}

TEST(Exponential, MomentsAndTail) {
  RNG rng(4);
  ExpectMoments(SampleMoments(kN, [&] { return rng.exponential(); }),
                kN, 1.0, 1.0);
  ExpectMoments(SampleMoments(kN, [&] { return rng.exponential(4.0); }),
                kN, 0.25, 1.0 / 16.0);
  int beyond_2 = 0;
  for (int i = 0; i < kN; ++i) {
    const double x = rng.exponential();
    ASSERT_GE(x, 0.0);
    beyond_2 += (x > 2.0);
  }
  const double p = std::exp(-2.0);
  EXPECT_NEAR(beyond_2, kN * p, 5.0 * std::sqrt(kN * p * (1 - p)));
}

TEST(Gamma, Moments) {
  RNG rng(5);
  // The second argument is a rate: mean = shape / rate.
  const double cases[][2] = {{0.3, 1.0}, {1.0, 1.0}, {2.5, 1.0}, {7.0, 3.0},
                             {0.8, 0.25}};
  for (const auto& c : cases) {
    const double shape = c[0], rate = c[1];
    ExpectMoments(SampleMoments(kN, [&] { return rng.gamma(shape, rate); }),
                  kN, shape / rate, shape / (rate * rate), 30.0);
  }
}

TEST(Gamma, InvalidParametersGiveNaN) {
  RNG rng(6);
  const double nan = std::numeric_limits<double>::quiet_NaN();
  EXPECT_TRUE(std::isnan(rng.gamma(0.0, 1.0)));
  EXPECT_TRUE(std::isnan(rng.gamma(-1.0, 1.0)));
  EXPECT_TRUE(std::isnan(rng.gamma(1.0, 0.0)));
  EXPECT_TRUE(std::isnan(rng.gamma(nan, 1.0)));  // Used to loop forever.
  EXPECT_TRUE(std::isnan(rng.gamma(1.0, nan)));
}

TEST(ChiSquare, Moments) {
  RNG rng(7);
  ExpectMoments(SampleMoments(kN, [&] { return rng.chi_square(5.0); }),
                kN, 5.0, 10.0, 30.0);
}

TEST(Beta, Moments) {
  RNG rng(8);
  const double a = 2.0, b = 5.0;
  ExpectMoments(SampleMoments(kN, [&] { return rng.beta(a, b); }), kN,
                a / (a + b), a * b / ((a + b) * (a + b) * (a + b + 1)));
}

TEST(Vectors, FillAndSample) {
  RNG a(9), b(9);
  std::vector<double> filled(50);
  a.normal(&filled, 1.0, 2.0);
  for (double x : filled) ASSERT_EQ(x, b.normal(1.0, 2.0));
  EXPECT_EQ(25u, a.normal_sample(25).size());
  EXPECT_EQ(25u, a.exponential_sample(25).size());
  std::vector<int> counts(50, -1);
  a.poisson(&counts, 4.0);
  for (int x : counts) ASSERT_GE(x, 0);
}

// ___________________________________________________________________________
// Discrete distributions

TEST(Poisson, Moments) {
  RNG rng(10);
  // Means below 10 use table lookup; larger means use the normal-based
  // algorithm.
  for (double mu : {0.2, 3.5, 9.99, 10.0, 50.0, 1234.5}) {
    ExpectMoments(SampleMoments(kN, [&] { return rng.poisson(mu); }),
                  kN, mu, mu, 3.0 + 1.0 / mu);
  }
}

TEST(Poisson, SmallMeanMatchesPmf) {
  RNG rng(11);
  const double mu = 2.0;
  int counts[8] = {};
  for (int i = 0; i < kN; ++i) {
    const int x = rng.poisson(mu);
    ASSERT_GE(x, 0);
    if (x < 8) ++counts[x];
  }
  double p = std::exp(-mu);
  for (int k = 0; k < 8; ++k) {
    EXPECT_NEAR(counts[k], kN * p, 5.0 * std::sqrt(kN * p * (1 - p)))
        << "k = " << k;
    p *= mu / (k + 1);
  }
}

TEST(Poisson, EdgeCases) {
  RNG rng(12);
  EXPECT_EQ(0, rng.poisson(0.0));
  EXPECT_EQ(0, rng.poisson(-1.0));
  // NaN used to loop forever.
  EXPECT_EQ(0, rng.poisson(std::numeric_limits<double>::quiet_NaN()));
  // Very large means used to overflow an int.
  const double huge = 2e9;
  for (int i = 0; i < 1000; ++i) {
    EXPECT_NEAR(rng.poisson(huge), huge, 6.0 * std::sqrt(huge));
  }
  EXPECT_EQ(INT_MAX, rng.poisson(1e300));
  EXPECT_EQ(INT_MAX, rng.poisson(std::numeric_limits<double>::infinity()));
}

TEST(Binomial, Moments) {
  RNG rng(13);
  // n * min(p, 1 - p) < 30 uses inversion; otherwise the BTPE algorithm.
  const struct { double p; int n; } cases[] = {
      {0.1, 10}, {0.5, 20}, {0.9, 25}, {0.3, 1000}, {0.5, 100},
      {0.97, 5000}, {0.001, 100000}, {0.25, 1000000}};
  for (const auto& c : cases) {
    const double mean = c.n * c.p, var = c.n * c.p * (1 - c.p);
    int min = INT_MAX, max = INT_MIN;
    const Moments m = SampleMoments(kN, [&] {
      const int x = rng.binomial(c.p, c.n);
      if (x < min) min = x;
      if (x > max) max = x;
      return x;
    });
    ExpectMoments(m, kN, mean, var, 3.0 + 1.0 / var);
    EXPECT_GE(min, 0);
    EXPECT_LE(max, c.n);
  }
}

TEST(Binomial, EdgeCases) {
  RNG rng(14);
  EXPECT_EQ(0, rng.binomial(0.0, 10));
  EXPECT_EQ(0, rng.binomial(-0.5, 10));
  EXPECT_EQ(10, rng.binomial(1.0, 10));
  EXPECT_EQ(10, rng.binomial(1.5, 10));
  EXPECT_EQ(0, rng.binomial(0.5, 0));
  EXPECT_EQ(0, rng.binomial(0.5, -3));
  EXPECT_EQ(0, rng.binomial(std::numeric_limits<double>::quiet_NaN(), 10));
}

// With n near INT_MAX an intermediate product used to overflow an int.
TEST(Binomial, VeryLargeN) {
  RNG rng(15);
  const int n = 2000000000;
  const double p = 0.5;
  const int draws = 100000;
  const Moments m = SampleMoments(draws, [&] {
    const int x = rng.binomial(p, n);
    EXPECT_TRUE(x >= 0 && x <= n);
    return x - n * p;  // Centered to keep the variance numerically stable.
  });
  ExpectMoments(m, draws, 0.0, n * p * (1 - p), 3.0);
}

// The cached setup must be refreshed whenever the parameters change.
TEST(Discrete, InterleavedParameters) {
  RNG rng(16);
  double sums[4] = {};
  const int n = 100000;
  for (int i = 0; i < n; ++i) {
    sums[0] += rng.poisson(3.0);
    sums[1] += rng.poisson(40.0);
    sums[2] += rng.binomial(0.2, 15);
    sums[3] += rng.binomial(0.7, 400);
  }
  EXPECT_NEAR(sums[0] / n, 3.0, 5.0 * std::sqrt(3.0 / n));
  EXPECT_NEAR(sums[1] / n, 40.0, 5.0 * std::sqrt(40.0 / n));
  EXPECT_NEAR(sums[2] / n, 3.0, 5.0 * std::sqrt(2.4 / n));
  EXPECT_NEAR(sums[3] / n, 280.0, 5.0 * std::sqrt(84.0 / n));
}

TEST(Multinomial, CountsSumAndMeans) {
  RNG rng(17);
  const std::vector<double> probs = {1.0, 2.0, 3.0, 4.0};  // Unnormalized.
  const uint n = 100;
  const int reps = 50000;
  std::vector<double> totals(probs.size(), 0.0);
  std::vector<uint> samp;
  for (int r = 0; r < reps; ++r) {
    rng.multinom(n, probs, samp);
    ASSERT_EQ(probs.size(), samp.size());
    uint sum = 0;
    for (size_t i = 0; i < samp.size(); ++i) {
      sum += samp[i];
      totals[i] += samp[i];
    }
    ASSERT_EQ(n, sum);
  }
  for (size_t i = 0; i < probs.size(); ++i) {
    const double p = probs[i] / 10.0;
    EXPECT_NEAR(totals[i] / reps, n * p, 5.0 * std::sqrt(n * p * (1 - p) / reps));
  }
}

// Rounding in the running total (0.6 - 0.1 - 0.2 != 0.3) must never let
// draws leak into a zero-probability category.
TEST(Multinomial, ZeroProbabilityCategoriesStayEmpty) {
  RNG rng(18);
  const std::vector<double> probs = {0.1, 0.0, 0.2, 0.3, 0.0};
  std::vector<uint> samp;
  for (int r = 0; r < 20000; ++r) {
    rng.multinom(1000000, probs, samp);
    ASSERT_EQ(0u, samp[1]);
    ASSERT_EQ(0u, samp[4]);
    ASSERT_EQ(1000000u, samp[0] + samp[2] + samp[3]);
  }
}

TEST(Multinomial, EdgeCases) {
  RNG rng(19);
  std::vector<uint> samp = {7, 7};
  rng.multinom(10, std::vector<double>(), samp);
  EXPECT_TRUE(samp.empty());

  rng.multinom(0, {0.5, 0.5}, samp);
  EXPECT_EQ((std::vector<uint>{0, 0}), samp);

  rng.multinom(10, {0.0, 0.0}, samp);
  EXPECT_EQ((std::vector<uint>{0, 0}), samp);

  rng.multinom(10, {0.0, 1.0, 0.0}, samp);
  EXPECT_EQ((std::vector<uint>{0, 10, 0}), samp);

  // Negative, NaN and infinite probabilities are ignored.
  const double nan = std::numeric_limits<double>::quiet_NaN();
  const double inf = std::numeric_limits<double>::infinity();
  rng.multinom(10, {-1.0, nan, 2.0, inf}, samp);
  EXPECT_EQ((std::vector<uint>{0, 0, 10, 0}), samp);
}

// ___________________________________________________________________________
// Threads

// Separate generators used on separate threads must not interfere with each
// other: each thread has to reproduce the single-threaded result exactly.
TEST(Threads, IndependentGeneratorsDoNotInterfere) {
  const auto work = [](uint64_t seed) {
    RNG rng(seed);
    double sum = 0.0;
    for (int i = 0; i < 20000; ++i) {
      // Different parameters per seed, to stress the cached setup.
      sum += rng.poisson(2.0 + seed) + rng.poisson(20.0 + seed);
      sum += rng.binomial(0.1 + 0.05 * seed, 50) +
             rng.binomial(0.3 + 0.05 * seed, 2000);
      sum += rng.normal() + rng.exponential() + rng.gamma(0.5 + seed);
    }
    return sum;
  };

  const int num_threads = 8;
  std::vector<double> expected(num_threads), actual(num_threads);
  for (int t = 0; t < num_threads; ++t) expected[t] = work(t);

  std::vector<std::thread> threads;
  for (int t = 0; t < num_threads; ++t) {
    threads.emplace_back([&, t] { actual[t] = work(t); });
  }
  for (auto& thread : threads) thread.join();
  for (int t = 0; t < num_threads; ++t) {
    EXPECT_EQ(expected[t], actual[t]) << "thread " << t;
  }
}

TEST(Threads, ConcurrentConstruction) {
  std::vector<std::thread> threads;
  std::vector<double> results(8);
  for (int t = 0; t < 8; ++t) {
    threads.emplace_back([&, t] {
      double sum = 0.0;
      for (int i = 0; i < 2000; ++i) {
        RNG rng;
        sum += rng.normal();
      }
      results[t] = sum;
    });
  }
  for (auto& thread : threads) thread.join();
  for (double r : results) EXPECT_TRUE(std::isfinite(r));
}

}  // namespace
