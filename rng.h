#ifndef RNG_H_
#define RNG_H_

#include <atomic>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <vector>

namespace rng_h {

constexpr double PI = 3.1415926535897932;

typedef unsigned int uint;

// A KISS ("Keep It Simple Stupid") random number generator, after
// Marsaglia (1999), plus samplers for several common distributions.
//
// Thread safety: distinct RNG objects may be used concurrently from
// different threads.  A single RNG object must not be shared between
// threads without external synchronization.
class RNG {
 private:
  // Generator state.  All arithmetic is done on exactly 32 bits so that a
  // given set of seeds yields the same sequence on every platform.
  uint32_t z, w, jsr, jcong;

  // Mixed into the seed of default-constructed RNGs so they all differ.
  static std::atomic<uint64_t> tm;

  // Ziggurat tables, filled once by zigset() and read-only afterwards.
  static uint32_t kn[128], ke[256];
  static double wn[128], fn[128], we[256], fe[256];

  // Per-instance cached setup for poisson(), reused while mu is unchanged.
  struct PoissonState {
    int l = 0, m = 0;
    double b1 = 0, b2 = 0, c = 0, c0 = 0, c1 = 0, c2 = 0, c3 = 0;
    double pp[36] = {}, p0 = 0, p = 0, q = 0, s = 0, d = 0, omega = 0;
    double big_l = 0;
    double muprev = 0.0, muprev2 = 0.0;
  } poisson_state_;

  // Per-instance cached setup for binomial(), reused while (p, n) is
  // unchanged.
  struct BinomialState {
    double c = 0, fm = 0, npq = 0, p1 = 0, p2 = 0, p3 = 0, p4 = 0;
    double qn = 0, xl = 0, xll = 0, xlr = 0, xm = 0, xr = 0;
    double psave = -1.0;
    int nsave = -1, m = 0;
  } binomial_state_;

  double nfix(int32_t h, uint32_t i);
  double efix(uint32_t j, uint32_t i);
  static void zigset();

 public:
  // Seeds from the system entropy source, the clock and a global counter.
  RNG() { zigset(); init(); }
  // Seeds deterministically: equal seeds give equal sequences everywhere.
  explicit RNG(uint64_t seed) { zigset(); init(seed); }
  // Sets the four KISS state words directly.  Values for which a component
  // generator would be stuck (e.g. 0) are replaced by fixed defaults.
  RNG(uint32_t z_, uint32_t w_, uint32_t jsr_, uint32_t jcong_) {
    zigset();
    init(z_, w_, jsr_, jcong_);
  }

  void init();
  void init(uint64_t seed);
  void init(uint32_t z_, uint32_t w_, uint32_t jsr_, uint32_t jcong_);

  // The component generators.  Each returns a value in [0,2^32-1].
  uint32_t znew()
    { return (z = 36969 * (z & 0xffffu) + (z >> 16)); }
  uint32_t wnew()
    { return (w = 18000 * (w & 0xffffu) + (w >> 16)); }
  uint32_t MWC()
    { return ((znew() << 16) + wnew()); }
  // Marsaglia's 1999 posting used the shifts (17, 13, 5), which do not give
  // the full period; (13, 17, 5) is the corrected triple with period 2^32-1.
  uint32_t SHR3()
    { jsr ^= (jsr << 13); jsr ^= (jsr >> 17); return (jsr ^= (jsr << 5)); }
  uint32_t CONG()
    { return (jcong = 69069 * jcong + 1234567); }

  uint32_t rand_int32()      // [0,2^32-1]
    { return ((MWC() ^ CONG()) + SHR3()); }
  uint64_t rand_int() {      // [0,2^64-1]
    const uint64_t lo = rand_int32();
    return (lo | (static_cast<uint64_t>(rand_int32()) << 32));
  }

  // Standard normal and standard exponential variates (ziggurat method).
  double RNOR() {
    const uint32_t u = rand_int32();
    const int32_t h = static_cast<int32_t>(u);
    const uint32_t i = u & 127;
    const uint32_t abs_h = (h < 0) ? (0u - u) : u;
    return ((abs_h < kn[i]) ? h * wn[i] : nfix(h, i));
  }
  double REXP() {
    const uint32_t j = rand_int32(), i = j & 255;
    return ((j < ke[i]) ? j * we[i] : efix(j, i));
  }

  // For a faster but lower quality RNG, uncomment the following
  // line, and comment out the original definition of rand_int32 above.
  // In practice, the faster RNG is fine for simulations
  // that do not simulate more than a few billion random numbers.
  // uint32_t rand_int32() { return SHR3(); }
  int32_t rand_int31()        // [0,2^31-1]
    { return static_cast<int32_t>(rand_int32() >> 1); }
  // The uniform generators below use the top 52 or 53 bits of rand_int(),
  // so every value is exactly representable and the stated endpoints are
  // honoured (no rounding up to 1.0 or down to 0.0).
  double rand_closed01()      // [0,1]
    { return (static_cast<double>(rand_int() >> 11) / 9007199254740991.0); }
  double rand_open01()        // (0,1)
    { return ((static_cast<double>(rand_int() >> 12) + 0.5) *
        (1.0 / 4503599627370496.0)); }
  double rand_halfclosed01()  // [0,1)
    { return (static_cast<double>(rand_int() >> 11) *
        (1.0 / 9007199254740992.0)); }
  double rand_halfopen01()    // (0,1]
    { return ((static_cast<double>(rand_int() >> 11) + 1.0) *
        (1.0 / 9007199254740992.0)); }

  // Continuous Distributions
  double uniform(double x = 0.0, double y = 1.0)
    { return rand_closed01() * (y - x) + x; }
  double normal(double mu = 0.0, double sd = 1.0)
    { return RNOR() * sd + mu; }
  double exponential(double lambda = 1)
    { return REXP() / lambda; }
  // Note: the second parameter is a rate (inverse scale); the mean is
  // shape / rate.  Returns NaN unless shape > 0 and rate > 0.
  double gamma(double shape = 1, double rate = 1);
  double chi_square(double df)
    { return gamma(df / 2.0, 0.5); }
  double beta(double a1, double a2)
    { const double x1 = gamma(a1, 1); return (x1 / (x1 + gamma(a2, 1))); }

  void uniform(std::vector<double>* res, double x = 0.0, double y = 1.0) {
    for (auto i = res->begin(); i != res->end(); ++i)
      *i = uniform(x, y);
  }
  void normal(std::vector<double>* res, double mu = 0.0, double sd = 1.0) {
    for (auto i = res->begin(); i != res->end(); ++i)
      *i = normal(mu, sd);
  }
  std::vector<double> normal_sample(
      size_t size, double mu = 0.0, double sd = 1.0) {
    std::vector<double> res(size);
    normal(&res, mu, sd);
    return res;
  }
  void exponential(std::vector<double>* res, double lambda = 1) {
    for (auto i = res->begin(); i != res->end(); ++i)
      *i = exponential(lambda);
  }
  std::vector<double> exponential_sample(size_t size, double lambda = 1) {
    std::vector<double> res(size);
    exponential(&res, lambda);
    return res;
  }
  void gamma(std::vector<double>* res, double shape = 1, double rate = 1) {
    for (auto i = res->begin(); i != res->end(); ++i)
      *i = gamma(shape, rate);
  }
  void chi_square(std::vector<double>* res, double df) {
    for (auto i = res->begin(); i != res->end(); ++i)
      *i = chi_square(df);
  }
  void beta(std::vector<double>* res, double a1, double a2) {
    for (auto i = res->begin(); i != res->end(); ++i)
      *i = beta(a1, a2);
  }

  // Discrete Distributions
  // poisson() returns 0 unless mu > 0.  binomial() returns 0 unless p > 0
  // and n > 0.  multinom() ignores negative and non-finite probabilities,
  // normalizes the rest, and requires n <= INT_MAX.
  int poisson(double mu);
  int binomial(double p, int n);
  void multinom(unsigned int n, const std::vector<double>& probs,
                std::vector<uint>& samp);
  void multinom(unsigned int n, const double* prob, uint K, uint* samp);

  void poisson(std::vector<int>* res, double lambda) {
    for (auto i = res->begin(); i != res->end(); ++i)
      *i = poisson(lambda);
  }
  void binomial(std::vector<int>* res, double p, int n) {
    for (auto i = res->begin(); i != res->end(); ++i)
      *i = binomial(p, n);
  }
};  // class RNG

}  // namespace rng_h

#endif  // RNG_H_
