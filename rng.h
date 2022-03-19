#ifndef RNG_H_
#define RNG_H_

#include <climits>
#include <cmath>
#include <cstdlib>
#include <cstdint>
#include <ctime>
#include <vector>

namespace rng_h {

#if UINT_FAST32_MAX == 4294967295ul
inline uint_fast32_t ULONG32(int_fast32_t x) {
  return static_cast<uint_fast32_t>(x);
}
inline uint_fast32_t ULONG32(uint_fast32_t x) {
  return x;
}
inline uint_fast32_t ULONG32(double x) {
  return static_cast<uint_fast32_t>(x);
}
inline int_fast32_t UL32toSL32(uint_fast32_t x) {
  return static_cast<int_fast32_t>(x);
}
#elif UINT_FAST32_MAX == 18446744073709551615ul
inline uint_fast32_t ULONG32(int_fast32_t x) {
  return static_cast<uint_fast32_t>(x) & 0xfffffffful;
}
inline uint_fast32_t ULONG32(uint_fast32_t x) {
  return x & 0xfffffffful;
}
inline uint_fast32_t ULONG32(double x) {
  return static_cast<uint_fast32_t>(x) & 0xfffffffful;
}
inline int_fast32_t UL32toSL32(uint_fast32_t x) {
  return (x < 0x80000000ul) ?  static_cast<int_fast32_t>(x) :
    -1 * (0x80000000ul - (x & 0x7ffffffful));
}
#else
#error "Unsupported uint_fast_32_t"
#endif

class RNG {
 private:
  uint_fast32_t z, w, jsr, jcong;  // Seeds

  static uint_fast32_t tm;  // Different seeds for differnet RNGs.
  static uint_fast32_t kn[128], ke[256];
  static double wn[128], fn[128], we[256], fe[256];

 public:
  RNG() { init(); zigset(); }
  explicit RNG(uint_fast32_t x_) :
    z(x_), w(x_), jsr(x_), jcong(x_) { zigset(); }
  RNG(uint_fast32_t z_, uint_fast32_t w_,
      uint_fast32_t jsr_, uint_fast32_t jcong_) :
    z(z_), w(w_), jsr(jsr_), jcong(jcong_) { zigset(); }
  ~RNG() { }

#if UINT_FAST32_MAX == 4294967295ul
  // 32 bit unsigned longs
  uint_fast32_t znew()
    { return (z = 36969 * (z & 0xfffful) + (z >> 16)); }
  uint_fast32_t wnew()
    { return (w = 18000 * (w & 0xfffful) + (w >> 16)); }
  uint_fast32_t MWC()
    { return (znew() << 16) + wnew(); }
  uint_fast32_t SHR3()
    { jsr ^= (jsr << 17); jsr ^= (jsr >> 13); return jsr ^= (jsr << 5); }
  uint_fast32_t CONG()
    { return (jcong = 69069 * jcong + 1234567); }
  uint_fast32_t rand_int32()       // [0,2^32-1]
    { return (MWC() ^ CONG()) + SHR3(); }
  uint_fast32_t rand_int()         // [0,2^32-1]
    { return (MWC() ^ CONG()) + SHR3(); }
#elif UINT_FAST32_MAX == 18446744073709551615ul
#ifdef RNG_C
#warning "Compiling RNG class for 64-bit architecture"
#endif
  // 64-bit unsigned longs
  // This is not as elegant and fast as it could be, but it works.
  uint_fast32_t znew()
    { return z = ((36969 * (z & 0xfffful) + (z >> 16)) & 0xfffffffful); }
  uint_fast32_t wnew()
    { return w = ((18000 * (w & 0xfffful) + (w >> 16)) & 0xfffffffful); }
  uint_fast32_t MWC()
    { return ((znew() << 16) + wnew()) & 0xfffffffful; }
  uint_fast32_t SHR3()
    { jsr ^= (jsr << 17); jsr ^= (jsr >> 13);
      return (jsr = ((jsr ^= (jsr << 5)) & 0xfffffffful)); }
  uint_fast32_t CONG()
    { return jcong = ((69069 * jcong + 1234567) & 0xfffffffful); }
  uint_fast32_t rand_int32()         // [0,2^32-1]
    { return ((MWC() ^ CONG()) + SHR3()) & 0xfffffffful; }
  uint_fast32_t rand_int()           // [0,2^64-1]
    { return rand_int32() | (rand_int32() << 32); }
#endif
  double RNOR() {
    int_fast32_t h = UL32toSL32(rand_int32()), i = h & 127;
    return ((uint_fast32_t)std::abs(h) < kn[i]) ? h * wn[i] : nfix(h, i);
  }
  double REXP() {
    uint_fast32_t j = rand_int32(), i = j & 255;
    return (j < ke[i]) ? j * we[i] : efix(j, i);
  }

  double nfix(int_fast32_t h, uint_fast32_t i);
  double efix(uint_fast32_t j, uint_fast32_t i);
  void zigset();

  void init()
    { z = w = jsr = jcong = static_cast<uint_fast32_t>(time(0)) + tm;
      tm += 123457; }
  void init(uint_fast32_t z_, uint_fast32_t w_,
      uint_fast32_t jsr_, uint_fast32_t jcong_ )
    { z = z_; w = w_; jsr = jsr_; jcong = jcong_; }

  // For a faster but lower quality RNG, uncomment the following
  // line, and comment out the original definition of rand_int above.
  // In practice, the faster RNG is fine for simulations
  // that do not simulate more than a few billion random numbers.
  // uint_fast32_t rand_int() { return SHR3(); }
  int_fast32_t rand_int31() {  // [0,2^31-1]
    return static_cast<int_fast32_t>(rand_int32()) >> 1;
  }
  double rand_closed01() {  // [0,1]
    return static_cast<double>(rand_int()) /
      static_cast<double>(UINT_FAST32_MAX);
  }
  double rand_open01() {  // (0,1)
    return (static_cast<double>(rand_int()) + 1.0) / (UINT_FAST32_MAX + 2.0);
  }
  double rand_halfclosed01() {  // [0,1)
    return static_cast<double>(rand_int()) / (UINT_FAST32_MAX + 1.0);
  }
  double rand_halfopen01() {  // (0,1]
    return (static_cast<double>(rand_int()) + 1.0) / (UINT_FAST32_MAX + 1.0);
  }

  // Continuous Distributions
  double uniform(double x = 0.0, double y = 1.0)
    { return rand_closed01() * (y - x) + x; }
  double normal(double mu = 0.0, double sd = 1.0)
    { return RNOR() * sd + mu; }
  double exponential(double lambda = 1)
    { return REXP() / lambda; }
  double gamma(double shape = 1, double scale = 1);
  double chi_square(double df)
    { return gamma(df / 2.0, 0.5); }
  double beta(double a1, double a2)
    { const double x1 = gamma(a1, 1); return x1 / (x1 + gamma(a2, 1)); }

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
  void gamma(std::vector<double>* res, double shape = 1, double scale = 1) {
    for (auto i = res->begin(); i != res->end(); ++i)
      *i = gamma(shape, scale);
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
  int poisson(double mu);
  int binomial(double p, int n);
  void multinom(unsigned int n, const std::vector<double>* probs,
                std::vector<unsigned int>* samp);
  void multinom(unsigned int n, const double* prob,
      unsigned int K, unsigned int* samp);

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
