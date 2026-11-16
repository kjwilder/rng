# RNG
A C++ Random Number Generator Class

## Note
- Most of the code in this repo was written around the year 2000. One might
  want to look elsewhere, perhaps the 
  [GNU Scientific Library](https://www.gnu.org/software/gsl/),
  for a modern C++ RNG.

## About
- This repo provides a fast and easy-to-use C++ random number
  generator class that generates high-quality random variates for
  several common distributions.
  - Continuous: Uniform, Normal, Gamma, Exponential, Chi-Square, Beta
  - Discrete: Binomial, Multinomial, Poisson
- The C++ RNGs in this repo use the simple and fast "KISS" (Keep It
  Simple Stupid) random number generator suggested by George Marsaglia
  in a Usenet posting from 1999.  He described it as "one of my
  favorite generators".  It generates high-quality random numbers by
  combining the results of three simple random number generators
  constructed from completely different algorithms, and it passes
  commonly used tests for randomness.  It does not have the ultra-long
  period of some other generators, but that is a "problem" that can be
  fixed easily.  The period is about 2^123.
- The generator state is four 32-bit words and all arithmetic is done
  on exactly 32 bits, so a given seed produces the same sequence on
  every platform.  One deliberate departure from Marsaglia's posting:
  the shift-register component uses the shifts (13, 17, 5) rather than
  the posted (17, 13, 5), which does not have the full period of
  2^32 - 1.
- The KISS algorithm is only used directly in the function
  `rand_int32`.  `rand_int32` is then used (directly or indirectly) by
  every other function that generates random numbers.  For faster
  random numbers, one can redefine `rand_int32` to return either
  MWC(), CONG(), or SHR3().  The speed will be two to three times
  faster, and the quality of the random numbers should be sufficient
  for many purposes.  The three alternatives are comparable to each
  other in terms of both speed and quality.
- The ziggurat method of Marsaglia is used to generate exponential and
  normal variates.  The method as well as source code can be found in
  the article "The Ziggurat Method for Generating Random Variables" by
  Marsaglia and Tsang, Journal of Statistical Software 5, 2000.
- The method for generating gamma variables appears in "A Simple
  Method for Generating Gamma Variables" by Marsaglia and Tsang, ACM
  Transactions on Mathematical Software, Vol. 26, No 3, Sep 2000,
  pages 363-372.

## Seeding
- `RNG x;` seeds from the system entropy source, the clock and a
  global counter, so every default-constructed generator is different.
- `RNG x(seed);` is deterministic: the 64-bit seed is hashed into the
  four state words, so nearby seeds give unrelated sequences.
- `RNG x(z, w, jsr, jcong);` sets the four state words directly.
  Values that would leave a component generator stuck (such as 0) are
  replaced by fixed defaults.

## Parameters
- `gamma(shape, rate)` takes a *rate* (inverse scale): the mean is
  `shape / rate`.  `exponential(lambda)` likewise has mean `1 / lambda`.
- Invalid parameters do not hang or crash: `gamma` returns NaN, and
  `poisson` and `binomial` return 0.  `poisson` results are capped at
  `INT_MAX`.

## Warnings
- Separate `RNG` objects can be used from separate threads.  A single
  `RNG` object must not be shared between threads without a lock.
- This is not a cryptographic generator.  Do not use it for anything
  security related.
- It is best NOT to use the generator as in the following line:
  ```
  for (int i = 0; i < 100; ++i) { RNG x; cout << x.uniform() << endl; }
  ```
  The concern is that each time through the loop, a new RNG 'x' is
  created, and that RNG is used to generate exactly one random number.
  The results may be satisfactory, but the class was designed to
  produce quality random numbers by calling a single RNG
  repeatedly.  The better way to write the above loop is:
  ```
  RNG x; for (int i = 0; i < 100; ++i) { cout << x.uniform() << endl; }
  ```

## Usage
- Include the header file `rng.h`
- Link C++ code with `rng.cc` (C++11 or later)
  ```
  #include "rng.h"

  rng_h::RNG x(12345);
  double u = x.uniform();          // [0, 1]
  double n = x.normal(10.0, 2.0);  // mean 10, standard deviation 2
  int k = x.poisson(3.5);
  ```

## Tests
- `bazel test :all`
