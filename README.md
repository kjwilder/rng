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
- The KISS algorithm is only used directly in the function
  `rand_int32`.  `rand_int32` is then used (directly or indirectly) by
  every other function that generates random numbers.  For faster
  random numbers, one can redefine `rand_int32` to return either
  WMC(), CONG(), or SHR3().  The speed will be two to three times
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

## Warnings
- My original intent was to use a single instance of a random number
  generator class on a platform with 32 bit ints. I made some effort
  to support 64 bit ints and multiple instances of the class, but they
  are not well tested and may have issues.
- This class has a potential race condition if used in a multi-threaded
  application. This can result in suspect random samples.
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
- Link C++ code with `rng.cc`
