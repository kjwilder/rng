load("@rules_cc//cc:defs.bzl", "cc_library", "cc_test")

cc_library(
    name = "rng",
    srcs = ["rng.cc"],
    hdrs = ["rng.h"],
)

cc_test(
  name = "rng_test",
  size = "small",
  srcs = ["rng.h", "tests/rng_test.cc"],
  deps = [":rng", "@googletest//:gtest_main"],
)
