# Tests

## Overview

Unit tests live in the `tests/` directory.
They use **GoogleTest (gtest)**, a lightweight C++ testing framework.
Tests are built and executed by running:

``` sh
make test
```


This compiles all test files (`tests/*.cpp`), links them against project libraries and gtest, and runs them automatically.
Passing tests will print `[ PASSED ]` lines in the terminal.

## Getting GoogleTest Source

This project vendors GoogleTest to make the test build self-contained.
Ensure the submodules are cloned so that the GoogleTest code is included:

```sh
git submodule update --init --recursive
```

This will populate `third_party/googletest/`, which contains the GoogleTest source built automatically when running `make test`.

## Cleaning Up

To remove test binaries:

``` sh
make clean-tests
```

## Writing a Test

Each test file is a small C++ program that includes both gtest and your project
headers.

You want to define test cases using the `TEST` macro provided by gtest, and
ensure that your functions behave as expected using assertions like `EXPECT_EQ`, `ASSERT_TRUE`, etc.

Example:

``` cpp
#include "util_math.h
#include <gtest/gtest.h>

TEST(Basic, Add) { EXPECT_EQ(add(2, 2), 4); }
```

## Running Tests

Run all tests:

``` sh
make test
```

Run a specific test suite or test case:

``` sh
make test_[pattern]
```

E.g., run all tests marked `UtilMath`, run:

``` sh
make test_UtilMath

```

Run a specific test case within a suite:

``` sh
make test_UtilMath.Twiddles
```
