#include "util_math.h"
#include <gtest/gtest.h>

int add(int a, int b) { return a + b; }

TEST(Math, Add) { EXPECT_EQ(add(2, 2), 4); }
