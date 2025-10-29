#include <gtest/gtest.h>

int add(int a, int b) { return a + b; }

TEST(Basic, Add) { EXPECT_EQ(add(2, 2), 4); }
