#include "util_math.h"
#include <gtest/gtest.h>

int add(int a, int b) { return a + b; }

TEST(Math, Add) { EXPECT_EQ(add(2, 2), 4); }

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
