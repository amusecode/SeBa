#include "util_math.h"
#include <gtest/gtest.h>

// basic numerical tolerance for comparisons
static constexpr real EPS = 1e-12;

TEST(UtilMath, Twiddles) {
    EXPECT_TRUE(twiddles(1.0, 1.0 + 1e-13, EPS));
    EXPECT_FALSE(twiddles(1.0, 1.1, EPS));
}

TEST(UtilMath, Sign) {
    EXPECT_EQ(sign(3.5), 1);
    EXPECT_EQ(sign(-2.1), -1);
    EXPECT_EQ(sign(0.0), 0);
}

TEST(UtilMath, AdjustNumberToPower) {
    // result should not exceed max_step_size
    real result = adjust_number_to_power(5.5, 4.0);
    EXPECT_LE(result, 4.0);
}

TEST(UtilMath, AsinhAcoshConsistency) {
    real x = 1.5;
    real s = sinh(x);
    real c = cosh(x);
    EXPECT_NEAR(asinh(s), x, 1e-12);
    EXPECT_NEAR(acosh(c), x, 1e-12);
}
