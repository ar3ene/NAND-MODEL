#include <math.h>

// Scaled complementary error function.
// Simple implementation with asymptotic handling for large x.
double erfcx(double x) {
    if (x > 25.0) {
        // Asymptotic expansion for large x.
        double inv = 1.0 / x;
        double inv2 = inv * inv;
        return inv / sqrt(M_PI) * (1.0 + 0.5 * inv2);
    }
    return exp(x * x) * erfc(x);
}

