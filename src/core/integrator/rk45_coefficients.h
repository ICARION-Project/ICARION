// ICARION: Ion Collision And Reaction IntegratiON
// MIT License - Copyright (c) 2025 ICARION Project Contributors

#pragma once

namespace ICARION {
namespace integrator {
namespace rk45 {

// Dormand-Prince coefficients exposed for reuse (RK4/5(4)7M)
struct Coefficients {
    // Time fractions c2..c7
    static constexpr double c2 = 1.0 / 5.0;
    static constexpr double c3 = 3.0 / 10.0;
    static constexpr double c4 = 4.0 / 5.0;
    static constexpr double c5 = 8.0 / 9.0;
    static constexpr double c6 = 1.0;
    static constexpr double c7 = 1.0;

    // a coefficients (lower triangular)
    static constexpr double a21 = 1.0 / 5.0;

    static constexpr double a31 = 3.0 / 40.0;
    static constexpr double a32 = 9.0 / 40.0;

    static constexpr double a41 = 44.0 / 45.0;
    static constexpr double a42 = -56.0 / 15.0;
    static constexpr double a43 = 32.0 / 9.0;

    static constexpr double a51 = 19372.0 / 6561.0;
    static constexpr double a52 = -25360.0 / 2187.0;
    static constexpr double a53 = 64448.0 / 6561.0;
    static constexpr double a54 = -212.0 / 729.0;

    static constexpr double a61 = 9017.0 / 3168.0;
    static constexpr double a62 = -355.0 / 33.0;
    static constexpr double a63 = 46732.0 / 5247.0;
    static constexpr double a64 = 49.0 / 176.0;
    static constexpr double a65 = -5103.0 / 18656.0;

    static constexpr double a71 = 35.0 / 384.0;
    static constexpr double a72 = 0.0;
    static constexpr double a73 = 500.0 / 1113.0;
    static constexpr double a74 = 125.0 / 192.0;
    static constexpr double a75 = -2187.0 / 6784.0;
    static constexpr double a76 = 11.0 / 84.0;

    // b4 (4th order) and b5 (5th order) weights
    static constexpr double b4_1 = 5179.0 / 57600.0;
    static constexpr double b4_2 = 0.0;
    static constexpr double b4_3 = 7571.0 / 16695.0;
    static constexpr double b4_4 = 393.0 / 640.0;
    static constexpr double b4_5 = -92097.0 / 339200.0;
    static constexpr double b4_6 = 187.0 / 2100.0;
    static constexpr double b4_7 = 1.0 / 40.0;

    static constexpr double b5_1 = 35.0 / 384.0;
    static constexpr double b5_2 = 0.0;
    static constexpr double b5_3 = 500.0 / 1113.0;
    static constexpr double b5_4 = 125.0 / 192.0;
    static constexpr double b5_5 = -2187.0 / 6784.0;
    static constexpr double b5_6 = 11.0 / 84.0;
    static constexpr double b5_7 = 0.0;
};

}  // namespace rk45
}  // namespace integrator
}  // namespace ICARION
