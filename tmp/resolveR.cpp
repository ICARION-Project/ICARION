namespace {

// Residual of the Orbitrap hyperlog surface at (r, z):
// z^2 = 0.5 * (r^2 - R^2) + Rm^2 * ln(R / r)
inline double orbitrap_surface_residual(double r, double z, double R, double Rm) {
    const double z2 = z * z;
    return 0.5 * (r * r - R * R) + Rm * Rm * std::log(R / r) - z2;
}

/**
 * @brief Compute the radial coordinate r(z) on an Orbitrap hyperbolic electrode.
 *
 * Implicit equation:
 *   z² = 0.5·(r² - R²) + R_m² · ln(R/r)
 *
 * Solves for r given z, R (electrode radius), and R_m (characteristic radius)
 * using bisection with adaptive bracketing.
 */
double orbitrap_r_for_z(double z, double R, double Rm) {
    const double z_abs = std::fabs(z);
    const double eps   = 1e-10;
    const int max_iter = 80;

    // For z=0, the solution is exactly r=R.
    if (z_abs < eps) {
        return R;
    }

    // Initial bracket around R. In practice the solution is close to R
    // for typical Orbitrap parameters.
    double r_lo = 0.3 * R;
    double r_hi = 3.0 * R;

    double f_lo = orbitrap_surface_residual(r_lo, z_abs, R, Rm);
    double f_hi = orbitrap_surface_residual(r_hi, z_abs, R, Rm);

    // Try to ensure bracketing. Expand a bit if needed.
    int expand_iter = 0;
    while (f_lo * f_hi > 0.0 && expand_iter < 10) {
        r_lo *= 0.5;
        r_hi *= 1.5;
        f_lo = orbitrap_surface_residual(r_lo, z_abs, R, Rm);
        f_hi = orbitrap_surface_residual(r_hi, z_abs, R, Rm);
        ++expand_iter;
    }

    // If still no sign change, fall back to R (safe-ish) to avoid NaNs.
    if (f_lo * f_hi > 0.0) {
        return R;
    }

    double r_mid = 0.5 * (r_lo + r_hi);
    for (int i = 0; i < max_iter; ++i) {
        r_mid = 0.5 * (r_lo + r_hi);
        double f_mid = orbitrap_surface_residual(r_mid, z_abs, R, Rm);
        if (std::fabs(f_mid) < eps) {
            break;
        }
        if (f_mid * f_lo > 0.0) {
            r_lo = r_mid;
            f_lo = f_mid;
        } else {
            r_hi = r_mid;
            f_hi = f_mid;
        }
    }

    return r_mid;
}

// Residual of the surface at a full 3D point (local coordinates):
inline double orbitrap_surface_residual(const Vec3& p_local, double R, double Rm) {
    const double r = std::sqrt(p_local.x * p_local.x + p_local.y * p_local.y);
    const double z = p_local.z;
    if (r <= 0.0) {
        // On axis: treat as outside (no physical electrode there)
        return 1.0;
    }
    return orbitrap_surface_residual(r, z, R, Rm);
}

} // anonymous namespace
