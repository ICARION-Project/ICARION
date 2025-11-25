void DomainManager::terminate_ion_at_boundary(IonState& ion, int domain_idx,
                                              const Vec3& pos_before_local,
                                              const Vec3& pos_after_local) const {
    const auto& dom = get_domain(domain_idx);

    // ─────────────────────────────────────────────────────────────
    // Orbitrap: Hyperlogarithmic electrode boundary (inner & outer)
    // ─────────────────────────────────────────────────────────────
    if (dom.instrument == config::Instrument::Orbitrap) {
        Vec3 dir = pos_after_local - pos_before_local;
        double step_len = std::sqrt(dir.x * dir.x + dir.y * dir.y + dir.z * dir.z);
        if (step_len < NUMERICAL_ZERO) {
            // Degenerate step: no movement
            ion.pos = local_to_global_pos(pos_before_local, domain_idx);
            ion.vel = {0.0, 0.0, 0.0};
            ion.active = false;
            return;
        }
        dir = dir * (1.0 / step_len); // normalize

        const double Rin = dom.geometry.radius_in_m;
        const double Rout = dom.geometry.radius_out_m;
        const double Rm   = dom.geometry.radius_char_m;
        const double z_max = 0.5 * dom.geometry.length_m;

        double t_min = 1.0;                 // param in [0,1] along the segment
        Vec3  intersection = pos_after_local;
        bool  hit = false;

        auto point_on_segment = [&](double t) {
            return pos_before_local + dir * (t * step_len);
        };

        // Helper: find intersection with one electrode surface defined by (R, Rm)
        auto intersect_surface = [&](double R) {
            const double f0 = orbitrap_surface_residual(pos_before_local, R, Rm);
            const double f1 = orbitrap_surface_residual(pos_after_local,  R, Rm);

            // If both on same side and not exactly on surface, no crossing.
            if (f0 * f1 > 0.0) {
                return;
            }

            double t_lo = 0.0;
            double t_hi = 1.0;
            double f_lo = f0;
            double f_hi = f1;
            const double eps = 1e-8;
            const int max_iter = 40;

            // Bisection along the segment in t
            for (int i = 0; i < max_iter; ++i) {
                double t_mid = 0.5 * (t_lo + t_hi);
                Vec3 p_mid = point_on_segment(t_mid);
                double f_mid = orbitrap_surface_residual(p_mid, R, Rm);

                if (std::fabs(f_mid) < eps) {
                    // Accept solution
                    if (t_mid < t_min) {
                        t_min = t_mid;
                        intersection = p_mid;
                        hit = true;
                    }
                    return;
                }

                if (f_mid * f_lo > 0.0) {
                    t_lo = t_mid;
                    f_lo = f_mid;
                } else {
                    t_hi = t_mid;
                    f_hi = f_mid;
                }
            }

            // After max_iter: take midpoint as best guess
            double t_mid = 0.5 * (t_lo + t_hi);
            if (t_mid < t_min) {
                t_min = t_mid;
                intersection = point_on_segment(t_mid);
                hit = true;
            }
        };

        // Check both inner and outer electrodes
        intersect_surface(Rin);
        intersect_surface(Rout);

        // Optional: axial clipping at |z| = z_max
        auto intersect_axial_plane = [&](double z_plane) {
            if (std::fabs(dir.z) < NUMERICAL_ZERO) {
                return;
            }
            double t_plane = (z_plane - pos_before_local.z) / (dir.z * step_len);
            if (t_plane <= 0.0 || t_plane >= 1.0) {
                return;
            }
            Vec3 p = point_on_segment(t_plane);
            // At that z, we also require that the point is between electrodes;
            // otherwise this "intersection" is outside physical domain.
            double r = std::sqrt(p.x * p.x + p.y * p.y);
            double r_in_allowed  = orbitrap_r_for_z(std::fabs(z_plane), Rin,  Rm);
            double r_out_allowed = orbitrap_r_for_z(std::fabs(z_plane), Rout, Rm);
            if (r >= r_in_allowed - DOMAIN_BOUNDARY_EPSILON &&
                r <= r_out_allowed + DOMAIN_BOUNDARY_EPSILON) {
                if (t_plane < t_min) {
                    t_min = t_plane;
                    intersection = p;
                    hit = true;
                }
            }
        };

        // Clip at axial edges of simulated domain
        intersect_axial_plane(+z_max);
        intersect_axial_plane(-z_max);

        if (!hit) {
            // No intersection found along the segment: just stop at end.
            intersection = pos_after_local;
        }

        ion.pos = local_to_global_pos(intersection, domain_idx);
        ion.vel = {0.0, 0.0, 0.0};
        ion.active = false;
        return;
    }

    // ─────────────────────────────────────────────────────────────
    // Default cylindrical ray-tracing branch (unchanged)
    // ─────────────────────────────────────────────────────────────
    Vec3 dir = pos_after_local - pos_before_local;
    double dir_len = std::sqrt(dir.x*dir.x + dir.y*dir.y + dir.z*dir.z);
    if (dir_len < NUMERICAL_ZERO) {
        ion.pos = local_to_global_pos(pos_before_local, domain_idx);
        ion.vel = {0.0, 0.0, 0.0};
        ion.active = false;
        return;
    }
    dir = dir * (1.0 / dir_len);

    double t_min = dir_len;
    Vec3 intersection = pos_after_local;

    // ... (dein bestehender zylindrischer Code wie gehabt)
}
