bool DomainManager::is_inside_domain(const config::DomainConfig& dom,
                                     const Vec3& globalPos) const {
    // Transform to local coordinates
    Vec3 local = dom.rotation_global_to_local * (globalPos - dom.geometry.origin_m);
    double r   = std::sqrt(local.x * local.x + local.y * local.y);

    // Non-Orbitrap: standard cylindrical domain
    if (dom.instrument != config::Instrument::Orbitrap) {
        return (local.z >= -DOMAIN_BOUNDARY_EPSILON &&
                local.z <  dom.geometry.length_m) &&
               (r < dom.geometry.radius_m);
    }

    // ─────────────────────────────────────────────────────────────
    // Orbitrap domain: hyperlogarithmic electrodes
    // ─────────────────────────────────────────────────────────────
    const double Rin = dom.geometry.radius_in_m;      // inner electrode radius
    const double Rout = dom.geometry.radius_out_m;    // outer electrode radius
    const double Rm  = dom.geometry.radius_char_m;    // characteristic radius

    // Optional axial window: finite simulation length around z=0
    const double z_max = 0.5 * dom.geometry.length_m;
    const double z_abs = std::fabs(local.z);

    // If we are beyond the simulated axial range, we are outside.
    if (z_abs > z_max + DOMAIN_BOUNDARY_EPSILON) {
        return false;
    }

    // Compute allowed radial corridor at this |z| from the analytical geometry.
    const double r_in_allowed  = orbitrap_r_for_z(z_abs, Rin,  Rm);
    const double r_out_allowed = orbitrap_r_for_z(z_abs, Rout, Rm);

    // Safety: if helper returned nonsense, treat as outside
    if (!(r_in_allowed > 0.0 && r_out_allowed > r_in_allowed)) {
        return false;
    }

    // Check if ion is between inner and outer electrode at this z.
    return (r >= r_in_allowed - DOMAIN_BOUNDARY_EPSILON) &&
           (r <= r_out_allowed + DOMAIN_BOUNDARY_EPSILON);
}