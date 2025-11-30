// ICARION: Ion Collision And Reaction IntegratiON
// MIT License - Copyright (c) 2025 ICARION Project Contributors

#include "fieldSampling.h"
#include <cmath>
#include <algorithm>

namespace ICARION {
namespace fields {

Vec3d sample_field(const FieldSnapshot& snapshot, const Vec3d& pos) {
    const auto& G = snapshot.grid;
    
    // Convert position to grid coordinates
    double fx = (pos.x - G.origin.x) / G.dx;
    double fy = (pos.y - G.origin.y) / G.dy;
    double fz = (pos.z - G.origin.z) / G.dz;
    
    // Clamp to grid bounds
    int i0 = std::max(0, std::min(static_cast<int>(std::floor(fx)), G.nx - 2));
    int j0 = std::max(0, std::min(static_cast<int>(std::floor(fy)), G.ny - 2));
    int k0 = std::max(0, std::min(static_cast<int>(std::floor(fz)), G.nz - 2));
    
    int i1 = i0 + 1;
    int j1 = j0 + 1;
    int k1 = k0 + 1;
    
    // Fractional parts for interpolation
    double tx = fx - i0;
    double ty = fy - j0;
    double tz = fz - k0;
    
    // Clamp to [0,1]
    tx = std::max(0.0, std::min(1.0, tx));
    ty = std::max(0.0, std::min(1.0, ty));
    tz = std::max(0.0, std::min(1.0, tz));
    
    // Trilinear interpolation weights
    double w000 = (1.0 - tx) * (1.0 - ty) * (1.0 - tz);
    double w100 = tx * (1.0 - ty) * (1.0 - tz);
    double w010 = (1.0 - tx) * ty * (1.0 - tz);
    double w110 = tx * ty * (1.0 - tz);
    double w001 = (1.0 - tx) * (1.0 - ty) * tz;
    double w101 = tx * (1.0 - ty) * tz;
    double w011 = (1.0 - tx) * ty * tz;
    double w111 = tx * ty * tz;
    
    // Helper to get flat index
    auto idx = [&G](int i, int j, int k) -> size_t {
        return static_cast<size_t>((k * G.ny + j) * G.nx + i);
    };
    
    // Sample all 8 corners
    size_t id000 = idx(i0, j0, k0);
    size_t id100 = idx(i1, j0, k0);
    size_t id010 = idx(i0, j1, k0);
    size_t id110 = idx(i1, j1, k0);
    size_t id001 = idx(i0, j0, k1);
    size_t id101 = idx(i1, j0, k1);
    size_t id011 = idx(i0, j1, k1);
    size_t id111 = idx(i1, j1, k1);
    
    // Interpolate each component
    Vec3d E_interp;
    E_interp.x = w000 * G.Ex[id000] + w100 * G.Ex[id100] + w010 * G.Ex[id010] + w110 * G.Ex[id110]
               + w001 * G.Ex[id001] + w101 * G.Ex[id101] + w011 * G.Ex[id011] + w111 * G.Ex[id111];
    
    E_interp.y = w000 * G.Ey[id000] + w100 * G.Ey[id100] + w010 * G.Ey[id010] + w110 * G.Ey[id110]
               + w001 * G.Ey[id001] + w101 * G.Ey[id101] + w011 * G.Ey[id011] + w111 * G.Ey[id111];
    
    E_interp.z = w000 * G.Ez[id000] + w100 * G.Ez[id100] + w010 * G.Ez[id010] + w110 * G.Ez[id110]
               + w001 * G.Ez[id001] + w101 * G.Ez[id001] + w011 * G.Ez[id011] + w111 * G.Ez[id111];
    
    return E_interp;
}

}  // namespace fields
}  // namespace ICARION
