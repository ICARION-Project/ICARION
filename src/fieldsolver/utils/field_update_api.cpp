// ICARION: Ion Collision And Reaction IntegratiON
// MIT License - Copyright (c) 2025 ICARION Project Contributors

// field_update_api.cpp
#include "field_update_api.h"
#include <cmath>

FieldServer::FieldServer() {}

void FieldServer::set_geometry(const GeometryConfig& g){
    std::lock_guard<std::mutex> lk(mtx_);
    geo_ = g;
    dirty_.store(true, std::memory_order_relaxed);
}

void FieldServer::set_grid(int nx, int ny, int nz, double dx, double dy, double dz, Vec3d origin){
    std::lock_guard<std::mutex> lk(mtx_);
    snapshot_.grid.nx = nx; snapshot_.grid.ny = ny; snapshot_.grid.nz = nz;
    snapshot_.grid.dx = dx; snapshot_.grid.dy = dy; snapshot_.grid.dz = dz;
    snapshot_.grid.origin = origin;
    snapshot_.grid.Ex.assign(nx*ny*nz, 0.0);
    snapshot_.grid.Ey.assign(nx*ny*nz, 0.0);
    snapshot_.grid.Ez.assign(nx*ny*nz, 0.0);
    dirty_.store(true, std::memory_order_relaxed);
}

void FieldServer::compute_locked_(){
    // Placeholder: synthesize a simple quadrupole-like field into grid
    auto& G = snapshot_.grid;
    const double R = std::max(1e-6, geo_.radius_m);
    for (int k=0; k<G.nz; ++k){
        for (int j=0; j<G.ny; ++j){
            for (int i=0; i<G.nx; ++i){
                const int id = (k*G.ny + j)*G.nx + i;
                double x = G.origin.x + i*G.dx;
                double y = G.origin.y + j*G.dy;
                // very rough toy model: Ex = +2(Vdc)x/R^2, Ey = -2(Vdc)y/R^2, Ez = 0
                double Vdc = geo_.quad_dc_V;
                double fac = (R>0.0) ? (2.0*Vdc/(R*R)) : 0.0;
                G.Ex[id] =  fac * x;
                G.Ey[id] = -fac * y;
                G.Ez[id] =  0.0;
            }
        }
    }
    snapshot_.timestamp = std::chrono::steady_clock::now();
    snapshot_.version = ++version_;
    dirty_.store(false, std::memory_order_release);
}

FieldSnapshot FieldServer::recompute_now(){
    std::lock_guard<std::mutex> lk(mtx_);
    compute_locked_();
    return snapshot_;
}

FieldSnapshot FieldServer::get_snapshot() const {
    std::lock_guard<std::mutex> lk(mtx_);
    return snapshot_;
}

bool FieldServer::has_newer(int known_version) const {
    return version_.load(std::memory_order_acquire) > known_version;
}
