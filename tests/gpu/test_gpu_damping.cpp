// SPDX-License-Identifier: GPL-3.0-only
// Copyright (c) 2026 Christoph Schaefer

#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>
#include <cmath>
#include <vector>
#include "core/gpu/core/GPUContext.h"
#include "core/gpu/kernels/damping.cuh"
#include "core/types/gpu/IonState_GPU.h"

using icarion::gpu::DeviceDamping;
using icarion::gpu::GPUContext;
using icarion::gpu::IonStateGPU;
using icarion::gpu::ion_state_conversion::upload_ions;
using icarion::gpu::ion_state_conversion::download_ions;
using icarion::gpu::launch_damping_kernel;
using ICARION::core::IonState;
using Catch::Approx;

TEST_CASE("GPU damping kernel matches analytic solution", "[gpu][damping]") {
    auto ctx = GPUContext::create(0);
    REQUIRE(ctx);

    IonState ion{};
    ion.pos = {0.0, 0.0, 0.0};
    ion.vel = {1.0, 0.0, 0.0};
    ion.mass_kg = 1.0;
    ion.ion_charge_C = 1.0;
    ion.active = true;
    ion.history_index = 0;
    ion.current_domain_index = 0;

    std::vector<IonState> ions = {ion};

    IonStateGPU ions_gpu;
    ions_gpu.allocate(ions.size());
    upload_ions(ions, ions_gpu, ctx->get_stream());

    DeviceDamping damping;
    damping.enabled = 1;
    damping.nu_const = 2.0f;  // 1/s
    damping.nu_per_ion = nullptr;

    double dt = 0.1;  // s
    launch_damping_kernel(ions_gpu, damping, dt, ctx->get_stream());

    download_ions(ions_gpu, ions, ctx->get_stream());
    ctx->synchronize();

    double expected = std::exp(-damping.nu_const * dt);
    REQUIRE(ions[0].vel.x == Approx(expected).margin(1e-6));
    REQUIRE(ions[0].vel.y == Approx(0.0).margin(1e-12));
    REQUIRE(ions[0].vel.z == Approx(0.0).margin(1e-12));

    ions_gpu.free();
}
