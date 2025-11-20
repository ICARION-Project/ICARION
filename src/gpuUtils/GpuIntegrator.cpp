// SPDX-License-Identifier: MIT
// SPDX-FileCopyrightText: 2025 ICARION Project Contributors

/**
 * =====================================================================
 *
 *   Ion Collision And Reaction IntegratiON (ICARION)
 *   ------------------------------------------------
 *   Modular framework for simulating ion trajectories in custom
 *   electric fields and background gas environments.
 *
 *   @file       GpuIntegrator.cpp
 *   @brief      GPU-accelerated RK4 integrator implementation
 *
 *   @details
 *   Implements a GPU-accelerated RK4 integrator that mirrors the CPU integrator interface.
 *
 *   @date       2025-11-14
 *   @version    1.0.0
 *   @authors    ICARION Development Team
 *
 * =====================================================================
 */
 
 #include "gpuUtils/GpuIntegrator.h"

#ifdef USE_GPU_ACCEL

#include <cuda_runtime.h>
#include <cuda_runtime_api.h>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <stdexcept>
#include <string>

#include "core/constants/constants.h"
#include "instrument/NeutralGasUtils.h"
#include "gpuUtils/integrate_rk4_kernel.cuh"
#include "gpuUtils/integrate_rk4_optimized.cuh"
// DISABLED: Experimental parity kernel
// #include "gpuUtils/integrate_rk4_parity_kernel.cuh"
#include "paramUtils/InstrumentEnums.h"
#include "core/debug/Debug.h"

// Optional debug harness (C linkage, defined in gpu_debug_single.cu)
extern "C" int gpu_debug_single_accel(const IonStateGPU*, const DomainGPU*, ::ICARION::core::Vec3*, ::ICARION::core::Vec3*);


namespace ICARION {
namespace gpu {

using core::Vec3;



bool gpuDebugEnabled() {
    static bool enabled = [] {
        const char* env = std::getenv("ICARION_DEBUG_GPU");
        if (!env || env[0] == '\0') return false;
        std::string s(env);
        for (auto &c : s) c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
        if (s == "1" || s == "true" || s == "yes" || s == "on") return true;
        // allow numeric > 0
        char* endptr = nullptr;
        long v = std::strtol(env, &endptr, 10);
        if (endptr != env) return v > 0;
        return false;
    }();
    return enabled;
}

namespace {

InstrumentGPU parseInstrument(const Json::Value& config) {
    std::string instrument = "IMS";
    if (config.isMember("instrument")) {
        if (config["instrument"].isString()) {
            instrument = config["instrument"].asString();
        } else if (config["instrument"].isObject() && config["instrument"].isMember("type")) {
            instrument = config["instrument"]["type"].asString();
        }
    }

    std::string key = instrument;
    std::transform(key.begin(), key.end(), key.begin(),
                   [](unsigned char c) { return static_cast<char>(std::toupper(c)); });

    if (key == "LQIT") return InstrumentGPU::LQIT;
    if (key == "IMS" || key == "SIFDT_MS" || key == "SIFDT" || key == "SIFDT-MS") return InstrumentGPU::IMS;
    if (key == "ORBITRAP") return InstrumentGPU::Orbitrap;
    if (key == "QUADRUPOLE" || key == "QUADRUPOLE_RF") return InstrumentGPU::QuadrupoleRF;
    if (key == "TOF") return InstrumentGPU::TOF;
    if (key == "FT_ICR" || key == "FTICR") return InstrumentGPU::FT_ICR;
    return InstrumentGPU::NoFixedInstrument;
}

Vec3 readVec3(const Json::Value& node, const Vec3& fallback = Vec3{0.0, 0.0, 0.0}) {
    if (node.isArray() && node.size() >= 3) {
        return Vec3(node[0].asDouble(), node[1].asDouble(), node[2].asDouble());
    }
    return fallback;
}

double readDouble(const Json::Value& obj, const char* key, double fallback) {
    if (obj.isObject() && obj.isMember(key)) {
        return obj[key].asDouble();
    }
    return fallback;
}

bool hasUniformField(const Json::Value& root) {
    return root.isMember("fields") &&
           root["fields"].isMember("electric") &&
           root["fields"]["electric"].isMember("uniform") &&
           root["fields"]["electric"]["uniform"].isMember("field");
}

Vec3 readUniformField(const Json::Value& root) {
    if (!hasUniformField(root)) {
        return Vec3{0.0, 0.0, 0.0};
    }
    return readVec3(root["fields"]["electric"]["uniform"]["field"]);
}

Vec3 readUniformMagnetic(const Json::Value& root) {
    if (!(root.isMember("fields") &&
          root["fields"].isMember("magnetic") &&
          root["fields"]["magnetic"].isMember("uniform") &&
          root["fields"]["magnetic"]["uniform"].isMember("field"))) {
        return Vec3{0.0, 0.0, 0.0};
    }
    return readVec3(root["fields"]["magnetic"]["uniform"]["field"]);
}

inline void checkCuda(cudaError_t err, const char* context) {
    if (err != cudaSuccess) {
        std::string message = std::string("[GpuIntegrator] ") + context + ": " + cudaGetErrorString(err);
        throw std::runtime_error(message);
    }
}

std::string readCollisionModelString(const Json::Value& node) {
    if (node.isString()) {
        return node.asString();
    }
    if (node.isObject()) {
        if (node.isMember("model") && node["model"].isString()) {
            return node["model"].asString();
        }
        if (node.isMember("type") && node["type"].isString()) {
            return node["type"].asString();
        }
    }
    return {};
}

int parseCollisionModel(const Json::Value& config) {
    std::string collision = "EHSS";
    bool collision_set = false;
    if (config.isMember("physics") && config["physics"].isObject()) {
        const auto& phys = config["physics"];
        if (phys.isMember("collision_model") && phys["collision_model"].isString()) {
            collision = phys["collision_model"].asString();
            collision_set = true;
        } else if (phys.isMember("collisions")) {
            const auto candidate = readCollisionModelString(phys["collisions"]);
            if (!candidate.empty()) {
                collision = candidate;
                collision_set = true;
            }
        }
    }
    if (!collision_set && config.isMember("collisions")) {
        const auto candidate = readCollisionModelString(config["collisions"]);
        if (!candidate.empty()) {
            collision = candidate;
            collision_set = true;
        }
    }
    if (!collision_set && config.isMember("collision_model") && config["collision_model"].isString()) {
        collision = config["collision_model"].asString();
        collision_set = true;
    }

    std::transform(collision.begin(), collision.end(), collision.begin(),
                   [](unsigned char c) { return static_cast<char>(std::tolower(c)); });

    if (collision == "hardsphere" || collision == "hard_sphere") return 0;
    if (collision == "langevin") return 1;
    if (collision == "friction") return 2;
    if (collision == "ehss") return 3;
    if (collision == "hsmc") return 4;
    if (collision == "nocollisions" || collision == "none") return 5;
    return 3;  // default to EHSS for consistency with CPU defaults
}

}  // namespace

GpuIntegrator::GpuIntegrator(const Json::Value& config) {
    gpu_available_ = detectGpuAvailability();
    if (gpu_available_) {
        // create two non-blocking streams for double-buffer overlap
        for (int i = 0; i < 2; ++i) {
            checkCuda(cudaStreamCreateWithFlags(&streams_[i], cudaStreamNonBlocking), "cudaStreamCreateWithFlags");
            stream_created_[i] = true;
        }
    } else {
    }
    // Read optional override for sync period from environment
    const char* sp = std::getenv("ICARION_GPU_SYNC_PERIOD");
    if (sp && sp[0] != '\0') {
        int v = std::atoi(sp);
        if (v >= 0) sync_period_ = v; // 0 means no periodic sync
        if (gpuDebugEnabled()) {
            std::cerr << "[GpuIntegrator] Using sync_period=" << sync_period_ << " (from ICARION_GPU_SYNC_PERIOD)" << std::endl;
        }
    }
    // Read optional override for D2H period (how many steps between device->host downloads)
    const char* dp = std::getenv("ICARION_GPU_D2H_PERIOD");
    if (dp && dp[0] != '\0') {
        int v = std::atoi(dp);
        if (v >= 1) d2h_period_ = v; // must be >=1
        if (gpuDebugEnabled()) {
            std::cerr << "[GpuIntegrator] Using d2h_period=" << d2h_period_ << " (from ICARION_GPU_D2H_PERIOD)" << std::endl;
        }
    }
    initializeFromConfig(config);
}

GpuIntegrator::~GpuIntegrator() {
    releaseResources();
}

void GpuIntegrator::releaseResources() {
    for (int i = 0; i < 2; ++i) {
        if (device_buffer_[i]) {
            cudaFree(device_buffer_[i]);
            device_buffer_[i] = nullptr;
        }
        if (host_pinned_buffer_[i]) {
            cudaFreeHost(host_pinned_buffer_[i]);
            host_pinned_buffer_[i] = nullptr;
        }
        // destroy events if created
        if (evt_h2d_start_[i]) {
            cudaEventDestroy(evt_h2d_start_[i]);
            evt_h2d_start_[i] = nullptr;
        }
        if (evt_h2d_end_[i]) {
            cudaEventDestroy(evt_h2d_end_[i]);
            evt_h2d_end_[i] = nullptr;
        }
        if (evt_kernel_end_[i]) {
            cudaEventDestroy(evt_kernel_end_[i]);
            evt_kernel_end_[i] = nullptr;
        }
        if (evt_d2h_end_[i]) {
            cudaEventDestroy(evt_d2h_end_[i]);
            evt_d2h_end_[i] = nullptr;
        }
        buffer_pending_[i] = false;
    }
    buffer_capacity_ = 0;

    if (device_domains_) {
        cudaFree(device_domains_);
        device_domains_ = nullptr;
        device_domain_count_ = 0;
    }

    for (int i = 0; i < 2; ++i) {
        if (stream_created_[i] && streams_[i]) {
            cudaStreamDestroy(streams_[i]);
            stream_created_[i] = false;
            streams_[i] = nullptr;
        }
    }
}

std::string GpuIntegrator::getName() const {
    return "GPU-RK4";
}

void GpuIntegrator::initializeFromConfig(const Json::Value& config) {
    domains_host_.clear();
    // Don't create generic domain from config - wait for instrument-specific domains
    // Domain creation is deferred to overrideDomains() call from instrument
    // if (instrument_domains_gpu_.empty()) {
    //     instrument_domains_gpu_.push_back(buildDomainFromConfig(config));
    // }
    
    // Allow initialization to continue even without domains initially
    // Domains will be provided later via overrideDomains() from the instrument

    domains_host_ = instrument_domains_gpu_;
    global_params_.collisionModel = parseCollisionModel(config);
    // collision model debug disabled in production build
    initialized_ = true;

    // Respect an explicit debug.validate_gpu flag: when present and true,
    // prefer the baseline kernel for parity/validation. If the debug block
    // is not present, preserve the default configured above (baseline).
    parity_mode_ = false;
    if (config.isMember("debug") && config["debug"].isObject() && config["debug"].isMember("validate_gpu")) {
        parity_mode_ = config["debug"].get("validate_gpu", false).asBool();
        use_optimized_kernel_ = !parity_mode_;
    }
    // Allow env var ICARION_DEBUG_VALIDATE_GPU to opt into validation without
    // changing test configs. If present and truthy, prefer the baseline kernel
    // for numerical parity runs.
    const char* env_validate = std::getenv("ICARION_DEBUG_VALIDATE_GPU");
    if (env_validate && env_validate[0] != '\0') {
        std::string s(env_validate);
        for (auto &c : s) c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
        bool env_val = false;
        if (s == "1" || s == "true" || s == "yes" || s == "on") {
            env_val = true;
        } else {
            char* endptr = nullptr;
            long v = std::strtol(env_validate, &endptr, 10);
            if (endptr != env_validate) env_val = (v > 0);
        }
        if (env_val) {
            parity_mode_ = true;
            use_optimized_kernel_ = false;
            if (gpuDebugEnabled()) {
                ICARION_DEBUG_LOGF("gpu", "ICARION_DEBUG_VALIDATE_GPU set - forcing parity kernel for validation");
            }
        }
    }

    // If validation mode requested, map stochastic collision models (EHSS/HSMC)
    // to NoCollisions so GPU parity runs match CPU semantics (CPU uses no
    // stochastic collisions and returns zero collision force for these models).
    bool validate_gpu_flag2 = false;
    if (config.isMember("debug") && config["debug"].isObject()) {
        validate_gpu_flag2 = config["debug"].get("validate_gpu", false).asBool();
    }
    const char* env_validate2 = std::getenv("ICARION_DEBUG_VALIDATE_GPU");
    if (env_validate2 && env_validate2[0] != '\0') {
        std::string s2(env_validate2);
        for (auto &c : s2) c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
        if (s2 == "1" || s2 == "true" || s2 == "yes" || s2 == "on") {
            validate_gpu_flag2 = true;
        } else {
            char* endptr2 = nullptr;
            long v2 = std::strtol(env_validate2, &endptr2, 10);
            if (endptr2 != env_validate2) validate_gpu_flag2 = (v2 > 0);
        }
    }
    if (validate_gpu_flag2) {
        if (global_params_.collisionModel == 3 || global_params_.collisionModel == 4) {
            if (gpuDebugEnabled()) {
                ICARION_DEBUG_LOGF("gpu", "Validation mode active - disabling stochastic collisions on GPU (collisionModel=%d -> NoCollisions)", global_params_.collisionModel);
            }
            global_params_.collisionModel = 5; // NoCollisions
        }
    }

    if (gpu_available_) {
        uploadDomainsToGpu();
    }
}

DomainGPU GpuIntegrator::buildDomainFromConfig(const Json::Value& config) const {
    DomainGPU dom{};
    dom.instrument = parseInstrument(config);

    const Json::Value empty_obj(Json::objectValue);
    const Json::Value& instrument_block = (config.isMember("instrument") && config["instrument"].isObject())
                                              ? config["instrument"]
                                              : empty_obj;
    const Json::Value& geometry = (config.isMember("geometry") && config["geometry"].isObject()) ? config["geometry"] : empty_obj;
    dom.geom.length_m = readDouble(geometry, "length_m", 1.0);
    dom.geom.radius_m = readDouble(geometry, "radius_m", 0.01);
    // Boundary fields: keep in sync with geometry defaults so GPU clipping
    // uses a canonical pair of fields dedicated for boundary tests.
    dom.boundary_radius_m = dom.geom.radius_m;
    dom.boundary_length_m = dom.geom.length_m;
    dom.geom.radius_char_m = readDouble(geometry, "radius_char_m", dom.geom.radius_m);
    dom.geom.radius_out_m = readDouble(geometry, "radius_out_m", dom.geom.radius_m * 1.05);
    dom.geom.radius_in_m = readDouble(geometry, "radius_in_m", dom.geom.radius_m * 0.95);
    dom.geom.radius_char_m = std::max(dom.geom.radius_char_m, 1e-6);
    dom.geom.acc_length_m = readDouble(geometry, "acc_length_m", dom.geom.length_m);
    dom.geom.end_aperture_m = readDouble(geometry, "end_aperture_m", dom.geom.radius_m);
    dom.geom.origin_m = readVec3(geometry["origin_m"], Vec3{0.0, 0.0, 0.0});
    dom.geom.rot_row0 = Vec3{1.0, 0.0, 0.0};
    dom.geom.rot_row1 = Vec3{0.0, 1.0, 0.0};
    dom.geom.rot_row2 = Vec3{0.0, 0.0, 1.0};

    const Json::Value& env = config["environment"];
    const double explicit_pressure = readDouble(env, "pressure_Pa", -1.0);
    const double explicit_temperature = readDouble(env, "temperature_K", -1.0);
    const auto gas = instrument::resolveNeutralGasProperties(config, instrument_block);
    const auto conditions = instrument::resolveGasConditions(
        config, instrument_block, explicit_pressure, explicit_temperature);
    const Vec3 flow = instrument::resolveFlowVelocity(config, instrument_block);
    instrument::populateEnvironment(dom.env, gas, conditions, flow);

    const Vec3 E = readUniformField(config);
    const double length = std::max(dom.geom.length_m, 1e-6);
    dom.DC.axial_V = E.z * length;
    dom.DC.radial_V = 0.0;
    dom.DC.quad_V = 0.0;
    dom.DC.enable_radial_voltage_sweep = 0;
    dom.DC.radial_slope_V_s = 0.0;
    dom.DC.radial_start_time_s = 0.0;
    dom.DC.radial_rise_time_s = 0.0;

    dom.RF.voltage_V = 0.0;
    dom.RF.omega_rad_s = 0.0;
    dom.RF.phase_rad = 0.0;

    dom.AC.voltage_V = 0.0;
    dom.AC.omega_rad_s = 0.0;
    dom.AC.enable_voltage_sweep = 0;
    dom.AC.enable_frequency_sweep = 0;
    dom.AC.ac_start_freq_Hz = 0.0;
    dom.AC.ac_sweep_slope_Hz_per_s = 0.0;
    dom.AC.start_time_s = 0.0;
    dom.AC.rise_time_s = 0.0;
    dom.AC.amplitude_slope_V_s = 0.0;
    dom.AC.lqit_lock_enable = 0;
    dom.AC.lqit_lock_phase_rad = 0.0;
    dom.AC.lqit_lock_bandwidth_Hz = 0.0;

    const Vec3 B = readUniformMagnetic(config);
    dom.B.enabled = (std::fabs(B.x) + std::fabs(B.y) + std::fabs(B.z)) > 0.0;
    dom.B.Bxyz = B;
    dom.B.B_gradient = Vec3{0.0, 0.0, 0.0};

    dom.FA.loaded = 0;

    if ((!config.isMember("geometry") || !config["geometry"].isObject()) &&
        config.isMember("instrument") && config["instrument"].isObject()) {
        const Json::Value& instr = config["instrument"];
        if (instr.isMember("length_m")) {
            dom.geom.length_m = instr["length_m"].asDouble();
            dom.boundary_length_m = dom.geom.length_m;
        }
        if (instr.isMember("acceleration_length_m")) {
            dom.geom.acc_length_m = instr["acceleration_length_m"].asDouble();
        } else {
            dom.geom.acc_length_m = dom.geom.length_m;
        }
        if (instr.isMember("radius_m")) {
            dom.geom.radius_m = instr["radius_m"].asDouble();
            dom.geom.radius_out_m = dom.geom.radius_m * 1.05;
            dom.geom.radius_in_m = dom.geom.radius_m * 0.95;
            dom.geom.radius_char_m = std::max(dom.geom.radius_m, 1e-6);
            dom.boundary_radius_m = dom.geom.radius_m;
        }
    }

    return dom;
}

void GpuIntegrator::uploadDomainsToGpu() {
    if (!gpu_available_) {
        return;
    }

    if (device_domains_) {
        cudaFree(device_domains_);
        device_domains_ = nullptr;
        device_domain_count_ = 0;
    }

    if (domains_host_.empty()) {
        return;
    }

    device_domain_count_ = static_cast<int>(domains_host_.size());
    checkCuda(cudaMalloc(reinterpret_cast<void**>(&device_domains_),
                         domains_host_.size() * sizeof(DomainGPU)),
              "cudaMalloc device_domains");
    // Optional debug: print DomainGPU contents on the host immediately before
    // the host->device memcpy so we can verify what is being uploaded.
    const char* dbg_dom = std::getenv("ICARION_DEBUG_GPU_DOMAIN_UPLOAD");
    if (dbg_dom && dbg_dom[0] != '\0') {
        try {
            for (size_t di = 0; di < domains_host_.size(); ++di) {
                const auto& d = domains_host_[di];
                std::cerr << "[GPU-DOMAIN-UPLOAD] idx=" << di
                          << " inst=" << static_cast<int>(d.instrument)
                          << " axial_V=" << d.DC.axial_V
                          << " quad_V=" << d.DC.quad_V
                          << " radial_V=" << d.DC.radial_V  // ADD THIS!
                          << " E_field=" << (d.DC.axial_V / (d.geom.length_m > 0.0 ? d.geom.length_m : 1.0))
                          << " neutral_mass_kg=" << d.env.neutral_mass_kg
                          << " density_m3=" << d.env.particle_density_m_3
                          << " temp_K=" << d.env.temperature_K
                          << " geom.length_m=" << d.geom.length_m
                          << std::endl;
            }
        } catch (...) {}
    }
#ifdef ICARION_DEBUG_ORBITRAP_FIELDS
    // Always print for Orbitrap debugging
    for (size_t di = 0; di < domains_host_.size(); ++di) {
        const auto& d = domains_host_[di];
        if (d.instrument == InstrumentGPU::Orbitrap) {
            std::cerr << "[GPU-DOMAIN-UPLOAD Orbitrap] host inst=" << static_cast<int>(d.instrument)
                      << " radial_V=" << d.DC.radial_V
                      << " r_in=" << d.geom.radius_in_m << " r_out=" << d.geom.radius_out_m
                      << " r_char=" << d.geom.radius_char_m << std::endl;
        }
    }
#endif
    // restored from old GPU pipeline: async device upload to avoid blocking the host
    // Ensure domains are available on device before any kernels run. Use a
    // synchronous memcpy here to avoid subtle stream-ordering races where a
    // kernel launched on a different stream may execute before the async copy
    // from stream[0] completes. The slight host-side stall is acceptable since
    // domains are small and uploaded only during initialization.
    checkCuda(cudaMemcpy(device_domains_,
                         domains_host_.data(),
                         domains_host_.size() * sizeof(DomainGPU),
                         cudaMemcpyHostToDevice),
              "cudaMemcpy device_domains");
#ifdef ICARION_DEBUG_ORBITRAP_FIELDS
    // Verify the device copy by downloading the first domain for Orbitrap runs
    if (device_domain_count_ > 0) {
        DomainGPU verify{};
        checkCuda(cudaMemcpy(&verify, device_domains_, sizeof(DomainGPU), cudaMemcpyDeviceToHost),
                  "cudaMemcpy verify_device_domain");
        std::cerr << "[GPU-DOMAIN-UPLOAD Orbitrap] device inst=" << static_cast<int>(verify.instrument)
                  << " radial_V=" << verify.DC.radial_V
                  << " r_in=" << verify.geom.radius_in_m
                  << " r_out=" << verify.geom.radius_out_m
                  << " r_char=" << verify.geom.radius_char_m
                  << std::endl;
    }
#endif
    for (const auto& dom : domains_host_) {
        ICARION_DEBUG_LOGF("gpu", "[GPU-Domain] mass=%e density=%e vth=%e", dom.env.neutral_mass_kg, dom.env.particle_density_m_3, dom.env.mean_thermal_velocity_m_s);
    }
}

void GpuIntegrator::overrideDomains(const std::vector<DomainGPU>& domains) {
    std::cerr << "[GpuIntegrator::overrideDomains] Called with " << domains.size() << " domains" << std::endl;
    
    if (domains.empty()) {
        return;
    }
    instrument_domains_gpu_ = domains;
    domains_host_ = domains;
    initialized_ = true;
    
    // Multi-domain simulations require the baseline kernel because the optimized kernel
    // only supports instrument_type==0 (LQIT) force evaluation. Disable optimized kernel
    // if multiple domains are present.
    if (false) { // TEMPORARILY DISABLE - domains.size() > 1) {
        use_optimized_kernel_ = false;
        std::cerr << "[GpuIntegrator] Multi-domain simulation detected (" << domains.size() 
                  << " domains) - using baseline kernel for full instrument support" << std::endl;
    }
    
    if (gpu_available_) {
        uploadDomainsToGpu();
    }
}

int GpuIntegrator::resolveSpeciesId(const std::string& species_id) {
    auto it = species_lut_.find(species_id);
    if (it != species_lut_.end()) {
        return it->second;
    }
    int new_id = static_cast<int>(species_lut_.size());
    species_lut_[species_id] = new_id;
    return new_id;
}

bool GpuIntegrator::detectGpuAvailability() const {
    // Quick reject: check that CUDA Runtime library is loaded
    cudaError_t err = cudaFree(0);  // triggers lazy initialization
    if (err != cudaSuccess) {
        ICARION_DEBUG_LOGF("gpu",
            "[GpuIntegrator] cudaFree(0) failed: %s (code=%d)",
            cudaGetErrorString(err), (int)err);
        return false;
    }
    int device_count = 0;
    err = cudaGetDeviceCount(&device_count);

    if (err == cudaSuccess) {
        if (device_count > 0) {
            ICARION_DEBUG_LOGF("gpu",
                "[GpuIntegrator] cudaGetDeviceCount reports %d device(s) - GPU available",
                device_count);
            return true;
        }
        // device_count == 0: check for a spurious error state to avoid false negatives
        cudaError_t last = cudaGetLastError();
        ICARION_DEBUG_LOGF("gpu",
            "[GpuIntegrator] cudaGetDeviceCount returned 0 devices; cudaGetLastError()=%s (code=%d)",
            cudaGetErrorString(last), (int)last);
        // fall through to runtime-version fallback
    } else {
        // Non-success from cudaGetDeviceCount — report and treat explicit NO_DEVICE as false
        ICARION_DEBUG_LOGF("gpu",
            "[GpuIntegrator] cudaGetDeviceCount failed: %s (code=%d), device_count=%d",
            cudaGetErrorString(err), (int)err, device_count);
        if (err == cudaErrorNoDevice) {
            return false;
        }
        // otherwise fall through to runtime-version fallback
    }

    // Fallback: if the runtime reports a version, treat that as evidence the CUDA runtime is present
    int runtime_ver = 0;
    cudaError_t ver_err = cudaRuntimeGetVersion(&runtime_ver);
    if (ver_err == cudaSuccess) {
        ICARION_DEBUG_LOGF("gpu",
            "[GpuIntegrator] cudaRuntimeGetVersion succeeded (ver=%d) - treating GPU as available",
            runtime_ver);
        return true;
    }

    ICARION_DEBUG_LOGF("gpu",
        "[GpuIntegrator] No CUDA devices found and cudaRuntimeGetVersion failed: %s (code=%d)",
        cudaGetErrorString(ver_err), (int)ver_err);
    return false;
}



void GpuIntegrator::ensureCapacity(size_t ion_count) {
    if (ion_count <= buffer_capacity_) {
        return;
    }

    // free any existing buffers and events
    for (int i = 0; i < 2; ++i) {
        if (device_buffer_[i]) { cudaFree(device_buffer_[i]); device_buffer_[i] = nullptr; }
        if (host_pinned_buffer_[i]) { cudaFreeHost(host_pinned_buffer_[i]); host_pinned_buffer_[i] = nullptr; }
        if (evt_h2d_start_[i]) { cudaEventDestroy(evt_h2d_start_[i]); evt_h2d_start_[i] = nullptr; }
        if (evt_h2d_end_[i]) { cudaEventDestroy(evt_h2d_end_[i]); evt_h2d_end_[i] = nullptr; }
        if (evt_kernel_end_[i]) { cudaEventDestroy(evt_kernel_end_[i]); evt_kernel_end_[i] = nullptr; }
        if (evt_d2h_end_[i]) { cudaEventDestroy(evt_d2h_end_[i]); evt_d2h_end_[i] = nullptr; }
        buffer_pending_[i] = false;
    }

    buffer_capacity_ = ion_count;
    for (int i = 0; i < 2; ++i) {
        checkCuda(cudaMalloc(reinterpret_cast<void**>(&device_buffer_[i]), buffer_capacity_ * sizeof(IonStateGPU)),
                  "cudaMalloc device_buffer");
    // restored from old GPU pipeline: allocate portable pinned host memory
    checkCuda(cudaHostAlloc(reinterpret_cast<void**>(&host_pinned_buffer_[i]), buffer_capacity_ * sizeof(IonStateGPU), cudaHostAllocPortable),
          "cudaHostAlloc (portable) host_buffer");

        // create timing/progress events for each buffer
        checkCuda(cudaEventCreate(&evt_h2d_start_[i]), "cudaEventCreate h2d_start");
        checkCuda(cudaEventCreate(&evt_h2d_end_[i]), "cudaEventCreate h2d_end");
        checkCuda(cudaEventCreate(&evt_kernel_end_[i]), "cudaEventCreate kernel_end");
        checkCuda(cudaEventCreate(&evt_d2h_end_[i]), "cudaEventCreate d2h_end");
    }
}

void GpuIntegrator::uploadIons(const std::vector<ICARION::core::IonState>& ions,
                               IonStateGPU* buffer,
                               double current_time,
                               double dt,
                               const DomainGPU& domain) {
    // Read debug counts once to avoid repeated getenv cost inside hot loop.
    const char* dbg = std::getenv("ICARION_DEBUG_GPU_UPLOAD");
    int upload_dbg_count = 0;
    if (dbg && dbg[0] != '\0') {
        const char* ccount = std::getenv("ICARION_DEBUG_GPU_UPLOAD_COUNT");
        if (ccount && ccount[0] != '\0') upload_dbg_count = std::max(1, std::atoi(ccount));
        else upload_dbg_count = 16; // default: print first 16 ions
    }

    const char* accel_dbg = std::getenv("ICARION_DEBUG_GPU_ACCEL_CHECK");
    int accel_dbg_count = 0;
    if (accel_dbg && accel_dbg[0] != '\0') {
        const char* acc_count = std::getenv("ICARION_DEBUG_GPU_ACCEL_CHECK_COUNT");
        if (acc_count && acc_count[0] != '\0') accel_dbg_count = std::max(1, std::atoi(acc_count));
        else accel_dbg_count = 8; // default: check first 8 ions
    }

    for (size_t i = 0; i < ions.size(); ++i) {
        const auto& ion = ions[i];
        IonStateGPU& dst = buffer[i];
        dst.pos = ion.pos;
        dst.vel = ion.vel;
        dst.mass_kg = ion.mass_kg;
        dst.reduced_mobility_cm2_Vs = ion.reduced_mobility_cm2_Vs;
        dst.ion_charge_C = ion.ion_charge_C;
        dst.CCS_m2 = ion.CCS_m2;
        dst.polarizability_m3 = ion.polarizability_m3;
        dst.birth_time_s = ion.birth_time_s;
        dst.history_index = ion.history_index;
        dst.active = ion.active ? 1 : 0;
        dst.born = ion.born ? 1 : 0;
        dst.current_domain_index = ion.current_domain_index;
        dst.species_id_int = resolveSpeciesId(ion.species_id);

        dst.domain_neutral_mass_kg =
            (ion.domain_neutral_mass_kg > 0.0) ? ion.domain_neutral_mass_kg
                                               : domain.env.neutral_mass_kg;
        dst.domain_temperature_K =
            (ion.domain_temperature_K > 0.0) ? ion.domain_temperature_K
                                             : domain.env.temperature_K;
        dst.domain_particle_density_m3 =
            (ion.domain_particle_density_m3 > 0.0)
                ? ion.domain_particle_density_m3
                : domain.env.particle_density_m_3;
        dst.domain_neutral_polarizability_m3 =
            (ion.domain_neutral_polarizability_m3 > 0.0)
                ? ion.domain_neutral_polarizability_m3
                : domain.env.neutral_polarizability_m3;
        dst.domain_gas_velocity_m_s =
            (ion.domain_gas_velocity_m_s.x != 0.0 ||
             ion.domain_gas_velocity_m_s.y != 0.0 ||
             ion.domain_gas_velocity_m_s.z != 0.0)
                ? ion.domain_gas_velocity_m_s
                : domain.env.flow_velocity_m_s;

        dst.t = current_time;
        dst.dt = dt;

        // Optional runtime diagnostic: dump the first ion and domain values when
        // ICARION_DEBUG_GPU_UPLOAD is set in the environment. This is a
        // short-lived debugging aid to verify host->device upload correctness
        // and should be inexpensive when disabled.
        if (upload_dbg_count > 0 && static_cast<int>(i) < upload_dbg_count) {
            try {
                std::cerr << "[GPU-UPLOAD] ion[" << i << "] pos=(" << dst.pos.x << "," << dst.pos.y << "," << dst.pos.z << ")"
                          << " vel=(" << dst.vel.x << "," << dst.vel.y << "," << dst.vel.z << ")"
                          << " mass_kg=" << dst.mass_kg
                          << " charge_C=" << dst.ion_charge_C
                          << " mobility_cm2/Vs=" << dst.reduced_mobility_cm2_Vs
                          << " domain_axial_V=" << domain.DC.axial_V
                          << " domain_E_field=" << (domain.DC.axial_V / (domain.geom.length_m > 0.0 ? domain.geom.length_m : 1.0))
                          << " domain_neutral_mass_kg=" << domain.env.neutral_mass_kg
                          << " domain_density_m3=" << domain.env.particle_density_m_3
                          << " domain_temperature_K=" << domain.env.temperature_K
                          << " domain_gas_velocity=(" << domain.env.flow_velocity_m_s.x << "," << domain.env.flow_velocity_m_s.y << "," << domain.env.flow_velocity_m_s.z << ")"
                          << std::endl;
            } catch (...) {}
        }

        // Optional: run a device-side single-ion acceleration check and print
        if (accel_dbg_count > 0 && static_cast<int>(i) < accel_dbg_count) {
                try {
                    // Ensure the debug harness uses the same collision model as the
                    // integrator unless the user explicitly overrides via env var.
                    const char* cm_env_local = std::getenv("ICARION_DEBUG_GPU_COLLISION_MODEL");
                    if (!cm_env_local || cm_env_local[0] == '\0') {
                        // set env from global_params_.collisionModel so the harness
                        // computes with the same deterministic collision model.
                        char buf[16];
                        snprintf(buf, sizeof(buf), "%d", global_params_.collisionModel);
                        setenv("ICARION_DEBUG_GPU_COLLISION_MODEL", buf, 1);
                    }

                    // call the debug harness (declared in global scope above)
                    Vec3 dev_vel, dev_acc;
                    int rc = gpu_debug_single_accel(&dst, &domain, &dev_vel, &dev_acc);
                if (rc == 0) {
                    // Compute CPU-expected acceleration for the identical IonState + Domain.
                    // For parity debugging we compute the deterministic collision-damping
                    // contribution on the host (mirroring the device helper) and print it
                    // alongside the device-returned acceleration. This helps identify
                    // differences in friction/hard-sphere/langevin formulas or uploaded params.
                    try {
                        Vec3 cpu_acc{0.0, 0.0, 0.0};
                        // compute collision-force (N) depending on collisionModel,
                        // then convert to acceleration by dividing with mass.
                        int cm = global_params_.collisionModel;
                        if (cm == 2) { // Friction
                            const double LOSCHMIDT_CONSTANT = 2.6867811e25;
                            double density = std::max(domain.env.particle_density_m_3, 1e-30);
                            double ion_mobility = (dst.reduced_mobility_cm2_Vs * 1e-4) * LOSCHMIDT_CONSTANT / density;
                            if (std::fabs(ion_mobility) > 1e-30 && std::fabs(dst.mass_kg) > 0.0) {
                                double damping = -dst.ion_charge_C / ion_mobility; // units: C / (m2/Vs) -> C*s*V/m2
                                Vec3 collision_force = dst.vel * damping; // N
                                cpu_acc = collision_force * (1.0 / dst.mass_kg);
                            }
                        } else if (cm == 0) { // HardSphere
                            double denom = dst.mass_kg + dst.domain_neutral_mass_kg;
                            if (dst.CCS_m2 > 0.0 && std::fabs(denom) > 1e-30 && std::fabs(dst.mass_kg) > 0.0) {
                                double momentum_transfer_rate = domain.env.particle_density_m_3 * dst.CCS_m2 * domain.env.mean_thermal_velocity_m_s * dst.mass_kg / denom;
                                double damping = -dst.mass_kg * momentum_transfer_rate;
                                Vec3 collision_force = dst.vel * damping;
                                cpu_acc = collision_force * (1.0 / dst.mass_kg);
                            }
                        } else if (cm == 1) { // Langevin
                            const double PI = 3.14159265358979323846;
                            const double EPSILON_0 = 8.8541878128e-12;
                            double vel_sq = dst.vel.x*dst.vel.x + dst.vel.y*dst.vel.y + dst.vel.z*dst.vel.z;
                            double v_mag = std::sqrt(std::max(vel_sq, 1e-18));
                            double denom_mass = dst.mass_kg + dst.domain_neutral_mass_kg;
                            if (denom_mass > 1e-30 && std::fabs(dst.mass_kg) > 0.0) {
                                double reduced_mass = (dst.mass_kg * dst.domain_neutral_mass_kg) / denom_mass;
                                double polarization_denom = 4.0 * PI * EPSILON_0 * std::max(reduced_mass, 1e-30);
                                if (polarization_denom > 1e-30) {
                                    double inside = dst.domain_neutral_polarizability_m3 / polarization_denom;
                                    if (inside > 0.0) {
                                        double cs = PI * std::fabs(dst.ion_charge_C) * std::sqrt(inside) / v_mag;
                                        double freq = domain.env.particle_density_m_3 * cs * domain.env.mean_thermal_velocity_m_s * dst.domain_neutral_mass_kg / denom_mass;
                                        double damping = -dst.mass_kg * freq;
                                        Vec3 collision_force = dst.vel * damping;
                                        cpu_acc = collision_force * (1.0 / dst.mass_kg);
                                    }
                                }
                            }
                        } else {
                            // EHSS/HSMC/no-collision: cpu_acc stays zero (stochastic)
                            cpu_acc = Vec3{0.0,0.0,0.0};
                        }

                        // Compute host-side instrument/field acceleration for parity
                        // debugging. Device compute_accelerations_device composes an
                        // instrument-derived acceleration (IMS/LQIT/TOF/etc) and then
                        // adds the collision contribution. To mirror that, compute a
                        // cpu_field_acc here and include it in the printed cpu_acc.
                        Vec3 cpu_field_acc{0.0,0.0,0.0};
                        
                        // Parity computation for CPU-GPU validation
                        
                        try {
                            // instrument mapping: use domain.instrument to pick the
                            // right host-side formula. We only implement the simple
                            // helpers used by device: IMS (ims_acceleration), LQIT,
                            // Quadrupole/TOF/FTICR minimal forms used in tests.
                            const int inst = domain.instrument;
                            
                            // Instrument-specific parity calculations
                            // IMS: target drift velocity K * E, relax toward it over tau
                            if (inst == static_cast<int>(InstrumentGPU::IMS)) {
                                const double length = std::max(domain.geom.length_m, 1e-12);
                                const double field_strength = (length > 0.0) ? domain.DC.axial_V / length : 0.0;
                                const double K = dst.reduced_mobility_cm2_Vs * 1e-4;
                                Vec3 v_target = Vec3{0.0, 0.0, field_strength} * K;
                                const double kMinRelaxationTime = 1e-4;
                                double tau = kMinRelaxationTime;
                                if (dst.dt > 0.0) tau = std::max(dst.dt, kMinRelaxationTime);
                                Vec3 dv = v_target - dst.vel;
                                cpu_field_acc = dv * (1.0 / tau);
                            } else if (inst == static_cast<int>(InstrumentGPU::QuadrupoleRF)) {
                                // simple axial + RF combined as in device helpers
                                const double length = std::max(domain.geom.length_m, 1e-12);
                                Vec3 axial{0.0,0.0,0.0};
                                axial.z = domain.DC.axial_V / length;
                                // RF quadrupole simplified as 2*x*(Vquad + Vrf*cos)/r0^2
                                double r0 = std::max(domain.geom.radius_m, 1e-12);
                                double rf_phase = domain.RF.omega_rad_s * dst.t + domain.RF.phase_rad;
                                double voltage = domain.DC.quad_V + domain.RF.voltage_V * std::cos(rf_phase);
                                double fac = 2.0 * voltage / (r0 * r0);
                                Vec3 rf_field = Vec3{ fac * dst.pos.x, -fac * dst.pos.y, 0.0 };
                                double q_over_m = (std::fabs(dst.mass_kg) > 0.0) ? (dst.ion_charge_C / dst.mass_kg) : 0.0;
                                cpu_field_acc = (rf_field + axial) * q_over_m;
                                
                                // Quadrupole RF field calculation completed for parity validation
                            } else if (inst == static_cast<int>(InstrumentGPU::TOF)) {
                                if (dst.pos.z < domain.geom.acc_length_m) {
                                    const double acc_length = std::max(domain.geom.acc_length_m, 1e-12);
                                    const double field_strength = domain.DC.axial_V / acc_length;
                                    double q_over_m = (std::fabs(dst.mass_kg) > 0.0) ? (dst.ion_charge_C / dst.mass_kg) : 0.0;
                                    cpu_field_acc = Vec3{0.0,0.0,field_strength} * q_over_m;
                                }
                            } else if (inst == static_cast<int>(InstrumentGPU::FT_ICR)) {
                                const double characteristic = std::max(domain.geom.radius_m, 1e-12);
                                const double factor = domain.DC.quad_V / (characteristic * characteristic);
                                double z_center = dst.pos.z - 0.5 * domain.geom.length_m;
                                Vec3 electric_field{ dst.pos.x * factor, dst.pos.y * factor, -2.0 * z_center * factor };
                                Vec3 magnetic_force{ dst.vel.y * domain.B.Bxyz.z - dst.vel.z * domain.B.Bxyz.y,
                                                     dst.vel.z * domain.B.Bxyz.x - dst.vel.x * domain.B.Bxyz.z,
                                                     dst.vel.x * domain.B.Bxyz.y - dst.vel.y * domain.B.Bxyz.x };
                                double q_over_m = (std::fabs(dst.mass_kg) > 0.0) ? (dst.ion_charge_C / dst.mass_kg) : 0.0;
                                cpu_field_acc = (electric_field + magnetic_force) * q_over_m;
                            }
                        } catch (...) {}

                        Vec3 cpu_total_acc = cpu_acc + cpu_field_acc;

                        // Print collision and field components separately for easier
                        // parity debugging (cpu_acc == collision contribution,
                        // cpu_field_acc == instrument/field contribution).
                        std::cerr << "[GPU-ACCEL-CHECK] dev_vel=(" << dev_vel.x << "," << dev_vel.y << "," << dev_vel.z << ")"
                                  << " dev_acc=(" << dev_acc.x << "," << dev_acc.y << "," << dev_acc.z << ")"
                                  << " cpu_collision_acc=(" << cpu_acc.x << "," << cpu_acc.y << "," << cpu_acc.z << ")"
                                  << " cpu_field_acc=(" << cpu_field_acc.x << "," << cpu_field_acc.y << "," << cpu_field_acc.z << ")"
                                  << " cpu_total_acc=(" << cpu_total_acc.x << "," << cpu_total_acc.y << "," << cpu_total_acc.z << ")"
                                  << " cm=" << cm
                                  << std::endl;
                    } catch (...) {
                        std::cerr << "[GPU-ACCEL-CHECK] dev_vel=(" << dev_vel.x << "," << dev_vel.y << "," << dev_vel.z << ")"
                                  << " dev_acc=(" << dev_acc.x << "," << dev_acc.y << "," << dev_acc.z << ")"
                                  << std::endl;
                    }
                } else {
                    std::cerr << "[GPU-ACCEL-CHECK] kernel returned error code=" << rc << std::endl;
                }
            } catch (...) {}
        }
    }
}

void GpuIntegrator::downloadIons(std::vector<ICARION::core::IonState>& ions,
                                 const IonStateGPU* buffer,
                                 size_t ion_count) {
    for (size_t i = 0; i < ion_count; ++i) {
        ions[i].pos = buffer[i].pos;
        ions[i].vel = buffer[i].vel;
        ions[i].active = buffer[i].active != 0;
        ions[i].born = buffer[i].born != 0;
    }

}

void GpuIntegrator::setDownloadInterval(int N) {
    if (N < 1) return;
    // If environment override exists, keep the env value (env takes precedence)
    const char* dp = std::getenv("ICARION_GPU_D2H_PERIOD");
    if (dp && dp[0] != '\0') {
        if (gpuDebugEnabled()) {
            ICARION_DEBUG_LOGF("gpuio", "[GpuIntegrator] ICARION_GPU_D2H_PERIOD present (%s) - honoring env over config", dp);
        }
        return;
    }
    gpu_download_interval_ = N;
    if (gpuDebugEnabled()) {
        ICARION_DEBUG_LOGF("gpuio", "[GpuIntegrator] Using gpu_download_interval=%d (from config)", gpu_download_interval_);
    }
}

void GpuIntegrator::forceDownload() {
    if (!gpu_available_) return;
    // Synchronize both streams to wait for kernel + D2H to complete
    for (int si = 0; si < 2; ++si) {
        if (stream_created_[si]) {
            cudaError_t e = cudaStreamSynchronize(streams_[si]);
            if (e != cudaSuccess) {
                // propagate as runtime error
                checkCuda(e, "cudaStreamSynchronize (forceDownload)");
            }
        }
    }
    // If we have a remembered host ion vector, download pending buffers into it
    if (last_host_ions_ptr_ && last_host_ion_count_ > 0) {
        for (int i = 0; i < 2; ++i) {
            if (buffer_pending_[i]) {
                // Ensure D2H has completed (events recorded earlier)
                if (evt_d2h_end_[i]) {
                    float unused = 0.0f;
                    checkCuda(cudaEventElapsedTime(&unused, evt_kernel_end_[i], evt_d2h_end_[i]), "cudaEventElapsedTime (forceDownload)");
                }
                // Copy from pinned host buffer into host vector
                downloadIons(*last_host_ions_ptr_, host_pinned_buffer_[i], last_host_ion_count_);
                buffer_pending_[i] = false;
            }
        }
    }
}

void GpuIntegrator::step(std::vector<ICARION::core::IonState>& ions,
                         double dt,
                         double current_time,
                         std::function<Vec3(const ICARION::core::IonState&, size_t)> cpu_force_fn) {
#if 0
    std::cout << "[GPU-INTEGRATOR-DEBUG] GpuIntegrator::step() called with " << ions.size() << " ions, dt=" << dt << std::endl;
#endif
    if (!initialized_) {
        throw std::runtime_error("GpuIntegrator used before initialization");
    }
    if (ions.empty()) {
        return;
    }
    if (domains_host_.empty()) {
        throw std::runtime_error("GpuIntegrator has no active GPU domains");
    }
    if (!gpu_available_) {
        throw std::runtime_error("GpuIntegrator requested without an available CUDA device");
    }

    ensureCapacity(ions.size());
    
    // Physically rigorous approach: check if CPU forces are zero and create
    // a zero-field domain for exact GPU/CPU parity. This feature is currently
    // disabled to avoid overwriting instrument-specific GPU domains (e.g.
    // Orbitrap radial_V) that were uploaded via overrideDomains()/uploadDomainsToGpu().
    DomainGPU active_domain = domains_host_.front();
#ifdef ICARION_DEBUG_ORBITRAP_FIELDS
    if (active_domain.instrument == InstrumentGPU::Orbitrap) {
        std::cerr << "[GpuIntegrator::step] Host active_domain radial_V=" << active_domain.DC.radial_V
                  << " r_in=" << active_domain.geom.radius_in_m
                  << " r_out=" << active_domain.geom.radius_out_m
                  << " r_char=" << active_domain.geom.radius_char_m
                  << std::endl;
    }
#endif
    (void)cpu_force_fn; // unused while zero-domain path is disabled

    // Always ensure instrument domains are uploaded once to the GPU.
    if (device_domains_ == nullptr) {
        uploadDomainsToGpu();
    }

#ifdef ICARION_DEBUG_ORBITRAP_FIELDS
    if (device_domain_count_ > 0) {
        DomainGPU verify{};
        checkCuda(cudaMemcpy(&verify, device_domains_, sizeof(DomainGPU), cudaMemcpyDeviceToHost),
                  "cudaMemcpy verify_device_domain_step");
        if (verify.instrument == InstrumentGPU::Orbitrap) {
            std::cerr << "[GpuIntegrator::step] Device domain radial_V=" << verify.DC.radial_V
                      << " r_in=" << verify.geom.radius_in_m
                      << " r_out=" << verify.geom.radius_out_m
                      << " r_char=" << verify.geom.radius_char_m
                      << std::endl;
        }
    }
#endif

    // select buffer indices
    const int cur = static_cast<int>(step_counter_ % 2);
    const int other = 1 - cur;

    uploadIons(ions, host_pinned_buffer_[cur], current_time, dt, active_domain);

    const size_t byte_count = ions.size() * sizeof(IonStateGPU);

    // Record pointer to latest host ions container so forceDownload() can
    // copy device state back into it on demand.
    last_host_ions_ptr_ = &ions;
    last_host_ion_count_ = ions.size();

    // Record event before H2D (marks start of H2D interval)
    checkCuda(cudaEventRecord(evt_h2d_start_[cur], streams_[cur]), "cudaEventRecord h2d_start");
    checkCuda(cudaMemcpyAsync(device_buffer_[cur], host_pinned_buffer_[cur], byte_count,
                              cudaMemcpyHostToDevice, streams_[cur]),
              "cudaMemcpyAsync H2D");
    // mark end of H2D
    checkCuda(cudaEventRecord(evt_h2d_end_[cur], streams_[cur]), "cudaEventRecord h2d_end");

    // Read optional debug flags from environment to enable a very narrow
    // quadrupole device-side diagnostic. Values default to off.
    int debug_quad_flag = 0;
    int debug_quad_index = 0;
    int debug_quad_max_steps = 10;
    const char* dbg_q = std::getenv("ICARION_DEBUG_GPU_QUAD");
    if (dbg_q && dbg_q[0] != '\0') debug_quad_flag = 1;
    const char* dbg_q_idx = std::getenv("ICARION_DEBUG_GPU_QUAD_INDEX");
    if (dbg_q_idx && dbg_q_idx[0] != '\0') debug_quad_index = std::atoi(dbg_q_idx);
    const char* dbg_q_steps = std::getenv("ICARION_DEBUG_GPU_QUAD_STEPS");
    if (dbg_q_steps && dbg_q_steps[0] != '\0') debug_quad_max_steps = std::atoi(dbg_q_steps);

    // Select kernel variant: prefer the optimized kernel for throughput,
    // but allow the integrator to run the baseline kernel when a parity
    // / validation run is requested (debug.validate_gpu=true). This keeps
    // numerical behaviour closer to the CPU integrator for tests.
    if (debug_quad_flag || step_counter_==0) {
        std::cout << "[GPU-INTEGRATOR-DT] current_time_s=" << current_time << " dt_s=" << dt << std::endl;
    }
    #if 1
        // Set parity_mode_ from debug flag or env
        parity_mode_ = false;
        if (std::getenv("ICARION_DEBUG_VALIDATE_GPU")) parity_mode_ = true;
        // ...existing code...
    #endif
        if (use_optimized_kernel_) {
    #if 1
        std::cout << "[GPU-INTEGRATOR-DEBUG] Calling integrate_rk4_step_gpu_stream_optimized with " << ions.size() << " ions" << std::endl;
    #endif
        integrate_rk4_step_gpu_stream_optimized(device_buffer_[cur],
                            static_cast<int>(ions.size()),
                            global_params_,
                            device_domains_,
                            device_domain_count_,
                            current_time,
                            dt,
                            current_time + dt,
                            nullptr,
                            nullptr,
                            debug_quad_flag, debug_quad_index, debug_quad_max_steps,
                            streams_[cur]);
        // DISABLED: Parity kernel experimental code - requires C++ wrapper or .cu file compilation
        // } else if (parity_mode_) {
        //     std::cout << "[GPU-INTEGRATOR-DEBUG] Calling integrate_rk4_parity_kernel with " << ions.size() << " ions" << std::endl;
        //     int threads = 128;
        //     int blocks = (ions.size() + threads - 1) / threads;
        //     integrate_rk4_parity_kernel<<<blocks, threads, 0, streams_[cur]>>>(device_buffer_[cur], global_params_, device_domains_, static_cast<int>(ions.size()), current_time, dt);
        //     checkCuda(cudaGetLastError(), "integrate_rk4_parity_kernel launch");
        } else {
    #if 1
        std::cout << "[GPU-INTEGRATOR-DEBUG] Calling integrate_rk4_step_gpu_stream baseline with " << ions.size() << " ions" << std::endl;
    #endif
        integrate_rk4_step_gpu_stream(device_buffer_[cur],
                          static_cast<int>(ions.size()),
                          global_params_,
                          device_domains_,
                          device_domain_count_,
                          current_time,
                          dt,
                          current_time + dt,
                          nullptr,
                          nullptr,
                          debug_quad_flag, debug_quad_index, debug_quad_max_steps,
                          streams_[cur]);
        }

    // mark kernel end (acts as kernel completion marker)
    checkCuda(cudaEventRecord(evt_kernel_end_[cur], streams_[cur]), "cudaEventRecord kernel_end");

    // enqueue device->host async copy and mark D2H end event according to
    // the configurable gpu_download_interval_. The environment variable
    // ICARION_GPU_D2H_PERIOD (read in ctor) takes precedence and will have
    // already populated d2h_period_. If an env var is present we honor it,
    // otherwise use gpu_download_interval_. This allows config-driven tuning
    // while letting users override via env var.
    const char* dp_env = std::getenv("ICARION_GPU_D2H_PERIOD");
    bool env_override = (dp_env && dp_env[0] != '\0');
    int effective_period = env_override ? d2h_period_ : gpu_download_interval_;
    bool do_d2h = (effective_period <= 1) || ((step_counter_ % effective_period) == 0);
    if (do_d2h) {
        checkCuda(cudaMemcpyAsync(host_pinned_buffer_[cur], device_buffer_[cur], byte_count,
                                  cudaMemcpyDeviceToHost, streams_[cur]),
                  "cudaMemcpyAsync D2H");
        checkCuda(cudaEventRecord(evt_d2h_end_[cur], streams_[cur]), "cudaEventRecord d2h_end");
        // Instrumentation: record that a D2H was scheduled for diagnostics/tests
        ++download_count_;
        // mark this buffer as pending (D2H+kernel in-flight)
        buffer_pending_[cur] = true;
    } else {
        // not downloading from device this step; leave buffer_pending false so
        // no host download or event query occurs for this buffer.
        buffer_pending_[cur] = false;
        if (gpuDebugEnabled()) {
            ICARION_DEBUG_LOGF("gpuio", "[GpuIntegrator] Skipping D2H this step (%zu), interval=%d", step_counter_, effective_period);
        }
    }

    // Periodic sync to avoid unbounded latency; synchronize every sync_period_ steps
    if (sync_period_ > 0 && (step_counter_ % sync_period_ == 0)) {
        // synchronize both streams periodically
        for (int si = 0; si < 2; ++si) {
            if (stream_created_[si]) checkCuda(cudaStreamSynchronize(streams_[si]), "cudaStreamSynchronize (periodic)");
        }
        // after sync, any pending buffer is done - download and clear
        for (int i = 0; i < 2; ++i) {
            if (buffer_pending_[i]) {
                // compute timing via events
                float h2d_ms = 0.0f, kernel_ms = 0.0f, d2h_ms = 0.0f, total_ms = 0.0f;
                checkCuda(cudaEventElapsedTime(&h2d_ms, evt_h2d_start_[i], evt_h2d_end_[i]), "cudaEventElapsedTime h2d");
                checkCuda(cudaEventElapsedTime(&kernel_ms, evt_h2d_end_[i], evt_kernel_end_[i]), "cudaEventElapsedTime kernel");
                checkCuda(cudaEventElapsedTime(&d2h_ms, evt_kernel_end_[i], evt_d2h_end_[i]), "cudaEventElapsedTime d2h");
                checkCuda(cudaEventElapsedTime(&total_ms, evt_h2d_start_[i], evt_d2h_end_[i]), "cudaEventElapsedTime total");
                if (gpuDebugEnabled()) {
                    std::cerr << "[GpuIntegrator-TIMINGS] H2D_ms=" << h2d_ms
                              << " kernel_ms=" << kernel_ms
                              << " D2H+sync_ms=" << d2h_ms
                              << " total_ms=" << total_ms
                              << " ions=" << ions.size() << std::endl;
                }
                downloadIons(ions, host_pinned_buffer_[i], ions.size());
                buffer_pending_[i] = false;

                // Optional parity dump: write host/device snapshot + cpu/device accel for each ion
                const char* p_dump = std::getenv("ICARION_DEBUG_PARITY_DUMP");
                if (p_dump && p_dump[0] != '\0') {
                    try {
                        std::string fname = "/tmp/icarion_parity_dump_step" + std::to_string(step_counter_) + "_buf" + std::to_string(i) + ".json";
                        Json::Value root;
                        root["step"] = static_cast<Json::UInt64>(step_counter_);
                        root["buffer_index"] = i;
                        Json::Value ions_arr(Json::arrayValue);
                        for (size_t ji = 0; ji < ions.size(); ++ji) {
                            const IonStateGPU& s = host_pinned_buffer_[i][ji];
                            Json::Value jion;
                            jion["index"] = static_cast<Json::UInt64>(ji);
                            jion["pos"] = Json::Value(Json::arrayValue);
                            jion["pos"].append(s.pos.x);
                            jion["pos"].append(s.pos.y);
                            jion["pos"].append(s.pos.z);
                            jion["vel"] = Json::Value(Json::arrayValue);
                            jion["vel"].append(s.vel.x);
                            jion["vel"].append(s.vel.y);
                            jion["vel"].append(s.vel.z);
                            jion["dt"] = s.dt;
                            jion["mass_kg"] = s.mass_kg;

                            // compute cpu-side components (re-use accel-check logic)
                            Vec3 cpu_collision_acc{0.0,0.0,0.0};
                            int cm = global_params_.collisionModel;
                            if (cm == 2) {
                                const double LOSCHMIDT_CONSTANT = 2.6867811e25;
                                double density = std::max(active_domain.env.particle_density_m_3, 1e-30);
                                double ion_mobility = (s.reduced_mobility_cm2_Vs * 1e-4) * LOSCHMIDT_CONSTANT / density;
                                if (std::fabs(ion_mobility) > 1e-30 && std::fabs(s.mass_kg) > 0.0) {
                                    double damping = -s.ion_charge_C / ion_mobility;
                                    Vec3 collision_force = s.vel * damping;
                                    cpu_collision_acc = collision_force * (1.0 / s.mass_kg);
                                }
                            } else if (cm == 0) {
                                double denom = s.mass_kg + s.domain_neutral_mass_kg;
                                if (s.CCS_m2 > 0.0 && std::fabs(denom) > 1e-30 && std::fabs(s.mass_kg) > 0.0) {
                                    double momentum_transfer_rate = active_domain.env.particle_density_m_3 * s.CCS_m2 * active_domain.env.mean_thermal_velocity_m_s * s.mass_kg / denom;
                                    double damping = -s.mass_kg * momentum_transfer_rate;
                                    Vec3 collision_force = s.vel * damping;
                                    cpu_collision_acc = collision_force * (1.0 / s.mass_kg);
                                }
                            } else if (cm == 1) {
                                const double PI = 3.14159265358979323846;
                                const double EPSILON_0 = 8.8541878128e-12;
                                double vel_sq = s.vel.x*s.vel.x + s.vel.y*s.vel.y + s.vel.z*s.vel.z;
                                double v_mag = std::sqrt(std::max(vel_sq, 1e-18));
                                double denom_mass = s.mass_kg + s.domain_neutral_mass_kg;
                                if (denom_mass > 1e-30 && std::fabs(s.mass_kg) > 0.0) {
                                    double reduced_mass = (s.mass_kg * s.domain_neutral_mass_kg) / denom_mass;
                                    double polarization_denom = 4.0 * PI * EPSILON_0 * std::max(reduced_mass, 1e-30);
                                    if (polarization_denom > 1e-30) {
                                        double inside = s.domain_neutral_polarizability_m3 / polarization_denom;
                                        if (inside > 0.0) {
                                            double cs = PI * std::fabs(s.ion_charge_C) * std::sqrt(inside) / v_mag;
                                            double freq = active_domain.env.particle_density_m_3 * cs * active_domain.env.mean_thermal_velocity_m_s * s.domain_neutral_mass_kg / denom_mass;
                                            double damping = -s.mass_kg * freq;
                                            Vec3 collision_force = s.vel * damping;
                                            cpu_collision_acc = collision_force * (1.0 / s.mass_kg);
                                        }
                                    }
                                }
                            }

                            // cpu field accel (IMS/TOF/etc) simplified
                            Vec3 cpu_field_acc{0.0,0.0,0.0};
                            int inst = active_domain.instrument;
                            if (inst == static_cast<int>(InstrumentGPU::IMS)) {
                                const double length = std::max(active_domain.geom.length_m, 1e-12);
                                const double field_strength = (length > 0.0) ? active_domain.DC.axial_V / length : 0.0;
                                const double K = s.reduced_mobility_cm2_Vs * 1e-4;
                                Vec3 v_target = Vec3{0.0,0.0,field_strength} * K;
                                const double kMinRelaxationTime = 1e-4;
                                double tau = kMinRelaxationTime;
                                if (s.dt > 0.0) tau = std::max(s.dt, kMinRelaxationTime);
                                Vec3 dv = v_target - s.vel;
                                cpu_field_acc = dv * (1.0 / tau);
                            }

                            jion["cpu_collision_acc"] = Json::Value(Json::arrayValue);
                            jion["cpu_collision_acc"].append(cpu_collision_acc.x);
                            jion["cpu_collision_acc"].append(cpu_collision_acc.y);
                            jion["cpu_collision_acc"].append(cpu_collision_acc.z);
                            jion["cpu_field_acc"] = Json::Value(Json::arrayValue);
                            jion["cpu_field_acc"].append(cpu_field_acc.x);
                            jion["cpu_field_acc"].append(cpu_field_acc.y);
                            jion["cpu_field_acc"].append(cpu_field_acc.z);

                            // recompute device-side accel via harness for comparision
                            Vec3 dev_vel{0.0,0.0,0.0}, dev_acc{0.0,0.0,0.0};
                            int rc = gpu_debug_single_accel(&s, &active_domain, &dev_vel, &dev_acc);
                            if (rc == 0) {
                                jion["dev_acc"] = Json::Value(Json::arrayValue);
                                jion["dev_acc"].append(dev_acc.x);
                                jion["dev_acc"].append(dev_acc.y);
                                jion["dev_acc"].append(dev_acc.z);
                            } else {
                                jion["dev_acc"] = Json::nullValue;
                            }

                            ions_arr.append(jion);
                        }
                        root["ions"] = ions_arr;
                        std::ofstream ofs(fname);
                        if (ofs) {
                            ofs << root;
                            ofs.close();
                            std::cerr << "[GPU-PARITY-DUMP] Wrote " << fname << std::endl;
                        }
                    } catch (...) {}
                }
            }
        }
    } else {
        // non-blocking path: check if the other buffer finished D2H and download if ready
        if (buffer_pending_[other]) {
            cudaError_t q = cudaEventQuery(evt_d2h_end_[other]);
            if (q == cudaSuccess) {
                float h2d_ms = 0.0f, kernel_ms = 0.0f, d2h_ms = 0.0f, total_ms = 0.0f;
                checkCuda(cudaEventElapsedTime(&h2d_ms, evt_h2d_start_[other], evt_h2d_end_[other]), "cudaEventElapsedTime h2d");
                checkCuda(cudaEventElapsedTime(&kernel_ms, evt_h2d_end_[other], evt_kernel_end_[other]), "cudaEventElapsedTime kernel");
                checkCuda(cudaEventElapsedTime(&d2h_ms, evt_kernel_end_[other], evt_d2h_end_[other]), "cudaEventElapsedTime d2h");
                checkCuda(cudaEventElapsedTime(&total_ms, evt_h2d_start_[other], evt_d2h_end_[other]), "cudaEventElapsedTime total");
                if (gpuDebugEnabled()) {
                    std::cerr << "[GpuIntegrator-TIMINGS] H2D_ms=" << h2d_ms
                              << " kernel_ms=" << kernel_ms
                              << " D2H+sync_ms=" << d2h_ms
                              << " total_ms=" << total_ms
                              << " ions=" << ions.size() << std::endl;
                }
                downloadIons(ions, host_pinned_buffer_[other], ions.size());
                buffer_pending_[other] = false;
                const char* p_dump_other = std::getenv("ICARION_DEBUG_PARITY_DUMP");
                if (p_dump_other && p_dump_other[0] != '\0') {
                    try {
                        std::string fname = "/tmp/icarion_parity_dump_step" + std::to_string(step_counter_) + "_buf" + std::to_string(other) + ".json";
                        Json::Value root;
                        root["step"] = static_cast<Json::UInt64>(step_counter_);
                        root["buffer_index"] = other;
                        Json::Value ions_arr(Json::arrayValue);
                        for (size_t ji = 0; ji < ions.size(); ++ji) {
                            const IonStateGPU& s = host_pinned_buffer_[other][ji];
                            Json::Value jion;
                            jion["index"] = static_cast<Json::UInt64>(ji);
                            jion["pos"] = Json::Value(Json::arrayValue);
                            jion["pos"].append(s.pos.x);
                            jion["pos"].append(s.pos.y);
                            jion["pos"].append(s.pos.z);
                            jion["vel"] = Json::Value(Json::arrayValue);
                            jion["vel"].append(s.vel.x);
                            jion["vel"].append(s.vel.y);
                            jion["vel"].append(s.vel.z);
                            jion["dt"] = s.dt;
                            jion["mass_kg"] = s.mass_kg;

                            // compute cpu collision + field (simplified reuse)
                            Vec3 cpu_collision_acc{0.0,0.0,0.0};
                            int cm = global_params_.collisionModel;
                            if (cm == 2) {
                                const double LOSCHMIDT_CONSTANT = 2.6867811e25;
                                double density = std::max(active_domain.env.particle_density_m_3, 1e-30);
                                double ion_mobility = (s.reduced_mobility_cm2_Vs * 1e-4) * LOSCHMIDT_CONSTANT / density;
                                if (std::fabs(ion_mobility) > 1e-30 && std::fabs(s.mass_kg) > 0.0) {
                                    double damping = -s.ion_charge_C / ion_mobility;
                                    Vec3 collision_force = s.vel * damping;
                                    cpu_collision_acc = collision_force * (1.0 / s.mass_kg);
                                }
                            }
                            Vec3 cpu_field_acc{0.0,0.0,0.0};
                            int inst = active_domain.instrument;
                            if (inst == static_cast<int>(InstrumentGPU::IMS)) {
                                const double length = std::max(active_domain.geom.length_m, 1e-12);
                                const double field_strength = (length > 0.0) ? active_domain.DC.axial_V / length : 0.0;
                                const double K = s.reduced_mobility_cm2_Vs * 1e-4;
                                Vec3 v_target = Vec3{0.0,0.0,field_strength} * K;
                                const double kMinRelaxationTime = 1e-4;
                                double tau = kMinRelaxationTime;
                                if (s.dt > 0.0) tau = std::max(s.dt, kMinRelaxationTime);
                                Vec3 dv = v_target - s.vel;
                                cpu_field_acc = dv * (1.0 / tau);
                            }
                            jion["cpu_collision_acc"] = Json::Value(Json::arrayValue);
                            jion["cpu_collision_acc"].append(cpu_collision_acc.x);
                            jion["cpu_collision_acc"].append(cpu_collision_acc.y);
                            jion["cpu_collision_acc"].append(cpu_collision_acc.z);
                            jion["cpu_field_acc"] = Json::Value(Json::arrayValue);
                            jion["cpu_field_acc"].append(cpu_field_acc.x);
                            jion["cpu_field_acc"].append(cpu_field_acc.y);
                            jion["cpu_field_acc"].append(cpu_field_acc.z);

                            Vec3 dev_vel{0.0,0.0,0.0}, dev_acc{0.0,0.0,0.0};
                            int rc = gpu_debug_single_accel(&s, &active_domain, &dev_vel, &dev_acc);
                            if (rc == 0) {
                                jion["dev_acc"] = Json::Value(Json::arrayValue);
                                jion["dev_acc"].append(dev_acc.x);
                                jion["dev_acc"].append(dev_acc.y);
                                jion["dev_acc"].append(dev_acc.z);
                            } else {
                                jion["dev_acc"] = Json::nullValue;
                            }
                            ions_arr.append(jion);
                        }
                        root["ions"] = ions_arr;
                        std::ofstream ofs(fname);
                        if (ofs) {
                            ofs << root;
                            ofs.close();
                            std::cerr << "[GPU-PARITY-DUMP] Wrote " << fname << std::endl;
                        }
                    } catch (...) {}
                }
            } else if (q != cudaErrorNotReady) {
                checkCuda(q, "cudaEventQuery d2h_end");
            }
        }
    }

    // optional D2H diagnostic
    const char* dbg_d2h = std::getenv("ICARION_DEBUG_GPU_D2H");
    if (dbg_d2h && dbg_d2h[0] != '\0') {
        try {
            if (host_pinned_buffer_[cur]) {
                const IonStateGPU& s = host_pinned_buffer_[cur][0];
                std::cout << "[D2H-DIAG] t_copied=" << s.t << " dt=" << s.dt
                          << " i=0 pos=(" << s.pos.x << "," << s.pos.y << "," << s.pos.z << ")"
                          << " vel=(" << s.vel.x << "," << s.vel.y << "," << s.vel.z << ")"
                          << " active=" << static_cast<int>(s.active)
                          << " born=" << static_cast<int>(s.born)
                          << std::endl;
            }
            cudaError_t last = cudaGetLastError();
            std::cout << "[D2H-DIAG] cudaGetLastError=" << cudaGetErrorString(last)
                      << " (" << static_cast<int>(last) << ")" << std::endl;
        } catch (...) {}
    }

    current_time_ = current_time + dt;
    ++step_counter_;
}

}  // namespace gpu
}  // namespace ICARION

#endif  // USE_GPU_ACCEL
