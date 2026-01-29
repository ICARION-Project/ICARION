// ICARION: Ion Collision And Reaction IntegratiON
// MIT License - Copyright (c) 2025 ICARION Project Contributors

/**
 * @file hdf5Writer.cpp
 * @brief Implementation of modern HDF5 writer using FullConfig
 * 
 * **v1.0 Improvements (Nov 2025):**
 * - Consistent metadata hierarchy: All under /metadata/ (no root attributes)
 * - Smart filtering: Only write species/reactions actually used in simulation
 *   - Reduces file size for large databases
 *   - Filters based on ion.species_id references
 *   - Reactions filtered by reactant species presence
 */

#include "hdf5Writer.h"
#include "core/log/Logger.h"
#include "core/config/conversion/EnumMapper.h"
#include "core/config/types/WaveformConfig.h"
#include "core/utils/hash.h"
#include <limits>
#include "core/types/IonEnsemble.h"
#include "core/integrator/strategies/RK45Strategy.h"
#include <chrono>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <set>

#ifdef _OPENMP
#include <omp.h>
#endif

#ifdef ICARION_ENABLE_CUDA
#include <cuda_runtime.h>
#endif

#include <sys/utsname.h>
#include <unistd.h>

namespace ICARION::io {

// Forward helpers
static void write_string_vector(H5::Group& group, const std::string& name, const std::vector<std::string>& data);
static void write_array_int(H5::Group& group, const std::string& name, const std::vector<int>& data);
static std::vector<std::string> collect_field_array_paths(const std::vector<config::DomainConfig>& domains);
static std::string read_file_if_small(const std::string& path, size_t max_bytes);
static std::vector<std::string> collect_field_array_paths(const std::vector<config::DomainConfig>& domains);
static std::vector<std::string> collect_field_array_paths(const std::vector<config::DomainConfig>& domains);

// ====================================================================
// Public API
// ====================================================================

void HDF5Writer::create_file(
    const std::string& filename,
    const config::FullConfig& config,
    const core::IonEnsemble& ions,
    const std::string& git_hash,
    const std::string& build_info
) {
    log::Logger::hdf5()->info("Creating HDF5 file: {}", filename);
    
    try {
        H5::H5File file(filename, H5F_ACC_TRUNC);
        
        write_config_metadata(file, config);
        write_reproducibility_metadata(file, config, git_hash, build_info);
        write_system_metadata(file);
        
        write_species_metadata(file, config.species_db, ions);
        
        if (config.physics.enable_reactions && !config.reaction_db.reactions.empty()) {
            write_reactions_metadata(file, config.reaction_db, ions, config.species_db);
        }
        
        write_domains(file, config.domains);
        
        write_ion_metadata(file, ions);
        
        file.createGroup("/trajectory");
        
        file.close();
        log::Logger::hdf5()->info("HDF5 file created successfully");
        
    } catch (const H5::Exception& e) {
        log::Logger::hdf5()->error("Failed to create HDF5 file: {}", e.getCDetailMsg());
        throw;
    }
}

void HDF5Writer::append_trajectory(
    const std::string& filename,
    double time,
    const core::IonEnsemble& ions
) {
    if (ions.empty()) {
        log::Logger::hdf5()->warn("Skipping trajectory append - no ions");
        return;
    }

    const size_t n_ions = ions.size();
    std::vector<double> positions;
    std::vector<double> velocities;
    std::vector<int> domain_indices;
    std::vector<const char*> species_ids;

    positions.reserve(n_ions * 3);
    velocities.reserve(n_ions * 3);
    domain_indices.reserve(n_ions);
    species_ids.reserve(n_ions);

    const auto* pos_x = ions.pos_x_data();
    const auto* pos_y = ions.pos_y_data();
    const auto* pos_z = ions.pos_z_data();
    const auto* vel_x = ions.vel_x_data();
    const auto* vel_y = ions.vel_y_data();
    const auto* vel_z = ions.vel_z_data();

    for (size_t i = 0; i < n_ions; ++i) {
        positions.push_back(pos_x[i]);
        positions.push_back(pos_y[i]);
        positions.push_back(pos_z[i]);

        velocities.push_back(vel_x[i]);
        velocities.push_back(vel_y[i]);
        velocities.push_back(vel_z[i]);

        domain_indices.push_back(ions.domain_index(i));
        species_ids.push_back(ions.species_id(i).c_str());
    }

    // Reuse AoS append implementation with flattened buffers
    try {
        H5::H5File file(filename, H5F_ACC_RDWR);
        H5::Group traj_group = file.openGroup("/trajectory");

        auto append_1d = [&](const std::string& name, const void* data, H5::PredType type, size_t count) {
            H5::DataSet dataset;
            if (traj_group.nameExists(name)) {
                H5::DataSet ds = traj_group.openDataSet(name);
                H5::DataSpace filespace = ds.getSpace();
                hsize_t old_dim = filespace.getSimpleExtentNpoints();
                hsize_t new_dim = old_dim + count;
                ds.extend(&new_dim);

                filespace = ds.getSpace();
                hsize_t start[1] = {old_dim};
                hsize_t write_count[1] = {count};
                H5::DataSpace memspace(1, write_count);
                filespace.selectHyperslab(H5S_SELECT_SET, write_count, start);
                ds.write(data, type, memspace, filespace);
            } else {
                hsize_t dims[1] = {count};
                hsize_t max_dims[1] = {H5S_UNLIMITED};
                H5::DataSpace space(1, dims, max_dims);
                H5::DSetCreatPropList plist;
                hsize_t chunk[1] = {std::min(count, size_t(1000))};
                plist.setChunk(1, chunk);
                plist.setDeflate(6);
                H5::DataSet ds = traj_group.createDataSet(name, type, space, plist);
                ds.write(data, type);
            }
        };

        auto append_2d = [&](const std::string& name, const void* data, H5::PredType type, 
                             size_t dim0, size_t dim1) {
            H5::DataSet dataset;
            if (traj_group.nameExists(name)) {
                dataset = traj_group.openDataSet(name);
                H5::DataSpace filespace = dataset.getSpace();
                hsize_t old_dims[2];
                filespace.getSimpleExtentDims(old_dims);

                hsize_t new_dims[2] = {old_dims[0] + dim0, old_dims[1]};
                dataset.extend(new_dims);

                filespace = dataset.getSpace();
                hsize_t start[2] = {old_dims[0], 0};
                hsize_t count[2] = {dim0, dim1};
                H5::DataSpace memspace(2, count);
                filespace.selectHyperslab(H5S_SELECT_SET, count, start);
                dataset.write(data, type, memspace, filespace);

            } else {
                hsize_t dims[2] = {dim0, dim1};
                hsize_t max_dims[2] = {H5S_UNLIMITED, dim1};
                H5::DataSpace space(2, dims, max_dims);
                H5::DSetCreatPropList plist;
                hsize_t chunk[2] = {std::min(dim0, size_t(100)), dim1};
                plist.setChunk(2, chunk);
                plist.setDeflate(6);
                dataset = traj_group.createDataSet(name, type, space, plist);
                dataset.write(data, type);
            }
        };

        auto append_3d = [&](const std::string& name, const std::vector<double>& data, 
                             size_t dim0, size_t dim1, size_t dim2) {
            H5::DataSet dataset;
            if (traj_group.nameExists(name)) {
                dataset = traj_group.openDataSet(name);
                H5::DataSpace filespace = dataset.getSpace();
                hsize_t old_dims[3];
                filespace.getSimpleExtentDims(old_dims);

                hsize_t new_dims[3] = {old_dims[0] + dim0, old_dims[1], old_dims[2]};
                dataset.extend(new_dims);

                filespace = dataset.getSpace();
                hsize_t start[3] = {old_dims[0], 0, 0};
                hsize_t count[3] = {dim0, dim1, dim2};
                H5::DataSpace memspace(3, count);
                filespace.selectHyperslab(H5S_SELECT_SET, count, start);
                dataset.write(data.data(), H5::PredType::NATIVE_DOUBLE, memspace, filespace);

            } else {
                hsize_t dims[3] = {dim0, dim1, dim2};
                hsize_t max_dims[3] = {H5S_UNLIMITED, dim1, dim2};
                H5::DataSpace space(3, dims, max_dims);
                H5::DSetCreatPropList plist;
                hsize_t chunk[3] = {std::min(dim0, size_t(100)), dim1, dim2};
                plist.setChunk(3, chunk);
                plist.setDeflate(6);
                dataset = traj_group.createDataSet(name, H5::PredType::NATIVE_DOUBLE, space, plist);
                dataset.write(data.data(), H5::PredType::NATIVE_DOUBLE);
            }
        };

        append_1d("time", &time, H5::PredType::NATIVE_DOUBLE, 1);
        append_3d("positions", positions, 1, n_ions, 3);
        append_3d("velocities", velocities, 1, n_ions, 3);
        append_2d("domain_indices", domain_indices.data(), H5::PredType::NATIVE_INT, 1, n_ions);

        H5::StrType str_type(H5::PredType::C_S1, H5T_VARIABLE);
        H5::DataSet dataset_species;
        if (traj_group.nameExists("species_ids")) {
            dataset_species = traj_group.openDataSet("species_ids");
            H5::DataSpace filespace = dataset_species.getSpace();
            hsize_t old_dims[2];
            filespace.getSimpleExtentDims(old_dims);

            hsize_t new_dims[2] = {old_dims[0] + 1, old_dims[1]};
            dataset_species.extend(new_dims);

            filespace = dataset_species.getSpace();
            hsize_t start[2] = {old_dims[0], 0};
            hsize_t count[2] = {1, n_ions};
            H5::DataSpace memspace(2, count);
            filespace.selectHyperslab(H5S_SELECT_SET, count, start);
            dataset_species.write(species_ids.data(), str_type, memspace, filespace);
        } else {
            hsize_t dims[2] = {1, n_ions};
            hsize_t max_dims[2] = {H5S_UNLIMITED, n_ions};
            H5::DataSpace space(2, dims, max_dims);
            H5::DSetCreatPropList plist;
            hsize_t chunk[2] = {1, n_ions};
            plist.setChunk(2, chunk);
            plist.setDeflate(6);
            dataset_species = traj_group.createDataSet("species_ids", str_type, space, plist);
            dataset_species.write(species_ids.data(), str_type);
        }

        H5Fflush(file.getId(), H5F_SCOPE_GLOBAL);
        file.close();
    } catch (const H5::Exception& e) {
        log::Logger::hdf5()->error("Failed to append trajectory (SoA): {}", e.getCDetailMsg());
        throw;
    }
}

void HDF5Writer::append_trajectory_batch(
    const std::string& filename,
    const std::vector<double>& times,
    const std::vector<core::IonEnsemble>& trajectories
) {
    if (times.empty() || trajectories.empty()) {
        log::Logger::hdf5()->warn("Skipping trajectory batch append - no data");
        return;
    }
    if (times.size() != trajectories.size()) {
        throw std::invalid_argument("HDF5Writer::append_trajectory_batch (SoA): times and trajectories size mismatch");
    }
    const size_t n_timesteps = times.size();
    const size_t n_ions = trajectories[0].size();
    for (size_t t = 0; t < n_timesteps; ++t) {
        if (trajectories[t].size() != n_ions) {
            throw std::invalid_argument("HDF5Writer::append_trajectory_batch (SoA): inconsistent ion count across timesteps");
        }
    }

    try {
        H5::H5File file(filename, H5F_ACC_RDWR);
        H5::Group traj_group = file.openGroup("/trajectory");

        std::vector<double> all_times;
        std::vector<double> all_positions, all_velocities;
        std::vector<int> all_domain_indices;
        std::vector<const char*> all_species_ids;

        all_times.reserve(n_timesteps);
        all_positions.reserve(n_timesteps * n_ions * 3);
        all_velocities.reserve(n_timesteps * n_ions * 3);
        all_domain_indices.reserve(n_timesteps * n_ions);
        all_species_ids.reserve(n_timesteps * n_ions);

        for (size_t t = 0; t < n_timesteps; ++t) {
            all_times.push_back(times[t]);
            const auto& ens = trajectories[t];

            const auto* pos_x = ens.pos_x_data();
            const auto* pos_y = ens.pos_y_data();
            const auto* pos_z = ens.pos_z_data();
            const auto* vel_x = ens.vel_x_data();
            const auto* vel_y = ens.vel_y_data();
            const auto* vel_z = ens.vel_z_data();

            for (size_t i = 0; i < n_ions; ++i) {
                all_positions.push_back(pos_x[i]);
                all_positions.push_back(pos_y[i]);
                all_positions.push_back(pos_z[i]);

                all_velocities.push_back(vel_x[i]);
                all_velocities.push_back(vel_y[i]);
                all_velocities.push_back(vel_z[i]);

                all_domain_indices.push_back(ens.domain_index(i));
                all_species_ids.push_back(ens.species_id(i).c_str());
            }
        }

        auto append_1d_batch = [&](const std::string& name, const void* data, H5::PredType type, size_t count) {
            H5::DataSet dataset;
            if (traj_group.nameExists(name)) {
                dataset = traj_group.openDataSet(name);
                H5::DataSpace filespace = dataset.getSpace();
                hsize_t old_dim = filespace.getSimpleExtentNpoints();
                hsize_t new_dim = old_dim + count;
                dataset.extend(&new_dim);

                filespace = dataset.getSpace();
                hsize_t start[1] = {old_dim};
                hsize_t write_count[1] = {count};
                H5::DataSpace memspace(1, write_count);
                filespace.selectHyperslab(H5S_SELECT_SET, write_count, start);
                dataset.write(data, type, memspace, filespace);
            } else {
                hsize_t dims[1] = {count};
                hsize_t max_dims[1] = {H5S_UNLIMITED};
                H5::DataSpace space(1, dims, max_dims);
                H5::DSetCreatPropList plist;
                hsize_t chunk[1] = {std::min(count, size_t(1000))};
                plist.setChunk(1, chunk);
                plist.setDeflate(6);
                dataset = traj_group.createDataSet(name, type, space, plist);
                dataset.write(data, type);
            }
        };

        auto append_2d_batch = [&](const std::string& name, const void* data, H5::PredType type, 
                                    size_t dim0, size_t dim1) {
            H5::DataSet dataset;
            if (traj_group.nameExists(name)) {
                dataset = traj_group.openDataSet(name);
                H5::DataSpace filespace = dataset.getSpace();
                hsize_t old_dims[2];
                filespace.getSimpleExtentDims(old_dims);

                hsize_t new_dims[2] = {old_dims[0] + dim0, old_dims[1]};
                dataset.extend(new_dims);

                filespace = dataset.getSpace();
                hsize_t start[2] = {old_dims[0], 0};
                hsize_t count[2] = {dim0, dim1};
                H5::DataSpace memspace(2, count);
                filespace.selectHyperslab(H5S_SELECT_SET, count, start);
                dataset.write(data, type, memspace, filespace);
            } else {
                hsize_t dims[2] = {dim0, dim1};
                hsize_t max_dims[2] = {H5S_UNLIMITED, dim1};
                H5::DataSpace space(2, dims, max_dims);
                H5::DSetCreatPropList plist;
                hsize_t chunk[2] = {std::min(dim0, size_t(100)), dim1};
                plist.setChunk(2, chunk);
                plist.setDeflate(6);
                dataset = traj_group.createDataSet(name, type, space, plist);
                dataset.write(data, type);
            }
        };

        auto append_3d_batch = [&](const std::string& name, const std::vector<double>& data, 
                                    size_t dim0, size_t dim1, size_t dim2) {
            H5::DataSet dataset;
            if (traj_group.nameExists(name)) {
                dataset = traj_group.openDataSet(name);
                H5::DataSpace filespace = dataset.getSpace();
                hsize_t old_dims[3];
                filespace.getSimpleExtentDims(old_dims);

                hsize_t new_dims[3] = {old_dims[0] + dim0, old_dims[1], old_dims[2]};
                dataset.extend(new_dims);

                filespace = dataset.getSpace();
                hsize_t start[3] = {old_dims[0], 0, 0};
                hsize_t count[3] = {dim0, dim1, dim2};
                H5::DataSpace memspace(3, count);
                filespace.selectHyperslab(H5S_SELECT_SET, count, start);
                dataset.write(data.data(), H5::PredType::NATIVE_DOUBLE, memspace, filespace);
            } else {
                hsize_t dims[3] = {dim0, dim1, dim2};
                hsize_t max_dims[3] = {H5S_UNLIMITED, dim1, dim2};
                H5::DataSpace space(3, dims, max_dims);
                H5::DSetCreatPropList plist;
                hsize_t chunk[3] = {std::min(dim0, size_t(100)), dim1, dim2};
                plist.setChunk(3, chunk);
                plist.setDeflate(6);
                dataset = traj_group.createDataSet(name, H5::PredType::NATIVE_DOUBLE, space, plist);
                dataset.write(data.data(), H5::PredType::NATIVE_DOUBLE);
            }
        };

        append_1d_batch("time", all_times.data(), H5::PredType::NATIVE_DOUBLE, n_timesteps);
        append_3d_batch("positions", all_positions, n_timesteps, n_ions, 3);
        append_3d_batch("velocities", all_velocities, n_timesteps, n_ions, 3);
        append_2d_batch("domain_indices", all_domain_indices.data(), H5::PredType::NATIVE_INT, n_timesteps, n_ions);

        H5::StrType str_type(H5::PredType::C_S1, H5T_VARIABLE);
        H5::DataSet dataset_species;
        if (traj_group.nameExists("species_ids")) {
            dataset_species = traj_group.openDataSet("species_ids");
            H5::DataSpace filespace = dataset_species.getSpace();
            hsize_t old_dims[2];
            filespace.getSimpleExtentDims(old_dims);

            hsize_t new_dims[2] = {old_dims[0] + n_timesteps, old_dims[1]};
            dataset_species.extend(new_dims);

            filespace = dataset_species.getSpace();
            hsize_t start[2] = {old_dims[0], 0};
            hsize_t count[2] = {n_timesteps, n_ions};
            H5::DataSpace memspace(2, count);
            filespace.selectHyperslab(H5S_SELECT_SET, count, start);
            dataset_species.write(all_species_ids.data(), str_type, memspace, filespace);
        } else {
            hsize_t dims[2] = {n_timesteps, n_ions};
            hsize_t max_dims[2] = {H5S_UNLIMITED, n_ions};
            H5::DataSpace space(2, dims, max_dims);
            H5::DSetCreatPropList plist;
            hsize_t chunk[2] = {100, n_ions};
            plist.setChunk(2, chunk);
            plist.setDeflate(6);
            dataset_species = traj_group.createDataSet("species_ids", str_type, space, plist);
            dataset_species.write(all_species_ids.data(), str_type);
        }

        H5Fflush(file.getId(), H5F_SCOPE_GLOBAL);
        file.close();
    } catch (const H5::Exception& e) {
        log::Logger::hdf5()->error("Failed to append trajectory batch (SoA): {}", e.getCDetailMsg());
        throw;
    }
}

void HDF5Writer::append_trajectory_batch_flat(
    const std::string& filename,
    const std::vector<double>& times,
    size_t n_ions,
    const std::vector<double>& positions,
    const std::vector<double>& velocities,
    const std::vector<int>& domain_indices,
    const std::vector<uint32_t>& species_indices,
    const std::vector<std::string>* species_pool,
    const std::vector<double>& per_ion_times
) {
    if (times.empty()) {
        log::Logger::hdf5()->warn("Skipping trajectory batch append - no data");
        return;
    }
    const size_t n_steps = times.size();
    const size_t expected_vec = n_steps * n_ions;
    const size_t expected_vec3 = expected_vec * 3;
    if (positions.size() != expected_vec3 ||
        velocities.size() != expected_vec3 ||
        domain_indices.size() != expected_vec ||
        species_indices.size() != expected_vec) {
        throw std::invalid_argument("HDF5Writer::append_trajectory_batch_flat: buffer size mismatch");
    }
    if (!per_ion_times.empty() && per_ion_times.size() != expected_vec) {
        throw std::invalid_argument("HDF5Writer::append_trajectory_batch_flat: per_ion_times size mismatch");
    }

    try {
        H5::H5File file(filename, H5F_ACC_RDWR);
        H5::Group traj_group = file.openGroup("/trajectory");

        auto append_1d_batch = [&](const std::string& name, const void* data, H5::PredType type, size_t count) {
            H5::DataSet dataset;
            if (traj_group.nameExists(name)) {
                dataset = traj_group.openDataSet(name);
                H5::DataSpace filespace = dataset.getSpace();
                hsize_t old_dim = filespace.getSimpleExtentNpoints();
                hsize_t new_dim = old_dim + count;
                dataset.extend(&new_dim);

                filespace = dataset.getSpace();
                hsize_t start[1] = {old_dim};
                hsize_t write_count[1] = {count};
                H5::DataSpace memspace(1, write_count);
                filespace.selectHyperslab(H5S_SELECT_SET, write_count, start);
                dataset.write(data, type, memspace, filespace);
            } else {
                hsize_t dims[1] = {count};
                hsize_t max_dims[1] = {H5S_UNLIMITED};
                H5::DataSpace space(1, dims, max_dims);
                H5::DSetCreatPropList plist;
                hsize_t chunk[1] = {std::min(count, size_t(1000))};
                plist.setChunk(1, chunk);
                plist.setDeflate(2);
                dataset = traj_group.createDataSet(name, type, space, plist);
                dataset.write(data, type);
            }
        };

        auto append_3d_batch = [&](const std::string& name, const std::vector<double>& data,
                                   size_t dim0, size_t dim1, size_t dim2) {
            H5::DataSet dataset;
            if (traj_group.nameExists(name)) {
                dataset = traj_group.openDataSet(name);
                H5::DataSpace filespace = dataset.getSpace();
                hsize_t old_dims[3];
                filespace.getSimpleExtentDims(old_dims);

                hsize_t new_dims[3] = {old_dims[0] + dim0, old_dims[1], old_dims[2]};
                dataset.extend(new_dims);

                filespace = dataset.getSpace();
                hsize_t start[3] = {old_dims[0], 0, 0};
                hsize_t count[3] = {dim0, dim1, dim2};
                H5::DataSpace memspace(3, count);
                filespace.selectHyperslab(H5S_SELECT_SET, count, start);
                dataset.write(data.data(), H5::PredType::NATIVE_DOUBLE, memspace, filespace);
            } else {
                hsize_t dims[3] = {dim0, dim1, dim2};
                hsize_t max_dims[3] = {H5S_UNLIMITED, dim1, dim2};
                H5::DataSpace space(3, dims, max_dims);
                H5::DSetCreatPropList plist;
                hsize_t chunk[3] = {std::min(dim0, size_t(100)), dim1, dim2};
                plist.setChunk(3, chunk);
                plist.setDeflate(2);
                dataset = traj_group.createDataSet(name, H5::PredType::NATIVE_DOUBLE, space, plist);
                dataset.write(data.data(), H5::PredType::NATIVE_DOUBLE);
            }
        };

        auto append_2d_batch_int = [&](const std::string& name, const void* data, H5::PredType type,
                                       size_t dim0, size_t dim1) {
            H5::DataSet dataset;
            if (traj_group.nameExists(name)) {
                dataset = traj_group.openDataSet(name);
                H5::DataSpace filespace = dataset.getSpace();
                hsize_t old_dims[2];
                filespace.getSimpleExtentDims(old_dims);

                hsize_t new_dims[2] = {old_dims[0] + dim0, old_dims[1]};
                dataset.extend(new_dims);

                filespace = dataset.getSpace();
                hsize_t start[2] = {old_dims[0], 0};
                hsize_t count[2] = {dim0, dim1};
                H5::DataSpace memspace(2, count);
                filespace.selectHyperslab(H5S_SELECT_SET, count, start);
                dataset.write(data, type, memspace, filespace);
            } else {
                hsize_t dims[2] = {dim0, dim1};
                hsize_t max_dims[2] = {H5S_UNLIMITED, dim1};
                H5::DataSpace space(2, dims, max_dims);
                H5::DSetCreatPropList plist;
                hsize_t chunk[2] = {std::min(dim0, size_t(100)), dim1};
                plist.setChunk(2, chunk);
                plist.setDeflate(2);
                dataset = traj_group.createDataSet(name, type, space, plist);
                dataset.write(data, type);
            }
        };

        append_1d_batch("time", times.data(), H5::PredType::NATIVE_DOUBLE, n_steps);
        append_3d_batch("positions", positions, n_steps, n_ions, 3);
        append_3d_batch("velocities", velocities, n_steps, n_ions, 3);
        append_2d_batch_int("domain_indices", domain_indices.data(), H5::PredType::NATIVE_INT, n_steps, n_ions);
        append_2d_batch_int("species_id_indices", species_indices.data(), H5::PredType::NATIVE_UINT32, n_steps, n_ions);
        if (!per_ion_times.empty()) {
            append_2d_batch_int("time_per_ion", per_ion_times.data(), H5::PredType::NATIVE_DOUBLE, n_steps, n_ions);
        }

        // Compatibility: also write species_ids (varlen strings) if pool available
        if (species_pool && !species_pool->empty()) {
            std::vector<const char*> species_ids;
            species_ids.reserve(expected_vec);
            for (size_t idx : species_indices) {
                if (idx < species_pool->size()) {
                    species_ids.push_back(species_pool->at(idx).c_str());
                } else {
                    species_ids.push_back("");
                }
            }

            H5::StrType str_type(H5::PredType::C_S1, H5T_VARIABLE);
            H5::DataSet dataset_species;
            if (traj_group.nameExists("species_ids")) {
                dataset_species = traj_group.openDataSet("species_ids");
                H5::DataSpace filespace = dataset_species.getSpace();
                hsize_t old_dims[2];
                filespace.getSimpleExtentDims(old_dims);

                hsize_t new_dims[2] = {old_dims[0] + n_steps, old_dims[1]};
                dataset_species.extend(new_dims);

                filespace = dataset_species.getSpace();
                hsize_t start[2] = {old_dims[0], 0};
                hsize_t count[2] = {n_steps, n_ions};
                H5::DataSpace memspace(2, count);
                filespace.selectHyperslab(H5S_SELECT_SET, count, start);
                dataset_species.write(species_ids.data(), str_type, memspace, filespace);
            } else {
                hsize_t dims[2] = {n_steps, n_ions};
                hsize_t max_dims[2] = {H5S_UNLIMITED, n_ions};
                H5::DataSpace space(2, dims, max_dims);
                H5::DSetCreatPropList plist;
                hsize_t chunk[2] = {std::min(n_steps, size_t(100)), n_ions};
                plist.setChunk(2, chunk);
                plist.setDeflate(2);
                dataset_species = traj_group.createDataSet("species_ids", str_type, space, plist);
                dataset_species.write(species_ids.data(), str_type);
            }
        }

        H5Fflush(file.getId(), H5F_SCOPE_GLOBAL);
        file.close();
    } catch (const H5::Exception& e) {
        log::Logger::hdf5()->error("Failed to append trajectory batch (flat SoA): {}", e.getCDetailMsg());
        throw;
    }
}

void HDF5Writer::update_death_times(
    const std::string& filename,
    const core::IonEnsemble& final_ensemble
) {
    std::vector<double> death_times;
    death_times.reserve(final_ensemble.size());
    const auto* death_ptr = final_ensemble.death_time_data();
    for (size_t i = 0; i < final_ensemble.size(); ++i) {
        death_times.push_back(death_ptr ? death_ptr[i] : -1.0);
    }
    try {
        H5::H5File file(filename, H5F_ACC_RDWR);
        H5::Group ion_group = file.openGroup("/ions");
        hsize_t dims[1] = {death_times.size()};
        H5::DataSpace space(1, dims);
        H5::DataSet dataset = ion_group.openDataSet("death_time_s");
        dataset.write(death_times.data(), H5::PredType::NATIVE_DOUBLE);
        file.close();
        log::Logger::hdf5()->debug("Updated death_time_s for {} ions (SoA)", death_times.size());
    } catch (const H5::Exception& e) {
        log::Logger::hdf5()->error("Failed to update death times (SoA): {}", e.getCDetailMsg());
        throw;
    }
}

void HDF5Writer::finalize(
    const std::string& filename,
    bool success,
    double final_time,
    size_t active_ions
) {
    try {
        H5::H5File file(filename, H5F_ACC_RDWR);
        H5::Group meta = file.openGroup("/metadata");
        
        // Create completion subgroup
        // Make finalize idempotent: if a previous finalize attempt created the group,
        // remove it so we can recreate it cleanly.
        if (H5Lexists(meta.getId(), "completion", H5P_DEFAULT) > 0) {
            (void)H5Ldelete(meta.getId(), "completion", H5P_DEFAULT);
        }
        H5::Group completion = meta.createGroup("completion");
        
        write_scalar(completion, "success", success ? 1 : 0);
        write_scalar(completion, "final_time_s", final_time);
        write_scalar(completion, "active_ions", static_cast<int>(active_ions));
        
        // Timestamp
        auto now = std::chrono::system_clock::now();
        time_t now_time = std::chrono::system_clock::to_time_t(now);
        char timestamp_buf[64];
        strftime(timestamp_buf, sizeof(timestamp_buf), "%Y-%m-%dT%H:%M:%SZ", gmtime(&now_time));
        write_string(completion, "completion_timestamp", timestamp_buf);
        
        // Note: Root-level attributes removed in v1.0 for consistency
        // All metadata now under /metadata/ hierarchy
        
        file.close();
        log::Logger::hdf5()->info("HDF5 file finalized");
        
    } catch (const H5::Exception& e) {
        log::Logger::hdf5()->error("Failed to finalize HDF5 file: {}", e.getCDetailMsg());
        throw;
    }
}

// ====================================================================
// Metadata Writers
// ====================================================================

void HDF5Writer::write_config_metadata(
    H5::H5File& file,
    const config::FullConfig& config
) {
    H5::Group meta = file.createGroup("/metadata");
    H5::Group cfg_group = meta.createGroup("config");
    H5::Group physics_group = meta.createGroup("physics");
    
    // === Write format version ===
    write_string(cfg_group, "format_version", "2.0.0");
    
    // === Write key parameters as datasets (for quick access) ===
    write_scalar(cfg_group, "dt_s", config.simulation.dt_s);
    write_scalar(cfg_group, "total_time_s", config.simulation.total_time_s);
    write_scalar(cfg_group, "total_steps", static_cast<int>(config.simulation.total_steps));
    write_scalar(cfg_group, "write_interval", config.simulation.write_interval);
    
    write_string(cfg_group, "integrator", config.simulation.integrator);  // Already a string
    write_string(cfg_group, "collision_model", config::EnumMapper::collision_model_to_string(config.physics.collision_model));
    
    write_scalar(cfg_group, "enable_reactions", config.physics.enable_reactions);
    write_scalar(cfg_group, "enable_space_charge", config.physics.enable_space_charge);
    write_scalar(cfg_group, "enable_space_charge_gpu", config.physics.enable_space_charge_gpu);
    write_scalar(cfg_group, "enable_gpu", config.simulation.enable_gpu);
    
    write_string(cfg_group, "output_file", config.output.trajectory_file);

    // Integrator metadata (useful for reproducibility)
    H5::Group integrator_group = cfg_group.createGroup("integrator_params");
    write_string(integrator_group, "name", config.simulation.integrator);
    if (config.simulation.integrator == "RK45") {
        // Prefer runtime settings captured from RK45Strategy; fallback to defaults
        integrator::RK45Strategy::AdaptiveConfig defaults;
        const auto& rt = config.simulation.rk45_runtime_settings;
        write_scalar(integrator_group, "rk45_atol", rt ? rt->atol : defaults.atol);
        write_scalar(integrator_group, "rk45_rtol", rt ? rt->rtol : defaults.rtol);
        write_scalar(integrator_group, "rk45_safety_factor", rt ? rt->safety_factor : defaults.safety_factor);
        write_scalar(integrator_group, "rk45_min_step_factor", rt ? rt->min_step_factor : defaults.min_step_factor);
        write_scalar(integrator_group, "rk45_max_step_factor", rt ? rt->max_step_factor : defaults.max_step_factor);
        write_scalar(integrator_group, "rk45_min_step_s", rt ? rt->absolute_min_step_s : config.simulation.rk45_min_step_s);
        write_scalar(integrator_group, "rk45_max_step_increase", rt ? rt->max_step_increase : defaults.max_step_increase);
        write_scalar(integrator_group, "rk45_max_step_decrease", rt ? rt->max_step_decrease : defaults.max_step_decrease);
        write_scalar(integrator_group, "rk45_absolute_min_step_s", rt ? rt->absolute_min_step_s : defaults.absolute_min_step_s);
    }
    write_scalar(integrator_group, "openmp_enabled", config.simulation.enable_openmp);
    write_scalar(integrator_group, "gpu_collision_threshold", 5000);
#ifdef ICARION_USE_GPU
    write_scalar(integrator_group, "gpu_space_charge_threshold", 1000);
#else
    write_scalar(integrator_group, "gpu_space_charge_threshold", 0);
#endif

    // Derived summary (helps reproducibility without full serialization)
    H5::Group derived_group = cfg_group.createGroup("derived_summary");
    write_scalar(derived_group, "num_domains", static_cast<int>(config.domains.size()));
    write_scalar(derived_group, "num_species_db", static_cast<int>(config.species_db.size()));
    write_scalar(derived_group, "num_reactions_db", static_cast<int>(config.reaction_db.reactions.size()));

    // Physics handlers
    write_string(physics_group, "collision_handler", config::EnumMapper::collision_model_to_string(config.physics.collision_model));
    write_string(physics_group, "reaction_handler", config.physics.enable_reactions ? "StochasticReactionHandler" : "None");
    write_scalar(physics_group, "reaction_gpu_threshold", 2000);
    write_scalar(physics_group, "collision_gpu_threshold", 5000);
    write_string(physics_group, "collision_mixture_limit", "8 (GPU helper)");

    // Embed resolved config JSON if available (prefer in-memory snapshot from FullConfig)
    try {
        if (!config.resolved_config_json.empty()) {
            write_string(cfg_group, "config_json", config.resolved_config_json);
            log::Logger::hdf5()->debug("Embedded config snapshot from memory");
        } else {
            std::filesystem::path out_dir = config.output.folder;
            std::filesystem::path traj_file = config.output.trajectory_file;
            std::string base = traj_file.stem().string();
            if (base.empty()) {
                base = "config_snapshot";
            }
            std::filesystem::path snapshot_path = out_dir / (base + ".config.json");
            std::ifstream in(snapshot_path);
            if (!in && !config.config_file_path.empty()) {
                in.open(config.config_file_path);
            }
            if (!in) {
                if (config.config_file_path.empty()) {
                    // Tests/alternate entrypoints may not set config_file_path; embed empty with warning
                    write_string(cfg_group, "config_json", "{}");
                    log::Logger::hdf5()->warn("Config snapshot not found and no config_file_path set; embedding empty object");
                    return;
                }
                throw std::runtime_error("Config snapshot not found: " + snapshot_path.string());
            }

            std::ostringstream buffer;
            buffer << in.rdbuf();
            write_string(cfg_group, "config_json", buffer.str());
            log::Logger::hdf5()->debug("Embedded config snapshot: {}", snapshot_path.string());
        }
    } catch (const std::exception& e) {
        log::Logger::hdf5()->error("Failed to embed config snapshot: {}", e.what());
        throw;
    }
    
    log::Logger::hdf5()->debug("Wrote config metadata");
}

void HDF5Writer::write_reproducibility_metadata(
    H5::H5File& file,
    const config::FullConfig& config,
    const std::string& git_hash,
    const std::string& build_info
) {
    H5::Group repro = file.openGroup("/metadata").createGroup("reproducibility");
    
    // === RNG ===
    write_scalar(repro, "global_seed", config.simulation.rng_seed);
    write_string(repro, "rng_algorithm", "std::mt19937_64");
    write_string(repro, "seed_scheme", "global_seed + ion_index");
    write_string(repro, "per_ion_rng_scope", "collisions,reactions,stochastic_forces");
    
    // === Git info ===
    write_string(repro, "git_hash", git_hash);
    
    #ifdef GIT_DIRTY
        write_scalar(repro, "git_dirty", true);
    #else
        write_scalar(repro, "git_dirty", false);
    #endif
    
    #ifdef ICARION_VERSION
        write_string(repro, "code_version", ICARION_VERSION);
    #else
        write_string(repro, "code_version", "unknown");
    #endif
    
    // === Build info ===
    #ifdef NDEBUG
        write_string(repro, "build_type", "Release");
    #else
        write_string(repro, "build_type", "Debug");
    #endif
    
    write_string(repro, "compiler_cxx", __VERSION__);
    write_string(repro, "build_info", build_info);
    
    #ifdef ICARION_ENABLE_CUDA
        write_string(repro, "cuda_version", 
                     std::to_string(CUDART_VERSION / 1000) + "." + 
                     std::to_string((CUDART_VERSION % 100) / 10));
    #endif
    
    // === Execution ===
    bool openmp_enabled = false;
    int openmp_threads = 1;
    #ifdef _OPENMP
        openmp_enabled = true;
        openmp_threads = omp_get_max_threads();
    #endif
    write_scalar(repro, "openmp_enabled", openmp_enabled);
    write_scalar(repro, "openmp_threads", openmp_threads);
    
    // === Input file hashes (SHA256) ===
    H5::Group hash_group = repro.createGroup("input_hash");
    
    // Config file hash
    if (!config.config_file_path.empty()) {
        std::string config_hash = utils::sha256_file_safe(config.config_file_path, "N/A");
        write_string(hash_group, "config_sha256", config_hash);
        if (config_hash != "N/A" && config_hash != "ERROR") {
            log::Logger::hdf5()->debug("Config SHA256: {}", config_hash);
        } else {
            log::Logger::hdf5()->warn("Could not hash config file: {}", config.config_file_path);
        }
    } else {
        write_string(hash_group, "config_sha256", "N/A");
    }
    
    // Species database hash (if external file exists)
    if (!config.species_database_path.empty()) {
        std::string species_hash = utils::sha256_file_safe(config.species_database_path, "N/A");
        write_string(hash_group, "species_db_sha256", species_hash);
    } else {
        write_string(hash_group, "species_db_sha256", "N/A");
    }
    
    // Reaction database hash (if external file exists)
    if (!config.reaction_database_path.empty()) {
        std::string reaction_hash = utils::sha256_file_safe(config.reaction_database_path, "N/A");
        write_string(hash_group, "reaction_db_sha256", reaction_hash);
    } else {
        write_string(hash_group, "reaction_db_sha256", "N/A");
    }

    // Field array hashes (paths from domains)
    auto field_paths = collect_field_array_paths(config.domains);
    std::vector<std::string> field_hashes;
    field_hashes.reserve(field_paths.size());
    for (const auto& p : field_paths) {
        field_hashes.push_back(utils::sha256_file_safe(p, "N/A"));
    }
    H5::Group field_hash_group = hash_group.createGroup("field_arrays");
    write_string_vector(field_hash_group, "files", field_paths);
    write_string_vector(field_hash_group, "sha256", field_hashes);

    // Optional: embed external inputs for maximal reproducibility
    H5::Group blobs = repro.createGroup("input_blobs");
    const size_t MAX_EMBED_BYTES = std::numeric_limits<size_t>::max(); // Embed regardless of size

    // Config JSON
    if (!config.config_file_path.empty()) {
        auto blob = read_file_if_small(config.config_file_path, MAX_EMBED_BYTES);
        if (!blob.empty()) {
            write_string(blobs, "config_json", blob);
        } else {
            write_string(blobs, "config_json", "{}");
        }
    } else {
        write_string(blobs, "config_json", "{}");
    }

    // Species / reaction DBs
    if (!config.species_database_path.empty()) {
        auto blob = read_file_if_small(config.species_database_path, MAX_EMBED_BYTES);
        if (!blob.empty()) {
            write_string(blobs, "species_db_json", blob);
        } else {
            write_string(blobs, "species_db_json", "{}");
        }
    } else {
        write_string(blobs, "species_db_json", "{}");
    }
    if (!config.reaction_database_path.empty()) {
        auto blob = read_file_if_small(config.reaction_database_path, MAX_EMBED_BYTES);
        if (!blob.empty()) {
            write_string(blobs, "reaction_db_json", blob);
        } else {
            write_string(blobs, "reaction_db_json", "{}");
        }
    } else {
        write_string(blobs, "reaction_db_json", "{}");
    }

    // Embed field array binaries when small enough
    H5::Group field_blobs = blobs.createGroup("field_arrays");
    const auto field_paths_embed = collect_field_array_paths(config.domains);
    for (size_t i = 0; i < field_paths_embed.size(); ++i) {
        const auto& path = field_paths_embed[i];
        auto blob = read_file_if_small(path, MAX_EMBED_BYTES);
        if (!blob.empty()) {
            std::string dataset = "blob_" + std::to_string(i);
            write_blob(field_blobs, dataset, blob);
            write_string(field_blobs, dataset + "_filename", path);
        } else {
            std::string dataset = "blob_" + std::to_string(i);
            write_string(field_blobs, dataset, "");
            write_string(field_blobs, dataset + "_filename", path);
        }
    }
    
    log::Logger::hdf5()->debug("Wrote reproducibility metadata with file hashes");
}

void HDF5Writer::write_system_metadata(H5::H5File& file) {
    H5::Group sys = file.openGroup("/metadata").createGroup("system");
    
    // === Hostname ===
    char hostname[256];
    gethostname(hostname, sizeof(hostname));
    write_string(sys, "hostname", hostname);
    
    // === Username ===
    char* username = getenv("USER");
    if (username) {
        write_string(sys, "username", username);
    }
    
    // === OS ===
    struct utsname sys_info;
    uname(&sys_info);
    write_string(sys, "os", std::string(sys_info.sysname) + " " + sys_info.release);
    write_string(sys, "kernel", sys_info.version);
    
    // === CPU ===
    std::ifstream cpuinfo("/proc/cpuinfo");
    if (cpuinfo) {
        std::string line, cpu_model;
        int cpu_cores = 0;
        while (std::getline(cpuinfo, line)) {
            if (line.find("model name") != std::string::npos && cpu_model.empty()) {
                cpu_model = line.substr(line.find(":") + 2);
            }
            if (line.find("processor") != std::string::npos) {
                cpu_cores++;
            }
        }
        write_string(sys, "cpu_model", cpu_model);
        write_scalar(sys, "cpu_cores", cpu_cores);
    }
    
    // === Memory ===
    std::ifstream meminfo("/proc/meminfo");
    if (meminfo) {
        std::string line;
        long mem_total_kb = 0;
        while (std::getline(meminfo, line)) {
            if (line.find("MemTotal:") == 0) {
                std::istringstream(line.substr(10)) >> mem_total_kb;
                break;
            }
        }
        if (mem_total_kb > 0) {
            write_scalar(sys, "memory_gb", mem_total_kb / 1024.0 / 1024.0);
        }
    }
    
    // === GPU ===
    #ifdef ICARION_ENABLE_CUDA
        int device_count;
        cudaGetDeviceCount(&device_count);
        if (device_count > 0) {
            cudaDeviceProp prop;
            cudaGetDeviceProperties(&prop, 0);
            write_string(sys, "gpu_model", prop.name);
            write_scalar(sys, "gpu_memory_gb", prop.totalGlobalMem / 1024.0 / 1024.0 / 1024.0);
            
            int driver_version;
            cudaDriverGetVersion(&driver_version);
            write_string(sys, "driver_version", 
                         std::to_string(driver_version / 1000) + "." + 
                         std::to_string((driver_version % 100) / 10));
        }
    #endif
    
    // === Timestamp (ISO 8601) ===
    auto now = std::chrono::system_clock::now();
    time_t now_time = std::chrono::system_clock::to_time_t(now);
    char timestamp_buf[64];
    strftime(timestamp_buf, sizeof(timestamp_buf), "%Y-%m-%dT%H:%M:%SZ", gmtime(&now_time));
    write_string(sys, "timestamp", timestamp_buf);
    
    log::Logger::hdf5()->debug("Wrote system metadata");
}

void HDF5Writer::write_species_metadata(
    H5::H5File& file,
    const config::SpeciesDatabase& species_db,
    const core::IonEnsemble& ions
) {
    if (species_db.size() == 0) {
        log::Logger::hdf5()->warn("No species in database - skipping species metadata");
        return;
    }

    std::set<std::string> used_species;
    const auto* species_pool = ions.species_pool();
    const auto* idx = ions.species_id_indices();
    for (size_t i = 0; i < ions.size(); ++i) {
        used_species.insert((*species_pool)[idx[i]]);
    }

    if (used_species.empty()) {
        log::Logger::hdf5()->warn("No ions provided - skipping species metadata");
        return;
    }

    H5::Group species_group = file.openGroup("/metadata").createGroup("species");

    std::vector<std::string> names;
    std::vector<double> masses_kg;
    std::vector<double> charges_C;
    std::vector<double> mobilities_m2Vs;
    std::vector<double> ccs_m2;

    for (const auto& species_name : used_species) {
        auto it = species_db.species.find(species_name);
        if (it == species_db.species.end()) {
            log::Logger::hdf5()->warn("Ion references unknown species '{}' - skipping", species_name);
            continue;
        }

        const auto& species = it->second;
        names.push_back(species_name);
        masses_kg.push_back(species.mass_kg);
        charges_C.push_back(species.charge_C);
        mobilities_m2Vs.push_back(species.mobility_m2Vs);
        ccs_m2.push_back(species.CCS_m2);
    }

    size_t n = names.size();
    if (n == 0) {
        log::Logger::hdf5()->warn("No valid species found - skipping species metadata");
        return;
    }

    hsize_t dims[1] = {n};
    H5::DataSpace space(1, dims);

    H5::StrType str_type(H5::PredType::C_S1, H5T_VARIABLE);
    H5::DataSet ds_names = species_group.createDataSet("names", str_type, space);
    std::vector<const char*> name_ptrs;
    for (const auto& name : names) {
        name_ptrs.push_back(name.c_str());
    }
    ds_names.write(name_ptrs.data(), str_type);

    write_array(species_group, "mass_kg", masses_kg);
    write_array(species_group, "charge_C", charges_C);
    write_array(species_group, "mobility_m2Vs", mobilities_m2Vs);
    write_array(species_group, "ccs_m2", ccs_m2);

    log::Logger::hdf5()->info("Wrote {} species to metadata (filtered from {} total)", n, species_db.size());
}

void HDF5Writer::write_reactions_metadata(
    H5::H5File& file,
    const config::ReactionDatabase& reaction_db,
    const core::IonEnsemble& ions,
    const config::SpeciesDatabase& species_db
) {
    if (reaction_db.reactions.empty()) {
        return;
    }
    (void)species_db; // currently unused; reactions metadata derived from reaction_db and ions only

    std::set<std::string> used_species;
    const auto* species_pool = ions.species_pool();
    const auto* idx = ions.species_id_indices();
    for (size_t i = 0; i < ions.size(); ++i) {
        used_species.insert((*species_pool)[idx[i]]);
    }

    std::vector<std::string> ids;
    std::vector<std::string> reactant_1;
    std::vector<std::string> reactant_2;
    std::vector<std::string> product_1;
    std::vector<double> rate_constants;
    std::vector<int> types;

    for (const auto& rxn : reaction_db.reactions) {
        if (used_species.count(rxn.reactant) == 0) {
            continue;
        }

        ids.push_back(rxn.id);
        reactant_1.push_back(rxn.reactant);
        reactant_2.push_back("");
        product_1.push_back(rxn.product);
        rate_constants.push_back(rxn.rate_constant);
        types.push_back(2);
    }

    size_t n = ids.size();
    if (n == 0) {
        log::Logger::hdf5()->warn("No reactions matched used species");
        return;
    }

    H5::Group rxn_group = file.openGroup("/metadata").createGroup("reactions");
    write_string_vector(rxn_group, "id", ids);
    write_string_vector(rxn_group, "reactant_1", reactant_1);
    write_string_vector(rxn_group, "reactant_2", reactant_2);
    write_string_vector(rxn_group, "product_1", product_1);
    write_array_int(rxn_group, "type", types);

    write_array(rxn_group, "rate_constant_m3s", rate_constants);

    log::Logger::hdf5()->info("Wrote {} reactions to metadata (filtered from {} total)", n, reaction_db.reactions.size());
}

// ====================================================================
// Domain Writers
// ====================================================================

void HDF5Writer::write_domains(
    H5::H5File& file,
    const std::vector<config::DomainConfig>& domains
) {
    H5::Group domains_group = file.createGroup("/domains");
    
    for (size_t i = 0; i < domains.size(); ++i) {
        write_domain(domains_group, domains[i], i);
    }
    
    log::Logger::hdf5()->debug("Wrote {} domains", domains.size());
}

void HDF5Writer::write_domain(
    H5::Group& parent,
    const config::DomainConfig& domain,
    size_t index
) {
    std::string group_name = "domain_" + std::to_string(index);
    H5::Group dom_group = parent.createGroup(group_name);
    
    // === Identification ===
    write_string(dom_group, "name", domain.name);
    write_string(dom_group, "instrument", config::EnumMapper::instrument_to_string(domain.instrument));
    write_string(dom_group, "solver", config::EnumMapper::solver_to_string(domain.solver));
    
    // === Geometry ===
    H5::Group geom = dom_group.createGroup("geometry");
    write_scalar(geom, "length_m", domain.geometry.length_m);
    write_scalar(geom, "radius_m", domain.geometry.radius_m);
    write_scalar(geom, "radius_in_m", domain.geometry.radius_in_m);
    write_scalar(geom, "radius_out_m", domain.geometry.radius_out_m);
    write_vec3(geom, "origin_m", domain.geometry.origin_m);
    
    // === Environment ===
    H5::Group env = dom_group.createGroup("environment");
    write_scalar(env, "pressure_Pa", domain.environment.pressure_Pa);
    write_scalar(env, "temperature_K", domain.environment.temperature_K);
    write_string(env, "gas_species", domain.environment.gas_species);
    write_scalar(env, "particle_density_m3", domain.environment.particle_density_m_3);
    write_scalar(env, "mean_thermal_velocity_ms", domain.environment.mean_thermal_velocity_m_s);
    write_vec3(env, "gas_velocity_ms", domain.environment.gas_velocity_m_s);
    
    // === Fields ===
    H5::Group fields = dom_group.createGroup("fields");
    
    // Waveform library (v1.1: store for reproducibility)
    // Collect both named waveforms AND inline waveforms
    std::map<std::string, config::Waveform> all_waveforms = domain.fields.waveform_library;
    
    // Helper to add inline waveform if present
    auto add_inline = [&](const std::string& field_name, const config::ValueOrWaveform& vow) {
        if (vow.waveform.has_value()) {
            all_waveforms[field_name + "_inline"] = vow.waveform.value();
        }
    };
    
    // Collect all inline waveforms from fields
    add_inline("dc_axial_V", domain.fields.dc.axial_V);
    add_inline("dc_radial_V", domain.fields.dc.radial_V);
    add_inline("dc_EN_Td", domain.fields.dc.EN_Td);
    add_inline("dc_quad_V", domain.fields.dc.quad_V);
    add_inline("rf_voltage_V", domain.fields.rf.voltage_V);
    add_inline("rf_frequency_Hz", domain.fields.rf.frequency_Hz);
    add_inline("ac_voltage_V", domain.fields.ac.voltage_V);
    add_inline("ac_frequency_Hz", domain.fields.ac.frequency_Hz);
    
    // Write all waveforms (both named and inline)
    if (!all_waveforms.empty()) {
        H5::Group waveform_group = fields.createGroup("waveforms");
        write_waveform_library(waveform_group, all_waveforms);
    }
    
    // Helper lambda to safely extract t=0 value (returns 0.0 if not set)
    auto get_t0_value = [&](const config::ValueOrWaveform& vow) -> double {
        if (vow.constant_value.has_value()) {
            return vow.constant_value.value();
        }
        try {
            return vow.evaluate(0.0, domain.fields.waveform_library);
        } catch (...) {
            return 0.0;  // Field not configured
        }
    };
    
    // DC (v1.0: write static value or t=0 evaluation)
    H5::Group dc = fields.createGroup("dc");
    write_scalar(dc, "axial_V", get_t0_value(domain.fields.dc.axial_V));
    write_scalar(dc, "EN_Td", get_t0_value(domain.fields.dc.EN_Td));
    write_scalar(dc, "quad_V", get_t0_value(domain.fields.dc.quad_V));
    write_scalar(dc, "radial_V", get_t0_value(domain.fields.dc.radial_V));
    
    // RF (v1.0: write static value or t=0 evaluation)
    H5::Group rf = fields.createGroup("rf");
    write_scalar(rf, "voltage_V", get_t0_value(domain.fields.rf.voltage_V));
    write_scalar(rf, "frequency_Hz", get_t0_value(domain.fields.rf.frequency_Hz));
    write_scalar(rf, "phase_rad", domain.fields.rf.phase_rad);
    
    // AC (v1.0: write static value or t=0 evaluation)
    H5::Group ac = fields.createGroup("ac");
    write_scalar(ac, "voltage_V", get_t0_value(domain.fields.ac.voltage_V));
    write_scalar(ac, "frequency_Hz", get_t0_value(domain.fields.ac.frequency_Hz));
}

// ====================================================================
// Ion Metadata
// ====================================================================

void HDF5Writer::write_ion_metadata(
    H5::H5File& file,
    const core::IonEnsemble& ions
) {
    if (ions.empty()) {
        log::Logger::hdf5()->warn("No ions - skipping ion metadata");
        return;
    }

    H5::Group ion_group = file.createGroup("/ions");
    size_t n = ions.size();

    std::vector<std::string> species_ids;
    std::vector<double> initial_pos_x, initial_pos_y, initial_pos_z;
    std::vector<double> initial_vel_x, initial_vel_y, initial_vel_z;
    std::vector<double> birth_times;
    std::vector<double> death_times;
    std::vector<double> charges;

    species_ids.reserve(n);
    initial_pos_x.reserve(n);
    initial_pos_y.reserve(n);
    initial_pos_z.reserve(n);
    initial_vel_x.reserve(n);
    initial_vel_y.reserve(n);
    initial_vel_z.reserve(n);
    birth_times.reserve(n);
    death_times.reserve(n);
    charges.reserve(n);

    const auto* pos_x = ions.pos_x_data();
    const auto* pos_y = ions.pos_y_data();
    const auto* pos_z = ions.pos_z_data();
    const auto* vel_x = ions.vel_x_data();
    const auto* vel_y = ions.vel_y_data();
    const auto* vel_z = ions.vel_z_data();
    const auto* birth = ions.birth_time_data();
    const auto* death = ions.death_time_data();
    const auto* charge = ions.charge_data();
    const auto* species_pool = ions.species_pool();
    const auto* species_idx = ions.species_id_indices();

    for (size_t i = 0; i < n; ++i) {
        species_ids.push_back((*species_pool)[species_idx[i]]);
        initial_pos_x.push_back(pos_x[i]);
        initial_pos_y.push_back(pos_y[i]);
        initial_pos_z.push_back(pos_z[i]);
        initial_vel_x.push_back(vel_x[i]);
        initial_vel_y.push_back(vel_y[i]);
        initial_vel_z.push_back(vel_z[i]);
        birth_times.push_back(birth ? birth[i] : 0.0);
        death_times.push_back(death ? death[i] : -1.0);
        charges.push_back(charge[i]);
    }

    hsize_t dims[1] = {n};
    H5::DataSpace space(1, dims);

    H5::StrType str_type(H5::PredType::C_S1, H5T_VARIABLE);
    H5::DataSet ds_species = ion_group.createDataSet("initial_species_id", str_type, space);
    std::vector<const char*> species_ptrs;
    for (const auto& s : species_ids) species_ptrs.push_back(s.c_str());
    ds_species.write(species_ptrs.data(), str_type);

    write_array(ion_group, "initial_pos_x", initial_pos_x);
    write_array(ion_group, "initial_pos_y", initial_pos_y);
    write_array(ion_group, "initial_pos_z", initial_pos_z);
    write_array(ion_group, "initial_vel_x", initial_vel_x);
    write_array(ion_group, "initial_vel_y", initial_vel_y);
    write_array(ion_group, "initial_vel_z", initial_vel_z);
    write_array(ion_group, "birth_time_s", birth_times);
    write_array(ion_group, "death_time_s", death_times);
    write_array(ion_group, "charge_C", charges);

    log::Logger::hdf5()->debug("Wrote metadata for {} ions (SoA)", n);
}

// ====================================================================
// Helper Functions
// ====================================================================

void HDF5Writer::write_scalar(H5::Group& group, const std::string& name, double value) {
    H5::DataSpace space(H5S_SCALAR);
    H5::DataSet dataset = group.createDataSet(name, H5::PredType::NATIVE_DOUBLE, space);
    dataset.write(&value, H5::PredType::NATIVE_DOUBLE);
}

void HDF5Writer::write_scalar(H5::Group& group, const std::string& name, int value) {
    H5::DataSpace space(H5S_SCALAR);
    H5::DataSet dataset = group.createDataSet(name, H5::PredType::NATIVE_INT, space);
    dataset.write(&value, H5::PredType::NATIVE_INT);
}

void HDF5Writer::write_scalar(H5::Group& group, const std::string& name, unsigned int value) {
    H5::DataSpace space(H5S_SCALAR);
    H5::DataSet dataset = group.createDataSet(name, H5::PredType::NATIVE_UINT, space);
    dataset.write(&value, H5::PredType::NATIVE_UINT);
}

void HDF5Writer::write_scalar(H5::Group& group, const std::string& name, bool value) {
    H5::DataSpace space(H5S_SCALAR);
    H5::DataSet dataset = group.createDataSet(name, H5::PredType::NATIVE_HBOOL, space);
    hbool_t hval = value ? 1 : 0;
    dataset.write(&hval, H5::PredType::NATIVE_HBOOL);
}

static void write_string_vector(H5::Group& group, const std::string& name, const std::vector<std::string>& data) {
    if (data.empty()) return;
    hsize_t dims[2] = {data.size(), 1};
    H5::DataSpace space(2, dims);
    H5::StrType str_type(H5::PredType::C_S1, H5T_VARIABLE);
    std::vector<const char*> ptrs;
    ptrs.reserve(data.size());
    for (const auto& s : data) ptrs.push_back(s.c_str());
    H5::DataSet ds = group.createDataSet(name, str_type, space);
    ds.write(ptrs.data(), str_type);
}

static void write_array_int(H5::Group& group, const std::string& name, const std::vector<int>& data) {
    if (data.empty()) return;
    hsize_t dims[1] = {data.size()};
    H5::DataSpace space(1, dims);
    H5::DataSet ds = group.createDataSet(name, H5::PredType::NATIVE_INT, space);
    ds.write(data.data(), H5::PredType::NATIVE_INT);
}

void HDF5Writer::write_string(H5::Group& group, const std::string& name, const std::string& value) {
    H5::StrType str_type(H5::PredType::C_S1, value.size() + 1);
    H5::DataSpace space(H5S_SCALAR);
    H5::DataSet dataset = group.createDataSet(name, str_type, space);
    dataset.write(value.c_str(), str_type);
}

void HDF5Writer::write_blob(H5::Group& group, const std::string& name, const std::string& data) {
    if (data.empty()) {
        write_string(group, name, "");
        return;
    }
    hsize_t dims[1] = {data.size()};
    H5::DataSpace space(1, dims);
    H5::DataSet dataset = group.createDataSet(name, H5::PredType::NATIVE_UCHAR, space);
    dataset.write(data.data(), H5::PredType::NATIVE_UCHAR);
}

void HDF5Writer::write_array(H5::Group& group, const std::string& name, const std::vector<double>& data) {
    if (data.empty()) return;
    
    hsize_t dims[1] = {data.size()};
    H5::DataSpace space(1, dims);
    H5::DataSet dataset = group.createDataSet(name, H5::PredType::NATIVE_DOUBLE, space);
    dataset.write(data.data(), H5::PredType::NATIVE_DOUBLE);
}

static std::vector<std::string> collect_field_array_paths(const std::vector<config::DomainConfig>& domains) {
    std::vector<std::string> paths;
    for (const auto& dom : domains) {
        const auto& fields = dom.fields;
        if (!fields.legacy_field_array_file.empty()) {
            paths.push_back(fields.legacy_field_array_file);
        }
        for (const auto& term : fields.field_array_terms) {
            if (!term.file.empty()) {
                paths.push_back(term.file);
            }
        }
    }
    return paths;
}

static std::string read_file_if_small(const std::string& path, size_t max_bytes) {
    std::ifstream in(path, std::ios::binary | std::ios::ate);
    if (!in) return {};
    std::streamsize size = in.tellg();
    if (size < 0 || static_cast<size_t>(size) > max_bytes) {
        return {};
    }
    in.seekg(0, std::ios::beg);
    std::string buffer(static_cast<size_t>(size), '\0');
    if (in.read(&buffer[0], size)) {
        return buffer;
    }
    return {};
}

void HDF5Writer::write_vec3(H5::Group& group, const std::string& name, const Vec3& vec) {
    std::vector<double> data = {vec.x, vec.y, vec.z};
    write_array(group, name, data);
}

void HDF5Writer::write_waveform_library(
    H5::Group& parent,
    const std::map<std::string, config::Waveform>& library
) {
    for (const auto& [name, waveform] : library) {
        write_waveform(parent, name, waveform);
    }
    log::Logger::hdf5()->debug("Wrote {} waveforms to library", library.size());
}

void HDF5Writer::write_waveform(
    H5::Group& parent,
    const std::string& name,
    const config::Waveform& waveform
) {
    H5::Group wf_group = parent.createGroup(name);
    
    std::visit([&](const auto& w) {
        using T = std::decay_t<decltype(w)>;
        
        if constexpr (std::is_same_v<T, config::ConstantWaveform>) {
            write_string(wf_group, "type", "constant");
            write_scalar(wf_group, "value", w.value);
        }
        else if constexpr (std::is_same_v<T, config::LinearWaveform>) {
            write_string(wf_group, "type", "linear");
            write_scalar(wf_group, "start_value", w.start_value);
            write_scalar(wf_group, "end_value", w.end_value);
            write_scalar(wf_group, "start_time_s", w.start_time_s);
            write_scalar(wf_group, "end_time_s", w.end_time_s);
            write_scalar(wf_group, "clamp", w.clamp);
        }
        else if constexpr (std::is_same_v<T, config::QuadraticWaveform>) {
            write_string(wf_group, "type", "quadratic");
            write_scalar(wf_group, "a", w.a);
            write_scalar(wf_group, "b", w.b);
            write_scalar(wf_group, "c", w.c);
            write_scalar(wf_group, "start_time_s", w.start_time_s);
            write_scalar(wf_group, "end_time_s", w.end_time_s);
        }
        else if constexpr (std::is_same_v<T, config::SinusoidalWaveform>) {
            write_string(wf_group, "type", "sinusoidal");
            write_scalar(wf_group, "offset", w.offset);
            write_scalar(wf_group, "amplitude", w.amplitude);
            write_scalar(wf_group, "frequency_Hz", w.frequency_Hz);
            write_scalar(wf_group, "phase_rad", w.phase_rad);
        }
        else if constexpr (std::is_same_v<T, config::PulsedWaveform>) {
            write_string(wf_group, "type", "pulsed");
            write_scalar(wf_group, "low_value", w.low_value);
            write_scalar(wf_group, "high_value", w.high_value);
            write_scalar(wf_group, "pulse_start_s", w.pulse_start_s);
            write_scalar(wf_group, "pulse_width_s", w.pulse_width_s);
        }
        else if constexpr (std::is_same_v<T, config::ArbitraryWaveform>) {
            write_string(wf_group, "type", "arbitrary");
            write_array(wf_group, "times_s", w.times_s);
            write_array(wf_group, "values", w.values);
            std::string interp_str;
            switch (w.interp) {
                case config::ArbitraryWaveform::Interpolation::Linear: interp_str = "linear"; break;
                case config::ArbitraryWaveform::Interpolation::Step: interp_str = "step"; break;
                case config::ArbitraryWaveform::Interpolation::Cubic: interp_str = "cubic"; break;
            }
            write_string(wf_group, "interpolation", interp_str);
        }
    }, waveform.data);
}

} // namespace ICARION::io
