/**
 * =====================================================================
 *
 *   ICARION: A Modular Framework for Ion Collision and Reaction Integration
 * =====================================================================
 *   A modular C++ framework for simulating ion trajectories 
 *   in user-defined electric fields and background gas environments.
 *
 *   @file        hdf5Writer.cpp
 *   @brief       Writes simulation parameters to HDF5 output files.
 *
 *   @details
 *   Provides functions to serialize global simulation parameters and
 *   instrument domain configurations into HDF5 files. Supports scalar,
 *   string, and array data types. Each instrument domain (geometry,
 *   fields, environment) is saved in a hierarchical structure.
 *
 *
 *   @date        2025-10-06
 *   @version     0.1
 *   @author      Christoph Schäfer
 *   @license     MIT License
 *
 * =====================================================================
 */

#include "core/io/hdf5Writer.h"
#include "utils/constants.h"
#include <iostream>
#include <string>
#include <vector>
#include <map>
#include <algorithm>
#include <type_traits>
#include <chrono>
#include <ctime>

namespace ICARION {
namespace io {

// -------------------------
// Type mapping specializations for HDF5
// -------------------------

template <>
const H5::PredType& getH5Type<int>() { return H5::PredType::NATIVE_INT; }

template <>
const H5::PredType& getH5Type<double>() { return H5::PredType::NATIVE_DOUBLE; }

template <>
const H5::PredType& getH5Type<float>() { return H5::PredType::NATIVE_FLOAT; }

template <>
const H5::PredType& getH5Type<bool>() { return H5::PredType::NATIVE_HBOOL; }

// -------------------------
// Scalar write
// -------------------------
/**
 * @brief Writes a scalar value of type T to an HDF5 group.
 *
 * @tparam T Data type (int, float, double, bool)
 * @param group HDF5 group to write to
 * @param name  Dataset name
 * @param value Value to write
 */
template <typename T>
void writeScalar(H5::Group& group, const std::string& name, const T& value) {
    H5::DataSpace space(H5S_SCALAR);
    H5::DataSet dataset = group.createDataSet(name, getH5Type<T>(), space);
    dataset.write(&value, getH5Type<T>());
}

/**
 * @brief Convenience wrapper to write a boolean.
 */
void writeBool(H5::Group& group, const std::string& name, bool value) {
    writeScalar(group, name, value);
}

/**
 * @brief Writes a string scalar to an HDF5 group.
 *
 * @param group HDF5 group
 * @param name Dataset name
 * @param value String value
 */
void writeString(H5::Group& group, const std::string& name, const std::string& value) {
    H5::StrType strType(0, H5T_VARIABLE);
    H5::DataSpace space(H5S_SCALAR);
    H5::DataSet dataset = group.createDataSet(name, strType, space);
    dataset.write(value, strType);
}

/**
 * @brief Writes a 1D array to HDF5.
 *
 * @tparam T Data type
 * @param group HDF5 group
 * @param name Dataset name
 * @param values Vector of values
 */
template <typename T>
void writeArray(H5::Group& group, const std::string& name, const std::vector<T>& values) {
    hsize_t dims[1] = { values.size() };
    H5::DataSpace space(1, dims);
    H5::DataSet dataset = group.createDataSet(name, getH5Type<T>(), space);
    dataset.write(values.data(), getH5Type<T>());
}

// --- Explicit template instantiations ---
template void writeScalar<int>(H5::Group&, const std::string&, const int&);
template void writeScalar<float>(H5::Group&, const std::string&, const float&);
template void writeScalar<double>(H5::Group&, const std::string&, const double&);
template void writeArray<int>(H5::Group&, const std::string&, const std::vector<int>&);
template void writeArray<float>(H5::Group&, const std::string&, const std::vector<float>&);
template void writeArray<double>(H5::Group&, const std::string&, const std::vector<double>&);

/**
 * @brief Writes GlobalParams and a vector of InstrumentDomain objects to HDF5.
 *
 * @param file Open HDF5 file
 * @param gParams Global simulation parameters
 * @param domains Vector of instrument domains
 *
 * @note Saves GlobalParams under "global_params" and domains under "domains".
 */
void write_params_to_HDF5(H5::H5File& file, const GlobalParams& gParams,
                          const std::vector<InstrumentDomain>& domains)
{
    // -----------------------------
    // GlobalParams
    // -----------------------------
    H5::CompType gType(sizeof(GlobalParams));
    gType.insertMember("num_ions", HOFFSET(GlobalParams, num_ions), H5::PredType::NATIVE_INT);
    gType.insertMember("sim_time_steps", HOFFSET(GlobalParams, sim_time_steps), H5::PredType::NATIVE_INT);
    gType.insertMember("write_interval", HOFFSET(GlobalParams, write_interval), H5::PredType::NATIVE_INT);
    gType.insertMember("dt_s", HOFFSET(GlobalParams, dt_s), H5::PredType::NATIVE_DOUBLE);
    gType.insertMember("parallelization", HOFFSET(GlobalParams, parallelization), H5::PredType::NATIVE_HBOOL);
    gType.insertMember("enable_reactions", HOFFSET(GlobalParams, enable_reactions), H5::PredType::NATIVE_HBOOL);
    gType.insertMember("print_results", HOFFSET(GlobalParams, print_results), H5::PredType::NATIVE_HBOOL);
    gType.insertMember("collisionModel", HOFFSET(GlobalParams, collisionModel), H5::PredType::NATIVE_INT);
    gType.insertMember("enable_space_charge", HOFFSET(GlobalParams, enable_space_charge), H5::PredType::NATIVE_HBOOL);
    gType.insertMember("enable_gpu", HOFFSET(GlobalParams, enable_gpu), H5::PredType::NATIVE_HBOOL);
    gType.insertMember("enable_ou_thermalization", HOFFSET(GlobalParams, enable_ou_thermalization), H5::PredType::NATIVE_HBOOL);

    hsize_t one = 1;
    H5::DataSpace gSpace(1, &one);
    H5::DataSet gDataset = file.createDataSet("global_params", gType, gSpace);
    gDataset.write(&gParams, gType);

    // -----------------------------
    // InstrumentDomain(s)
    // -----------------------------
    H5::Group domainsGroup = file.createGroup("domains");

for (size_t i = 0; i < domains.size(); ++i) {
    const auto& dom = domains[i];
    std::string domName = "domain_" + std::to_string(i);
    H5::Group gDom = domainsGroup.createGroup(domName);

    // Instrument + Solver
    writeScalar(gDom, "instrument", static_cast<int>(dom.instrument));
    writeScalar(gDom, "solver_type", static_cast<int>(dom.solver_type));
    
    // Geometry group
    H5::Group gGeom = gDom.createGroup("geom");
    writeScalar(gGeom, "length_m", dom.geom.length_m);
    writeScalar(gGeom, "radius_m", dom.geom.radius_m);
    writeScalar(gGeom, "radius_in_m", dom.geom.radius_in_m);
    writeScalar(gGeom, "radius_out_m", dom.geom.radius_out_m);
    writeScalar(gGeom, "radius_char_m", dom.geom.radius_char_m);
    writeScalar(gGeom, "acc_length_m", dom.geom.acc_length_m);
    std::vector<double> origin = { dom.geom.origin_m.x,
                                dom.geom.origin_m.y,
                                dom.geom.origin_m.z };
    writeArray(gGeom, "origin_m", origin);

    // Environment group
    H5::Group gEnv = gDom.createGroup("env");
    writeScalar(gEnv, "pressure_Pa", dom.env.pressure_Pa);
    writeScalar(gEnv, "temperature_K", dom.env.temperature_K);
    writeScalar(gEnv, "particle_density_m_3", dom.env.particle_density_m_3);
    writeScalar(gEnv, "mean_thermal_velocity_m_s", dom.env.mean_thermal_velocity_m_s);
    writeScalar(gEnv, "neutral_mass_kg", dom.env.neutral_mass_kg);
    writeScalar(gEnv, "neutral_polarizability_m3", dom.env.neutral_polarizability_m3);
    std::vector<double> gas_vel = { dom.env.gas_velocity_m_s.x,
                                dom.env.gas_velocity_m_s.y,
                                dom.env.gas_velocity_m_s.z };
    writeArray(gEnv, "gas_velocity_m_s", gas_vel);
    writeString(gEnv, "neutral_species_id", dom.env.neutral_species_id);

    // DC, RF, AC group
    H5::Group gDC = gDom.createGroup("DC");
    writeScalar(gDC, "axial_V", dom.DC.axial_V);
    writeScalar(gDC, "EN_Td", dom.DC.EN_Td);
    writeScalar(gDC, "EN_Vm2", dom.DC.EN_Vm2);
    writeScalar(gDC, "radial_V", dom.DC.radial_V);
    writeScalar(gDC, "quad_V", dom.DC.quad_V);
    writeBool(gDC, "enable_radial_voltage_sweep", dom.DC.enable_radial_voltage_sweep);
    writeScalar(gDC, "radial_slope_V_s", dom.DC.radial_slope_V_s);
    writeScalar(gDC, "radial_start_time_s", dom.DC.radial_start_time_s);
    writeScalar(gDC, "radial_rise_time_s", dom.DC.radial_rise_time_s);

    H5::Group gRF = gDom.createGroup("RF");
    writeScalar(gRF, "voltage_V", dom.RF.voltage_V);
    writeScalar(gRF, "frequency_Hz", dom.RF.frequency_Hz);
    writeScalar(gRF, "angular_frequency_rad_s", dom.RF.angular_frequency_rad_s);
    writeScalar(gRF, "phase_rad", dom.RF.phase_rad);

    H5::Group gAC = gDom.createGroup("AC");
    writeScalar(gAC, "voltage_V", dom.AC.voltage_V);
    writeScalar(gAC, "frequency_Hz", dom.AC.frequency_Hz);
    writeScalar(gAC, "angular_frequency_rad_s", dom.AC.angular_frequency_rad_s);
    writeBool(gAC, "enable_voltage_sweep", dom.AC.enable_voltage_sweep);
    writeScalar(gAC, "amplitude_slope_V_s", dom.AC.amplitude_slope_V_s);
    writeScalar(gAC, "start_time_s", dom.AC.start_time_s);
    writeScalar(gAC, "rise_time_s", dom.AC.rise_time_s);
    // AC sweep parameters (start frequency and linear slope). These were
    // previously present in the parameter loader and runtime but not written
    // to HDF5 outputs; include them here so output files capture the AC
    // frequency sweep configuration.
    writeScalar(gAC, "ac_start_freq_Hz", dom.AC.ac_start_freq_Hz);
    writeScalar(gAC, "ac_sweep_slope_Hz_per_s", dom.AC.ac_sweep_slope_Hz_per_s);

    // Magnetic field group
    H5::Group gB = gDom.createGroup("B");
    writeBool(gB, "enabled", dom.B.enabled);
    std::vector<double> B_field = { dom.B.field_strength_T.x,
                                   dom.B.field_strength_T.y,
                                   dom.B.field_strength_T.z };
    writeArray(gB, "field_strength_T", B_field);
    std::vector<double> B_grad = { dom.B.field_gradient_T_m.x,
                                   dom.B.field_gradient_T_m.y,
                                   dom.B.field_gradient_T_m.z };
    writeArray(gB, "field_gradient_T_m", B_grad);
}

    hsize_t nDomains = domains.size();
    H5::DataSpace dSpace(1, &nDomains);
}

/**
 * @brief Append buffered time steps and ion states to HDF5 file.
 *
 * @param[in] filename        HDF5 file name
 * @param[in] times_buffer    Vector of times to append
 * @param[in] trajectory_buffer Vector of IonState vectors per timestep
 */

void append_to_HDF5(const std::string& filename,
                    const std::vector<double>& times_buffer,
                    const std::vector<std::vector<IonState>>& trajectory_buffer) {

    H5::H5File file;
    bool file_exists = false;

    // Try to open existing file
    try {
        file.openFile(filename, H5F_ACC_RDWR);
        file_exists = true;
    } catch (...) {
        file.openFile(filename, H5F_ACC_TRUNC);
        file_exists = false;
    }

    const size_t n_new_times = times_buffer.size();
    const size_t n_ions      = trajectory_buffer[0].size();

    // Safety check: cannot write empty trajectory
    if (n_ions == 0 || n_new_times == 0) {
        std::cerr << "Warning: Skipping HDF5 write - empty trajectory buffer (n_ions=" 
                  << n_ions << ", n_times=" << n_new_times << ")\n";
        file.close();
        return;
    }

    // --- Write species metadata on first write ---
    if (!file_exists && !trajectory_buffer.empty() && !trajectory_buffer[0].empty()) {
        std::cout << "[DEBUG] Calling write_species_metadata with " << trajectory_buffer[0].size() << " ions\n";
        write_species_metadata(file, trajectory_buffer[0]);
    }

    // --- Time dataset ---
    {
        H5::DataSet dataset_time;
        if (file_exists && file.exists("time")) {
            dataset_time = file.openDataSet("time");
            H5::DataSpace filespace = dataset_time.getSpace();
            hsize_t old_dim = filespace.getSimpleExtentNpoints();
            hsize_t new_dim = old_dim + n_new_times;
            dataset_time.extend(&new_dim);

            filespace = dataset_time.getSpace();
            hsize_t start[1] = {old_dim};
            hsize_t count[1] = {n_new_times};
            H5::DataSpace memspace(1, count);
            filespace.selectHyperslab(H5S_SELECT_SET, count, start);
            dataset_time.write(times_buffer.data(), H5::PredType::NATIVE_DOUBLE, memspace, filespace);

        } else {
            hsize_t dims[1]     = {n_new_times};
            hsize_t max_dims[1] = {H5S_UNLIMITED};
            H5::DataSpace space(1, dims, max_dims);
            H5::DSetCreatPropList plist;
            hsize_t chunk[1] = {1};
            plist.setChunk(1, chunk);
            dataset_time = file.createDataSet("time", H5::PredType::NATIVE_DOUBLE, space, plist);
            dataset_time.write(times_buffer.data(), H5::PredType::NATIVE_DOUBLE);
        }
    }

    // --- Flatten ion data (SLIM: only dynamic properties) ---
    std::vector<double> flat_positions, flat_velocities;
    std::vector<int>    flat_born, flat_active, flat_ion_id;
    std::vector<const char*> flat_species;

    flat_positions.reserve(n_new_times * n_ions * 3);
    flat_velocities.reserve(n_new_times * n_ions * 3);
    flat_born.reserve(n_new_times * n_ions);
    flat_active.reserve(n_new_times * n_ions);
    flat_ion_id.reserve(n_new_times * n_ions);
    flat_species.reserve(n_new_times * n_ions);

    for (size_t t_idx = 0; t_idx < trajectory_buffer.size(); ++t_idx) {
        for (size_t ion_idx = 0; ion_idx < trajectory_buffer[t_idx].size(); ++ion_idx) {
            const auto& ion = trajectory_buffer[t_idx][ion_idx];
            
            flat_positions.push_back(ion.pos.x);
            flat_positions.push_back(ion.pos.y);
            flat_positions.push_back(ion.pos.z);

            flat_velocities.push_back(ion.vel.x);
            flat_velocities.push_back(ion.vel.y);
            flat_velocities.push_back(ion.vel.z);

            flat_born.push_back(ion.born ? 1 : 0);
            flat_active.push_back(ion.active ? 1 : 0);
            flat_ion_id.push_back(static_cast<int>(ion_idx));  // Ion index as ID
            flat_species.push_back(ion.species_id.c_str());
        }
    }

    // --- Helper for 2D datasets (double or int) ---
    auto create_or_extend_2D = [&](const std::string& name, auto& data, H5::PredType type) {
        H5::DataSet dataset;
        if (file_exists && file.exists(name)) {
            dataset = file.openDataSet(name);
            H5::DataSpace filespace = dataset.getSpace();
            hsize_t old_dims[2];
            filespace.getSimpleExtentDims(old_dims);

            hsize_t new_dims[2] = {old_dims[0] + n_new_times, old_dims[1]};
            dataset.extend(new_dims);

            filespace = dataset.getSpace();
            hsize_t start[2] = {old_dims[0], 0};
            hsize_t count[2] = {n_new_times, old_dims[1]};
            H5::DataSpace memspace(2, count);
            filespace.selectHyperslab(H5S_SELECT_SET, count, start);
            dataset.write(data.data(), type, memspace, filespace);

        } else {
            hsize_t dims[2]     = {n_new_times, n_ions};
            hsize_t max_dims[2] = {H5S_UNLIMITED, n_ions};
            H5::DataSpace space(2, dims, max_dims);
            H5::DSetCreatPropList plist;
            hsize_t chunk[2] = {1, n_ions};
            plist.setChunk(2, chunk);
            dataset = file.createDataSet(name, type, space, plist);
            dataset.write(data.data(), type);
        }
    };

    // Write SLIM datasets (only dynamic properties per timestep)
    create_or_extend_2D("born", flat_born, H5::PredType::NATIVE_INT);
    create_or_extend_2D("active", flat_active, H5::PredType::NATIVE_INT);
    create_or_extend_2D("ion_id", flat_ion_id, H5::PredType::NATIVE_INT);

    // --- Strings ---
    {
        H5::StrType str_type(H5::PredType::C_S1, H5T_VARIABLE);
        H5::DataSet dataset;
        if (file_exists && file.exists("ion_species_id")) {
            dataset = file.openDataSet("ion_species_id");
            H5::DataSpace filespace = dataset.getSpace();
            hsize_t old_dims[2];
            filespace.getSimpleExtentDims(old_dims);

            hsize_t new_dims[2] = {old_dims[0] + n_new_times, old_dims[1]};
            dataset.extend(new_dims);

            filespace = dataset.getSpace();
            hsize_t start[2] = {old_dims[0], 0};
            hsize_t count[2] = {n_new_times, old_dims[1]};
            H5::DataSpace memspace(2, count);
            filespace.selectHyperslab(H5S_SELECT_SET, count, start);
            dataset.write(flat_species.data(), str_type, memspace, filespace);

        } else {
            hsize_t dims[2]     = {n_new_times, n_ions};
            hsize_t max_dims[2] = {H5S_UNLIMITED, n_ions};
            H5::DataSpace space(2, dims, max_dims);
            H5::DSetCreatPropList plist;
            hsize_t chunk[2] = {1, n_ions};
            plist.setChunk(2, chunk);
            dataset = file.createDataSet("ion_species_id", str_type, space, plist);
            dataset.write(flat_species.data(), str_type);
        }
    }

    // --- Positions & Velocities (3D) ---
    auto create_or_extend_3D = [&](const std::string& name, std::vector<double>& data) {
        H5::DataSet dataset;
        if (file_exists && file.exists(name)) {
            dataset = file.openDataSet(name);
            H5::DataSpace filespace = dataset.getSpace();
            hsize_t old_dims[3];
            filespace.getSimpleExtentDims(old_dims);

            hsize_t new_dims[3] = {old_dims[0] + n_new_times, old_dims[1], old_dims[2]};
            dataset.extend(new_dims);

            filespace = dataset.getSpace();
            hsize_t start[3] = {old_dims[0], 0, 0};
            hsize_t count[3] = {n_new_times, old_dims[1], old_dims[2]};
            H5::DataSpace memspace(3, count);
            filespace.selectHyperslab(H5S_SELECT_SET, count, start);
            dataset.write(data.data(), H5::PredType::NATIVE_DOUBLE, memspace, filespace);

        } else {
            hsize_t dims[3]     = {n_new_times, n_ions, 3};
            hsize_t max_dims[3] = {H5S_UNLIMITED, n_ions, 3};
            H5::DataSpace space(3, dims, max_dims);
            H5::DSetCreatPropList plist;
            hsize_t chunk[3] = {1, n_ions, 3};
            plist.setChunk(3, chunk);
            dataset = file.createDataSet(name, H5::PredType::NATIVE_DOUBLE, space, plist);
            dataset.write(data.data(), H5::PredType::NATIVE_DOUBLE);
        }
    };

    create_or_extend_3D("positions", flat_positions);
    create_or_extend_3D("velocities", flat_velocities);

    // Ensure data is flushed to disk so readers (tests) see the datasets immediately
    // Use C API flush via file id
    H5Fflush(file.getId(), H5F_SCOPE_GLOBAL);

    file.close();
}

// -----------------------------
// Species metadata support
// -----------------------------

/**
 * @brief Write species-constant properties to /metadata/species group.
 * 
 * Extracts unique species from ion ensemble and writes their constant properties
 * (mass, charge, mobility, CCS) once. This significantly reduces HDF5 file size
 * by avoiding redundant storage of these values at every timestep.
 */
void write_species_metadata(H5::H5File& file, const std::vector<IonState>& ions) {
    if (ions.empty()) {
        std::cerr << "Warning: write_species_metadata called with empty ion vector\n";
        return;
    }
    
    // Extract unique species
    std::map<std::string, IonState> unique_species;
    for (const auto& ion : ions) {
        if (unique_species.find(ion.species_id) == unique_species.end()) {
            unique_species[ion.species_id] = ion;
        }
    }
    
    if (unique_species.empty()) {
        std::cerr << "Warning: No unique species found in ion vector\n";
        return;
    }
    
    try {
        // Create /metadata group if not exists
        H5::Group metadata_group;
        try {
            metadata_group = file.openGroup("/metadata");
        } catch (...) {
            metadata_group = file.createGroup("/metadata");
        }
        
        // Create /metadata/species group
        H5::Group species_group;
        try {
            species_group = metadata_group.openGroup("species");
        } catch (...) {
            species_group = metadata_group.createGroup("species");
        }
        
        // Prepare species data - use std::string to keep data valid
        const size_t n_species = unique_species.size();
        std::vector<std::string> species_name_strings;  // Keep strings alive
        std::vector<const char*> species_names;
        std::vector<double> masses, mobilities, ccs_values, charges;
        
        species_name_strings.reserve(n_species);
        species_names.reserve(n_species);
        masses.reserve(n_species);
        mobilities.reserve(n_species);
        ccs_values.reserve(n_species);
        charges.reserve(n_species);
        
        for (const auto& [name, ion] : unique_species) {
            species_name_strings.push_back(name);
            species_names.push_back(species_name_strings.back().c_str());
            masses.push_back(ion.mass_kg);
            mobilities.push_back(ion.reduced_mobility_cm2_Vs);
            ccs_values.push_back(ion.CCS_m2);
            charges.push_back(ion.ion_charge_C);
        }
        
        // Create 1D datasets for each property
        hsize_t dims[1] = {n_species};
        H5::DataSpace space(1, dims);
        
        // Species names (variable-length strings)
        H5::StrType str_type(H5::PredType::C_S1, H5T_VARIABLE);
        H5::DataSet ds_names = species_group.createDataSet("name", str_type, space);
        ds_names.write(species_names.data(), str_type);
        ds_names.close();
        
        // Mass
        H5::DataSet ds_mass = species_group.createDataSet("mass_kg", H5::PredType::NATIVE_DOUBLE, space);
        ds_mass.write(masses.data(), H5::PredType::NATIVE_DOUBLE);
        ds_mass.close();
        
        // Mobility
        H5::DataSet ds_mobility = species_group.createDataSet("reduced_mobility_cm2_Vs", 
                                                              H5::PredType::NATIVE_DOUBLE, space);
        ds_mobility.write(mobilities.data(), H5::PredType::NATIVE_DOUBLE);
        ds_mobility.close();
        
        // CCS
        H5::DataSet ds_ccs = species_group.createDataSet("CCS_m2", H5::PredType::NATIVE_DOUBLE, space);
        ds_ccs.write(ccs_values.data(), H5::PredType::NATIVE_DOUBLE);
        ds_ccs.close();
        
        // Charge
        H5::DataSet ds_charge = species_group.createDataSet("charge_C", H5::PredType::NATIVE_DOUBLE, space);
        ds_charge.write(charges.data(), H5::PredType::NATIVE_DOUBLE);
        ds_charge.close();
        
        species_group.close();
        metadata_group.close();
        
        std::cout << "✓ Wrote metadata for " << n_species << " species to /metadata/species\n";
    } catch (const H5::Exception& e) {
        std::cerr << "HDF5 Error in write_species_metadata: " << e.getDetailMsg() << "\n";
    } catch (const std::exception& e) {
        std::cerr << "Error in write_species_metadata: " << e.what() << "\n";
    }
}

// -----------------------------
// Continue mode support
// -----------------------------

/**
 * @brief Write simulation completion metadata to HDF5 file.
 * 
 * Stores attributes indicating whether simulation completed successfully,
 * how many ions were still active at the end, RNG seed, git commit hash,
 * and timestamp for full reproducibility.
 */
void write_simulation_metadata(H5::H5File& file, bool complete, int active_ions, double final_time_s,
                               unsigned int rng_seed, const std::string& git_hash, 
                               const std::string& config_json) {
    H5::DataSpace scalar(H5S_SCALAR);
    
    // Simulation complete flag
    H5::Attribute attr_complete = file.createAttribute("simulation_complete", 
                                                       H5::PredType::NATIVE_HBOOL, 
                                                       scalar);
    attr_complete.write(H5::PredType::NATIVE_HBOOL, &complete);
    
    // Active ions at end
    H5::Attribute attr_active = file.createAttribute("active_ions_at_end",
                                                      H5::PredType::NATIVE_INT,
                                                      scalar);
    attr_active.write(H5::PredType::NATIVE_INT, &active_ions);
    
    // Final simulation time
    H5::Attribute attr_time = file.createAttribute("final_time_s",
                                                    H5::PredType::NATIVE_DOUBLE,
                                                    scalar);
    attr_time.write(H5::PredType::NATIVE_DOUBLE, &final_time_s);
    
    // RNG seed for reproducibility
    if (rng_seed != 0) {
        H5::Attribute attr_seed = file.createAttribute("rng_seed",
                                                       H5::PredType::NATIVE_UINT,
                                                       scalar);
        attr_seed.write(H5::PredType::NATIVE_UINT, &rng_seed);
    }
    
    // Git commit hash (if available)
    if (!git_hash.empty()) {
        H5::StrType str_type(H5::PredType::C_S1, git_hash.length() + 1);
        H5::Attribute attr_git = file.createAttribute("git_commit_hash", str_type, scalar);
        attr_git.write(str_type, git_hash.c_str());
    }
    
    // ISO 8601 timestamp
    auto now = std::chrono::system_clock::now();
    time_t now_time = std::chrono::system_clock::to_time_t(now);
    char timestamp_buf[64];
    strftime(timestamp_buf, sizeof(timestamp_buf), "%Y-%m-%dT%H:%M:%SZ", gmtime(&now_time));
    std::string timestamp_str(timestamp_buf);
    
    H5::StrType timestamp_type(H5::PredType::C_S1, timestamp_str.length() + 1);
    H5::Attribute attr_timestamp = file.createAttribute("timestamp_utc", timestamp_type, scalar);
    attr_timestamp.write(timestamp_type, timestamp_str.c_str());
    
    // ICARION version
    #ifdef ICARION_VERSION
        std::string version_str = ICARION_VERSION;
        H5::StrType version_type(H5::PredType::C_S1, version_str.length() + 1);
        H5::Attribute attr_version = file.createAttribute("icarion_version", version_type, scalar);
        attr_version.write(version_type, version_str.c_str());
    #endif
    
    // Compiler information
    std::string compiler_info;
    #ifdef __GNUC__
        compiler_info = "GCC " + std::to_string(__GNUC__) + "." + 
                       std::to_string(__GNUC_MINOR__) + "." + 
                       std::to_string(__GNUC_PATCHLEVEL__);
    #elif defined(__clang__)
        compiler_info = "Clang " + std::to_string(__clang_major__) + "." + 
                       std::to_string(__clang_minor__) + "." + 
                       std::to_string(__clang_patchlevel__);
    #elif defined(_MSC_VER)
        compiler_info = "MSVC " + std::to_string(_MSC_VER);
    #else
        compiler_info = "Unknown";
    #endif
    
    H5::StrType compiler_type(H5::PredType::C_S1, compiler_info.length() + 1);
    H5::Attribute attr_compiler = file.createAttribute("compiler", compiler_type, scalar);
    attr_compiler.write(compiler_type, compiler_info.c_str());
    
    // Build type
    std::string build_type;
    #ifdef ICARION_BUILD_CORE_ONLY
        build_type = "Core-Only";
    #else
        build_type = "Full";
    #endif
    #ifdef USE_CUDA
        build_type += " + GPU";
    #endif
    
    H5::StrType build_type_t(H5::PredType::C_S1, build_type.length() + 1);
    H5::Attribute attr_build = file.createAttribute("build_type", build_type_t, scalar);
    attr_build.write(build_type_t, build_type.c_str());
    
    // Configuration JSON (optional)
    if (!config_json.empty()) {
        H5::StrType config_type(H5::PredType::C_S1, config_json.length() + 1);
        H5::Attribute attr_config = file.createAttribute("configuration_json", config_type, scalar);
        attr_config.write(config_type, config_json.c_str());
    }
}

/**
 * @brief Load final ion states from HDF5 file for continuation.
 * 
 * Reads the last time slice of positions, velocities, and other ion properties
 * to resume simulation from a checkpoint.
 */
std::vector<IonState> load_final_state_from_HDF5(const std::string& filename) {
    std::vector<IonState> ions;
    
    try {
        H5::H5File file(filename, H5F_ACC_RDONLY);
        
        // Read positions dataset shape: [T, N, 3]
        H5::DataSet ds_pos = file.openDataSet("positions");
        H5::DataSpace fs_pos = ds_pos.getSpace();
        hsize_t dims[3];
        fs_pos.getSimpleExtentDims(dims);
        const hsize_t T = dims[0];
        const hsize_t N = dims[1];
        
        if (T == 0 || N == 0) {
            throw std::runtime_error("Empty HDF5 dataset");
        }
        
        // Read last time slice [T-1, :, :]
        hsize_t start[3] = {T - 1, 0, 0};
        hsize_t count[3] = {1, N, 3};
        H5::DataSpace memspace(3, count);
        fs_pos.selectHyperslab(H5S_SELECT_SET, count, start);
        
        std::vector<double> positions(N * 3);
        ds_pos.read(positions.data(), H5::PredType::NATIVE_DOUBLE, memspace, fs_pos);
        
        // Read velocities
        H5::DataSet ds_vel = file.openDataSet("velocities");
        H5::DataSpace fs_vel = ds_vel.getSpace();
        fs_vel.selectHyperslab(H5S_SELECT_SET, count, start);
        
        std::vector<double> velocities(N * 3);
        ds_vel.read(velocities.data(), H5::PredType::NATIVE_DOUBLE, memspace, fs_vel);
        
        // Read 2D datasets [T-1, :]
        hsize_t start2d[2] = {T - 1, 0};
        hsize_t count2d[2] = {1, N};
        H5::DataSpace memspace2d(2, count2d);
        
        H5::DataSet ds_mass = file.openDataSet("ion_mass_kg");
        H5::DataSpace fs_mass = ds_mass.getSpace();
        fs_mass.selectHyperslab(H5S_SELECT_SET, count2d, start2d);
        std::vector<double> masses(N);
        ds_mass.read(masses.data(), H5::PredType::NATIVE_DOUBLE, memspace2d, fs_mass);
        
        H5::DataSet ds_mob = file.openDataSet("reduced_mobility_cm2_Vs");
        H5::DataSpace fs_mob = ds_mob.getSpace();
        fs_mob.selectHyperslab(H5S_SELECT_SET, count2d, start2d);
        std::vector<double> mobilities(N);
        ds_mob.read(mobilities.data(), H5::PredType::NATIVE_DOUBLE, memspace2d, fs_mob);
        
        H5::DataSet ds_ccs = file.openDataSet("CCS_m2");
        H5::DataSpace fs_ccs = ds_ccs.getSpace();
        fs_ccs.selectHyperslab(H5S_SELECT_SET, count2d, start2d);
        std::vector<double> ccs(N);
        ds_ccs.read(ccs.data(), H5::PredType::NATIVE_DOUBLE, memspace2d, fs_ccs);
        
        H5::DataSet ds_active = file.openDataSet("active");
        H5::DataSpace fs_active = ds_active.getSpace();
        fs_active.selectHyperslab(H5S_SELECT_SET, count2d, start2d);
        std::vector<int> active_flags(N);
        ds_active.read(active_flags.data(), H5::PredType::NATIVE_INT, memspace2d, fs_active);
        
        H5::DataSet ds_born = file.openDataSet("born");
        H5::DataSpace fs_born = ds_born.getSpace();
        fs_born.selectHyperslab(H5S_SELECT_SET, count2d, start2d);
        std::vector<int> born_flags(N);
        ds_born.read(born_flags.data(), H5::PredType::NATIVE_INT, memspace2d, fs_born);
        
        // Read species (strings)
        H5::StrType str_type(H5::PredType::C_S1, H5T_VARIABLE);
        H5::DataSet ds_species = file.openDataSet("ion_species_id");
        H5::DataSpace fs_species = ds_species.getSpace();
        fs_species.selectHyperslab(H5S_SELECT_SET, count2d, start2d);
        std::vector<const char*> species_ptrs(N);
        ds_species.read(species_ptrs.data(), str_type, memspace2d, fs_species);
        
        // Construct IonState objects
        ions.reserve(N);
        for (hsize_t i = 0; i < N; ++i) {
            IonState ion;
            ion.pos.x = positions[i * 3 + 0];
            ion.pos.y = positions[i * 3 + 1];
            ion.pos.z = positions[i * 3 + 2];
            
            ion.vel.x = velocities[i * 3 + 0];
            ion.vel.y = velocities[i * 3 + 1];
            ion.vel.z = velocities[i * 3 + 2];
            
            ion.mass_kg = masses[i];
            ion.reduced_mobility_cm2_Vs = mobilities[i];
            ion.CCS_m2 = ccs[i];
            ion.active = (active_flags[i] != 0);
            ion.born = (born_flags[i] != 0);
            ion.species_id = species_ptrs[i];
            
            // Charge from mass (simple heuristic, could be improved)
            ion.ion_charge_C = ELEM_CHARGE_C;
            
            ions.push_back(ion);
        }
        
        // Free string memory
        for (auto ptr : species_ptrs) {
            if (ptr) free(const_cast<char*>(ptr));
        }
        
        file.close();
        
        std::cout << "Loaded " << ions.size() << " ions from " << filename << "\n";
        size_t active_count = std::count_if(ions.begin(), ions.end(), 
                                           [](const IonState& i){ return i.active; });
        std::cout << "  Active ions: " << active_count << " ("
                  << (100.0 * active_count / ions.size()) << "%)\n";
        
    } catch (const std::exception& e) {
        throw std::runtime_error("Failed to load HDF5 state: " + std::string(e.what()));
    }
    
    return ions;
}

}  // namespace io
}  // namespace ICARION
