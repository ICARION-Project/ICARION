// ICARION: Ion Collision And Reaction IntegratiON
// SPDX-License-Identifier: GPL-3.0-only
// Copyright (c) 2026 Christoph Schaefer
#ifndef ICARION_CONFIG_UNIT_CONVERTER_H
#define ICARION_CONFIG_UNIT_CONVERTER_H

#include "utils/constants.h"
#include <cmath>

namespace ICARION::config {

/**
 * @brief Centralized unit conversion utilities
 * 
 * All unit conversions in one place for consistency and testability.
 * Uses constants from utils/constants.h.
 */
class UnitConverter {
public:
    // === Electric field ===
    
    /**
     * @brief Convert reduced field strength from Townsend to V·m² (EN -> EN) (SI)
     * 
     * @param EN_Td Field strength in Townsend [Td]
     * @return Reduced field strength in V·m² [V·m²]
     */
    static constexpr double townsend_to_Vm2(double EN_Td) {
        return EN_Td * 1e-21;  // 1 Td = 10^-21 V·m²
    }
    
    /**
     * @brief Convert E/N to voltage given particle density and length (EN -> U) (SI)
     * 
     * @param EN_Vm2 Reduced field [V·m²]
     * @param density_m3 Particle number density [m⁻³]
     * @param length_m Characteristic length [m]
     * @return Voltage [V]
     */
    static double Vm2_to_V(double EN_Vm2, double density_m3, double length_m) {
        return EN_Vm2 * density_m3 * length_m;
    }
    
    /**
     * @brief Convert Townsend to V/m given particle density (EN -> E) (SI)
     * 
     * @param EN_Td Reduced field [Td]
     * @param density_m3 Particle number density [m⁻³]
     * @return Electric field [V/m]
     */
    static double townsend_to_Vm(double EN_Td, double density_m3) {
        return townsend_to_Vm2(EN_Td) * density_m3;
    }
    
    // === Frequency ===
    
    /**
     * @brief Convert frequency from Hz to rad/s (SI)
     * 
     * @param freq_Hz Frequency [Hz]
     * @return Angular frequency [rad/s]
     */
    static constexpr double Hz_to_rad_s(double freq_Hz) {
        return freq_Hz * 2.0 * M_PI;
    }
    
    /**
     * @brief Convert angular frequency from rad/s to Hz (SI)
     * 
     * @param omega_rad_s Angular frequency [rad/s]
     * @return Frequency [Hz]
     */
    static constexpr double rad_s_to_Hz(double omega_rad_s) {
        return omega_rad_s / (2.0 * M_PI);
    }
    
    // === Mobility ===
    
    /**
     * @brief Convert mobility from m²/(V·s) to cm²/(V·s) (both SI, but cm²/(V·s) reported in literature)
     * 
     * @param mobility_SI Mobility [m²/(V·s)]
     * @return Mobility [cm²/(V·s)]
     */
    static constexpr double m2Vs_to_cm2Vs(double mobility_SI) {
        return mobility_SI * 1e4;
    }
    
    /**
     * @brief Convert mobility from cm²/(V·s) to m²/(V·s) (both SI, but cm²/(V·s) reported in literature)
     * 
     * @param mobility_cgs Mobility [cm²/(V·s)]
     * @return Mobility [m²/(V·s)]
     */
    static constexpr double cm2Vs_to_m2Vs(double mobility_cgs) {
        return mobility_cgs * 1e-4;
    }
    
    // === Energy ===
    
    /**
     * @brief Convert electron volts to Joules
     * 
     * @param eV Energy [eV]
     * @return Energy [J]
     */
    static constexpr double eV_to_J(double eV) {
        return eV * ELEM_CHARGE_C;
    }
    
    /**
     * @brief Convert Joules to electron volts
     * 
     * @param J Energy [J]
     * @return Energy [eV]
     */
    static constexpr double J_to_eV(double J) {
        return J / ELEM_CHARGE_C;
    }
    
    // === Temperature ===
    
    /**
     * @brief Convert temperature from Celsius to Kelvin
     * 
     * @param celsius Temperature [°C]
     * @return Temperature [K]
     */
    static constexpr double celsius_to_kelvin(double celsius) {
        return celsius + 273.15;
    }
    
    /**
     * @brief Convert temperature from Kelvin to Celsius
     * 
     * @param kelvin Temperature [K]
     * @return Temperature [°C]
     */
    static constexpr double kelvin_to_celsius(double kelvin) {
        return kelvin - 273.15;
    }
};

} // namespace ICARION::config

#endif // ICARION_CONFIG_UNIT_CONVERTER_H
