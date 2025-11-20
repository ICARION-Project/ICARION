/**
 * =====================================================================
 *
 *   ICARION: A Modular Framework for Ion Collision and Reaction Integration
 * =====================================================================
 *   A modular C++ framework for simulating ion trajectories 
 *   in user-defined electric fields and background gas environments.
 *
 *   @file        defineFields.cpp
 *   @brief       Calculates the electric fields.
 *
 *   @details
 *   Includes quadrupolar RF fields, axial DC fields, AC fields in
 *   arbitrary directions, and Orbitrap-like potentials.
 *
 *
 *   @date        2025-10-08
 *   @version     0.1
 *   @author      Christoph Schäfer
 *   @license     MIT License
 *
 * =====================================================================
 */
#include "defineFields.h"
#include <functional>
#include <cmath>
#include <algorithm>

/**
 * @brief Computes radial RF field in a quadrupole (x and y), zero along z.
 * 
 * @param[in] ion        Current ion state.
 * @param[in] RF_voltage RF amplitude [V].
 * @param[in] DC_voltage DC voltage offset [V].
 * @param[in] omega      RF angular frequency [rad/s].
 * @param[in] radius     Quadrupole characteristic radius [m].
 * @param[in] t          Current time [s].
 * @return Vec3 Electric field vector [V/m].
 */
Vec3 RFField(const IonState& ion, double RF_voltage, double DC_voltage, double omega, double radius, double t)
{
    // Radial RF field, axial component is zero
    double Ex =  2.0 * ion.pos.x * (DC_voltage + RF_voltage * cos(omega * t)) / (radius*radius);
    double Ey = -2.0 * ion.pos.y * (DC_voltage + RF_voltage * cos(omega * t)) / (radius*radius);
    double Ez = 0.0;

    return Vec3{Ex, Ey, Ez}; 
}

/**
 * @brief Computes axial DC field; x and y components are zero.
 * 
 * @param[in] ion     Current ion state.
 * @param[in] voltage Applied DC voltage [V].
 * @param[in] length  Characteristic axial length [m].
 * @return Vec3 Electric field vector [V/m].
 */

Vec3 DCField(const IonState& ion, double voltage, double length) 
{
    //no radial fields, only axial 
    double Ex = 0.0;
    double Ey = 0.0;
    double Ez = voltage / length;
    
    return Vec3{Ex, Ey, Ez};
}

/**
 * @brief Computes an AC field along a specified direction.
 *
 * The direction vector is normalized internally to ensure correct amplitude.
 *
 * @param[in] ion       Current ion state.
 * @param[in] voltage   AC amplitude [V].
 * @param[in] omega     AC angular frequency [rad/s].
 * @param[in] radius    Characteristic radius [m].
 * @param[in] t         Current time [s].
 * @param[in] direction Direction of the applied field (not normalized).
 * @return Vec3 Electric field vector [V/m].
 */
Vec3 ACField(const IonState& ion, double voltage, double omega, double radius, double t,
                    const Vec3& direction)
{
    Vec3 dir_unit = normalize(direction);
    double mag = -1.0 / radius * (voltage * cos(omega * t));
    return dir_unit * mag;
}

/**
 * @brief Computes the electric field in an Orbitrap-like potential.
 *
 * Radial and axial components follow a harmonic/quadrupole approximation.
 * Safeguard avoids division by zero for ions at the origin.
 *
 * @param[in] ion      Current ion state.
 * @param[in] k        Force constant [N/C or appropriate units].
 * @param[in] r_char   Characteristic radial distance [m].
 * @param[in] length_m Trap length along z [m].
 * @return Vec3 Electric field vector [V/m].
 */
Vec3 OrbitrapField(const IonState& ion, double k, double r_char, double length_m)
{
    double r2 = std::max(1e-18, ion.pos.x*ion.pos.x + ion.pos.y*ion.pos.y);
    double C = 1.0 - r_char * r_char / r2;
    double z_center = ion.pos.z - 0.5 * length_m; 
    
    double Ex = 0.5 * k * ion.pos.x * C;
    double Ey = 0.5 * k * ion.pos.y * C;
    double Ez = - k * z_center;

    return Vec3{Ex, Ey, Ez};
}

/**
 * @brief Computes the electric field in an FT-ICR-like potential.
 *
 * Quadrupole field approximation.
 *
 * @param[in] ion      Current ion state.
 * @param[in] voltage  Applied voltage [V].
 * @param[in] characteristic_length   Characteristic length [m].
 * @return Vec3 Electric field vector [V/m].
 */
Vec3 FTICRField(const IonState& ion, double voltage, double characteristic_length, double instrument_length_axial_m) {

    double factor = voltage / (characteristic_length * characteristic_length);
    double Ex =  ion.pos.x * factor;
    double Ey = ion.pos.y * factor;
    double Ez = - 2.0 * (ion.pos.z - 0.5 * instrument_length_axial_m) * factor;

    return Vec3{Ex, Ey, Ez};
}

/**
 * @brief Computes the magnetic field vector.
 *
 * Sums a constant field and a linear gradient component.
 *
 * @param[in] ion      Current ion state.
 * @param[in] dom      Instrument domain parameters.
 * @return Vec3 Magnetic field vector [T].
 */
Vec3 MagneticFieldVec(const IonState& ion, const InstrumentDomain& dom) {
    if (!dom.B.enabled) return Vec3{0.0, 0.0, 0.0};

    //linear model: B = B0 + grad * r
    return dom.B.field_strength_T + Vec3{
        dom.B.field_gradient_T_m.x * ion.pos.x,
        dom.B.field_gradient_T_m.y * ion.pos.y,
        dom.B.field_gradient_T_m.z * ion.pos.z
    };
}

/**
 * @brief Computes the total electric field vector at the ion's position and time.
 *
 * Sums contributions from DC, RF, AC, and Orbitrap-like fields as applicable.
 * More modular approach allows easy extension with additional field types instead of
 * only using the hardcoded instrument types for more flexibility.
 *
 * @param[in] ion        Current ion state.
 * @param[in] dom        Instrument domain parameters.
 * @param[in] t_global   Current global time [s].
 * @return Vec3 Total electric field vector [V/m].
 */
Vec3 ElectricFieldVec(const IonState& ion, const InstrumentDomain& dom, double t_global) {
    Vec3 E_total{0.0, 0.0, 0.0};
    std::vector<std::function<Vec3(const IonState&, double)>> fieldModules;

    // --- Add field modules independent of instrument type ---
    // DC field (axial), if DC voltage is non-zero
    if (std::fabs(dom.DC.axial_V) > 1e-12) {
        fieldModules.push_back([&dom](const IonState& ion, double t) {
            return DCField(ion, dom.DC.axial_V, dom.geom.length_m); // in V/m
        });
    }

    // RF field (radial), if RF voltage is non-zero
    if (std::fabs(dom.RF.voltage_V) > 1e-12) {
        fieldModules.push_back([&dom](const IonState& ion, double t) {
            return RFField(ion, dom.RF.voltage_V, dom.DC.quad_V,
                                        dom.RF.angular_frequency_rad_s, dom.geom.radius_m, t);
        }); 
    }

    // AC field, if AC voltage is non-zero
    if (std::fabs(dom.AC.voltage_V) > 1e-12) {
        // Example direction: along x-axis; can be modified as needed
        fieldModules.push_back([&dom](const IonState& ion, double t) {
            return ACField(ion, dom.AC.voltage_V,
                                        dom.AC.angular_frequency_rad_s, dom.geom.radius_m, t, Vec3(1, 0, 0));
        });
    }   

    // --- Compute total electric field ---
    for (const auto& fieldModule : fieldModules) {
        E_total += fieldModule(ion, t_global);
    }

    return E_total;
}