// ICARION: Ion Collision And Reaction IntegratiON
// MIT License - Copyright (c) 2025 ICARION Project Contributors

#include "defineSpaceChargeField.h"
#include "utils/constants.h"
#include <cmath>
#include <omp.h>

/**
 * @brief Compute total Coulombic space-charge field on a single ion
 *        from all other active ions (direct N² method).
 *
 * @param ion   Ion for which the field is computed
 * @param ions  Vector of all ions in the simulation
 * @param eps0  Permittivity of free space
 * @return Vec3 Electric field [V/m] due to other ions
 */
Vec3 SpaceChargeField(const IonState& ion,
                      const std::vector<IonState>& ions,
                      double eps0)
{
    Vec3 E_total{0.0, 0.0, 0.0};

    // Parallelized accumulation
    double Ex_sum = 0.0, Ey_sum = 0.0, Ez_sum = 0.0;

    //#pragma omp parallel for reduction(+:Ex_sum,Ey_sum,Ez_sum)
    for (int j = 0; j < static_cast<int>(ions.size()); ++j)
    {
        const IonState& other = ions[j];
        if (&other == &ion || !other.active || !other.born)
            continue;

        Vec3 dr = ion.pos - other.pos;
        double r2 = dr.x*dr.x + dr.y*dr.y + dr.z*dr.z + 1e-18; // avoid div0
        double inv_r3 = 1.0 / (r2 * std::sqrt(r2));
        double coeff = other.ion_charge_C / (4.0 * M_PI * eps0);
        E_total.x += coeff * dr.x * inv_r3;
        E_total.y += coeff * dr.y * inv_r3;
        E_total.z += coeff * dr.z * inv_r3;
    }

    return E_total;
}
