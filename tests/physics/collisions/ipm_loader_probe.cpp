// ICARION: Ion Collision And Reaction IntegratiON
// SPDX-License-Identifier: GPL-3.0-only
#include "core/physics/collisions/InteractionPotentialOfflineSampleSet.h"
#include <iostream>

int main(int argc, char** argv) {
    if (argc != 2) {
        std::cerr << "usage: ipm_loader_probe FILE.h5\n";
        return 2;
    }
    ICARION::physics::InteractionPotentialOfflineSampleSet samples;
    std::string error;
    if (!ICARION::physics::load_interaction_potential_offline_sample_set_file(argv[1], samples, &error)) {
        std::cerr << error << '\n';
        return 1;
    }
    return 0;
}
