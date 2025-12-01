#include <iostream>
#include "src/core/io/moleculeLoader.h"
int main() {
    try {
        std::cout << "Loading H3O+ geometry..." << std::endl;
        auto molecule = ICARION::io::load_molecule("data/molecules/H3O+.json");
        std::cout << "Success! Molecule: " << molecule.name << " with " << molecule.atoms.size() << " atoms" << std::endl;
        return 0;
    } catch (const std::exception& e) {
        std::cout << "Error: " << e.what() << std::endl;
        return 1;
    }
}
