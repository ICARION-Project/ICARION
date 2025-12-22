// ICARION: Ion Collision And Reaction IntegratiON
// MIT License - Copyright (c) 2025 ICARION Project Contributors

#pragma once
#include <string>
#include "trajsim/fieldArrayLoader.h"

// Write a FieldArray to HDF5 file
void write_fieldarray_hdf5(const std::string& path, const FieldArray& fld);
