// ICARION: Ion Collision And Reaction IntegratiON
// SPDX-License-Identifier: GPL-3.0-only
// Copyright (c) 2025 ICARION Project Contributors
#pragma once
#include <string>

#include "core/io/fieldArrayLoader.h"

// Write a FieldArray to HDF5 file
void write_fieldarray_hdf5(const std::string& path, const FieldArray& fld);
