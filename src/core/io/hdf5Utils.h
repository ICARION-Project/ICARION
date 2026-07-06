// ICARION: Ion Collision And Reaction IntegratiON
// SPDX-License-Identifier: GPL-3.0-only
// Copyright (c) 2026 Christoph Schaefer
#pragma once

#include <H5Cpp.h>
#include <string>

namespace ICARION::io {

std::string read_hdf5_attr_string(H5::H5File& file, const std::string& name);
bool read_hdf5_attr_long_long(H5::H5File& file, const std::string& name, long long& value);

} // namespace ICARION::io
