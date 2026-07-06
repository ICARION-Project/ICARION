// ICARION: Ion Collision And Reaction IntegratiON
// SPDX-License-Identifier: GPL-3.0-only
// Copyright (c) 2026 Christoph Schaefer
#include "hdf5Utils.h"

namespace ICARION::io {

std::string read_hdf5_attr_string(H5::H5File& file, const std::string& name) {
    if (!file.attrExists(name)) {
        return {};
    }
    H5::Attribute attr = file.openAttribute(name);
    H5::StrType type = attr.getStrType();
    std::string value;
    attr.read(type, value);
    return value;
}

bool read_hdf5_attr_long_long(H5::H5File& file, const std::string& name, long long& value) {
    if (!file.attrExists(name)) {
        return false;
    }
    H5::Attribute attr = file.openAttribute(name);
    attr.read(H5::PredType::NATIVE_LLONG, &value);
    return true;
}

} // namespace ICARION::io
