// ICARION: Ion Collision And Reaction IntegratiON
// SPDX-License-Identifier: GPL-3.0-only
// Copyright (c) 2025 Christoph Schaefer
#include "FieldProviderModel.h"
#include "fieldsolver/utils/IFieldProvider.h"

namespace ICARION::config {

Vec3 FieldProviderModel::E(const Vec3& global_pos, double t) const {
    if (!provider_) {
        return Vec3{0.0, 0.0, 0.0};
    }
    return provider_->get_E(global_pos, t);
}

} // namespace ICARION::config
