// ICARION: Ion Collision And Reaction IntegratiON
// SPDX-License-Identifier: GPL-3.0-only
// Copyright (c) 2026 Christoph Schaefer
#ifndef ICARION_CONFIG_OUTPUT_ENUMS_H
#define ICARION_CONFIG_OUTPUT_ENUMS_H

namespace ICARION::config {

enum class DeepAnalysisMode {
    Off,
    Summary,
    SampledEvents,
    FullEvents
};

} // namespace ICARION::config

#endif // ICARION_CONFIG_OUTPUT_ENUMS_H
