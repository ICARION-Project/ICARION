// ICARION: Ion Collision And Reaction IntegratiON
// SPDX-License-Identifier: GPL-3.0-only
#pragma once

#include <H5Cpp.h>
#include <filesystem>
#include <string>

namespace ICARION::io {

struct UserAnnotation {
    bool present = false;
    std::string note;
    std::string source;
    std::string source_filename;
    std::string note_sha256;
};

constexpr std::size_t USER_ANNOTATION_MAX_BYTES = 65536;

UserAnnotation resolve_inline_annotation(std::string note);
UserAnnotation resolve_file_annotation(const std::filesystem::path& path,
                                       const std::filesystem::path& resolution_base = {});
void validate_annotation(const UserAnnotation& annotation);
void write_annotation_hdf5(H5::Group& metadata_group, const UserAnnotation& annotation);
UserAnnotation read_annotation_hdf5(H5::Group& metadata_group);

} // namespace ICARION::io
