// ICARION: Ion Collision And Reaction IntegratiON
// SPDX-License-Identifier: GPL-3.0-only
#include "UserAnnotation.h"
#include "core/utils/hash.h"

#include <cstdint>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <utility>

namespace ICARION::io {
namespace {
bool is_valid_utf8(const std::string& value) {
    const auto* bytes = reinterpret_cast<const unsigned char*>(value.data());
    std::size_t i = 0;
    while (i < value.size()) {
        const unsigned char first = bytes[i++];
        if (first <= 0x7f) continue;

        std::uint32_t codepoint = 0;
        std::size_t continuation_count = 0;
        std::uint32_t minimum = 0;
        if (first >= 0xc2 && first <= 0xdf) {
            codepoint = first & 0x1f;
            continuation_count = 1;
            minimum = 0x80;
        } else if (first >= 0xe0 && first <= 0xef) {
            codepoint = first & 0x0f;
            continuation_count = 2;
            minimum = 0x800;
        } else if (first >= 0xf0 && first <= 0xf4) {
            codepoint = first & 0x07;
            continuation_count = 3;
            minimum = 0x10000;
        } else {
            return false;
        }
        if (i + continuation_count > value.size()) return false;
        for (std::size_t j = 0; j < continuation_count; ++j) {
            const unsigned char continuation = bytes[i++];
            if ((continuation & 0xc0) != 0x80) return false;
            codepoint = (codepoint << 6) | (continuation & 0x3f);
        }
        if (codepoint < minimum || (codepoint >= 0xd800 && codepoint <= 0xdfff) ||
            codepoint > 0x10ffff) {
            return false;
        }
    }
    return true;
}

void validate_note_bytes(const std::string& note) {
    if (note.empty()) throw std::runtime_error("User annotation note must not be empty");
    if (note.size() > USER_ANNOTATION_MAX_BYTES)
        throw std::runtime_error("User annotation note exceeds 65,536 UTF-8 bytes");
    if (note.find('\0') != std::string::npos)
        throw std::runtime_error("User annotation note contains an embedded NUL byte");
    if (!is_valid_utf8(note))
        throw std::runtime_error("User annotation note is not valid UTF-8");
}

void write_string(H5::Group& group, const std::string& name, const std::string& value) {
    H5::StrType type(H5::PredType::C_S1, H5T_VARIABLE);
    type.setCset(H5T_CSET_UTF8);
    H5::DataSpace space(H5S_SCALAR);
    auto dataset = group.createDataSet(name, type, space);
    const char* ptr = value.c_str();
    dataset.write(&ptr, type);
}

std::string read_string(H5::Group& group, const std::string& name) {
    if (H5Lexists(group.getId(), name.c_str(), H5P_DEFAULT) <= 0)
        throw std::runtime_error("User annotation integrity error: missing dataset '" + name + "'");
    auto dataset = group.openDataSet(name);
    if (dataset.getTypeClass() != H5T_STRING)
        throw std::runtime_error("User annotation integrity error: dataset '" + name + "' is not a string");
    auto space = dataset.getSpace();
    if (space.getSimpleExtentType() != H5S_SCALAR)
        throw std::runtime_error("User annotation integrity error: dataset '" + name + "' is not scalar");
    auto type = dataset.getStrType();
    if (!type.isVariableStr())
        throw std::runtime_error("User annotation integrity error: dataset '" + name + "' is not variable-length");
    if (type.getCset() != H5T_CSET_UTF8)
        throw std::runtime_error("User annotation integrity error: dataset '" + name + "' is not declared UTF-8");
    char* ptr = nullptr;
    dataset.read(&ptr, type);
    std::string value = ptr ? ptr : "";
    if (ptr) H5free_memory(ptr);
    return value;
}
} // namespace

UserAnnotation resolve_inline_annotation(std::string note) {
    UserAnnotation annotation{true, std::move(note), "inline", "", ""};
    annotation.note_sha256 = utils::sha256_bytes(annotation.note);
    validate_annotation(annotation);
    return annotation;
}

UserAnnotation resolve_file_annotation(const std::filesystem::path& path,
                                       const std::filesystem::path& resolution_base) {
    const auto resolved = path.is_absolute() || resolution_base.empty() ? path : resolution_base / path;
    std::ifstream input(resolved, std::ios::binary);
    if (!input) throw std::runtime_error("Cannot open user annotation file: " + resolved.string());
    std::ostringstream bytes;
    bytes << input.rdbuf();
    if (input.bad()) throw std::runtime_error("Cannot read user annotation file: " + resolved.string());
    auto note = bytes.str();
    UserAnnotation annotation{true, std::move(note), "file", path.filename().string(), ""};
    annotation.note_sha256 = utils::sha256_bytes(annotation.note);
    validate_annotation(annotation);
    return annotation;
}

void validate_annotation(const UserAnnotation& annotation) {
    if (!annotation.present) return;
    validate_note_bytes(annotation.note);
    if (annotation.source != "inline" && annotation.source != "file")
        throw std::runtime_error("Invalid user annotation source");
    if (annotation.source == "inline" && !annotation.source_filename.empty())
        throw std::runtime_error("Inline user annotation has a source filename");
    if (annotation.source == "file" && annotation.source_filename.empty())
        throw std::runtime_error("File user annotation lacks a source filename");
    if (!annotation.source_filename.empty() && !is_valid_utf8(annotation.source_filename))
        throw std::runtime_error("User annotation source filename is not valid UTF-8");
    if (annotation.source == "file") {
        if (annotation.source_filename.find('/') != std::string::npos ||
            annotation.source_filename.find('\\') != std::string::npos ||
            annotation.source_filename == "." || annotation.source_filename == ".." ||
            std::filesystem::path(annotation.source_filename).filename().string() != annotation.source_filename) {
            throw std::runtime_error("User annotation source filename must be a basename");
        }
    }
    if (annotation.note_sha256 != utils::sha256_bytes(annotation.note))
        throw std::runtime_error("User annotation SHA-256 does not match note content");
}

void write_annotation_hdf5(H5::Group& metadata_group, const UserAnnotation& annotation) {
    if (!annotation.present) return;
    validate_annotation(annotation);
    auto group = metadata_group.createGroup("annotations");
    write_string(group, "note", annotation.note);
    write_string(group, "source", annotation.source);
    write_string(group, "source_filename", annotation.source_filename);
    write_string(group, "note_sha256", annotation.note_sha256);
}

UserAnnotation read_annotation_hdf5(H5::Group& metadata_group) {
    try {
        const htri_t annotation_status = H5Lexists(metadata_group.getId(), "annotations", H5P_DEFAULT);
        if (annotation_status < 0)
            throw std::runtime_error("User annotation integrity error: cannot inspect annotations group");
        if (annotation_status == 0) return {};
        auto group = metadata_group.openGroup("annotations");
        UserAnnotation result{true, read_string(group, "note"), read_string(group, "source"),
                              read_string(group, "source_filename"), read_string(group, "note_sha256")};
        validate_annotation(result);
        return result;
    } catch (const H5::Exception& e) {
        throw std::runtime_error("User annotation integrity error: " + e.getDetailMsg());
    }
}
} // namespace ICARION::io
