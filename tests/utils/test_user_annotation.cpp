// ICARION: Ion Collision And Reaction IntegratiON
// SPDX-License-Identifier: GPL-3.0-only
#include <catch2/catch_test_macros.hpp>
#include "core/io/UserAnnotation.h"

#include <filesystem>
#include <fstream>
#include <chrono>
#include <functional>

using namespace ICARION::io;

namespace {
class TemporaryDirectory {
public:
    TemporaryDirectory() : path_(std::filesystem::temp_directory_path() /
        ("icarion_annotation_" + std::to_string(
            std::chrono::steady_clock::now().time_since_epoch().count()))) {
        std::filesystem::create_directories(path_);
    }
    ~TemporaryDirectory() { std::error_code error; std::filesystem::remove_all(path_, error); }
    const std::filesystem::path& path() const { return path_; }
private:
    std::filesystem::path path_;
};

void replace_with_vlen_string(H5::Group& group, const char* name, const std::string& value,
                              H5T_cset_t charset = H5T_CSET_UTF8) {
    H5Ldelete(group.getId(), name, H5P_DEFAULT);
    H5::StrType type(H5::PredType::C_S1, H5T_VARIABLE);
    type.setCset(charset);
    auto dataset = group.createDataSet(name, type, H5::DataSpace(H5S_SCALAR));
    const char* pointer = value.c_str();
    dataset.write(&pointer, type);
}

void require_corrupt_annotation_rejected(const std::function<void(H5::Group&)>& corrupt) {
    TemporaryDirectory temporary;
    const auto path = temporary.path() / "corrupt.h5";
    {
        H5::H5File file(path.string(), H5F_ACC_TRUNC);
        auto metadata = file.createGroup("/metadata");
        write_annotation_hdf5(metadata, resolve_inline_annotation("valid note"));
        auto annotations = metadata.openGroup("annotations");
        corrupt(annotations);
    }
    H5::H5File file(path.string(), H5F_ACC_RDONLY);
    auto metadata = file.openGroup("/metadata");
    REQUIRE_THROWS_AS(read_annotation_hdf5(metadata), std::runtime_error);
}
} // namespace

TEST_CASE("User annotations preserve and validate exact bytes", "[annotation][io]") {
    SECTION("inline and multiline notes") {
        auto single = resolve_inline_annotation("known content");
        REQUIRE(single.present);
        REQUIRE(single.source == "inline");
        REQUIRE(single.source_filename.empty());
        REQUIRE(single.note_sha256 == "41277d8d0b0610e58f13bdc06b732c629a2fd3ff93c382f40af3f60cfe5e5c9e");
        auto multiline = resolve_inline_annotation(" line one \r\nline two\n");
        REQUIRE(multiline.note == " line one \r\nline two\n");
        const auto unicode = resolve_inline_annotation("Grüße – 東京 – 😀");
        REQUIRE(unicode.note == "Grüße – 東京 – 😀");
    }
    SECTION("file notes use basename and exact bytes") {
        TemporaryDirectory temporary;
        std::filesystem::create_directories(temporary.path() / "nested");
        const auto path = temporary.path() / "nested" / "note.md";
        std::ofstream(path, std::ios::binary) << "file\r\nnote\n";
        auto note = resolve_file_annotation("nested/note.md", temporary.path());
        REQUIRE(note.note == "file\r\nnote\n");
        REQUIRE(note.source == "file");
        REQUIRE(note.source_filename == "note.md");
    }
    SECTION("invalid sizes and NUL are rejected") {
        REQUIRE_THROWS(resolve_inline_annotation(""));
        REQUIRE_THROWS(resolve_inline_annotation(std::string("a\0b", 3)));
        REQUIRE_NOTHROW(resolve_inline_annotation(std::string(USER_ANNOTATION_MAX_BYTES, 'a')));
        REQUIRE_THROWS(resolve_inline_annotation(std::string(USER_ANNOTATION_MAX_BYTES + 1, 'a')));
        TemporaryDirectory temporary;
        const auto empty = temporary.path() / "empty.txt";
        std::ofstream(empty, std::ios::binary);
        REQUIRE_THROWS(resolve_file_annotation(empty));
    }
    SECTION("invalid UTF-8 is rejected") {
        REQUIRE_THROWS(resolve_inline_annotation(std::string("\x80", 1)));
        REQUIRE_THROWS(resolve_inline_annotation(std::string("\xc2", 1)));
        REQUIRE_THROWS(resolve_inline_annotation(std::string("\xc0\xaf", 2)));
        REQUIRE_THROWS(resolve_inline_annotation(std::string("\xed\xa0\x80", 3)));
        REQUIRE_THROWS(resolve_inline_annotation(std::string("\xf4\x90\x80\x80", 4)));
        REQUIRE_THROWS(resolve_inline_annotation(std::string("\xe2\x28\xa1", 3)));
    }
    SECTION("file provenance must be a UTF-8 basename") {
        const auto valid = resolve_inline_annotation("content");
        for (const auto& filename : {"../note.txt", "/tmp/note.txt", "dir/note.txt",
                                     "dir\\note.txt", ".", ".."}) {
            auto annotation = valid;
            annotation.source = "file";
            annotation.source_filename = filename;
            REQUIRE_THROWS(validate_annotation(annotation));
        }
        auto invalid_utf8 = valid;
        invalid_utf8.source = "file";
        invalid_utf8.source_filename = std::string("bad\xff.txt", 8);
        REQUIRE_THROWS(validate_annotation(invalid_utf8));
    }
    SECTION("common HDF5 layout round trips") {
        TemporaryDirectory temporary;
        const auto path = temporary.path() / "annotation.h5";
        {
            H5::H5File file(path.string(), H5F_ACC_TRUNC);
            auto metadata = file.createGroup("/metadata");
            write_annotation_hdf5(metadata, resolve_inline_annotation("hello\nworld"));
        }
        {
            H5::H5File file(path.string(), H5F_ACC_RDONLY);
            auto metadata = file.openGroup("/metadata");
            auto restored = read_annotation_hdf5(metadata);
            REQUIRE(restored.note == "hello\nworld");
            REQUIRE(restored.source == "inline");
        }
    }
}

TEST_CASE("Malformed HDF5 annotations are rejected safely", "[annotation][io][corruption]") {
    SECTION("missing dataset") {
        require_corrupt_annotation_rejected([](H5::Group& group) {
            H5Ldelete(group.getId(), "source", H5P_DEFAULT);
        });
    }
    SECTION("numeric note") {
        require_corrupt_annotation_rejected([](H5::Group& group) {
            H5Ldelete(group.getId(), "note", H5P_DEFAULT);
            auto dataset = group.createDataSet("note", H5::PredType::NATIVE_INT, H5::DataSpace(H5S_SCALAR));
            int value = 1;
            dataset.write(&value, H5::PredType::NATIVE_INT);
        });
    }
    SECTION("nonscalar note") {
        require_corrupt_annotation_rejected([](H5::Group& group) {
            H5Ldelete(group.getId(), "note", H5P_DEFAULT);
            H5::StrType type(H5::PredType::C_S1, H5T_VARIABLE);
            type.setCset(H5T_CSET_UTF8);
            hsize_t dimension = 1;
            auto dataset = group.createDataSet("note", type, H5::DataSpace(1, &dimension));
            const char* value = "valid note";
            dataset.write(&value, type);
        });
    }
    SECTION("fixed-length note") {
        require_corrupt_annotation_rejected([](H5::Group& group) {
            H5Ldelete(group.getId(), "note", H5P_DEFAULT);
            H5::StrType type(H5::PredType::C_S1, 11);
            type.setCset(H5T_CSET_UTF8);
            auto dataset = group.createDataSet("note", type, H5::DataSpace(H5S_SCALAR));
            const char value[11] = "valid note";
            dataset.write(value, type);
        });
    }
    SECTION("non-UTF-8 declaration") {
        require_corrupt_annotation_rejected([](H5::Group& group) {
            replace_with_vlen_string(group, "note", "valid note", H5T_CSET_ASCII);
        });
    }
    SECTION("invalid UTF-8 note bytes") {
        require_corrupt_annotation_rejected([](H5::Group& group) {
            replace_with_vlen_string(group, "note", std::string("bad\xff", 4));
        });
    }
    SECTION("invalid hash") {
        require_corrupt_annotation_rejected([](H5::Group& group) {
            replace_with_vlen_string(group, "note_sha256", std::string(64, '0'));
        });
    }
}
