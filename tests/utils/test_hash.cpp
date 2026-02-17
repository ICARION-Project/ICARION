// ICARION: Ion Collision And Reaction IntegratiON
// SPDX-License-Identifier: GPL-3.0-only
// Copyright (c) 2025 ICARION Project Contributors
/**
 * @file test_hash.cpp
 * @brief Unit tests for SHA256 file hashing utility
 */

#include <catch2/catch_test_macros.hpp>
#include "core/utils/hash.h"
#include <fstream>
#include <filesystem>

using namespace ICARION::utils;
namespace fs = std::filesystem;

TEST_CASE("SHA256 hashing works correctly", "[hash][utils]") {
    // Create temporary test file
    std::string test_file = "/tmp/test_hash_file.txt";
    
    SECTION("Hash of known content") {
        // Write known content
        std::ofstream ofs(test_file);
        ofs << "Hello, World!";
        ofs.close();
        
        // Compute hash
        std::string hash = sha256_file(test_file);
        
        // Verify hash format (64 hex characters)
        REQUIRE(hash.length() == 64);
        REQUIRE(hash.find_first_not_of("0123456789abcdef") == std::string::npos);
        
        // Verify against known SHA256 of "Hello, World!"
        // echo -n "Hello, World!" | sha256sum
        // dffd6021bb2bd5b0af676290809ec3a53191dd81c7f70a4b28688a362182986f
        REQUIRE(hash == "dffd6021bb2bd5b0af676290809ec3a53191dd81c7f70a4b28688a362182986f");
        
        // Cleanup
        fs::remove(test_file);
    }
    
    SECTION("Hash is deterministic") {
        // Write content
        std::ofstream ofs(test_file);
        ofs << "Test content for determinism";
        ofs.close();
        
        // Compute hash twice
        std::string hash1 = sha256_file(test_file);
        std::string hash2 = sha256_file(test_file);
        
        // Should be identical
        REQUIRE(hash1 == hash2);
        
        fs::remove(test_file);
    }
    
    SECTION("Different content produces different hash") {
        // Write first content
        {
            std::ofstream ofs(test_file);
            ofs << "Content A";
        }
        std::string hash_a = sha256_file(test_file);
        
        // Write different content
        {
            std::ofstream ofs(test_file);
            ofs << "Content B";
        }
        std::string hash_b = sha256_file(test_file);
        
        // Hashes should differ
        REQUIRE(hash_a != hash_b);
        
        fs::remove(test_file);
    }
    
    SECTION("Empty file hash") {
        // Create empty file
        std::ofstream(test_file).close();
        
        std::string hash = sha256_file(test_file);
        
        // SHA256 of empty string is well-known
        // echo -n "" | sha256sum
        REQUIRE(hash == "e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855");
        
        fs::remove(test_file);
    }
    
    SECTION("Large file handling") {
        // Create 10MB file
        {
            std::ofstream ofs(test_file, std::ios::binary);
            std::vector<char> data(10 * 1024 * 1024, 'X');
            ofs.write(data.data(), data.size());
        }
        
        // Should not crash or timeout
        std::string hash = sha256_file(test_file);
        REQUIRE(hash.length() == 64);
        
        fs::remove(test_file);
    }
}

TEST_CASE("SHA256 error handling", "[hash][utils][error]") {
    SECTION("Non-existent file throws exception") {
        REQUIRE_THROWS_AS(sha256_file("/nonexistent/file.txt"), std::runtime_error);
    }
    
    SECTION("Safe wrapper returns default on error") {
        std::string hash = sha256_file_safe("/nonexistent/file.txt", "N/A");
        REQUIRE(hash == "N/A");
    }
    
    SECTION("Safe wrapper with custom default") {
        std::string hash = sha256_file_safe("/nonexistent/file.txt", "FILE_NOT_FOUND");
        REQUIRE(hash == "FILE_NOT_FOUND");
    }
    
    SECTION("Safe wrapper returns ERROR by default") {
        std::string hash = sha256_file_safe("/nonexistent/file.txt");
        REQUIRE(hash == "ERROR");
    }
}

TEST_CASE("SHA256 integration with HDF5Writer", "[hash][integration]") {
    // This tests the actual usage pattern in HDF5Writer v2
    std::string config_file = "/tmp/test_config.json";
    
    // Create mock config
    std::ofstream ofs(config_file);
    ofs << R"({"simulation": {"dt_s": 1e-9}})";
    ofs.close();
    
    // Compute hash
    std::string hash = sha256_file_safe(config_file);
    
    // Should be valid hash, not error
    REQUIRE(hash != "ERROR");
    REQUIRE(hash != "N/A");
    REQUIRE(hash.length() == 64);
    
    // Verify reproducibility
    std::string hash2 = sha256_file_safe(config_file);
    REQUIRE(hash == hash2);
    
    fs::remove(config_file);
}
