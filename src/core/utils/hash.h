// ICARION: Ion Collision And Reaction IntegratiON
// MIT License - Copyright (c) 2025 ICARION Project Contributors

/**
 * @file hash.h
 * @brief File hashing utilities for reproducibility verification
 * @date 2025-11-21
 */

#pragma once

#include <string>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <stdexcept>
#include <openssl/sha.h>

namespace ICARION::utils {

/**
 * @brief Compute SHA256 hash of a file
 * 
 * This function is critical for ensuring reproducibility by verifying
 * that input files (config, species DB, reactions DB) have not been modified.
 * 
 * @param filepath Path to the file to hash
 * @return Hex-encoded SHA256 hash (64 characters)
 * @throws std::runtime_error If file cannot be opened or read
 * 
 * @example
 * ```cpp
 * std::string hash = utils::sha256_file("config.json");
 * // Returns: "e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855"
 * ```
 */
inline std::string sha256_file(const std::string& filepath) {
    // Open file in binary mode
    std::ifstream file(filepath, std::ios::binary);
    if (!file.is_open()) {
        throw std::runtime_error("Cannot open file for SHA256 hashing: " + filepath);
    }
    
    // Initialize SHA256 context
    SHA256_CTX ctx;
    if (SHA256_Init(&ctx) != 1) {
        throw std::runtime_error("SHA256_Init failed");
    }
    
    // Read file in 8KB chunks and update hash
    constexpr size_t BUFFER_SIZE = 8192;
    char buffer[BUFFER_SIZE];
    
    while (file.read(buffer, BUFFER_SIZE) || file.gcount() > 0) {
        if (SHA256_Update(&ctx, buffer, file.gcount()) != 1) {
            throw std::runtime_error("SHA256_Update failed");
        }
    }
    
    if (file.bad()) {
        throw std::runtime_error("Error reading file for SHA256 hashing: " + filepath);
    }
    
    // Finalize hash computation
    unsigned char hash[SHA256_DIGEST_LENGTH];
    if (SHA256_Final(hash, &ctx) != 1) {
        throw std::runtime_error("SHA256_Final failed");
    }
    
    // Convert binary hash to hex string
    std::ostringstream oss;
    for (int i = 0; i < SHA256_DIGEST_LENGTH; ++i) {
        oss << std::hex << std::setw(2) << std::setfill('0') << static_cast<int>(hash[i]);
    }
    
    return oss.str();
}

/**
 * @brief Safe wrapper for sha256_file that returns error message instead of throwing
 * 
 * @param filepath Path to file
 * @param default_value Value to return if hashing fails (default: "ERROR")
 * @return SHA256 hash or default_value on error
 */
inline std::string sha256_file_safe(const std::string& filepath, const std::string& default_value = "ERROR") {
    try {
        return sha256_file(filepath);
    } catch (const std::exception& e) {
        return default_value;
    }
}

} // namespace ICARION::utils
