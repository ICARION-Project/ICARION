// ICARION: Ion Collision And Reaction IntegratiON
// SPDX-License-Identifier: GPL-3.0-only
// Copyright (c) 2025 Christoph Schaefer
#ifndef ICARION_CONFIG_VALIDATION_RESULT_H
#define ICARION_CONFIG_VALIDATION_RESULT_H

#include <string>
#include <vector>
#include <iostream>

namespace ICARION::config {

/**
 * @brief Result of configuration validation
 * 
 * Collects errors (blocking) and warnings (non-blocking) during validation.
 */
struct ValidationResult {
    bool valid = true;
    std::vector<std::string> errors;
    std::vector<std::string> warnings;
    
    /**
     * @brief Add a blocking error
     */
    void add_error(const std::string& msg) {
        errors.push_back(msg);
        valid = false;
    }
    
    /**
     * @brief Add a non-blocking warning
     */
    void add_warning(const std::string& msg) {
        warnings.push_back(msg);
    }
    
    /**
     * @brief Merge another validation result
     */
    void merge(const ValidationResult& other) {
        errors.insert(errors.end(), other.errors.begin(), other.errors.end());
        warnings.insert(warnings.end(), other.warnings.begin(), other.warnings.end());
        valid = valid && other.valid;
    }
    
    /**
     * @brief Print all errors and warnings
     */
    void print() const {
        if (!errors.empty()) {
            std::cerr << "=== VALIDATION ERRORS ===" << std::endl;
            for (const auto& err : errors) {
                std::cerr << "  [ERROR] " << err << std::endl;
            }
        }
        if (!warnings.empty()) {
            std::cout << "=== VALIDATION WARNINGS ===" << std::endl;
            for (const auto& warn : warnings) {
                std::cout << "  [WARN]  " << warn << std::endl;
            }
        }
    }
};

} // namespace ICARION::config

#endif // ICARION_CONFIG_VALIDATION_RESULT_H
