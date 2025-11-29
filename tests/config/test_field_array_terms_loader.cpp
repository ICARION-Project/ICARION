/**
 * @file test_field_array_terms_loader.cpp
 * @brief Unit tests for field_array_terms JSON loading
 * 
 * Tests DomainConfigLoader's ability to parse field_array_terms from JSON.
 */

#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>

#include "core/config/loader/DomainConfigLoader.h"
#include <json/json.h>

using Catch::Approx;
using namespace ICARION::config;

// Helper to create minimal domain config with field_array_terms
Json::Value create_domain_with_field_terms(const Json::Value& field_array_terms) {
    Json::Value json;
    json["name"] = "test_domain";
    json["instrument"] = "IMS";
    
    // Minimal geometry
    Json::Value geometry;
    geometry["origin_m"] = Json::arrayValue;
    geometry["origin_m"].append(0.0);
    geometry["origin_m"].append(0.0);
    geometry["origin_m"].append(0.0);
    geometry["length_m"] = 0.1;
    geometry["radius_m"] = 0.05;
    json["geometry"] = geometry;
    
    // Minimal environment
    Json::Value env;
    env["temperature_K"] = 300.0;
    env["pressure_Pa"] = 101325.0;
    json["env"] = env;
    
    // Fields with field_array_terms
    Json::Value fields;
    fields["field_array_terms"] = field_array_terms;
    json["fields"] = fields;
    
    return json;
}

using namespace ICARION::config;

TEST_CASE("FieldArrayTerms: Load single constant term", "[config][field_array]") {
    Json::Value terms = Json::arrayValue;
    
    Json::Value term;
    term["file"] = "test_field.h5";
    term["scale_type"] = "Constant";
    term["constant_V"] = 100.0;
    terms.append(term);
    
    Json::Value domain_json = create_domain_with_field_terms(terms);
    DomainConfig domain = DomainConfigLoader::load(domain_json);
    
    REQUIRE(domain.fields.field_array_terms.size() == 1);
    REQUIRE(domain.fields.field_array_terms[0].file == "test_field.h5");
    REQUIRE(domain.fields.field_array_terms[0].kind == FieldsConfig::FieldArrayTerm::ScaleKind::Constant);
    REQUIRE(domain.fields.field_array_terms[0].constant == 100.0);
}

TEST_CASE("FieldArrayTerms: Load DC_Axial term", "[config][field_array]") {
    Json::Value terms = Json::arrayValue;
    
    Json::Value term;
    term["file"] = "bem_field.h5";
    term["scale_type"] = "DC_Axial";
    terms.append(term);
    
    Json::Value domain_json = create_domain_with_field_terms(terms);
    DomainConfig domain = DomainConfigLoader::load(domain_json);
    
    REQUIRE(domain.fields.field_array_terms.size() == 1);
    REQUIRE(domain.fields.field_array_terms[0].file == "bem_field.h5");
    REQUIRE(domain.fields.field_array_terms[0].kind == FieldsConfig::FieldArrayTerm::ScaleKind::DC_Axial);
}

TEST_CASE("FieldArrayTerms: Load RF term with phase", "[config][field_array]") {
    Json::Value terms = Json::arrayValue;
    
    Json::Value term;
    term["file"] = "rf_field.h5";
    term["scale_type"] = "RF";
    term["phase_rad"] = 1.57;  // π/2
    term["frequency_Hz"] = 1e6;  // 1 MHz
    terms.append(term);
    
    Json::Value domain_json = create_domain_with_field_terms(terms);
    DomainConfig domain = DomainConfigLoader::load(domain_json);
    
    REQUIRE(domain.fields.field_array_terms.size() == 1);
    REQUIRE(domain.fields.field_array_terms[0].file == "rf_field.h5");
    REQUIRE(domain.fields.field_array_terms[0].kind == FieldsConfig::FieldArrayTerm::ScaleKind::RF);
    REQUIRE(domain.fields.field_array_terms[0].phase_rad == Approx(1.57));
    REQUIRE(domain.fields.field_array_terms[0].frequency_Hz == Approx(1e6));
}

TEST_CASE("FieldArrayTerms: Load multiple terms", "[config][field_array]") {
    Json::Value terms = Json::arrayValue;
    
    // Term 1: DC Axial
    Json::Value term1;
    term1["file"] = "dc_axial.h5";
    term1["scale_type"] = "DC_Axial";
    terms.append(term1);
    
    // Term 2: DC Quad
    Json::Value term2;
    term2["file"] = "dc_quad.h5";
    term2["scale_type"] = "DC_Quad";
    terms.append(term2);
    
    // Term 3: RF
    Json::Value term3;
    term3["file"] = "rf_trap.h5";
    term3["scale_type"] = "RF";
    term3["phase_rad"] = 0.0;
    terms.append(term3);
    
    Json::Value domain_json = create_domain_with_field_terms(terms);
    DomainConfig domain = DomainConfigLoader::load(domain_json);
    
    REQUIRE(domain.fields.field_array_terms.size() == 3);
    
    REQUIRE(domain.fields.field_array_terms[0].kind == FieldsConfig::FieldArrayTerm::ScaleKind::DC_Axial);
    REQUIRE(domain.fields.field_array_terms[1].kind == FieldsConfig::FieldArrayTerm::ScaleKind::DC_Quad);
    REQUIRE(domain.fields.field_array_terms[2].kind == FieldsConfig::FieldArrayTerm::ScaleKind::RF);
}

TEST_CASE("FieldArrayTerms: Default values", "[config][field_array]") {
    Json::Value terms = Json::arrayValue;
    
    Json::Value term;
    term["file"] = "field.h5";
    // No scale_type → defaults to Constant
    // No constant_V → defaults to 1.0
    // No phase_rad → defaults to 0.0
    // No frequency_Hz → defaults to 0.0
    terms.append(term);
    
    Json::Value domain_json = create_domain_with_field_terms(terms);
    DomainConfig domain = DomainConfigLoader::load(domain_json);
    
    REQUIRE(domain.fields.field_array_terms.size() == 1);
    REQUIRE(domain.fields.field_array_terms[0].kind == FieldsConfig::FieldArrayTerm::ScaleKind::Constant);
    REQUIRE(domain.fields.field_array_terms[0].constant == Approx(1.0));
    REQUIRE(domain.fields.field_array_terms[0].phase_rad == Approx(0.0));
    REQUIRE(domain.fields.field_array_terms[0].frequency_Hz == Approx(0.0));
}

TEST_CASE("FieldArrayTerms: All scale types", "[config][field_array]") {
    std::vector<std::string> scale_types = {
        "Constant", "DC_Axial", "DC_Quad", "DC_Radial", "RF"
    };
    
    std::vector<FieldsConfig::FieldArrayTerm::ScaleKind> expected_kinds = {
        FieldsConfig::FieldArrayTerm::ScaleKind::Constant,
        FieldsConfig::FieldArrayTerm::ScaleKind::DC_Axial,
        FieldsConfig::FieldArrayTerm::ScaleKind::DC_Quad,
        FieldsConfig::FieldArrayTerm::ScaleKind::DC_Radial,
        FieldsConfig::FieldArrayTerm::ScaleKind::RF
    };
    
    for (size_t i = 0; i < scale_types.size(); ++i) {
        Json::Value terms = Json::arrayValue;
        
        Json::Value term;
        term["file"] = "test.h5";
        term["scale_type"] = scale_types[i];
        terms.append(term);
        
        Json::Value domain_json = create_domain_with_field_terms(terms);
        DomainConfig domain = DomainConfigLoader::load(domain_json);
        
        REQUIRE(domain.fields.field_array_terms[0].kind == expected_kinds[i]);
    }
}

TEST_CASE("FieldArrayTerms: Error - missing file", "[config][field_array]") {
    Json::Value terms = Json::arrayValue;
    
    Json::Value term;
    // No "file" field!
    term["scale_type"] = "Constant";
    terms.append(term);
    
    Json::Value domain_json = create_domain_with_field_terms(terms);
    REQUIRE_THROWS_AS(DomainConfigLoader::load(domain_json), std::runtime_error);
}

TEST_CASE("FieldArrayTerms: Error - unknown scale_type", "[config][field_array]") {
    Json::Value terms = Json::arrayValue;
    
    Json::Value term;
    term["file"] = "test.h5";
    term["scale_type"] = "InvalidType";
    terms.append(term);
    
    Json::Value domain_json = create_domain_with_field_terms(terms);
    REQUIRE_THROWS_AS(DomainConfigLoader::load(domain_json), std::runtime_error);
}

TEST_CASE("FieldArrayTerms: Empty array is valid", "[config][field_array]") {
    Json::Value terms = Json::arrayValue;  // Empty array
    
    Json::Value domain_json = create_domain_with_field_terms(terms);
    DomainConfig domain = DomainConfigLoader::load(domain_json);
    
    REQUIRE(domain.fields.field_array_terms.empty());
}

TEST_CASE("FieldArrayTerms: No field_array_terms key is valid", "[config][field_array]") {
    Json::Value json;
    json["name"] = "test_domain";
    json["instrument"] = "IMS";
    
    Json::Value geometry;
    geometry["origin_m"] = Json::arrayValue;
    geometry["origin_m"].append(0.0);
    geometry["origin_m"].append(0.0);
    geometry["origin_m"].append(0.0);
    geometry["length_m"] = 0.1;
    geometry["radius_m"] = 0.05;
    json["geometry"] = geometry;
    
    Json::Value env;
    env["temperature_K"] = 300.0;
    env["pressure_Pa"] = 101325.0;
    json["env"] = env;
    // No "fields" or "field_array_terms" key at all
    
    DomainConfig domain = DomainConfigLoader::load(json);
    
    REQUIRE(domain.fields.field_array_terms.empty());
}
