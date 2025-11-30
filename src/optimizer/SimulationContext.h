// ICARION: Ion Collision And Reaction IntegratiON
// MIT License - Copyright (c) 2025 ICARION Project Contributors

/*
 * Canonical SimulationContext header.
 *
 * This file contains the single, lightweight SimulationContext used by the
 * integrator core and new code paths. The previous split between a
 * namespaced "light" definition and a separate legacy global definition has
 * been removed: the legacy rich representation was deleted and the
 * lightweight type is now the single authoritative definition.
 *
 * To minimise breakage for existing code, this header also provides a
 * convenience alias in the global namespace:
 *   using ::SimulationContext = optimization::SimulationContext;
 * This alias does NOT recreate the old legacy rich struct; it simply makes
 * the new lightweight type available under the former global name. Please
 * prefer `optimization::SimulationContext` in new code.
 */

#pragma once

#include <vector>
#include <string>
#include <iostream>
#include <unordered_map>
#include <memory>

#include "fieldsolver/io/SolverConfig.h"
// Field provider abstraction (optional injected field source)
#include "fieldsolver/utils/IFieldProvider.h"

namespace optimization {

struct SimulationContext {
	GlobalParams gParams;
	std::vector<InstrumentDomain> domains;
	std::unordered_map<std::string, Species> speciesDB;
	// std::vector<ReactionEntry> reaction_list;  // DISABLED: ReactionEntry is legacy
	std::vector<IonState> ions;

	// Backward-compatible auxiliary fields (migrated from the old legacy
	// SimulationContext). These make it easier to populate the context from
	// adapters without needing a separate legacy type.
	std::string run_id;
	unsigned int rng_seed = 0;

	struct InjectionDescriptor {
		std::string type = "ring_seed";
		int count = 100;
		std::string species; // species name
		double energy_eV = 1.0;
		double radius_m = 0.0; // for ring_seed
		double z_m = 0.0;      // injection z coordinate
		double radial_jitter_m = 0.0; // positional jitter
		double angular_spread_rad = 2.0 * 3.14159265358979323846; // full circle by default
	} injectionOwned;
	InjectionDescriptor* injection = nullptr; // optional pointer to descriptor

	// Optional solver configuration populated by adapter (meshes, grid settings)
	std::unique_ptr<SolverConfig> solverConfigOwned;
	SolverConfig* solverConfig = nullptr;

	// Optional injected field provider. When set, integrator will prefer
	// this provider for evaluating E/phi instead of legacy domain-based logic.
	std::unique_ptr<IFieldProvider> fieldProviderOwned;
	const IFieldProvider* fieldProvider = nullptr;

	SimulationContext() = default;

	// Constructor without ReactionEntry (legacy type removed)
	SimulationContext(const GlobalParams& g,
					  const std::vector<InstrumentDomain>& d,
					  const std::unordered_map<std::string, Species>& s,
					  const std::vector<IonState>& i)
		: gParams(g), domains(d), speciesDB(s), ions(i) {}

	// Basic validation to catch common misconfigurations early.
	bool validate(std::string* out_msg = nullptr) const {
		if (gParams.t_eval.empty()) {
			if (out_msg) *out_msg = "GlobalParams.t_eval is empty";
			return false;
		}
		if (domains.empty()) {
			if (out_msg) *out_msg = "No instrument domains defined";
			return false;
		}
		if (ions.empty()) {
			if (out_msg) *out_msg = "Ion ensemble is empty";
			return false;
		}
		return true;
	}

	void print_summary(std::ostream& os = std::cout) const {
		os << "SimulationContext summary:\n";
		os << "  num_domains: " << domains.size() << "\n";
		os << "  num_species: " << speciesDB.size() << "\n";
		os << "  num_reactions: " << reaction_list.size() << "\n";
		os << "  num_ions: " << ions.size() << "\n";
		if (!gParams.t_eval.empty()) {
			os << "  t_start: " << gParams.t_eval.front() << " s\n";
			os << "  t_end:   " << gParams.t_eval.back()  << " s\n";
		}
	}
};

} // namespace optimization

// Note: the namespaced lightweight type is available as
// `optimization::SimulationContext`. For minimal breakage we also provide
// an unqualified alias in the global namespace so existing code that uses
// `SimulationContext` (without the namespace) keeps compiling.

using SimulationContext = optimization::SimulationContext;

