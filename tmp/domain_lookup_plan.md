Domain Lookup Performance Plan
==============================

Goal
- Replace O(N_ions × N_domains) per-step scans with indexed lookup while preserving boundary checks/actions.

Current Pain
- `SimulationEngine` calls DomainManager to find domains via linear scans; repeated twice per step.
- Multi-domain IMS/drift setups blow up OpenMP runtime; cache-unfriendly.

Plan

1) Choose Index Strategy
- Use a simple uniform grid or axis-aligned interval bins along z (common for drift/TOF/IMS) plus radius check for cylinders.
- Keep a fallback linear scan for complex geometries (Orbitrap).
- Index built once from DomainConfig (geometry origin/length/radius).

2) DomainManager Extension
- Add a spatial index structure (e.g., per-axis bins -> vector<vector<domain_id>>).
- API: `find_domain_index(pos)` uses index first; fallback to brute force if ambiguous/complex geometry.
- Maintain existing boundary actions; after index hit, still perform precise geometry check.

3) SimulationEngine Integration
- Replace O(domains) loops with indexed lookup helper (DomainManager) in the hot path.
- Keep behavior identical for edge cases (domain transitions, boundary actions).

4) Validation
- Unit tests: synthetic multi-domain configs (stacked cylinders) with random points; index vs brute-force results identical.
- Performance micro-benchmark: compare per-step domain lookup time before/after for 1e5 ions, 10 domains.

5) Rollout
- Default to indexed lookup; allow config flag to disable if needed (for debugging).
- Document limitation: Orbitrap/complex geometries fall back to linear scan.

Out of Scope
- Full 3D spatial acceleration (BVH/KD-tree) for arbitrary geometries.
- GPU-side domain lookup.
