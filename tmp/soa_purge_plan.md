## Goal
Make the integrator stack SoA-only (no AoS APIs or buffers). Remove mixed paths, wrappers, and duplicated logic.

## Scope
- Public APIs: SimulationEngine, OutputManager, PhysicsSetup/main entry.
- Integration strategies / GPU fallbacks: SoA only.
- Tests/benchmarks: call SoA API (convert once at setup if needed).
- Docs: reflect SoA-only.

## Steps
1) Identify AoS entry points
   - SimulationEngine: run(std::vector<IonState>&), initialize(std::vector<IonState>&), try_gpu_* AoS signatures, AoS allocations in GPU helpers.
   - OutputManager: initialize/log_step/finalize AoS overloads, AoS buffers.
   - PhysicsSetup/main: generate ions as AoS and pass through; space-charge selection uses AoS.
   - Tests/benchmarks that call AoS run (e.g., benchmark_soa_performance, main integration tests).
2) API cleanup
   - Remove AoS overloads from SimulationEngine/OutputManager headers and impl.
   - keep a single SoA API: run(core::IonEnsemble&), initialize_soa, log_step_soa, finalize_soa.
   - Provide explicit AoS→SoA conversion helpers at call sites (main/tests), not inside core.
3) Call-site updates
   - main.cpp: generate ions AoS -> convert to IonEnsemble once; PhysicsSetup takes SoA (adjust signature) and returns modules; OutputManager initialized via SoA path.
   - PhysicsSetup: accept IonEnsemble (for space charge sizing) or pass ion count; remove AoS includes.
   - Tests/benchmarks: where AoS is created, wrap with IonEnsemble::from_legacy before run; remove any AoS-engine usage.
4) Strategies/GPU
   - Remove AoS signatures from integration strategies (RK4/RK45/Boris) if present; ensure GPU helper calls are SoA-only.
   - Space charge GPU/direct already SoA-aware; drop any AoS fallbacks if still present.
5) Output buffers
   - OutputManager: delete AoS trajectory_buffer_, keep only SoA buffer; simplify flush().
6) Docs
   - Update DEVELOPERS_GUIDE / ARCHITECTURE to state SoA-only pipeline and removal of AoS APIs.
7) Build & tests
   - Full build; fix failing tests due to API changes.

## Risks / notes
- Main and PhysicsSetup need coordinated signature changes.
- Benchmarks that compare AoS vs SoA may need removal or rewording (only SoA path remains; can keep AoS input by converting once).
- Be careful with GPU tests expecting AoS inputs; convert at setup.
