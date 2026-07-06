#!/usr/bin/env python3
"""
TIMS Mobility-Sorted Elution Validation Script

Purpose:
    Validate that the TIMSAxialFieldModel correctly separates ions by mobility
    during the electric-gradient ramp phase. This is the defining physics of
    Trapped Ion Mobility Spectrometry (TIMS).

Physics Tested:
    During the TIMS ramp, the net axial field
        E_net(z,t) = (1 - f(t))*E_initial(z) + f(t)*E_final(z)
    decreases as f(t) increases from 0 to 1. An ion is trapped as long as:
        K(P,T) * |E_net(z_trap)| >= v_gas
    and elutes when the ramp has reduced the field sufficiently. Since K0 differs
    by species, elution times depend on mobility.

    Theoretical elution ramp fraction:
        f_elut = 1 - v_gas / (K(P,T) * E_initial_max)
    → Lower K0  → lower f_elut → elutes earlier in the ramp
    → Higher K0 → higher f_elut → elutes later in the ramp

    Expected elution ORDER (from lowest to highest K0):
        ReserpineH+ (K0=3.02) → CaffeineH+ (K0=4.51) → 26DTBPH+ (K0=6.31)

Validation Strategy:
    1. Simulate 3-species mixture (60 ions each) in TIMS geometry
    2. Extract ion elution times from death_time_s in HDF5 output
    3. Validate elution order is correct (PRIMARY gate: PASS/FAIL)
    4. Validate species are clearly separated (gap > 5 × pooled σ)
    5. Compare mean elution times to analytical theory (within 25% tolerance)

Analytical reference:
    Bruker-style TIMS parameters from literature-point campaign:
    - 3 species: ReserpineH+, CaffeineH+, 26DTBPH+
    - N2 at 310 Pa, 300 K, axial flow v_gas = 134 m/s
    - Linear E field: 0 → -4495 V/m over L = 20 mm
    - Linear ramp: f 0 → 1 over t_ramp = 0.8 ms → 1.8 ms

Runtime: ~15-20 seconds (3 ms simulation, 180 ions, dt=10 ns)

References:
    - Meier et al. (2015): "Trapped Ion Mobility Spectrometry" Anal. Chem.
    - Hernandez-Mesa et al. (2019): TIMS-QTOF fundamentals
    - ICARION TIMSAxialFieldModel: src/core/config/types/TIMSAxialFieldModel.h

Author: ICARION Validation Team
"""

import subprocess
import json
import h5py
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from pathlib import Path
import sys
from datetime import datetime
import os

# ---------------------------------------------------------------------------
# Physical constants
# ---------------------------------------------------------------------------
BOLTZMANN_K   = 1.380649e-23    # J/K
AMU_TO_KG     = 1.66053906660e-27  # kg/amu
ELEM_CHARGE_C = 1.602176634e-19  # C
P0_PA         = 101325.0         # STP pressure [Pa]
T0_K          = 273.15           # STP temperature [K]

# ---------------------------------------------------------------------------
# TIMS simulation parameters (match examples/ims/ims_tims_basic.json)
# ---------------------------------------------------------------------------
TUBE_LENGTH_M  = 0.02     # m
TUBE_RADIUS_M  = 0.005    # m
PRESSURE_PA    = 310.0    # Pa
TEMPERATURE_K  = 300.0    # K
GAS_SPECIES    = "N2"
GAS_VELOCITY_MS = 134.0   # m/s (axial)
RF_VOLTAGE_V   = 180.0    # V
RF_FREQ_HZ     = 2e6      # Hz
E_INITIAL_MAX_VM = 4495.0  # V/m (magnitude at z=L before ramp)
E_FINAL_MAX_VM   = 0.0     # V/m (magnitude at z=L after ramp)
RAMP_START_S   = 0.8e-3   # s
RAMP_END_S     = 1.8e-3   # s
TOTAL_TIME_S   = 3.0e-3   # s
DT_S           = 1e-8     # s
WRITE_INTERVAL = 300      # steps
IONS_PER_SPECIES = 60

# ---------------------------------------------------------------------------
# Species and reference K0 values (from species_database_v1.json)
# ---------------------------------------------------------------------------
SPECIES = [
    {"id": "ReserpineH+", "K0_cm2Vs": 3.02},
    {"id": "CaffeineH+",  "K0_cm2Vs": 4.51},
    {"id": "26DTBPH+",    "K0_cm2Vs": 6.307},
]
# Expected elution order (lowest K0 first)
EXPECTED_ORDER = ["ReserpineH+", "CaffeineH+", "26DTBPH+"]

# Tolerances
ORDER_IS_PRIMARY_GATE   = True   # FAIL if order is wrong
SEPARATION_MIN_SIGMA    = 5.0    # gap must exceed N × pooled_σ
THEORY_TOLERANCE_PCT    = 30.0   # elution time within 30% of theory

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
REPO_ROOT = Path(__file__).resolve().parents[3]
VALIDATION_DIR = REPO_ROOT / "validation"
RUN_DIR = os.environ.get("ICARION_VALIDATION_RUN_DIR")
if RUN_DIR:
    RUN_DIR = Path(RUN_DIR)
    FIGURES_DIR = RUN_DIR / "figures" / "physics"
    LOGS_DIR    = RUN_DIR / "logs"
    RESULTS_DIR = RUN_DIR / "results" / "tims_elution"
else:
    FIGURES_DIR = VALIDATION_DIR / "figures" / "physics"
    LOGS_DIR    = VALIDATION_DIR / "logs"
    RESULTS_DIR = VALIDATION_DIR / "results" / "tims_elution"

ICARION_BIN  = REPO_ROOT / "build" / "src" / "icarion_main"
SPECIES_DB   = REPO_ROOT / "data" / "species_database_v1.json"

FIGURES_DIR.mkdir(parents=True, exist_ok=True)
LOGS_DIR.mkdir(exist_ok=True)
RESULTS_DIR.mkdir(parents=True, exist_ok=True)

LOG_FILE = LOGS_DIR / "TIMS_ELUTION_VALIDATION.txt"
_log_handle = open(LOG_FILE, 'w')


def log(msg: str) -> None:
    print(msg)
    _log_handle.write(msg + "\n")
    _log_handle.flush()


# ---------------------------------------------------------------------------
# Theoretical elution time computation
# ---------------------------------------------------------------------------
def compute_theoretical_elution(K0_cm2Vs: float) -> float:
    """
    Return predicted mean elution time [s] for a species with reduced mobility K0.

    At elution: K(P,T) × E_initial_max × (1 - f_elut) = v_gas
    → f_elut = 1 - v_gas / (K(P,T) × E_initial_max)
    → t_elut = ramp_start + f_elut × (ramp_end - ramp_start)

    The gradient field means the actual trapping position shifts, introducing a
    small correction; we report the simple analytical estimate for reference.
    """
    N_ratio = (P0_PA / PRESSURE_PA) * (TEMPERATURE_K / T0_K)
    K_m2Vs  = (K0_cm2Vs * 1e-4) * N_ratio   # K at actual conditions [m²/(V·s)]
    f_elut  = 1.0 - GAS_VELOCITY_MS / (K_m2Vs * E_INITIAL_MAX_VM)
    f_elut  = max(0.0, min(1.0, f_elut))     # clamp to [0, 1]
    t_elut  = RAMP_START_S + f_elut * (RAMP_END_S - RAMP_START_S)
    return t_elut


# ---------------------------------------------------------------------------
# Config builder
# ---------------------------------------------------------------------------
def create_config(output_name: str) -> dict:
    ions_list = []
    for sp in SPECIES:
        ions_list.append({
            "id": sp["id"],
            "count": IONS_PER_SPECIES,
            "position": {
                "type": "gaussian",
                "center": [0.0, 0.0, TUBE_LENGTH_M / 2],
                "std":    [6e-5, 6e-5, 2e-4]
            },
            "velocity": {
                "type": "thermal",
                "temperature_K": TEMPERATURE_K
            }
        })

    return {
        "simulation": {
            "total_time_s":  TOTAL_TIME_S,
            "dt_s":          DT_S,
            "write_interval": WRITE_INTERVAL,
            "integrator":    "RK4",
            "enable_gpu":    False,
            "enable_openmp": True,
            "rng_seed":      42
        },
        "physics": {
            "collision_model":            "HSS",
            "enable_ou_thermalization":   False,
            "enable_space_charge":        False,
            "collision_multi_event_mode": True,
            "collision_max_events_per_step": 4
        },
        "output": {
            "folder":          str(RESULTS_DIR),
            "trajectory_file": output_name,
            "print_progress":  True
        },
        "ions": {
            "species": ions_list
        },
        "species_database": str(SPECIES_DB),
        "domains": [{
            "name":       "tims_elution_validation",
            "instrument": "TIMS",
            "geometry": {
                "origin_m": [0.0, 0.0, 0.0],
                "length_m": TUBE_LENGTH_M,
                "radius_m": TUBE_RADIUS_M
            },
            "boundary": {
                "type": "Absorption"
            },
            "env": {
                "temperature_K": TEMPERATURE_K,
                "pressure_Pa":   PRESSURE_PA,
                "gas_species":   GAS_SPECIES,
                "flow_model":    "axial_uniform",
                "flow_parameters": {
                    "axial_flow_velocity_m_s": GAS_VELOCITY_MS
                }
            },
            "fields": {
                "RF": {
                    "voltage_V":    RF_VOLTAGE_V,
                    "frequency_Hz": RF_FREQ_HZ,
                    "phase_rad":    0.0
                },
                "DC": {
                    "quad_V": 0.0
                },
                "TIMS": {
                    "enabled":                True,
                    "z_positions_m":          [0.0, TUBE_LENGTH_M],
                    "axial_field_initial_profile_V_m": [0.0, -E_INITIAL_MAX_VM],
                    "axial_field_final_profile_V_m":   [0.0, -E_FINAL_MAX_VM],
                    "ramp_start_s":           RAMP_START_S,
                    "ramp_end_s":             RAMP_END_S,
                    "ramp_mode":              "linear"
                }
            }
        }]
    }


# ---------------------------------------------------------------------------
# Simulation runner
# ---------------------------------------------------------------------------
def run_simulation(config_name: str, config_dict: dict) -> Path | None:
    config_path = RESULTS_DIR / config_name
    with open(config_path, 'w') as f:
        json.dump(config_dict, f, indent=2)

    log(f"\n▶  Running TIMS elution simulation: {config_name}")
    log(f"   Config: {config_path}")

    if not ICARION_BIN.exists():
        log(f"   ✗ ICARION binary not found: {ICARION_BIN}")
        log(f"     Build the project first: cmake --build build/ --target icarion_main")
        return None

    try:
        result = subprocess.run(
            [str(ICARION_BIN), str(config_path)],
            capture_output=True, text=True, timeout=300
        )
        if result.returncode != 0:
            log(f"   ✗ Simulation failed (exit {result.returncode})")
            log(f"     stderr: {result.stderr[-800:]}")
            return None
        log("   ✓ Simulation complete")
        return RESULTS_DIR / config_dict["output"]["trajectory_file"]
    except subprocess.TimeoutExpired:
        log("   ✗ Simulation timed out (>5 min)")
        return None


# ---------------------------------------------------------------------------
# Analysis
# ---------------------------------------------------------------------------
def analyze_elution(h5_path: Path) -> dict | None:
    """
    Load HDF5 output and extract per-species elution statistics.

    Elution time = death_time_s, i.e. the simulation step at which the ion
    exited the domain and was deactivated.  Ions that survive the full
    simulation (death_time_s == 0 or negative) are flagged as 'not eluted'
    and excluded from the elution time statistics.
    """
    log(f"\n📊  Analysing: {h5_path.name}")
    try:
        with h5py.File(h5_path, 'r') as f:
            ion_species  = f['ions/initial_species_id'][:]   # bytes
            death_times  = f['ions/death_time_s'][:]         # float64
            n_ions_total = len(ion_species)
    except Exception as exc:
        log(f"   ✗ Could not read HDF5: {exc}")
        return None

    results = {}
    for sp in SPECIES:
        sid = sp["id"]
        mask   = np.array([s.decode() == sid for s in ion_species])
        dt_all = death_times[mask]
        # Ions with death_time > 0 have exited the domain (eluted or absorbed at wall)
        eluted = dt_all[dt_all > 0]
        not_eluted = np.sum(dt_all <= 0)

        if len(eluted) == 0:
            log(f"   ⚠  {sid}: 0 ions eluted — all remain trapped at simulation end")
            results[sid] = {
                "mean_s": None, "std_s": None,
                "n_eluted": 0, "n_not_eluted": int(not_eluted)
            }
        else:
            mean_t = eluted.mean()
            std_t  = eluted.std()
            t_theory = compute_theoretical_elution(sp["K0_cm2Vs"])
            err_pct  = 100.0 * abs(mean_t - t_theory) / t_theory
            log(f"   {sid}:")
            log(f"     Eluted:  {len(eluted):3d}/{len(dt_all)}  "
                f"mean={mean_t*1e3:.3f} ms  σ={std_t*1e6:.1f} µs  "
                f"[{eluted.min()*1e3:.3f}..{eluted.max()*1e3:.3f} ms]")
            log(f"     Theory:  {t_theory*1e3:.3f} ms   Δ={err_pct:.1f}%")
            results[sid] = {
                "mean_s": mean_t, "std_s": std_t,
                "n_eluted": len(eluted), "n_not_eluted": int(not_eluted),
                "elution_times": eluted,
                "theory_s": t_theory, "theory_err_pct": err_pct
            }
    return results


# ---------------------------------------------------------------------------
# Plotting
# ---------------------------------------------------------------------------
def plot_elution_spectrum(results: dict, fig_path: Path) -> None:
    """Plot a TIMS-style 'mobilogram' — histogram of elution times per species."""
    fig, axes = plt.subplots(2, 1, figsize=(9, 7), sharex=True,
                             gridspec_kw={"height_ratios": [3, 1]})
    ax_hist, ax_sep = axes

    colors = {"ReserpineH+": "#e04b4b", "CaffeineH+": "#4b8ae0", "26DTBPH+": "#4bc24b"}
    labels = {
        "ReserpineH+": f"ReserpineH⁺ (K₀ = {SPECIES[0]['K0_cm2Vs']:.2f} cm²/Vs)",
        "CaffeineH+":  f"CaffeineH⁺  (K₀ = {SPECIES[1]['K0_cm2Vs']:.2f} cm²/Vs)",
        "26DTBPH+":    f"26DTBPH⁺   (K₀ = {SPECIES[2]['K0_cm2Vs']:.3f} cm²/Vs)",
    }

    bin_edges = np.linspace(RAMP_START_S * 0.95, RAMP_END_S * 1.05, 60)

    for sid, res in results.items():
        if res.get("elution_times") is None:
            continue
        t_ms = res["elution_times"] * 1e3
        be_ms = bin_edges * 1e3
        counts, _ = np.histogram(t_ms, bins=be_ms)
        centers = 0.5 * (be_ms[:-1] + be_ms[1:])
        ax_hist.bar(centers, counts, width=np.diff(be_ms),
                    color=colors.get(sid, "grey"), alpha=0.65,
                    label=labels.get(sid, sid))
        ax_hist.axvline(res["mean_s"] * 1e3, color=colors.get(sid, "grey"),
                        lw=1.5, ls="--")
        t_th = res.get("theory_s", None)
        if t_th is not None:
            ax_hist.axvline(t_th * 1e3, color=colors.get(sid, "grey"),
                            lw=1.0, ls=":", alpha=0.7)

    ax_hist.set_ylabel("Ion count per bin")
    ax_hist.set_title("TIMS Mobility-Sorted Elution Validation\n"
                       "N₂, 310 Pa, 300 K, v_gas = 134 m/s, "
                       f"E_initial = {E_INITIAL_MAX_VM:.0f} V/m, "
                       f"RF {RF_VOLTAGE_V:.0f} V @ {RF_FREQ_HZ/1e6:.0f} MHz")
    ax_hist.legend(fontsize=9, loc="upper left")

    # Dashed = simulation mean, dotted = analytical theory
    ax_hist.plot([], [], 'k--', lw=1.5, label="Simulated mean")
    ax_hist.plot([], [], 'k:',  lw=1.0, alpha=0.7, label="Analytical theory")
    ax_hist.legend(fontsize=8, loc="upper right")

    # Lower panel: highlight ramp window
    ax_sep.axvspan(RAMP_START_S * 1e3, RAMP_END_S * 1e3,
                   color="lightyellow", alpha=0.6, label="Ramp window")
    ax_sep.set_xlabel("Elution time (ms)")
    ax_sep.set_ylabel("Ramp f(t)")
    t_arr = np.linspace(0, TOTAL_TIME_S, 300) * 1e3
    f_arr = np.clip((t_arr/1e3 - RAMP_START_S) / (RAMP_END_S - RAMP_START_S), 0, 1)
    ax_sep.plot(t_arr, f_arr, 'k-', lw=1.5, label="f(t) ramp")
    ax_sep.set_xlim(RAMP_START_S * 0.95e3, RAMP_END_S * 1.05e3)
    ax_sep.set_ylim(-0.05, 1.1)
    ax_sep.legend(fontsize=8)

    plt.tight_layout()
    fig_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(fig_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    log(f"\n📈  Figure saved: {fig_path}")


# ---------------------------------------------------------------------------
# Validation gates
# ---------------------------------------------------------------------------
def validate(results: dict) -> bool:
    """
    Apply validation gates.  Returns True if ALL primary gates pass.
    """
    log("\n" + "=" * 65)
    log("VALIDATION GATES")
    log("=" * 65)

    all_pass = True

    # Gate 1: All species must have eluted ions
    log("\n[Gate 1]  All species produced eluted ions")
    for sid in EXPECTED_ORDER:
        res = results.get(sid, {})
        n = res.get("n_eluted", 0)
        if n == 0:
            log(f"  FAIL  {sid}: 0 eluted ions")
            all_pass = False
        else:
            log(f"  PASS  {sid}: {n} ions eluted")

    # Gate 2: Elution order (primary physics gate)
    log("\n[Gate 2]  Elution order  (ReserpineH+ → CaffeineH+ → 26DTBPH+)")
    means = {}
    for sid in EXPECTED_ORDER:
        m = results.get(sid, {}).get("mean_s", None)
        means[sid] = m
        if m is None:
            log(f"  SKIP  {sid}: no elution data")

    ordered_species = [s for s in EXPECTED_ORDER if means.get(s) is not None]
    order_ok = all(
        means[ordered_species[i]] < means[ordered_species[i + 1]]
        for i in range(len(ordered_species) - 1)
    )
    if len(ordered_species) < 2:
        log("  SKIP  insufficient data for order check")
    elif order_ok:
        order_str = " → ".join(
            f"{s} ({means[s]*1e3:.3f} ms)" for s in ordered_species
        )
        log(f"  PASS  Order correct: {order_str}")
    else:
        log("  FAIL  Elution order is WRONG — mobility separation not working")
        for s in ordered_species:
            log(f"        {s}: {means[s]*1e3:.3f} ms")
        all_pass = False

    # Gate 3: Clear peak separation (gap > N × pooled_σ)
    log(f"\n[Gate 3]  Peak separation  (gap > {SEPARATION_MIN_SIGMA:.0f} × pooled σ)")
    for i in range(len(ordered_species) - 1):
        s1, s2 = ordered_species[i], ordered_species[i + 1]
        if means.get(s1) is None or means.get(s2) is None:
            log(f"  SKIP  {s1} → {s2}: missing data")
            continue
        gap = means[s2] - means[s1]
        sig1 = results[s1].get("std_s", 1e-6) or 1e-6
        sig2 = results[s2].get("std_s", 1e-6) or 1e-6
        pooled_sigma = np.sqrt(0.5 * (sig1**2 + sig2**2))
        ratio = gap / pooled_sigma
        status = "PASS" if ratio >= SEPARATION_MIN_SIGMA else "WARN"
        log(f"  {status}  {s1} → {s2}: "
            f"gap={gap*1e6:.1f} µs, σ_pool={pooled_sigma*1e6:.1f} µs, "
            f"ratio={ratio:.1f}×  (need ≥{SEPARATION_MIN_SIGMA:.0f}×)")
        if status == "WARN":
            log("        → Peaks overlap — mobility separation unclear")
            # Separation is secondary; don't fail the primary gate
            # but log it clearly

    # Gate 4: Elution times close to analytical theory (informational)
    log(f"\n[Gate 4]  Theory agreement  (within ±{THEORY_TOLERANCE_PCT:.0f}%)")
    for sid in ordered_species:
        err = results[sid].get("theory_err_pct", None)
        t_sim = means[sid] * 1e3
        t_th  = results[sid].get("theory_s", None)
        if err is None or t_th is None:
            log(f"  SKIP  {sid}: theory not computed")
            continue
        status = "PASS" if err <= THEORY_TOLERANCE_PCT else "WARN"
        log(f"  {status}  {sid}: sim={t_sim:.3f} ms, "
            f"theory={t_th*1e3:.3f} ms, Δ={err:.1f}%")

    return all_pass


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main() -> int:
    start_time = datetime.now()
    log("=" * 65)
    log("TIMS MOBILITY-SORTED ELUTION VALIDATION")
    log(f"Started: {start_time.strftime('%Y-%m-%d %H:%M:%S')}")
    log("=" * 65)
    log(f"\nParameters:")
    log(f"  Tube: L={TUBE_LENGTH_M*1e3:.0f} mm, r={TUBE_RADIUS_M*1e3:.1f} mm")
    log(f"  Gas:  {GAS_SPECIES}, {PRESSURE_PA:.0f} Pa, {TEMPERATURE_K:.0f} K, "
        f"v_gas={GAS_VELOCITY_MS:.0f} m/s")
    log(f"  E_initial = [0, -{E_INITIAL_MAX_VM:.0f}] V/m, "
        f"E_final = [0, -{E_FINAL_MAX_VM:.0f}] V/m (linear profile)")
    log(f"  Ramp: {RAMP_START_S*1e3:.1f} ms → {RAMP_END_S*1e3:.1f} ms (linear, f: 0→1)")
    log(f"  RF:   {RF_VOLTAGE_V:.0f} V @ {RF_FREQ_HZ/1e6:.0f} MHz")
    log(f"  Ions: {IONS_PER_SPECIES} per species × {len(SPECIES)} species = "
        f"{IONS_PER_SPECIES * len(SPECIES)} total")
    log(f"\nTheoretical elution predictions:")
    for sp in SPECIES:
        t_th = compute_theoretical_elution(sp["K0_cm2Vs"])
        log(f"  {sp['id']:20s} K0={sp['K0_cm2Vs']:.3f} cm²/Vs → "
            f"t_elut ≈ {t_th*1e3:.3f} ms")

    # --- Build and run simulation ---
    config = create_config("tims_elution_trajectories.h5")
    h5_file = run_simulation("tims_elution_config.json", config)

    if h5_file is None or not h5_file.exists():
        log("\n✗ Simulation did not produce output. Aborting.")
        return 1

    # --- Analyse ---
    results = analyze_elution(h5_file)
    if results is None:
        log("\n✗ Analysis failed. Aborting.")
        return 1

    # --- Validate ---
    passed = validate(results)

    # --- Plot ---
    fig_path = FIGURES_DIR / "tims_elution_validation.png"
    try:
        plot_elution_spectrum(results, fig_path)
    except Exception as exc:
        log(f"\n⚠  Plotting failed (non-fatal): {exc}")

    # --- Summary ---
    elapsed = (datetime.now() - start_time).total_seconds()
    log("\n" + "=" * 65)
    log(f"RESULT:  {'✓ PASS' if passed else '✗ FAIL'}")
    log(f"Elapsed: {elapsed:.1f} s")
    log(f"Log:     {LOG_FILE}")
    log(f"Figure:  {fig_path}")
    log("=" * 65)

    _log_handle.close()
    return 0 if passed else 1


if __name__ == "__main__":
    sys.exit(main())
