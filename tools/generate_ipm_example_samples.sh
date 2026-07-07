#!/usr/bin/env bash
set -euo pipefail

# Regenerate the small IPM example tables used by examples/ims/ims_ipm_basic.json.
# Defaults are intentionally small so the script is useful as a reproducible
# example. Increase N_TRIALS/N_ORIENTATIONS for production-quality sample sets.

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
BIN="${BIN:-$ROOT/build/src/interaction_potential_precompute}"
OUT_DIR="${OUT_DIR:-$ROOT/data/molecules/precomputed_ipm}"
SPECIES_DB="${SPECIES_DB:-$ROOT/data/species_database_v1.json}"
SPECIES_LIST="${SPECIES_LIST:-H3O+ PentanalH+}"
GAS="${GAS:-He}"
N_ORIENTATIONS="${N_ORIENTATIONS:-12}"
N_TRIALS="${N_TRIALS:-1000}"
V_BINS="${V_BINS:-16}"
V_MIN="${V_MIN:-200}"
V_MAX="${V_MAX:-20000}"
THREADS="${THREADS:-0}"
POTENTIAL="${POTENTIAL:-lj1264}"
POLARIZATION="${POLARIZATION:-partial}"
MIXING_RULE="${MIXING_RULE:-lb}"
ORIENT_GRID="${ORIENT_GRID:-random}"
FORCEFIELD="${FORCEFIELD:-$OUT_DIR/ipm_example_forcefield.json}"

if [[ ! -x "$BIN" ]]; then
  echo "Missing precompute binary: $BIN" >&2
  echo "Build it with: cmake --build build --target interaction_potential_precompute" >&2
  exit 2
fi

mkdir -p "$OUT_DIR"

cat > "$FORCEFIELD" <<'JSON'
{
  "gases": {
    "He": {
      "sigma_A": 2.551,
      "epsilon_eV": 0.0008798,
      "alpha_A3": 0.204956
    }
  },
  "elements": {
    "H": { "sigma_A": 2.261, "epsilon_eV": 0.0989235 },
    "C": { "sigma_A": 3.0126, "epsilon_eV": 0.21252132 },
    "N": { "sigma_A": 3.0126, "epsilon_eV": 0.21252132 },
    "O": { "sigma_A": 2.4344, "epsilon_eV": 0.1717344 }
  }
}
JSON

echo "Writing IPM example samples to: $OUT_DIR"
echo "Settings: gas=$GAS orientations=$N_ORIENTATIONS trials=$N_TRIALS v_bins=$V_BINS v=[$V_MIN,$V_MAX]"

for SPECIES in $SPECIES_LIST; do
  OUT="$OUT_DIR/${SPECIES}_ipm_samples_${GAS}.h5"
  echo
  echo "==> $SPECIES / $GAS"
  "$BIN" \
    --input "$SPECIES_DB" \
    --species "$SPECIES" \
    --gas "$GAS" \
    --output "$OUT" \
    --gas-params "$FORCEFIELD" \
    --element-params "$FORCEFIELD" \
    --potential "$POTENTIAL" \
    --polarization "$POLARIZATION" \
    --mixing-rule "$MIXING_RULE" \
    --orient-grid "$ORIENT_GRID" \
    --n-orientations "$N_ORIENTATIONS" \
    --n-trials "$N_TRIALS" \
    --v-bins "$V_BINS" \
    --v-min "$V_MIN" \
    --v-max "$V_MAX" \
    --threads "$THREADS"
done

echo
echo "Generated:"
find "$OUT_DIR" -maxdepth 1 -type f -name '*_ipm_samples_*.h5' -printf '  %p\n'
