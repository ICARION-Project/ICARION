#!/bin/bash
# Sequentially execute the instrument validation runners.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}" )" && pwd)"
VALIDATION_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
RUNNER="$SCRIPT_DIR/run_instrument_tests.sh"

DEFAULT_INSTRUMENTS=(ims fticr lqit orbitrap tof quadrupole)
SELECTED_INSTRUMENTS=()
JOBS=""
THREADS=""
BINARY=""
CONFIG_ROOT=""
OUTPUT_ROOT=""

print_usage() {
  cat <<'EOF'
Usage: ./run_instrument_suite.sh [options] [instrument ...]

Options:
  -j, --jobs N          Parallel jobs passed to each run_instrument_tests invocation
  -t, --threads N       Threads per simulation forwarded to the runner
  -b, --binary PATH     Override path to icarion_main
  -c, --config-root DIR Use DIR/<instrument> for config overrides
  -o, --output-root DIR Use DIR/<instrument> for output overrides
  --list                Show normalized instrument keys and exit
  -h, --help            Show this help text

If no instruments are specified, the suite runs: ims, fticr, lqit, orbitrap, tof, quadrupole.
EOF
}

list_instruments() {
  cat <<'EOF'
Available instruments:
  ims          Ion Mobility Spectrometry (delegates to IMS runner)
  fticr        Fourier Transform Ion Cyclotron Resonance
  lqit         Linear Quadrupole Ion Trap
  orbitrap     Orbitrap axial frequency suite
  tof          Time-of-flight mass spectrometer
  quadrupole   Quadrupole stability map
EOF
}

abs_path() {
  local target="$1"
  if [[ "$target" = /* ]]; then
    printf '%s' "$target"
  else
    printf '%s/%s' "$(pwd)" "$target"
  fi
}

normalize_instrument() {
  local token="${1,,}"
  case "$token" in
    ims|ion-mobility|ionmobility) echo "ims" ;;
    fticr|ft-icr) echo "fticr" ;;
    lqit) echo "lqit" ;;
    orbitrap) echo "orbitrap" ;;
    tof|time-of-flight|timeofflight) echo "tof" ;;
    quadrupole|quad|quadrupole_stability) echo "quadrupole" ;;
    all) echo "all" ;;
    *) return 1 ;;
  esac
}

dedupe_array() {
  local -n arr=$1
  local -A seen=()
  local unique=()
  for item in "${arr[@]:-}"; do
    if [[ -n "$item" && -z "${seen[$item]:-}" ]]; then
      unique+=("$item")
      seen[$item]=1
    fi
  done
  arr=("${unique[@]:-}")
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    -j|--jobs)
      [[ $# -lt 2 ]] && print_usage && exit 1
      JOBS="$2"
      shift 2
      ;;
    -t|--threads)
      [[ $# -lt 2 ]] && print_usage && exit 1
      THREADS="$2"
      shift 2
      ;;
    -b|--binary)
      [[ $# -lt 2 ]] && print_usage && exit 1
      BINARY="$2"
      shift 2
      ;;
    -c|--config-root)
      [[ $# -lt 2 ]] && print_usage && exit 1
      CONFIG_ROOT="$2"
      shift 2
      ;;
    -o|--output-root)
      [[ $# -lt 2 ]] && print_usage && exit 1
      OUTPUT_ROOT="$2"
      shift 2
      ;;
    --list)
      list_instruments
      exit 0
      ;;
    -h|--help)
      print_usage
      exit 0
      ;;
    --)
      shift
      break
      ;;
    -*)
      echo "Unknown option: $1" >&2
      print_usage >&2
      exit 1
      ;;
    *)
      if ! normalized=$(normalize_instrument "$1"); then
        echo "Unknown instrument: $1" >&2
        exit 1
      fi
      if [[ "$normalized" == "all" ]]; then
        SELECTED_INSTRUMENTS+=("${DEFAULT_INSTRUMENTS[@]}")
      else
        SELECTED_INSTRUMENTS+=("$normalized")
      fi
      shift
      ;;
  esac
done

if [[ ${#SELECTED_INSTRUMENTS[@]} -eq 0 ]]; then
  SELECTED_INSTRUMENTS=("${DEFAULT_INSTRUMENTS[@]}")
fi

dedupe_array SELECTED_INSTRUMENTS

if [[ ! -x "$RUNNER" ]]; then
  echo "Error: Expected runner not found at $RUNNER" >&2
  exit 1
fi

if [[ -n "$CONFIG_ROOT" ]]; then
  CONFIG_ROOT="$(abs_path "$CONFIG_ROOT")"
fi
if [[ -n "$OUTPUT_ROOT" ]]; then
  OUTPUT_ROOT="$(abs_path "$OUTPUT_ROOT")"
fi

TOTAL=0
PASSED=0
FAILED=0
FAILED_LIST=()

printf '==============================================\n'
printf 'Instrument Validation Suite\n'
printf '==============================================\n'
printf 'Runner      : %s\n' "$RUNNER"
printf 'Jobs        : %s\n' "${JOBS:-default}"
printf 'Threads     : %s\n' "${THREADS:-default}"
printf 'Binary      : %s\n' "${BINARY:-auto}"
printf 'Config root : %s\n' "${CONFIG_ROOT:-validation/configs/instruments}" 
printf 'Output root : %s\n' "${OUTPUT_ROOT:-validation/results/instruments}"
printf 'Targets     : %s\n' "${SELECTED_INSTRUMENTS[*]}"
printf '==============================================\n\n'

run_single() {
  local instrument="$1"
  local cmd=("$RUNNER")
  if [[ -n "$JOBS" ]]; then
    cmd+=("-j" "$JOBS")
  fi
  if [[ -n "$THREADS" ]]; then
    cmd+=("-t" "$THREADS")
  fi
  if [[ -n "$BINARY" ]]; then
    cmd+=("-b" "$BINARY")
  fi
  if [[ -n "$CONFIG_ROOT" ]]; then
    cmd+=("-c" "$CONFIG_ROOT/$instrument")
  fi
  if [[ -n "$OUTPUT_ROOT" ]]; then
    cmd+=("-o" "$OUTPUT_ROOT/$instrument")
  fi
  cmd+=("$instrument")

  echo "[${instrument}] starting..."
  if "${cmd[@]}"; then
    echo "[${instrument}] ✅ completed"
    return 0
  else
    echo "[${instrument}] ❌ failed" >&2
    return 1
  fi
}

for instrument in "${SELECTED_INSTRUMENTS[@]}"; do
  TOTAL=$((TOTAL + 1))
  if run_single "$instrument"; then
    PASSED=$((PASSED + 1))
  else
    FAILED=$((FAILED + 1))
    FAILED_LIST+=("$instrument")
  fi
done

printf '\n==============================================\n'
printf 'Suite summary\n'
printf '==============================================\n'
printf 'Total instruments : %d\n' "$TOTAL"
printf 'Succeeded         : %d\n' "$PASSED"
printf 'Failed            : %d\n' "$FAILED"
if [[ $FAILED -gt 0 ]]; then
  printf 'Failures          : %s\n' "${FAILED_LIST[*]}"
fi
printf 'Finished          : %s\n' "$(date)"
printf '==============================================\n'

if [[ $FAILED -gt 0 ]]; then
  exit 1
fi

exit 0
