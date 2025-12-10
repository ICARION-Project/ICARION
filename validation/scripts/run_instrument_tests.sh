#!/bin/bash
# Unified entrypoint for ICARION instrument validation suites.
# Dispatches to instrument-specific runners when available, otherwise runs
# a generic batch executor over the curated configs under validation/configs.

set -euo pipefail

usage() {
  cat <<'EOF'
Usage: ./scripts/run_instrument_tests.sh [options] <instrument> [instrument-args]

Options:
  -j, --jobs <N>          Number of concurrent simulations (default: 2)
  -t, --threads <N>       Threads passed to icarion_main (default: 4)
  -b, --binary <path>     Path to icarion_main (default: ../build/src/icarion_main)
  -c, --config-dir <dir>  Override config directory (generic mode only)
  -o, --output-root <dir> Override log/output root for this run (generic mode only)
  -h, --help              Show this help message

Examples:
  ./scripts/run_instrument_tests.sh ims
  ./scripts/run_instrument_tests.sh -j 4 quadrupole

When an instrument has a dedicated runner (e.g., IMS or Quadrupole), this
wrapper will delegate to that script so existing workflows keep working.
EOF
}

normalize_instrument() {
  local raw="${1,,}"
  case "$raw" in
    ims|ion-mobility|ionmobility)
      printf '%s' "ims"
      ;;
    quadrupole|quad|quadrupole_stability)
      printf '%s' "quadrupole"
      ;;
    fticr|ft-icr)
      printf '%s' "fticr"
      ;;
    lqit)
      printf '%s' "lqit"
      ;;
    orbitrap)
      printf '%s' "orbitrap"
      ;;
    tof)
      printf '%s' "tof"
      ;;
    *)
      printf '%s' "$raw"
      ;;
  esac
}

require_positive_integer() {
  local value="$1"
  local label="$2"
  if ! [[ "$value" =~ ^[1-9][0-9]*$ ]]; then
    echo "ERROR: $label must be a positive integer (got '$value')." >&2
    exit 1
  fi
}

to_abs_path() {
  local path="$1"
  if [[ -z "$path" ]]; then
    return 1
  fi
  if [[ "$path" = /* ]]; then
    printf '%s' "$path"
  else
    printf '%s/%s' "$(pwd)" "$path"
  fi
}

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
VALIDATION_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
REPO_ROOT="$(cd "$VALIDATION_DIR/.." && pwd)"
DEFAULT_BINARY="$REPO_ROOT/build/src/icarion_main"
DEFAULT_CONFIG_ROOT="$VALIDATION_DIR/configs/instruments"
DEFAULT_OUTPUT_ROOT="$VALIDATION_DIR/results/v1.0_test/instruments"

JOBS=2
THREADS=4
CONFIG_OVERRIDE=""
OUTPUT_OVERRIDE=""
BINARY_OVERRIDE=""
CUSTOM_OPTIONS_USED=false

while [[ $# -gt 0 ]]; do
  case "$1" in
    -j|--jobs)
      [[ $# -lt 2 ]] && usage && exit 1
      JOBS="$2"
      CUSTOM_OPTIONS_USED=true
      shift 2
      ;;
    -t|--threads)
      [[ $# -lt 2 ]] && usage && exit 1
      THREADS="$2"
      CUSTOM_OPTIONS_USED=true
      shift 2
      ;;
    -b|--binary)
      [[ $# -lt 2 ]] && usage && exit 1
      BINARY_OVERRIDE="$2"
      CUSTOM_OPTIONS_USED=true
      shift 2
      ;;
    -c|--config-dir)
      [[ $# -lt 2 ]] && usage && exit 1
      CONFIG_OVERRIDE="$2"
      CUSTOM_OPTIONS_USED=true
      shift 2
      ;;
    -o|--output-root)
      [[ $# -lt 2 ]] && usage && exit 1
      OUTPUT_OVERRIDE="$2"
      CUSTOM_OPTIONS_USED=true
      shift 2
      ;;
    -h|--help)
      usage
      exit 0
      ;;
    --)
      shift
      break
      ;;
    -*)
      echo "Unknown option: $1" >&2
      usage
      exit 1
      ;;
    *)
      break
      ;;
  esac
done

if [[ $# -lt 1 ]]; then
  usage
  exit 1
fi

RAW_INSTRUMENT="$1"
shift
INSTRUMENT_ARGS=("$@")
INSTRUMENT="$(normalize_instrument "$RAW_INSTRUMENT")"

require_positive_integer "$JOBS" "jobs"
require_positive_integer "$THREADS" "threads"

DELEGATE=""
case "$INSTRUMENT" in
  ims)
    DELEGATE="$SCRIPT_DIR/instrumentation/ims/run_ims_tests.sh"
    ;;
  quadrupole)
    DELEGATE="$SCRIPT_DIR/instrumentation/quadrupole/run_quad_tests.sh"
    ;;
  *)
    DELEGATE=""
    ;;
esac

if [[ -n "$DELEGATE" && -x "$DELEGATE" ]]; then
  if [[ "$CUSTOM_OPTIONS_USED" == true ]]; then
    echo "NOTE: '$INSTRUMENT' uses dedicated runner $DELEGATE; generic options (-j/-t/-c/-o/-b) are ignored."
  fi
  echo "Delegating to instrument-specific runner: $DELEGATE"
  exec "$DELEGATE" "${INSTRUMENT_ARGS[@]}"
fi

if [[ ${#INSTRUMENT_ARGS[@]} -gt 0 ]]; then
  echo "ERROR: Extra arguments '${INSTRUMENT_ARGS[*]}' are not supported in generic mode." >&2
  exit 1
fi

if [[ -n "$BINARY_OVERRIDE" ]]; then
  ICARION_BIN="$(to_abs_path "$BINARY_OVERRIDE")"
else
  ICARION_BIN="$DEFAULT_BINARY"
fi
if [[ ! -x "$ICARION_BIN" ]]; then
  echo "ERROR: ICARION binary not found or not executable: $ICARION_BIN" >&2
  echo "Build it via: cd \"$REPO_ROOT\"/build && make -j\$(nproc)" >&2
  exit 1
fi

if [[ -n "$CONFIG_OVERRIDE" ]]; then
  CONFIG_DIR="$(to_abs_path "$CONFIG_OVERRIDE")"
else
  CONFIG_DIR="$DEFAULT_CONFIG_ROOT/$INSTRUMENT"
fi
if [[ ! -d "$CONFIG_DIR" ]]; then
  echo "ERROR: Config directory not found: $CONFIG_DIR" >&2
  exit 1
fi
CONFIG_DIR="$(cd "$CONFIG_DIR" && pwd)"

if [[ -n "$OUTPUT_OVERRIDE" ]]; then
  OUTPUT_ROOT="$(to_abs_path "$OUTPUT_OVERRIDE")"
else
  OUTPUT_ROOT="$DEFAULT_OUTPUT_ROOT/$INSTRUMENT"
fi
mkdir -p "$OUTPUT_ROOT"
OUTPUT_ROOT="$(cd "$OUTPUT_ROOT" && pwd)"
SESSION_STAMP="$(date +%Y%m%d_%H%M%S)"
SESSION_LOG_DIR="$OUTPUT_ROOT/run_logs/$SESSION_STAMP"
mkdir -p "$SESSION_LOG_DIR"

mapfile -t CONFIG_FILES < <(find "$CONFIG_DIR" -type f -name '*.json' | sort)
CONFIG_TOTAL=${#CONFIG_FILES[@]}
if [[ $CONFIG_TOTAL -eq 0 ]]; then
  echo "ERROR: No JSON configs found under $CONFIG_DIR" >&2
  exit 1
fi

printf '==============================================\n'
printf 'ICARION Instrument Validation (generic mode)\n'
printf '==============================================\n'
printf 'Instrument : %s\n' "$INSTRUMENT"
printf 'Configs    : %s\n' "$CONFIG_DIR"
printf 'Binary     : %s\n' "$ICARION_BIN"
printf 'Jobs       : %s\n' "$JOBS"
printf 'Threads    : %s\n' "$THREADS"
printf 'Log dir    : %s\n' "$SESSION_LOG_DIR"
printf 'Started    : %s\n' "$(date)"
printf 'Total runs : %d\n' "$CONFIG_TOTAL"
printf '==============================================\n\n'

PASSED=0
FAILED=0
declare -a FAILED_CONFIGS=()
declare -a PIDS=()
declare -A PID_TO_CONFIG=()

cleanup() {
  if [[ ${#PIDS[@]} -gt 0 ]]; then
    kill "${PIDS[@]}" 2>/dev/null || true
  fi
}
trap cleanup INT TERM

run_config() {
  local config="$1"
  local ordinal="$2"
  local total="$3"
  local name
  name="$(basename "$config" .json)"
  local stdout_log="$SESSION_LOG_DIR/${name}_stdout.log"
  local stderr_log="$SESSION_LOG_DIR/${name}_stderr.log"

  printf '[%s] [%d/%d] Starting %s\n' "$(date '+%H:%M:%S')" "$ordinal" "$total" "$name"
  if (cd "$REPO_ROOT" && "$ICARION_BIN" "$config" --threads "$THREADS" >"$stdout_log" 2>"$stderr_log"); then
    printf '[%s] [%d/%d] ✅ %s\n' "$(date '+%H:%M:%S')" "$ordinal" "$total" "$name"
    return 0
  else
    printf '[%s] [%d/%d] ❌ %s\n' "$(date '+%H:%M:%S')" "$ordinal" "$total" "$name"
    tail -n 20 "$stderr_log" >&2 || true
    return 1
  fi
}

wait_on_pid() {
  local pid="$1"
  local config="${PID_TO_CONFIG[$pid]}"
  if wait "$pid"; then
    ((PASSED++))
  else
    ((FAILED++))
    FAILED_CONFIGS+=("$config")
  fi
  unset PID_TO_CONFIG["$pid"]
}

for idx in "${!CONFIG_FILES[@]}"; do
  config="${CONFIG_FILES[$idx]}"
  ordinal=$((idx + 1))
  run_config "$config" "$ordinal" "$CONFIG_TOTAL" &
  pid=$!
  PIDS+=("$pid")
  PID_TO_CONFIG["$pid"]="$config"

  if [[ ${#PIDS[@]} -ge $JOBS ]]; then
    wait_on_pid "${PIDS[0]}"
    PIDS=("${PIDS[@]:1}")
  fi
done

for pid in "${PIDS[@]}"; do
  wait_on_pid "$pid"
done

SUCCESS_RATE=$(awk -v p="$PASSED" -v t="$CONFIG_TOTAL" 'BEGIN { if (t == 0) { print "0.0" } else { printf "%.1f", (p/t)*100 } }')

printf '\n==============================================\n'
printf 'Validation summary\n'
printf '==============================================\n'
printf 'Total configs : %d\n' "$CONFIG_TOTAL"
printf 'Passed        : %d\n' "$PASSED"
printf 'Failed        : %d\n' "$FAILED"
printf 'Success rate  : %s%%%%\n' "$SUCCESS_RATE"
printf 'Logs captured : %s\n' "$SESSION_LOG_DIR"
printf 'Finished      : %s\n' "$(date)"
printf '==============================================\n'

if [[ $FAILED -gt 0 ]]; then
  printf '\nFailed configurations:\n'
  for cfg in "${FAILED_CONFIGS[@]}"; do
    printf '  - %s\n' "$(basename "$cfg")"
  done
  printf '\nReview stderr logs under %s for details.\n' "$SESSION_LOG_DIR"
fi

printf "\nNote: Simulation outputs follow each config's 'output.folder'.\n"
if [[ $FAILED -eq 0 ]]; then
  STATUS="OK"
else
  STATUS="FAIL"
fi
printf 'Generic runner exit status: %s\n' "$STATUS"

if [[ $FAILED -gt 0 ]]; then
  exit 1
fi
