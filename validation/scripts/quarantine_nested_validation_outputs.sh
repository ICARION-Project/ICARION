#!/bin/bash
# Safely move accidental nested validation output trees out of the way.
# This does NOT delete anything.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
VALIDATION_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
REPO_ROOT="$(cd "$VALIDATION_DIR/.." && pwd)"

TS=$(date +%Y%m%d_%H%M%S)
QUARANTINE_DIR="$VALIDATION_DIR/_quarantine/$TS"

mkdir -p "$QUARANTINE_DIR"

moved_any=false
move_if_exists() {
  local src="$1"
  local name="$2"

  if [[ -d "$src" ]]; then
    echo "→ Moving '$src' → '$QUARANTINE_DIR/$name'"
    mkdir -p "$QUARANTINE_DIR"
    mv "$src" "$QUARANTINE_DIR/$name"
    moved_any=true
  fi
}

move_if_exists "$VALIDATION_DIR/validation" "validation"
move_if_exists "$VALIDATION_DIR/scripts/validation" "scripts_validation"
move_if_exists "$VALIDATION_DIR/scripts/results" "scripts_results"

if ! $moved_any; then
  echo "Nothing to quarantine."
  exit 0
fi

echo ""
echo "Quarantined under: $QUARANTINE_DIR"
echo "If everything looks good, you can delete that folder later."
