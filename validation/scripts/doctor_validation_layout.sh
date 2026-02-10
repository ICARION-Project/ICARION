#!/bin/bash
# Diagnose common validation output layout issues.
# In particular: accidental nested output trees like validation/validation/results
# caused by running icarion_main with a CWD inside validation/.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
VALIDATION_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
REPO_ROOT="$(cd "$VALIDATION_DIR/.." && pwd)"

print_dir_info() {
  local path="$1"
  if [[ -d "$path" ]]; then
    # show size summary without failing if dir is empty
    local size
    size=$(du -sh "$path" 2>/dev/null | awk '{print $1}' || true)
    echo "- Found: $path${size:+ (size: $size)}"
  fi
}

echo "============================================================"
echo "ICARION validation layout doctor"
echo "Repo root     : $REPO_ROOT"
echo "Validation dir: $VALIDATION_DIR"
echo "============================================================"
echo ""

FOUND=false
for p in \
  "$VALIDATION_DIR/validation" \
  "$VALIDATION_DIR/scripts/validation" \
  "$VALIDATION_DIR/scripts/results"; do
  if [[ -d "$p" ]]; then
    FOUND=true
    print_dir_info "$p"
  fi
done

if ! $FOUND; then
  echo "✅ No accidental nested validation output directories detected."
  exit 0
fi

echo ""
echo "Why this happens"
echo "- Many JSON configs set output.folder to a relative path like 'validation/results/...'."
echo "- If you run icarion_main with CWD inside 'validation/', it resolves to 'validation/validation/results/...'."
echo ""
echo "What to do"
echo "- Preferred: use suite runners; they run icarion_main from the repo root."
echo "- Safe cleanup: move these folders aside (no deletion)."
echo "  Run: validation/scripts/quarantine_nested_validation_outputs.sh"
echo ""
