#!/bin/bash
# Export a validation run into a git-friendly bundle (logs/figures/config snapshots),
# excluding heavy simulation artifacts like .h5.

set -euo pipefail

print_usage() {
  cat <<'EOF'
Usage: export_run_bundle.sh --run-dir PATH --out-dir PATH

Example:
  ./validation/scripts/export_run_bundle.sh \
    --run-dir validation/runs/my_run \
    --out-dir validation/published/my_run

What it copies:
  - manifest.json (if present)
  - logs/
  - figures/
  - any JSON config snapshots under results/** (but NOT *.h5)

It does NOT copy:
  - any *.h5
EOF
}

RUN_DIR=""
OUT_DIR=""

while [[ $# -gt 0 ]]; do
  case "$1" in
    --run-dir)
      RUN_DIR="$2"; shift 2 ;;
    --out-dir)
      OUT_DIR="$2"; shift 2 ;;
    -h|--help)
      print_usage; exit 0 ;;
    *)
      echo "Unknown argument: $1" >&2
      print_usage >&2
      exit 1
      ;;
  esac
done

if [[ -z "$RUN_DIR" || -z "$OUT_DIR" ]]; then
  print_usage >&2
  exit 1
fi

# Normalize to absolute paths
if [[ "$RUN_DIR" != /* ]]; then
  RUN_DIR="$(pwd)/$RUN_DIR"
fi
if [[ "$OUT_DIR" != /* ]]; then
  OUT_DIR="$(pwd)/$OUT_DIR"
fi

if [[ ! -d "$RUN_DIR" ]]; then
  echo "Error: run dir not found: $RUN_DIR" >&2
  exit 1
fi

mkdir -p "$OUT_DIR"

copy_tree_excluding_h5() {
  local src="$1"
  local dst="$2"
  [[ -d "$src" ]] || return 0

  # Use tar for portability and to support excludes without rsync.
  tar --exclude='*.h5' -C "$src" -cf - . | tar -C "$dst" -xf -
}

# manifests
mkdir -p "$OUT_DIR/manifests"
if [[ -f "$RUN_DIR/manifest.json" ]]; then
  cp "$RUN_DIR/manifest.json" "$OUT_DIR/manifest.json"
fi
for mf in "$RUN_DIR"/manifest.*.json; do
  [[ -f "$mf" ]] || continue
  cp "$mf" "$OUT_DIR/manifests/$(basename "$mf")"
done

# logs and figures
mkdir -p "$OUT_DIR/logs" "$OUT_DIR/figures"
copy_tree_excluding_h5 "$RUN_DIR/logs" "$OUT_DIR/logs"
copy_tree_excluding_h5 "$RUN_DIR/figures" "$OUT_DIR/figures"

# config snapshots (heuristic): copy json files from results tree
mkdir -p "$OUT_DIR/configs"
if [[ -d "$RUN_DIR/results" ]]; then
  while IFS= read -r -d '' f; do
    rel="${f#$RUN_DIR/results/}"
    mkdir -p "$OUT_DIR/configs/$(dirname "$rel")"
    cp "$f" "$OUT_DIR/configs/$rel"
  done < <(find "$RUN_DIR/results" -type f -name '*.json' -print0)
fi

# lightweight index
{
  echo "# Published validation bundle"
  echo ""
  echo "Source run dir: $RUN_DIR"
  echo "Generated: $(date -Iseconds)"
  echo ""
  echo "Contents:"
  [[ -f "$OUT_DIR/manifest.json" ]] && echo "- manifest.json"
  if [[ -d "$OUT_DIR/manifests" ]] && [[ -n "$(ls -A "$OUT_DIR/manifests" 2>/dev/null)" ]]; then
    echo "- manifests/ (suite-level manifests)"
  fi
  echo "- logs/"
  echo "- figures/"
  echo "- configs/ (snapshotted JSON configs found under results/)"
  echo ""
  echo "Notes:"
  echo "- Heavy artifacts (*.h5) are intentionally excluded."
} > "$OUT_DIR/README.md"

echo "Exported bundle to: $OUT_DIR"
