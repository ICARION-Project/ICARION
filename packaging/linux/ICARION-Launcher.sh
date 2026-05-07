#!/usr/bin/env bash
set -u

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"

if [[ -f "$SCRIPT_DIR/bin/icarion" ]]; then
  ROOT="$SCRIPT_DIR"
  ICARION_EXE="$SCRIPT_DIR/bin/icarion"
  if [[ ! -x "$ICARION_EXE" ]]; then
    chmod +x "$ICARION_EXE" 2>/dev/null || true
  fi
  if [[ -d "$SCRIPT_DIR/examples" ]]; then
    EXAMPLES_DIR="$SCRIPT_DIR/examples"
  else
    EXAMPLES_DIR="$SCRIPT_DIR/share/icarion/examples"
  fi
  LOG_DIR="$SCRIPT_DIR/launcher-logs"
else
  ROOT="$SCRIPT_DIR"
  ICARION_EXE="$(command -v icarion || true)"
  EXAMPLES_DIR="/usr/share/icarion/examples"
  LOG_DIR="${XDG_STATE_HOME:-$HOME/.local/state}/icarion/launcher-logs"
fi

has_cmd() {
  command -v "$1" >/dev/null 2>&1
}

show_error() {
  if has_cmd zenity; then
    zenity --error --title="ICARION Launcher" --width=520 --text="$1"
  else
    printf 'ICARION Launcher error: %s\n' "$1" >&2
  fi
}

show_info() {
  if has_cmd zenity; then
    zenity --info --title="ICARION Launcher" --width=520 --text="$1"
  else
    printf '%s\n' "$1"
  fi
}

pick_config_gui() {
  local start_dir="$EXAMPLES_DIR"
  [[ -d "$start_dir" ]] || start_dir="$ROOT"
  zenity \
    --file-selection \
    --title="Select ICARION JSON config" \
    --filename="$start_dir/" \
    --file-filter="JSON configs | *.json" \
    --file-filter="All files | *"
}

pick_config_terminal() {
  local default_config=""
  if [[ -f "$EXAMPLES_DIR/ims/ims_basic.json" ]]; then
    default_config="$EXAMPLES_DIR/ims/ims_basic.json"
  fi

  if [[ -n "$default_config" ]]; then
    printf 'Config file [%s]: ' "$default_config" >&2
  else
    printf 'Config file: ' >&2
  fi
  read -r config_path
  if [[ -z "$config_path" ]]; then
    config_path="$default_config"
  fi
  printf '%s\n' "$config_path"
}

find_terminal() {
  local candidate
  for candidate in gnome-terminal konsole xfce4-terminal mate-terminal lxterminal xterm; do
    if has_cmd "$candidate"; then
      printf '%s\n' "$candidate"
      return 0
    fi
  done
  return 1
}

run_in_terminal() {
  local terminal="$1"
  local command_text="$2"

  case "$terminal" in
    gnome-terminal|mate-terminal|xfce4-terminal)
      "$terminal" -- bash -lc "$command_text"
      ;;
    konsole|xterm|lxterminal)
      "$terminal" -e bash -lc "$command_text"
      ;;
    *)
      "$terminal" -e bash -lc "$command_text"
      ;;
  esac
}

if [[ -z "$ICARION_EXE" || ! -x "$ICARION_EXE" ]]; then
  show_error "Cannot find the ICARION executable. Expected bin/icarion next to this launcher, or icarion in PATH."
  exit 1
fi

if has_cmd zenity; then
  CONFIG_PATH="$(pick_config_gui)"
  status=$?
  if [[ $status -ne 0 || -z "$CONFIG_PATH" ]]; then
    exit 0
  fi
else
  CONFIG_PATH="$(pick_config_terminal)"
fi

if [[ ! -f "$CONFIG_PATH" ]]; then
  show_error "This is not a file:
$CONFIG_PATH"
  exit 1
fi

if [[ "${CONFIG_PATH##*.}" != "json" ]]; then
  show_error "ICARION configs must be .json files:
$CONFIG_PATH"
  exit 1
fi

mkdir -p "$LOG_DIR"
RUN_STAMP="$(date '+%Y%m%d-%H%M%S')"
LOG_PATH="$LOG_DIR/icarion-run-$RUN_STAMP.log"

{
  printf 'ICARION: %s\n' "$ICARION_EXE"
  printf 'Config:  %s\n' "$CONFIG_PATH"
  printf 'Log:     %s\n\n' "$LOG_PATH"
} > "$LOG_PATH"

RUN_COMMAND=$(printf 'cd %q; printf "ICARION: %%s\\nConfig:  %%s\\nLog:     %%s\\n\\n" %q %q %q; %q %q 2>&1 | tee -a %q; code=${PIPESTATUS[0]}; printf "\\nFinished with exit code %%s\\n" "$code" | tee -a %q; printf "\\nPress Enter to close..."; read -r _; exit "$code"' \
  "$ROOT" "$ICARION_EXE" "$CONFIG_PATH" "$LOG_PATH" "$ICARION_EXE" "$CONFIG_PATH" "$LOG_PATH" "$LOG_PATH")

if [[ -t 0 && -t 1 ]]; then
  show_info "Running ICARION in this terminal.

Log:
$LOG_PATH"
  bash -lc "$RUN_COMMAND"
elif terminal="$(find_terminal)"; then
  run_in_terminal "$terminal" "$RUN_COMMAND"
  show_info "ICARION started.

Log:
$LOG_PATH"
else
  show_info "No graphical terminal was found. Running in this terminal.

Log:
$LOG_PATH"
  bash -lc "$RUN_COMMAND"
fi
