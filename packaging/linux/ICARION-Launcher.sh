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
  if [[ -d "$SCRIPT_DIR/analysis" ]]; then
    ANALYSIS_DIR="$SCRIPT_DIR/analysis"
  else
    ANALYSIS_DIR="$SCRIPT_DIR/share/icarion/analysis"
  fi
else
  ROOT="$SCRIPT_DIR"
  ICARION_EXE="$(command -v icarion || true)"
  EXAMPLES_DIR="/usr/share/icarion/examples"
  LOG_DIR="${XDG_STATE_HOME:-$HOME/.local/state}/icarion/launcher-logs"
  ANALYSIS_DIR="/usr/share/icarion/analysis"
fi

has_cmd() {
  command -v "$1" >/dev/null 2>&1
}

has_gui_dialog() {
  has_cmd zenity && [[ -n "${DISPLAY:-}${WAYLAND_DISPLAY:-}" ]]
}

show_error() {
  if has_gui_dialog; then
    zenity --error --title="ICARION Launcher" --width=520 --text="$1"
  else
    printf 'ICARION Launcher error: %s\n' "$1" >&2
  fi
}

show_info() {
  if has_gui_dialog; then
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

pick_trajectory_gui() {
  local start_dir="$ROOT"
  zenity \
    --file-selection \
    --title="Select ICARION trajectory HDF5 file" \
    --filename="$start_dir/" \
    --file-filter="HDF5 trajectories | *.h5 *.hdf5" \
    --file-filter="All files | *"
}

pick_trajectory_terminal() {
  printf 'Trajectory HDF5 file: ' >&2
  read -r traj_path
  printf '%s\n' "$traj_path"
}

pick_action_gui() {
  zenity \
    --list \
    --title="ICARION Launcher" \
    --width=420 \
    --height=260 \
    --column="Action" \
    "Run simulation" \
    "Analyze IMS mobility" \
    "Plot arrival times" \
    "Plot trajectories" \
    "Plot mean positions" \
    "Plot eliminations" \
    "Animate trajectories"
}

pick_action_terminal() {
  printf '1) Run simulation\n2) Analyze IMS mobility\n3) Plot arrival times\n4) Plot trajectories\n5) Plot mean positions\n6) Plot eliminations\n7) Animate trajectories\nSelect action [1]: ' >&2
  read -r action
  case "$action" in
    2) printf 'Analyze IMS mobility\n' ;;
    3) printf 'Plot arrival times\n' ;;
    4) printf 'Plot trajectories\n' ;;
    5) printf 'Plot mean positions\n' ;;
    6) printf 'Plot eliminations\n' ;;
    7) printf 'Animate trajectories\n' ;;
    *) printf 'Run simulation\n' ;;
  esac
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

if has_gui_dialog; then
  ACTION="$(pick_action_gui)"
  status=$?
  if [[ $status -ne 0 || -z "$ACTION" ]]; then
    exit 0
  fi
else
  ACTION="$(pick_action_terminal)"
fi

if has_gui_dialog; then
  if [[ "$ACTION" == "Run simulation" ]]; then
    CONFIG_PATH="$(pick_config_gui)"
    status=$?
    if [[ $status -ne 0 || -z "$CONFIG_PATH" ]]; then
      exit 0
    fi
  else
    TRAJ_PATH="$(pick_trajectory_gui)"
    status=$?
    if [[ $status -ne 0 || -z "$TRAJ_PATH" ]]; then
      exit 0
    fi
  fi
else
  if [[ "$ACTION" == "Run simulation" ]]; then
    CONFIG_PATH="$(pick_config_terminal)"
  else
    TRAJ_PATH="$(pick_trajectory_terminal)"
  fi
fi

if [[ "$ACTION" == "Run simulation" ]]; then
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
else
  if [[ ! -f "$TRAJ_PATH" ]]; then
    show_error "This is not a file:
$TRAJ_PATH"
    exit 1
  fi
  case "$TRAJ_PATH" in
    *.h5|*.hdf5) ;;
    *)
      show_error "ICARION analysis expects a .h5 or .hdf5 trajectory file:
$TRAJ_PATH"
      exit 1
      ;;
  esac
  if [[ ! -d "$ANALYSIS_DIR" ]]; then
    show_error "Cannot find packaged analysis scripts:
$ANALYSIS_DIR"
    exit 1
  fi
  if ! has_cmd python3; then
    show_error "Python 3 was not found. Install Python 3 and packages from requirements-analysis.txt."
    exit 1
  fi
fi

mkdir -p "$LOG_DIR"
RUN_STAMP="$(date '+%Y%m%d-%H%M%S')"
if [[ "$ACTION" == "Run simulation" ]]; then
  LOG_PATH="$LOG_DIR/icarion-run-$RUN_STAMP.log"

  {
    printf 'ICARION: %s\n' "$ICARION_EXE"
    printf 'Config:  %s\n' "$CONFIG_PATH"
    printf 'Log:     %s\n\n' "$LOG_PATH"
  } > "$LOG_PATH"

  RUN_COMMAND=$(printf 'cd %q; printf "ICARION: %%s\\nConfig:  %%s\\nLog:     %%s\\n\\n" %q %q %q; %q %q 2>&1 | tee -a %q; code=${PIPESTATUS[0]}; printf "\\nFinished with exit code %%s\\n" "$code" | tee -a %q; printf "\\nPress Enter to close..."; read -r _; exit "$code"' \
    "$ROOT" "$ICARION_EXE" "$CONFIG_PATH" "$LOG_PATH" "$ICARION_EXE" "$CONFIG_PATH" "$LOG_PATH" "$LOG_PATH")
else
  OUT_DIR="$ROOT/analysis-output"
  mkdir -p "$OUT_DIR"
  case "$ACTION" in
    "Analyze IMS mobility")
      ANALYSIS_NAME="ims-mobility"
      SCRIPT_PATH="$ANALYSIS_DIR/ims_mobility.py"
      OUTPUT_EXT=".png"
      CSV_OPTION="--out-csv"
      EXTRA_ARGS=""
      ;;
    "Plot arrival times")
      ANALYSIS_NAME="arrival-times"
      SCRIPT_PATH="$ANALYSIS_DIR/arrival_time_distribution.py"
      OUTPUT_EXT=".png"
      CSV_OPTION="--out-csv"
      EXTRA_ARGS=""
      ;;
    "Plot trajectories")
      ANALYSIS_NAME="trajectories"
      SCRIPT_PATH="$ANALYSIS_DIR/plot_trajectories.py"
      OUTPUT_EXT=".png"
      CSV_OPTION=""
      EXTRA_ARGS="--time-stride 5 --max-ions 120 --max-per-species 40"
      ;;
    "Plot mean positions")
      ANALYSIS_NAME="mean-positions"
      SCRIPT_PATH="$ANALYSIS_DIR/mean_positions.py"
      OUTPUT_EXT=".png"
      CSV_OPTION="--csv"
      EXTRA_ARGS="--time-stride 5 --max-ions 400 --max-per-species 200"
      ;;
    "Plot eliminations")
      ANALYSIS_NAME="eliminations"
      SCRIPT_PATH="$ANALYSIS_DIR/elimination_histograms.py"
      OUTPUT_EXT=".png"
      CSV_OPTION=""
      EXTRA_ARGS="--max-ions 500 --max-per-species 200"
      ;;
    "Animate trajectories")
      ANALYSIS_NAME="animation"
      SCRIPT_PATH="$ANALYSIS_DIR/animate_trajectories.py"
      OUTPUT_EXT=".gif"
      CSV_OPTION=""
      EXTRA_ARGS="--projection xy --style trail --theme dark --max-ions 80 --max-per-species 25 --time-stride 4 --frame-step 2 --max-frames 300 --writer pillow --fps 15 --dpi 120"
      ;;
    *)
      show_error "Unknown action: $ACTION"
      exit 1
      ;;
  esac
  LOG_PATH="$LOG_DIR/$ANALYSIS_NAME-$RUN_STAMP.log"
  PLOT_PATH="$OUT_DIR/$ANALYSIS_NAME-$RUN_STAMP$OUTPUT_EXT"
  CSV_PATH=""
  if [[ -n "$CSV_OPTION" ]]; then
    CSV_PATH="$OUT_DIR/$ANALYSIS_NAME-$RUN_STAMP.csv"
  fi
  PER_SPECIES_PATH=""
  PLOT_ARGS=""
  if [[ "$ACTION" == "Plot arrival times" ]]; then
    PER_SPECIES_PATH="$OUT_DIR/$ANALYSIS_NAME-$RUN_STAMP-per-species.png"
    PLOT_ARGS="--out-per-species $(printf '%q' "$PER_SPECIES_PATH")"
  fi
  CSV_ARGS=""
  if [[ -n "$CSV_OPTION" ]]; then
    CSV_ARGS="$CSV_OPTION $(printf '%q' "$CSV_PATH")"
  fi

  {
    printf 'Analysis: %s\n' "$ACTION"
    printf 'Trajectory: %s\n' "$TRAJ_PATH"
    printf 'Plot: %s\n' "$PLOT_PATH"
    if [[ -n "$PER_SPECIES_PATH" ]]; then
      printf 'Per-species plot: %s\n' "$PER_SPECIES_PATH"
    fi
    printf 'CSV: %s\n' "$CSV_PATH"
    printf 'Log: %s\n\n' "$LOG_PATH"
  } > "$LOG_PATH"

  RUN_COMMAND=$(printf 'cd %q; export PYTHONPATH=%q:%q:${PYTHONPATH:-}; printf "Analysis: %%s\\nTrajectory: %%s\\nPlot: %%s\\nCSV: %%s\\nLog: %%s\\n\\n" %q %q %q %q %q; python3 %q --traj %q --out %q %s %s %s 2>&1 | tee -a %q; code=${PIPESTATUS[0]}; printf "\\nFinished with exit code %%s\\n" "$code" | tee -a %q; printf "\\nPress Enter to close..."; read -r _; exit "$code"' \
    "$ROOT" "$ROOT" "$ANALYSIS_DIR" "$ACTION" "$TRAJ_PATH" "$PLOT_PATH" "$CSV_PATH" "$LOG_PATH" "$SCRIPT_PATH" "$TRAJ_PATH" "$PLOT_PATH" "$CSV_ARGS" "$PLOT_ARGS" "$EXTRA_ARGS" "$LOG_PATH" "$LOG_PATH")
fi

if [[ -t 0 && -t 1 ]]; then
  show_info "Running selected ICARION action in this terminal.

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
