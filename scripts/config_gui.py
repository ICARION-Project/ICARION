#!/usr/bin/env python3
"""
ICARION JSON Config GUI (AC/RF editor)

Simple desktop GUI to load an ICARION JSON config, edit AC/RF frequency and voltage
for each domain, and save the result.

Run:
    python3 scripts/config_gui.py
    python3 scripts/config_gui.py --config examples/lqit/lqit_basic.json
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import tkinter as tk
from tkinter import filedialog, messagebox, scrolledtext, ttk
from typing import Any


class ConfigGui:
    def __init__(self, root: tk.Tk, initial_path: Path | None = None) -> None:
        self.root = root
        self.root.title("ICARION Config GUI (AC/RF)")
        self.root.geometry("980x620")
        self.root.minsize(900, 560)

        self.config_path: Path | None = None
        self.config: dict[str, Any] | None = None

        self.domain_index_var = tk.StringVar(value="")
        self.status_var = tk.StringVar(value="Ready")

        self.rf_voltage_var = tk.StringVar(value="")
        self.rf_frequency_var = tk.StringVar(value="")
        self.ac_voltage_var = tk.StringVar(value="")
        self.ac_frequency_var = tk.StringVar(value="")

        self.rf_voltage_info = tk.StringVar(value="")
        self.rf_frequency_info = tk.StringVar(value="")
        self.ac_voltage_info = tk.StringVar(value="")
        self.ac_frequency_info = tk.StringVar(value="")
        self.path_setter_var = tk.StringVar(value="")
        self.value_setter_var = tk.StringVar(value="")

        self._configure_style()
        self._build_ui()

        if initial_path is not None:
            self.load_file(initial_path)

    def _configure_style(self) -> None:
        style = ttk.Style(self.root)

        for theme in ("clam", "alt", "default"):
            try:
                style.theme_use(theme)
                break
            except tk.TclError:
                continue

        bg = "#111827"
        panel = "#1f2937"
        text = "#e5e7eb"
        muted = "#9ca3af"
        accent = "#2563eb"
        accent_active = "#1d4ed8"
        input_bg = "#0f172a"

        self.root.configure(bg=bg)

        style.configure("Root.TFrame", background=bg)
        style.configure("Panel.TFrame", background=panel)
        style.configure("HeaderTitle.TLabel", background=bg, foreground="#f9fafb", font=("Segoe UI", 18, "bold"))
        style.configure("HeaderSub.TLabel", background=bg, foreground=muted, font=("Segoe UI", 10))

        style.configure("TLabel", background=panel, foreground=text, font=("Segoe UI", 10))
        style.configure("Muted.TLabel", background=panel, foreground=muted, font=("Segoe UI", 9))
        style.configure("Status.TLabel", background=bg, foreground="#86efac", font=("Segoe UI", 10, "bold"))

        style.configure(
            "TButton",
            font=("Segoe UI", 10, "bold"),
            padding=(12, 8),
            background=accent,
            foreground="#ffffff",
            borderwidth=0,
        )
        style.map(
            "TButton",
            background=[("active", accent_active), ("pressed", accent_active)],
            foreground=[("disabled", "#6b7280")],
        )

        style.configure(
            "TEntry",
            fieldbackground=input_bg,
            foreground=text,
            insertcolor=text,
            borderwidth=1,
            padding=6,
        )
        style.configure(
            "TCombobox",
            fieldbackground=input_bg,
            foreground=text,
            borderwidth=1,
            padding=6,
        )
        style.map(
            "TCombobox",
            fieldbackground=[("readonly", input_bg)],
            foreground=[("readonly", text)],
            selectbackground=[("readonly", input_bg)],
            selectforeground=[("readonly", text)],
        )

        style.configure("Card.TLabelframe", background=panel, borderwidth=0, relief="flat")
        style.configure("Card.TLabelframe.Label", background=panel, foreground="#cbd5e1", font=("Segoe UI", 10, "bold"))

    def _build_ui(self) -> None:
        outer = ttk.Frame(self.root, style="Root.TFrame", padding=18)
        outer.pack(fill=tk.BOTH, expand=True)

        header = ttk.Frame(outer, style="Root.TFrame")
        header.pack(fill=tk.X, pady=(0, 12))
        ttk.Label(header, text="ICARION Config Editor", style="HeaderTitle.TLabel").pack(anchor="w")
        ttk.Label(
            header,
            text="Guided AC/RF editing + full JSON editor for all config settings",
            style="HeaderSub.TLabel",
        ).pack(anchor="w", pady=(2, 0))

        toolbar = ttk.Frame(outer, style="Panel.TFrame", padding=12)
        toolbar.pack(fill=tk.X, pady=(0, 12))
        ttk.Button(toolbar, text="Open JSON", command=self.open_dialog).pack(side=tk.LEFT)
        ttk.Button(toolbar, text="Save", command=self.save_current).pack(side=tk.LEFT, padx=(8, 0))
        ttk.Button(toolbar, text="Save As", command=self.save_as_dialog).pack(side=tk.LEFT, padx=(8, 0))
        ttk.Button(toolbar, text="Reload", command=self.reload_current).pack(side=tk.LEFT, padx=(8, 0))

        self.path_label = ttk.Label(toolbar, text="No file loaded", style="Muted.TLabel")
        self.path_label.pack(side=tk.LEFT, padx=(14, 0), fill=tk.X, expand=True)

        content = ttk.Frame(outer, style="Panel.TFrame", padding=14)
        content.pack(fill=tk.BOTH, expand=True)

        notebook = ttk.Notebook(content)
        notebook.pack(fill=tk.BOTH, expand=True)

        tab_guided = ttk.Frame(notebook, style="Panel.TFrame", padding=12)
        tab_json = ttk.Frame(notebook, style="Panel.TFrame", padding=12)
        notebook.add(tab_guided, text="Guided AC/RF")
        notebook.add(tab_json, text="Full JSON")

        select_card = ttk.LabelFrame(tab_guided, text="Domain Selection", style="Card.TLabelframe", padding=12)
        select_card.pack(fill=tk.X, pady=(0, 12))

        ttk.Label(select_card, text="Domain", style="TLabel").pack(anchor="w", pady=(0, 6))
        self.domain_combo = ttk.Combobox(
            select_card,
            state="readonly",
            textvariable=self.domain_index_var,
            width=80,
            values=[],
        )
        self.domain_combo.pack(fill=tk.X)
        self.domain_combo.bind("<<ComboboxSelected>>", lambda _e: self.populate_fields())

        edit_card = ttk.LabelFrame(tab_guided, text="AC / RF Parameters", style="Card.TLabelframe", padding=12)
        edit_card.pack(fill=tk.X)

        grid = ttk.Frame(edit_card, style="Panel.TFrame")
        grid.pack(fill=tk.X)

        ttk.Label(grid, text="Field").grid(row=0, column=0, sticky="w", padx=6, pady=(2, 8))
        ttk.Label(grid, text="voltage_V").grid(row=0, column=1, sticky="w", padx=6, pady=(2, 8))
        ttk.Label(grid, text="frequency_Hz").grid(row=0, column=2, sticky="w", padx=6, pady=(2, 8))

        ttk.Label(grid, text="RF").grid(row=1, column=0, sticky="w", padx=6, pady=6)
        ttk.Entry(grid, textvariable=self.rf_voltage_var).grid(row=1, column=1, sticky="we", padx=6, pady=6)
        ttk.Entry(grid, textvariable=self.rf_frequency_var).grid(row=1, column=2, sticky="we", padx=6, pady=6)

        ttk.Label(grid, text="AC").grid(row=2, column=0, sticky="w", padx=6, pady=6)
        ttk.Entry(grid, textvariable=self.ac_voltage_var).grid(row=2, column=1, sticky="we", padx=6, pady=6)
        ttk.Entry(grid, textvariable=self.ac_frequency_var).grid(row=2, column=2, sticky="we", padx=6, pady=6)

        ttk.Label(grid, textvariable=self.rf_voltage_info, style="Muted.TLabel").grid(row=3, column=1, sticky="w", padx=6)
        ttk.Label(grid, textvariable=self.rf_frequency_info, style="Muted.TLabel").grid(row=3, column=2, sticky="w", padx=6)
        ttk.Label(grid, textvariable=self.ac_voltage_info, style="Muted.TLabel").grid(row=4, column=1, sticky="w", padx=6)
        ttk.Label(grid, textvariable=self.ac_frequency_info, style="Muted.TLabel").grid(row=4, column=2, sticky="w", padx=6)

        for c in (1, 2):
            grid.columnconfigure(c, weight=1)

        hint = ttk.Label(
            tab_guided,
            style="Muted.TLabel",
            justify=tk.LEFT,
            text=(
                "Empty input = value unchanged. "
                "If a waveform object or @ref is currently set, only the current type is shown. "
                "Entering a number overwrites the target field with a static value."
            ),
            wraplength=900,
        )
        hint.pack(anchor="w", pady=(12, 10))

        ttk.Button(tab_guided, text="Apply to selected domain", command=self.apply_to_domain).pack(anchor="w")

        json_card = ttk.LabelFrame(tab_json, text="Full JSON Editor (all settings)", style="Card.TLabelframe", padding=12)
        json_card.pack(fill=tk.BOTH, expand=True)

        editor_toolbar = ttk.Frame(json_card, style="Panel.TFrame")
        editor_toolbar.pack(fill=tk.X, pady=(0, 8))

        ttk.Button(editor_toolbar, text="Sync from model", command=self.sync_json_text_from_config).pack(side=tk.LEFT)
        ttk.Button(editor_toolbar, text="Apply JSON", command=self.apply_json_text_to_config).pack(side=tk.LEFT, padx=(8, 0))
        ttk.Button(editor_toolbar, text="Validate", command=self.validate_json_text).pack(side=tk.LEFT, padx=(8, 0))
        ttk.Button(editor_toolbar, text="Format", command=self.format_json_text).pack(side=tk.LEFT, padx=(8, 0))

        path_row = ttk.Frame(json_card, style="Panel.TFrame")
        path_row.pack(fill=tk.X, pady=(0, 8))
        ttk.Label(path_row, text="Path", style="Muted.TLabel").pack(side=tk.LEFT)
        ttk.Entry(path_row, textvariable=self.path_setter_var, width=40).pack(side=tk.LEFT, padx=(6, 8), fill=tk.X, expand=True)
        ttk.Label(path_row, text="Value (JSON)", style="Muted.TLabel").pack(side=tk.LEFT)
        ttk.Entry(path_row, textvariable=self.value_setter_var, width=30).pack(side=tk.LEFT, padx=(6, 8), fill=tk.X, expand=True)
        ttk.Button(path_row, text="Set path", command=self.set_path_value).pack(side=tk.LEFT)

        self.json_text = scrolledtext.ScrolledText(
            json_card,
            wrap=tk.NONE,
            undo=True,
            height=16,
            bg="#0b1220",
            fg="#e5e7eb",
            insertbackground="#e5e7eb",
            relief=tk.FLAT,
            borderwidth=0,
            font=("Consolas", 10),
        )
        self.json_text.pack(fill=tk.BOTH, expand=True)

        status_frame = ttk.Frame(outer, style="Root.TFrame", padding=(2, 12, 2, 0))
        status_frame.pack(fill=tk.X)
        ttk.Label(status_frame, textvariable=self.status_var, style="Status.TLabel").pack(anchor="w")

    def open_dialog(self) -> None:
        fp = filedialog.askopenfilename(
            title="Open ICARION config",
            filetypes=[("JSON files", "*.json"), ("All files", "*.*")],
        )
        if fp:
            self.load_file(Path(fp))

    def load_file(self, path: Path) -> None:
        try:
            with path.open("r", encoding="utf-8") as f:
                cfg = json.load(f)
            if not isinstance(cfg, dict):
                raise ValueError("Top-level JSON must be an object")
            domains = cfg.get("domains", [])
            if not isinstance(domains, list):
                raise ValueError("'domains' must be an array")
        except Exception as exc:
            messagebox.showerror("Load error", f"Could not load JSON:\n{exc}")
            self.status_var.set("Load failed")
            return

        self.config = cfg
        self.config_path = path
        self.path_label.configure(text=str(path))
        self._refresh_domain_list()
        self.sync_json_text_from_config()
        self.status_var.set(f"Loaded: {path}")

    def _refresh_domain_list(self) -> None:
        if self.config is None:
            self.domain_combo["values"] = []
            self.domain_index_var.set("")
            return

        domains = self.config.get("domains", [])
        values: list[str] = []
        for i, dom in enumerate(domains):
            name = dom.get("name", f"domain_{i}") if isinstance(dom, dict) else f"domain_{i}"
            inst = dom.get("instrument", "?") if isinstance(dom, dict) else "?"
            values.append(f"{i}: {name} [{inst}]")

        self.domain_combo["values"] = values
        if values:
            self.domain_combo.current(0)
            self.populate_fields()
        else:
            self.domain_index_var.set("")
            self._clear_inputs()

    def _clear_inputs(self) -> None:
        for v in (self.rf_voltage_var, self.rf_frequency_var, self.ac_voltage_var, self.ac_frequency_var):
            v.set("")
        for info in (self.rf_voltage_info, self.rf_frequency_info, self.ac_voltage_info, self.ac_frequency_info):
            info.set("")

    def _selected_domain(self) -> dict[str, Any] | None:
        if self.config is None:
            return None
        sel = self.domain_index_var.get().strip()
        if not sel:
            return None
        try:
            idx = int(sel.split(":", 1)[0])
        except Exception:
            return None

        domains = self.config.get("domains", [])
        if not isinstance(domains, list) or idx < 0 or idx >= len(domains):
            return None
        dom = domains[idx]
        if not isinstance(dom, dict):
            return None
        return dom

    @staticmethod
    def _fmt_type_info(value: Any) -> str:
        if isinstance(value, (int, float)) and not isinstance(value, bool):
            return "current type: number"
        if isinstance(value, str):
            return "current type: string/ref"
        if isinstance(value, dict):
            wtype = value.get("type", "object")
            return f"current type: waveform/object ({wtype})"
        if value is None:
            return "current type: missing"
        return f"current type: {type(value).__name__}"

    def populate_fields(self) -> None:
        dom = self._selected_domain()
        self._clear_inputs()
        if dom is None:
            return

        fields = dom.get("fields", {}) if isinstance(dom.get("fields", {}), dict) else {}
        rf = fields.get("RF", {}) if isinstance(fields.get("RF", {}), dict) else {}
        ac = fields.get("AC", {}) if isinstance(fields.get("AC", {}), dict) else {}

        rf_v = rf.get("voltage_V")
        rf_f = rf.get("frequency_Hz")
        ac_v = ac.get("voltage_V")
        ac_f = ac.get("frequency_Hz")

        if isinstance(rf_v, (int, float)) and not isinstance(rf_v, bool):
            self.rf_voltage_var.set(str(rf_v))
        if isinstance(rf_f, (int, float)) and not isinstance(rf_f, bool):
            self.rf_frequency_var.set(str(rf_f))
        if isinstance(ac_v, (int, float)) and not isinstance(ac_v, bool):
            self.ac_voltage_var.set(str(ac_v))
        if isinstance(ac_f, (int, float)) and not isinstance(ac_f, bool):
            self.ac_frequency_var.set(str(ac_f))

        self.rf_voltage_info.set(self._fmt_type_info(rf_v))
        self.rf_frequency_info.set(self._fmt_type_info(rf_f))
        self.ac_voltage_info.set(self._fmt_type_info(ac_v))
        self.ac_frequency_info.set(self._fmt_type_info(ac_f))

    @staticmethod
    def _parse_optional_float(raw: str) -> float | None:
        txt = raw.strip()
        if txt == "":
            return None
        return float(txt)

    def apply_to_domain(self) -> None:
        dom = self._selected_domain()
        if dom is None:
            messagebox.showwarning("No domain", "Please select a domain first.")
            return

        try:
            rf_v = self._parse_optional_float(self.rf_voltage_var.get())
            rf_f = self._parse_optional_float(self.rf_frequency_var.get())
            ac_v = self._parse_optional_float(self.ac_voltage_var.get())
            ac_f = self._parse_optional_float(self.ac_frequency_var.get())
        except ValueError:
            messagebox.showerror("Invalid number", "Please enter valid numeric values (or leave empty).")
            return

        fields = dom.setdefault("fields", {})
        if not isinstance(fields, dict):
            messagebox.showerror("Invalid JSON", "Selected domain has non-object 'fields'.")
            return

        rf = fields.setdefault("RF", {})
        ac = fields.setdefault("AC", {})
        if not isinstance(rf, dict) or not isinstance(ac, dict):
            messagebox.showerror("Invalid JSON", "'RF'/'AC' is not an object.")
            return

        updates = 0
        if rf_v is not None:
            rf["voltage_V"] = rf_v
            updates += 1
        if rf_f is not None:
            rf["frequency_Hz"] = rf_f
            updates += 1
        if ac_v is not None:
            ac["voltage_V"] = ac_v
            updates += 1
        if ac_f is not None:
            ac["frequency_Hz"] = ac_f
            updates += 1

        self.populate_fields()
        self.sync_json_text_from_config()
        self.status_var.set(f"Applied {updates} value(s) to selected domain")

    def sync_json_text_from_config(self) -> None:
        if self.config is None:
            return
        text = json.dumps(self.config, indent=2, ensure_ascii=False) + "\n"
        self.json_text.delete("1.0", tk.END)
        self.json_text.insert("1.0", text)

    def _apply_json_text_to_config(self, show_message: bool = True) -> bool:
        if self.config is None:
            if show_message:
                messagebox.showwarning("No file", "No JSON loaded.")
            return False

        raw = self.json_text.get("1.0", tk.END).strip()
        if raw == "":
            if show_message:
                messagebox.showerror("Invalid JSON", "Editor is empty.")
            return False

        try:
            parsed = json.loads(raw)
        except json.JSONDecodeError as exc:
            if show_message:
                messagebox.showerror("Invalid JSON", f"JSON parse error:\n{exc}")
            return False

        if not isinstance(parsed, dict):
            if show_message:
                messagebox.showerror("Invalid JSON", "Top-level JSON must be an object.")
            return False

        domains = parsed.get("domains")
        if domains is not None and not isinstance(domains, list):
            if show_message:
                messagebox.showerror("Invalid JSON", "'domains' must be an array if present.")
            return False

        self.config = parsed
        self._refresh_domain_list()
        if show_message:
            self.status_var.set("Applied JSON editor changes to model")
        return True

    def apply_json_text_to_config(self) -> None:
        self._apply_json_text_to_config(show_message=True)

    def validate_json_text(self) -> None:
        ok = self._apply_json_text_to_config(show_message=False)
        if ok:
            self.status_var.set("JSON valid")
        else:
            messagebox.showerror("Validation failed", "JSON is invalid. Fix syntax/structure and try again.")

    def format_json_text(self) -> None:
        raw = self.json_text.get("1.0", tk.END).strip()
        if raw == "":
            return
        try:
            parsed = json.loads(raw)
        except json.JSONDecodeError as exc:
            messagebox.showerror("Format error", f"Cannot format invalid JSON:\n{exc}")
            return
        pretty = json.dumps(parsed, indent=2, ensure_ascii=False) + "\n"
        self.json_text.delete("1.0", tk.END)
        self.json_text.insert("1.0", pretty)
        self.status_var.set("JSON formatted")

    @staticmethod
    def _parse_path_tokens(path: str) -> list[str | int]:
        tokens: list[str | int] = []
        for part in path.split("."):
            p = part.strip()
            if p == "":
                continue
            if p.isdigit():
                tokens.append(int(p))
            else:
                tokens.append(p)
        return tokens

    @staticmethod
    def _set_by_tokens(root_obj: Any, tokens: list[str | int], value: Any) -> None:
        if not tokens:
            raise ValueError("Path is empty")

        cur = root_obj
        for i, tok in enumerate(tokens[:-1]):
            nxt = tokens[i + 1]
            if isinstance(tok, str):
                if not isinstance(cur, dict):
                    raise ValueError(f"Expected object at '{tok}'")
                if tok not in cur or cur[tok] is None:
                    cur[tok] = [] if isinstance(nxt, int) else {}
                cur = cur[tok]
            else:
                if not isinstance(cur, list):
                    raise ValueError(f"Expected array at index {tok}")
                if tok < 0:
                    raise ValueError("Negative indices are not supported")
                while len(cur) <= tok:
                    cur.append({} if not isinstance(nxt, int) else [])
                cur = cur[tok]

        last = tokens[-1]
        if isinstance(last, str):
            if not isinstance(cur, dict):
                raise ValueError("Final parent is not an object")
            cur[last] = value
        else:
            if not isinstance(cur, list):
                raise ValueError("Final parent is not an array")
            if last < 0:
                raise ValueError("Negative indices are not supported")
            while len(cur) <= last:
                cur.append(None)
            cur[last] = value

    def set_path_value(self) -> None:
        if self.config is None:
            messagebox.showwarning("No file", "No JSON loaded.")
            return

        path = self.path_setter_var.get().strip()
        raw_value = self.value_setter_var.get().strip()
        if path == "":
            messagebox.showerror("Missing path", "Please provide a path, e.g. domains.0.env.pressure_Pa")
            return
        if raw_value == "":
            messagebox.showerror("Missing value", "Please provide a JSON value, e.g. 0.5 or \"N2\"")
            return

        try:
            value = json.loads(raw_value)
        except json.JSONDecodeError:
            messagebox.showerror("Invalid value", "Value must be valid JSON (number/string/object/array/bool/null).")
            return

        try:
            tokens = self._parse_path_tokens(path)
            self._set_by_tokens(self.config, tokens, value)
        except ValueError as exc:
            messagebox.showerror("Path error", str(exc))
            return

        self._refresh_domain_list()
        self.sync_json_text_from_config()
        self.status_var.set(f"Set {path}")

    def save_current(self) -> None:
        if self.config is None:
            messagebox.showwarning("No file", "No JSON loaded.")
            return
        if not self._apply_json_text_to_config(show_message=False):
            messagebox.showerror("Save blocked", "JSON editor contains invalid content. Fix it or click 'Sync from model'.")
            return
        if self.config_path is None:
            self.save_as_dialog()
            return
        self._write_json(self.config_path)

    def save_as_dialog(self) -> None:
        if self.config is None:
            messagebox.showwarning("No file", "No JSON loaded.")
            return
        if not self._apply_json_text_to_config(show_message=False):
            messagebox.showerror("Save blocked", "JSON editor contains invalid content. Fix it or click 'Sync from model'.")
            return
        fp = filedialog.asksaveasfilename(
            title="Save ICARION config as",
            defaultextension=".json",
            filetypes=[("JSON files", "*.json"), ("All files", "*.*")],
        )
        if fp:
            self.config_path = Path(fp)
            self.path_label.configure(text=str(self.config_path))
            self._write_json(self.config_path)

    def _write_json(self, path: Path) -> None:
        assert self.config is not None
        try:
            with path.open("w", encoding="utf-8") as f:
                json.dump(self.config, f, indent=2, ensure_ascii=False)
                f.write("\n")
            self.status_var.set(f"Saved: {path}")
        except Exception as exc:
            messagebox.showerror("Save error", f"Could not save JSON:\n{exc}")
            self.status_var.set("Save failed")

    def reload_current(self) -> None:
        if self.config_path is None:
            messagebox.showwarning("No file", "No JSON loaded.")
            return
        self.load_file(self.config_path)


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="ICARION GUI editor for AC/RF config values")
    p.add_argument("--config", type=Path, help="Optional path to JSON config to open")
    return p.parse_args()


def main() -> int:
    args = parse_args()

    root = tk.Tk()
    app = ConfigGui(root, args.config)
    _ = app
    root.mainloop()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
