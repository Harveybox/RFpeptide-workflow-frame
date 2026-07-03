#!/usr/bin/env python3
import csv
import heapq
import json
import math
import shlex
import subprocess
import sys
import threading
import traceback
from copy import deepcopy
from datetime import datetime
from pathlib import Path
from tkinter import (
    BooleanVar,
    Button,
    Canvas,
    Checkbutton,
    END,
    Entry,
    Frame,
    Label,
    LabelFrame,
    Listbox,
    StringVar,
    Tk,
    Text,
    Toplevel,
    filedialog,
    messagebox,
    ttk,
)

SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

import generate_workflow
import cluster_ops
import credential_store
import pdb_preprocess


FIELD_GROUPS = [
    (
        "Project",
        [
            ("target", "Target name", str),
            ("pilot", "Pilot name", str),
            ("remarks", "Remarks / pilot notes", str),
            ("cluster_user", "Cluster user", str),
            ("scratch_date", "Scratch date", str),
            ("project_dir_name", "Scratch project dir", str),
            ("home_project_dir", "Cluster HOME project dir", str),
        ],
    ),
    (
        "RFDiffusion",
        [
            ("input_pdb", "Input PDB path", str),
            ("target_chain", "Target chain", str),
            ("binder_chain", "Binder chain", str),
            ("contigs", "RFdiffusion contigs", str),
            ("rfdiffusion.hotspot_res", "Hotspot residues (for example A326 or A67,A141)", str),
            ("binder_length", "Binder length note", str),
            ("rfdiffusion.queue", "GPU queue", str),
            ("rfdiffusion.gpu_req", "LSF GPU request", str),
            ("rfdiffusion.gpu_ncpu", "GPU job CPU cores", int),
            ("rfdiffusion.gpu_span", "GPU span request", str),
            ("rfdiffusion.env_name", "Environment", str),
            ("rfdiffusion.rfdiffusion_dir", "RFdiffusion dir", str),
            ("rfdiffusion.n_shards", "Shards", int),
            ("rfdiffusion.designs_per_shard", "Designs per shard", int),
            ("rfdiffusion.diffuser_T", "Diffuser T", int),
        ],
    ),
    (
        "ProteinMPNN",
        [
            ("proteinmpnn.queue", "CPU queue", str),
            ("proteinmpnn.env_name", "Environment", str),
            ("proteinmpnn.script", "MPNN script", str),
            ("proteinmpnn.n_shards", "Shards", int),
            ("proteinmpnn.ncpu", "LSF CPU cores (#BSUB -n)", int),
            ("proteinmpnn.resource_req", "LSF resource request (#BSUB -R, optional)", str),
            ("proteinmpnn.seqs_per_struct", "Seqs per structure", int),
            ("proteinmpnn.relax_cycles", "Relax cycles", int),
        ],
    ),
    (
        "AfCycDesign",
        [
            ("afcyc.gpu_queue", "GPU queue", str),
            ("afcyc.cpu_queue", "CPU queue", str),
            ("afcyc.gpu_req", "LSF GPU request", str),
            ("afcyc.gpu_ncpu", "GPU job CPU cores", int),
            ("afcyc.gpu_span", "GPU span request", str),
            ("afcyc.cpu_ncpu", "RMSD CPU cores", int),
            ("afcyc.cpu_span", "RMSD span request", str),
            ("afcyc.env_name", "Environment", str),
            ("afcyc.afcyc_script", "AfCyc script", str),
            ("afcyc.local_afcyc_script", "Local bundled AfCyc script", str),
            ("afcyc.rmsd_script", "RMSD script", str),
            ("afcyc.local_rmsd_script", "Local bundled RMSD script", str),
            ("afcyc.merge_script", "Merge script", str),
            ("afcyc.local_merge_script", "Local bundled merge script", str),
            ("afcyc.n_shards", "Shards", int),
            ("afcyc.num_recycles", "Recycles", int),
            ("afcyc.num_models", "Models", int),
            ("afcyc.af_params_dir", "AF params dir", str),
        ],
    ),
    (
        "PyRosetta",
        [
            ("pyrosetta.queue", "CPU queue", str),
            ("pyrosetta.env_name", "Environment", str),
            ("pyrosetta.script", "Scoring script", str),
            ("pyrosetta.local_script", "Local bundled scoring script", str),
            ("pyrosetta.merge_script", "Merge script", str),
            ("pyrosetta.local_merge_script", "Local bundled merge script", str),
            ("pyrosetta.n_shards", "Shards", int),
            ("pyrosetta.ncpu", "CPU cores", int),
            ("pyrosetta.resource_req", "LSF resource request (#BSUB -R, optional)", str),
            ("pyrosetta.pack_input", "Pack input", bool),
            ("pyrosetta.pack_separated", "Pack separated", bool),
            ("pyrosetta.packstat", "Packstat", bool),
        ],
    ),
]

FIELD_DEFAULTS = {
    "proteinmpnn.ncpu": 1,
    "proteinmpnn.resource_req": "",
    "proteinmpnn.relax_cycles": 4,
    "pyrosetta.resource_req": "",
}

AFCYC_SCORE_PRESETS = [
    ("i_pae_reported", "lower", "10"),
    ("i_pae_norm_31", "lower", "10"),
    ("i_pae", "lower", "10"),
    ("binder_ca_rmsd_target_align", "lower", "10"),
    ("plddt_binder", "higher", "90"),
]

PYROSETTA_SCORE_PRESETS = [
    ("interface_dG", "lower", "10"),
    ("separated_interface_energy", "lower", "10"),
    ("cms", "higher", "90"),
    ("interface_delta_sasa", "higher", "90"),
    ("sap_bound", "lower", "10"),
    ("sap_binder", "lower", "10"),
    ("interface_delta_hbond_unsat", "lower", "10"),
    ("interface_packstat", "higher", "90"),
]

HIT_SCREEN_DEFAULTS = {
    "afcyc_ipae_column": "i_pae_reported",
    "afcyc_ipae_max": "0.3",
    "afcyc_rmsd_column": "binder_ca_rmsd_target_align",
    "afcyc_rmsd_max": "1.5",
    "afcyc_plddt_column": "plddt_binder",
    "afcyc_plddt_min": "",
    "pyro_dg_column": "interface_dG",
    "pyro_dg_max": "-30",
    "pyro_sap_column": "sap_bound",
    "pyro_sap_max": "35",
    "pyro_cms_column": "cms",
    "pyro_cms_min": "300",
    "pyro_sasa_column": "interface_delta_sasa",
    "pyro_sasa_min": "",
    "pyro_packstat_column": "interface_packstat",
    "pyro_packstat_min": "",
}


HIT_FILTER_OUTPUT_COLUMNS = [
    "i_pae",
    "rmsd",
    "plddt",
    "interface_dG",
    "sap_bound",
    "cms",
    "interface_delta_sasa",
    "interface_packstat",
]


PYROSETTA_HIT_OUTPUT_COLUMNS = [
    "interface_dG",
    "separated_interface_energy",
    "interface_delta_sasa",
    "interface_packstat",
    "interface_delta_hbond_unsat",
    "num_interface_residues",
    "gly_interface_energy",
    "sap_bound",
    "sap_free",
    "sap_binder",
    "cms",
    "total_score_after_iam",
    "swi_score",
    "cyclization_seq_feasibility",
    "gravy",
    "aromaticity",
    "isoelectric_point",
    "instability_index",
    "charge_pH7",
    "frac_charged",
    "frac_aromatic",
    "n_term_bulky",
    "c_term_bulky",
    "junction_turn_count",
]


def get_nested(data: dict, dotted_key: str):
    current = data
    for part in dotted_key.split("."):
        if not isinstance(current, dict):
            return ""
        current = current.get(part, "")
    return current


def set_nested(data: dict, dotted_key: str, value):
    current = data
    parts = dotted_key.split(".")
    for part in parts[:-1]:
        current = current.setdefault(part, {})
    current[parts[-1]] = value


class WorkflowGui:
    def __init__(self, root: Tk):
        self.root = root
        self.root.title("RFpeptide Workflow")
        self.root.geometry("1220x860")
        self.config_path = generate_workflow.DEFAULT_CONFIG
        self.output_dir = SCRIPT_DIR / "generated_examples"
        self.config_data = generate_workflow.load_config(self.config_path)
        self.cluster_profile_path = cluster_ops.DEFAULT_PROFILE
        self.cluster_profile = self._load_initial_cluster_profile()
        self.variables = {}
        self.result_entries = []
        self.result_sort_column = "mtime"
        self.result_sort_reverse = True
        self.debug_log_path = SCRIPT_DIR / "logs" / "workflow_gui_debug.log"
        self.status_var = StringVar(value="Idle")
        self.log_lines = []
        self.raw_output_blocks = []
        self.log_text = None
        self.raw_output_text = None
        self.log_window = None
        self.raw_output_window = None

        self._build()
        self._load_values(self.config_data)
        self._log(f"Loaded config: {self.config_path}")

    def _configure_theme(self):
        self.colors = {
            "app_bg": "#eef3f8",
            "card_bg": "#ffffff",
            "nav_bg": "#102a43",
            "nav_hover": "#173c5e",
            "nav_active": "#1f6feb",
            "nav_text": "#dbeafe",
            "header_bg": "#0f2742",
            "header_text": "#f8fafc",
            "muted": "#64748b",
            "border": "#d7e0ea",
            "primary": "#1f6feb",
            "primary_dark": "#174ea6",
            "text": "#1e293b",
            "soft_blue": "#dbeafe",
            "soft_panel": "#f8fafc",
        }
        self.root.configure(bg=self.colors["app_bg"])
        self.root.option_add("*Font", "{Segoe UI} 9")
        self.root.option_add("*Button.Background", "#f8fafc")
        self.root.option_add("*Button.Foreground", self.colors["text"])
        self.root.option_add("*Button.ActiveBackground", "#e2e8f0")
        self.root.option_add("*Button.Relief", "flat")
        self.root.option_add("*Button.BorderWidth", 1)
        self.root.option_add("*Text.Background", "#fbfdff")
        self.root.option_add("*Text.Foreground", self.colors["text"])
        self.root.option_add("*Text.Relief", "solid")
        self.root.option_add("*Text.BorderWidth", 1)
        self.root.option_add("*Entry.Background", "#ffffff")
        self.root.option_add("*Entry.Relief", "solid")
        self.root.option_add("*Entry.BorderWidth", 1)
        self.root.option_add("*Label.Background", self.colors["app_bg"])
        self.root.option_add("*Frame.Background", self.colors["app_bg"])
        self.root.option_add("*Labelframe.Background", self.colors["card_bg"])
        self.root.option_add("*Labelframe.Foreground", self.colors["text"])
        style = ttk.Style()
        try:
            style.theme_use("clam")
        except Exception:
            pass
        style.configure("Treeview", rowheight=24, background="#ffffff", fieldbackground="#ffffff", foreground=self.colors["text"])
        style.configure("Treeview.Heading", font="{Segoe UI} 9 bold", background="#e8eef6", foreground=self.colors["text"])
        style.map("Treeview", background=[("selected", "#cfe3ff")], foreground=[("selected", "#0f172a")])
        style.configure("TCombobox", fieldbackground="#ffffff", background="#ffffff")
        style.configure("Vertical.TScrollbar", background="#d7e0ea")
        style.configure("Horizontal.TScrollbar", background="#d7e0ea")

    def _build(self):
        self._configure_theme()

        self.config_path_var = StringVar(value=str(self.config_path))
        self.output_dir_var = StringVar(value=str(self.output_dir))
        self.top_local_pdb_var = StringVar()
        self.current_project_var = StringVar(value="Current project: not loaded")

        header = Frame(self.root, bg=self.colors["header_bg"])
        header.pack(fill="x")
        Label(
            header,
            text="RFpeptide Workflow",
            bg=self.colors["header_bg"],
            fg=self.colors["header_text"],
            font="{Segoe UI} 15 bold",
        ).pack(side="left", padx=16, pady=(10, 2))
        Label(
            header,
            textvariable=self.current_project_var,
            bg=self.colors["header_bg"],
            fg="#bfdbfe",
            anchor="e",
        ).pack(side="right", fill="x", expand=True, padx=16, pady=(12, 2))

        top = Frame(self.root, bg=self.colors["card_bg"], bd=1, relief="solid")
        top.pack(fill="x", padx=12, pady=(10, 8))
        Label(top, text="Config", bg=self.colors["card_bg"], fg=self.colors["muted"]).grid(row=0, column=0, sticky="w", padx=(10, 4), pady=8)
        Entry(top, textvariable=self.config_path_var).grid(row=0, column=1, sticky="ew", padx=6, pady=8)
        Button(top, text="Open", command=self.open_config).grid(row=0, column=2, padx=3, pady=8)
        Button(top, text="Save", command=self.save_config).grid(row=0, column=3, padx=3, pady=8)
        Button(top, text="Save as", command=self.save_config_as).grid(row=0, column=4, padx=(3, 10), pady=8)
        Label(top, text="Local PDB", bg=self.colors["card_bg"], fg=self.colors["muted"]).grid(row=1, column=0, sticky="w", padx=(10, 4), pady=(0, 8))
        Entry(top, textvariable=self.top_local_pdb_var).grid(row=1, column=1, sticky="ew", padx=6, pady=(0, 8))
        Button(top, text="Choose PDB", command=self.choose_local_input_pdb).grid(row=1, column=2, padx=3, pady=(0, 8))
        Button(top, text="Preview/Clean", command=self.preview_top_local_pdb).grid(row=1, column=3, padx=3, pady=(0, 8))
        Label(top, text="Output", bg=self.colors["card_bg"], fg=self.colors["muted"]).grid(row=2, column=0, sticky="w", padx=(10, 4), pady=(0, 8))
        Entry(top, textvariable=self.output_dir_var).grid(row=2, column=1, sticky="ew", padx=6, pady=(0, 8))
        Button(top, text="Choose", command=self.choose_output_dir).grid(row=2, column=2, padx=3, pady=(0, 8))
        generate_button = Button(top, text="Generate workflow", command=self.generate, bg=self.colors["primary"], fg="#ffffff", activebackground=self.colors["primary_dark"])
        generate_button.grid(row=2, column=3, columnspan=2, sticky="ew", padx=(3, 10), pady=(0, 8))
        top.columnconfigure(1, weight=1)

        shell = Frame(self.root, bg=self.colors["app_bg"])
        shell.pack(fill="both", expand=True, padx=12, pady=(0, 10))
        shell.columnconfigure(1, weight=1)
        shell.rowconfigure(0, weight=1)

        self.nav_frame = Frame(shell, bg=self.colors["nav_bg"], width=210)
        self.nav_frame.grid(row=0, column=0, sticky="nsw")
        self.nav_frame.grid_propagate(False)
        self.page_container = Frame(shell, bg=self.colors["app_bg"])
        self.page_container.grid(row=0, column=1, sticky="nsew", padx=(10, 0))
        self.page_container.rowconfigure(0, weight=1)
        self.page_container.columnconfigure(0, weight=1)
        self.page_frames = {}
        self.nav_buttons = {}

        self.pdb_preprocess_frame = self._create_page("PDB Preprocess", "PDB Preprocess", "Fetch, inspect, select chains, and clean input structures.", self._build_pdb_preprocess, scrollable=True)
        for group_name, fields in FIELD_GROUPS:
            self._create_page(group_name, group_name, "Edit target-specific workflow parameters.", lambda parent, name=group_name, group_fields=fields: self._build_group(parent, name, group_fields), scrollable=True)
        for traced_key in ("target", "pilot", "scratch_date", "project_dir_name"):
            if traced_key in self.variables:
                self.variables[traced_key][0].trace_add("write", lambda *_: self._update_current_project_banner())
        cluster_settings_frame = self._create_page("Cluster Settings", "Cluster Settings", "Connection, authentication, scan, and local tool settings.", self._build_cluster_settings)
        self.cluster_settings_frame = cluster_settings_frame
        self._create_page("Cluster Dashboard", "Cluster Dashboard", "Submit, monitor, inspect, and stop cluster jobs.", self._build_cluster_dashboard, scrollable=True)
        self._create_page("Results Browser", "Results Browser", "Scan scratch, filter outputs, download results, and preview files.", self._build_results_browser, scrollable=True)
        self._create_page("Data Processing", "Data Processing", "Preview merged CSV files, plot metrics, and shortlist candidates.", self._build_data_processing, scrollable=True)
        self._create_page("Scratch Safety", "Scratch Safety", "Move scratch date folders safely to avoid cluster cleanup.", self._build_scratch_safety, scrollable=True)

        self._build_navigation([
            ("Run", ["Cluster Dashboard", "Results Browser", "Data Processing"]),
            ("Configure", ["Project", "RFDiffusion", "ProteinMPNN", "AfCycDesign", "PyRosetta", "PDB Preprocess", "Cluster Settings"]),
            ("Maintenance", ["Scratch Safety"]),
        ])
        self._show_page("Cluster Dashboard")

        self._build_output_bar()

    def _create_page(self, key: str, title: str, subtitle: str, builder, scrollable: bool = False):
        page = Frame(self.page_container, bg=self.colors["app_bg"])
        page.grid(row=0, column=0, sticky="nsew")
        page.columnconfigure(0, weight=1)
        page.rowconfigure(1, weight=1)
        header = Frame(page, bg=self.colors["app_bg"])
        header.grid(row=0, column=0, sticky="ew", pady=(0, 8))
        Label(header, text=title, bg=self.colors["app_bg"], fg=self.colors["text"], font="{Segoe UI} 14 bold").pack(anchor="w")
        Label(header, text=subtitle, bg=self.colors["app_bg"], fg=self.colors["muted"]).pack(anchor="w", pady=(2, 0))
        body = Frame(page, bg=self.colors["app_bg"])
        body.grid(row=1, column=0, sticky="nsew")
        self.page_frames[key] = page
        builder(self._scrollable_page_body(body) if scrollable else body)
        return page

    def _scrollable_page_body(self, parent: Frame) -> Frame:
        outer = Frame(parent, bg=self.colors["app_bg"])
        outer.pack(fill="both", expand=True)
        canvas = Canvas(outer, highlightthickness=0, bg=self.colors["app_bg"])
        y_scrollbar = ttk.Scrollbar(outer, orient="vertical", command=canvas.yview)
        canvas.configure(yscrollcommand=y_scrollbar.set)
        y_scrollbar.pack(side="right", fill="y")
        canvas.pack(side="left", fill="both", expand=True)
        body = Frame(canvas, bg=self.colors["app_bg"])
        body_window = canvas.create_window((0, 0), window=body, anchor="nw")

        def update_scroll_region(_event=None):
            viewport_width = max(canvas.winfo_width(), 1)
            canvas.itemconfigure(body_window, width=viewport_width)
            canvas.configure(scrollregion=(0, 0, viewport_width, body.winfo_reqheight()))

        def fit_body_width(event):
            canvas.itemconfigure(body_window, width=max(event.width, 1))
            canvas.configure(scrollregion=(0, 0, max(event.width, 1), body.winfo_reqheight()))

        def on_mousewheel(event):
            canvas.yview_scroll(int(-1 * (event.delta / 120)), "units")

        def enable_mousewheel(_event=None):
            canvas.bind_all("<MouseWheel>", on_mousewheel)

        def disable_mousewheel(_event=None):
            canvas.unbind_all("<MouseWheel>")

        body.bind("<Configure>", update_scroll_region)
        canvas.bind("<Configure>", fit_body_width)
        body.bind("<Enter>", enable_mousewheel)
        body.bind("<Leave>", disable_mousewheel)
        return body

    def _build_output_bar(self):
        bar = Frame(self.root, bg=self.colors["card_bg"], bd=1, relief="solid")
        bar.pack(fill="x", padx=12, pady=(0, 8))
        Label(bar, text="Status", bg=self.colors["card_bg"], fg=self.colors["muted"]).pack(side="left", padx=(10, 6), pady=6)
        Label(bar, textvariable=self.status_var, bg=self.colors["card_bg"], fg=self.colors["text"], anchor="w").pack(side="left", fill="x", expand=True, padx=(0, 8), pady=6)
        Button(bar, text="Open Log", command=self.open_log_window).pack(side="right", padx=(4, 10), pady=4)
        Button(bar, text="Raw SSH/LSF", command=self.open_raw_output_window).pack(side="right", padx=4, pady=4)
        Button(bar, text="Clear raw", command=self.clear_raw_output).pack(side="right", padx=4, pady=4)
        Button(bar, text="Debug folder", command=self.open_debug_log_folder).pack(side="right", padx=4, pady=4)

    def _open_output_window(self, title: str, content: str, wrap: str, on_close):
        window = Toplevel(self.root)
        window.title(title)
        window.geometry("1050x650")
        frame = Frame(window)
        frame.pack(fill="both", expand=True, padx=8, pady=8)
        text = Text(frame, wrap=wrap)
        y_scroll = ttk.Scrollbar(frame, orient="vertical", command=text.yview)
        x_scroll = ttk.Scrollbar(frame, orient="horizontal", command=text.xview)
        text.configure(yscrollcommand=y_scroll.set, xscrollcommand=x_scroll.set)
        text.grid(row=0, column=0, sticky="nsew")
        y_scroll.grid(row=0, column=1, sticky="ns")
        if wrap == "none":
            x_scroll.grid(row=1, column=0, sticky="ew")
        frame.rowconfigure(0, weight=1)
        frame.columnconfigure(0, weight=1)
        text.insert(END, content)
        text.see(END)

        def close_window():
            on_close()
            window.destroy()

        window.protocol("WM_DELETE_WINDOW", close_window)
        return window, text

    def open_log_window(self):
        if self.log_window is not None and self.log_window.winfo_exists():
            self.log_window.lift()
            return
        content = "\n".join(self.log_lines)
        if content:
            content += "\n"
        self.log_window, self.log_text = self._open_output_window(
            "Workflow Log",
            content,
            "word",
            lambda: self._clear_log_window_refs(),
        )

    def open_raw_output_window(self):
        if self.raw_output_window is not None and self.raw_output_window.winfo_exists():
            self.raw_output_window.lift()
            return
        self.raw_output_window, self.raw_output_text = self._open_output_window(
            "Raw SSH/LSF output",
            "".join(self.raw_output_blocks),
            "none",
            lambda: self._clear_raw_output_window_refs(),
        )

    def _clear_log_window_refs(self):
        self.log_window = None
        self.log_text = None

    def _clear_raw_output_window_refs(self):
        self.raw_output_window = None
        self.raw_output_text = None

    def _build_navigation(self, sections: list[tuple[str, list[str]]]):
        Label(
            self.nav_frame,
            text="WORKFLOW",
            bg=self.colors["nav_bg"],
            fg="#93c5fd",
            font="{Segoe UI} 9 bold",
        ).pack(anchor="w", padx=14, pady=(16, 8))
        for section_title, page_keys in sections:
            Label(
                self.nav_frame,
                text=section_title.upper(),
                bg=self.colors["nav_bg"],
                fg="#7dd3fc",
                font="{Segoe UI} 8 bold",
            ).pack(anchor="w", padx=14, pady=(12, 4))
            for page_key in page_keys:
                button = Button(
                    self.nav_frame,
                    text=page_key,
                    anchor="w",
                    bd=0,
                    relief="flat",
                    bg=self.colors["nav_bg"],
                    fg=self.colors["nav_text"],
                    activebackground=self.colors["nav_hover"],
                    activeforeground="#ffffff",
                    command=lambda key=page_key: self._show_page(key),
                )
                button.pack(fill="x", padx=8, pady=1, ipady=7)
                self.nav_buttons[page_key] = button

    def _show_page(self, key: str):
        frame = self.page_frames.get(key)
        if not frame:
            return
        frame.tkraise()
        for page_key, button in self.nav_buttons.items():
            if page_key == key:
                button.configure(bg=self.colors["nav_active"], fg="#ffffff")
            else:
                button.configure(bg=self.colors["nav_bg"], fg=self.colors["nav_text"])

    def _build_status_bar(self, parent: Frame, variable: StringVar, title: str = "Status") -> Frame:
        frame = Frame(parent, bg="#fff3cd", bd=1, relief="solid")
        badge = Label(frame, text=title, bg="#7a4f00", fg="#ffffff", width=11)
        badge.pack(side="left", fill="y")
        text_label = Label(frame, textvariable=variable, bg="#fff3cd", fg="#1f1f1f", anchor="w")
        text_label.pack(side="left", fill="x", expand=True, padx=8, pady=4)

        def refresh_colors(*_args):
            message = variable.get().lower()
            if any(token in message for token in ("failed", "error", "not downloaded", "失败")):
                badge_bg, body_bg, fg = "#8b1e1e", "#fde2e1", "#3c0b0b"
            elif any(token in message for token in ("downloading", "scanning", "running", "moving", "upload", "submit", "preview")):
                badge_bg, body_bg, fg = "#0b4f71", "#dff3ff", "#0b2f42"
            elif any(token in message for token in ("finished", "complete", "loaded", "found", "succeeded", "download finished", "moved")):
                badge_bg, body_bg, fg = "#2f6b2f", "#e4f6e4", "#123112"
            else:
                badge_bg, body_bg, fg = "#7a4f00", "#fff3cd", "#1f1f1f"
            frame.configure(bg=body_bg)
            badge.configure(bg=badge_bg)
            text_label.configure(bg=body_bg, fg=fg)

        variable.trace_add("write", refresh_colors)
        refresh_colors()
        return frame

    def _build_pdb_preprocess(self, parent: Frame):
        body = Frame(parent)
        body.pack(fill="both", expand=True, padx=12, pady=12)

        self.pdb_rcsb_id_var = StringVar()
        self.pdb_raw_path_var = StringVar()
        self.pdb_original_path_var = StringVar()
        self.pdb_clean_output_var = StringVar()
        self.pdb_keep_protein_only_var = BooleanVar(value=True)
        self.pdb_renumber_atoms_var = BooleanVar(value=True)
        self.pdb_renumber_residues_var = BooleanVar(value=False)
        self.pdb_renumber_chains_var = BooleanVar(value=False)
        self.pdb_set_as_input_var = BooleanVar(value=True)
        self.pdb_preprocess_status_var = StringVar(value="Select a PDB and preview it.")
        self.pdb_chain_ids = []
        self.pdb_chain_info = {}

        fetch_frame = LabelFrame(body, text="Experimental RCSB PDB fetch")
        fetch_frame.pack(fill="x", pady=(0, 8))
        Label(fetch_frame, text="PDB ID").grid(row=0, column=0, sticky="w", padx=6, pady=6)
        Entry(fetch_frame, textvariable=self.pdb_rcsb_id_var, width=12).grid(row=0, column=1, sticky="w", padx=6, pady=6)
        Button(fetch_frame, text="Fetch PDB", command=self.fetch_rcsb_pdb).grid(row=0, column=2, padx=3, pady=6)
        Label(
            fetch_frame,
            text="Downloads the legacy PDB-format structure from RCSB into inputs_preprocessed/rcsb.",
            fg=self.colors["muted"],
        ).grid(row=0, column=3, sticky="w", padx=8, pady=6)
        fetch_frame.columnconfigure(3, weight=1)

        input_frame = LabelFrame(body, text="Raw PDB")
        input_frame.pack(fill="x", pady=(0, 8))
        Label(input_frame, text="PDB file").grid(row=0, column=0, sticky="w", padx=6, pady=6)
        Entry(input_frame, textvariable=self.pdb_raw_path_var).grid(row=0, column=1, sticky="ew", padx=6, pady=6)
        Button(input_frame, text="Browse", command=self.choose_preprocess_pdb).grid(row=0, column=2, padx=3, pady=6)
        Button(input_frame, text="Preview structure", command=self.preview_preprocess_pdb).grid(row=0, column=3, padx=3, pady=6)
        Button(input_frame, text="Use as workflow input", command=self.use_preprocess_pdb_as_input).grid(row=0, column=4, padx=3, pady=6)
        Button(input_frame, text="Restore original PDB", command=self.restore_original_pdb).grid(row=0, column=5, padx=3, pady=6)
        Label(input_frame, text="Original snapshot").grid(row=1, column=0, sticky="w", padx=6, pady=(0, 6))
        Entry(input_frame, textvariable=self.pdb_original_path_var, state="readonly").grid(row=1, column=1, columnspan=5, sticky="ew", padx=6, pady=(0, 6))
        input_frame.columnconfigure(1, weight=1)

        chain_frame = LabelFrame(body, text="Chains to keep")
        chain_frame.pack(fill="x", pady=(0, 8))
        Label(
            chain_frame,
            text="Select one or more chains (Ctrl/Shift for multi-select). Cleanup keeps only selected chains.",
            fg=self.colors["muted"],
        ).pack(anchor="w", padx=6, pady=(6, 2))
        chain_content = Frame(chain_frame)
        chain_content.pack(fill="x", padx=6, pady=(0, 6))
        self.pdb_chain_listbox = Listbox(chain_content, selectmode="extended", exportselection=False, height=6)
        chain_scroll = ttk.Scrollbar(chain_content, orient="vertical", command=self.pdb_chain_listbox.yview)
        self.pdb_chain_listbox.configure(yscrollcommand=chain_scroll.set)
        self.pdb_chain_listbox.pack(side="left", fill="x", expand=True)
        chain_scroll.pack(side="left", fill="y")
        chain_actions = Frame(chain_content)
        chain_actions.pack(side="left", fill="y", padx=(8, 0))
        Button(chain_actions, text="Select all", command=self.select_all_pdb_chains).pack(fill="x", pady=(0, 4))
        Button(chain_actions, text="Protein chains", command=self.select_protein_pdb_chains).pack(fill="x", pady=4)
        Button(chain_actions, text="Clear", command=lambda: self.pdb_chain_listbox.selection_clear(0, END)).pack(fill="x", pady=4)

        options_frame = LabelFrame(body, text="Cleanup options")
        options_frame.pack(fill="x", pady=(0, 8))
        Label(options_frame, text="Cleaned output").grid(row=0, column=0, sticky="w", padx=6, pady=6)
        Entry(options_frame, textvariable=self.pdb_clean_output_var).grid(row=0, column=1, sticky="ew", padx=6, pady=6)
        Button(options_frame, text="Choose output", command=self.choose_clean_pdb_output).grid(row=0, column=2, padx=3, pady=6)
        Checkbutton(options_frame, text="Remove non-protein components", variable=self.pdb_keep_protein_only_var).grid(row=1, column=0, sticky="w", padx=6, pady=(0, 6))
        Checkbutton(options_frame, text="Renumber atom serials", variable=self.pdb_renumber_atoms_var).grid(row=1, column=1, sticky="w", padx=6, pady=(0, 6))
        Checkbutton(options_frame, text="Renumber residues continuously by chain", variable=self.pdb_renumber_residues_var).grid(row=1, column=2, sticky="w", padx=6, pady=(0, 6))
        Checkbutton(options_frame, text="Rename retained chains from A", variable=self.pdb_renumber_chains_var).grid(row=1, column=3, sticky="w", padx=6, pady=(0, 6))
        Checkbutton(options_frame, text="Set cleaned PDB as workflow input", variable=self.pdb_set_as_input_var).grid(row=2, column=0, columnspan=2, sticky="w", padx=6, pady=(0, 6))
        Button(options_frame, text="Apply cleanup / keep chains", command=self.clean_preprocess_pdb).grid(row=2, column=2, sticky="ew", padx=6, pady=(0, 6))
        self._build_status_bar(options_frame, self.pdb_preprocess_status_var, "PDB").grid(row=2, column=3, sticky="ew", padx=6, pady=(0, 6))
        options_frame.columnconfigure(1, weight=1)
        options_frame.columnconfigure(3, weight=1)

        component_frame = LabelFrame(body, text="Structure components")
        component_frame.pack(fill="both", expand=True, pady=(0, 8))
        self.pdb_component_text = Text(component_frame, height=16, wrap="word")
        component_y_scroll = ttk.Scrollbar(component_frame, orient="vertical", command=self.pdb_component_text.yview)
        self.pdb_component_text.configure(yscrollcommand=component_y_scroll.set)
        self.pdb_component_text.grid(row=0, column=0, sticky="nsew")
        component_y_scroll.grid(row=0, column=1, sticky="ns")
        component_frame.rowconfigure(0, weight=1)
        component_frame.columnconfigure(0, weight=1)

        preview_frame = Frame(body)
        preview_frame.pack(fill="both", expand=True, pady=(0, 8))

        summary_frame = LabelFrame(preview_frame, text="Sequence summary")
        summary_frame.pack(fill="both", expand=True, pady=(0, 6))
        self.pdb_summary_text = Text(summary_frame, height=16, wrap="none")
        summary_y_scroll = ttk.Scrollbar(summary_frame, orient="vertical", command=self.pdb_summary_text.yview)
        summary_x_scroll = ttk.Scrollbar(summary_frame, orient="horizontal", command=self.pdb_summary_text.xview)
        self.pdb_summary_text.configure(yscrollcommand=summary_y_scroll.set, xscrollcommand=summary_x_scroll.set)
        self.pdb_summary_text.grid(row=0, column=0, sticky="nsew")
        summary_y_scroll.grid(row=0, column=1, sticky="ns")
        summary_x_scroll.grid(row=1, column=0, sticky="ew")
        summary_frame.rowconfigure(0, weight=1)
        summary_frame.columnconfigure(0, weight=1)

        table_frame = LabelFrame(preview_frame, text="Residues")
        table_frame.pack(fill="both", expand=True)
        columns = ("chain", "resseq", "icode", "resname", "aa", "atoms")
        self.pdb_residue_tree = ttk.Treeview(table_frame, columns=columns, show="headings", height=16)
        residue_y_scroll = ttk.Scrollbar(table_frame, orient="vertical", command=self.pdb_residue_tree.yview)
        residue_x_scroll = ttk.Scrollbar(table_frame, orient="horizontal", command=self.pdb_residue_tree.xview)
        self.pdb_residue_tree.configure(yscrollcommand=residue_y_scroll.set, xscrollcommand=residue_x_scroll.set)
        for column in columns:
            self.pdb_residue_tree.heading(column, text=column)
            self.pdb_residue_tree.column(column, width=85, minwidth=50, anchor="w", stretch=False)
        self.pdb_residue_tree.grid(row=0, column=0, sticky="nsew")
        residue_y_scroll.grid(row=0, column=1, sticky="ns")
        residue_x_scroll.grid(row=1, column=0, sticky="ew")
        table_frame.rowconfigure(0, weight=1)
        table_frame.columnconfigure(0, weight=1)

    def _build_cluster_settings(self, parent: Frame):
        body = self._scrollable_page_body(parent)
        body.configure(padx=12, pady=12)

        self.cluster_profile_path_var = StringVar(value=str(self.cluster_profile_path))
        self.cluster_settings_status_var = StringVar(value="Cluster settings loaded.")
        self.cluster_host_var = StringVar()
        self.cluster_user_var = StringVar()
        self.cluster_port_var = StringVar()
        self.auth_method_var = StringVar()
        self.password_var = StringVar()
        self.remember_password_var = BooleanVar()
        self.cluster_ssh_key_var = StringVar()
        self.cluster_remote_upload_parent_var = StringVar()
        self.cluster_remote_workflow_dir_var = StringVar()
        self.cluster_local_download_dir_var = StringVar()
        self.cluster_pymol_var = StringVar()
        self.cluster_bjobs_user_var = StringVar()
        self.cluster_scan_depth_var = StringVar()
        self.cluster_log_depth_var = StringVar()
        self.cluster_result_patterns_var = StringVar()
        self.cluster_ssh_extra_args_var = StringVar()

        profile_frame = LabelFrame(body, text="Cluster profile file")
        profile_frame.pack(fill="x", pady=(0, 8))
        Label(profile_frame, text="Profile").grid(row=0, column=0, sticky="w", padx=6, pady=6)
        Entry(profile_frame, textvariable=self.cluster_profile_path_var).grid(row=0, column=1, sticky="ew", padx=6, pady=6)
        Button(profile_frame, text="Open", command=self.open_cluster_profile).grid(row=0, column=2, padx=3, pady=6)
        Button(profile_frame, text="Save", command=self.save_cluster_profile).grid(row=0, column=3, padx=3, pady=6)
        Button(profile_frame, text="Save as", command=self.save_cluster_profile_as).grid(row=0, column=4, padx=3, pady=6)
        profile_frame.columnconfigure(1, weight=1)

        connection_frame = LabelFrame(body, text="Connection")
        connection_frame.pack(fill="x", pady=(0, 8))
        for label, variable, row, column in [
            ("Host", self.cluster_host_var, 0, 0),
            ("User", self.cluster_user_var, 0, 2),
            ("Port", self.cluster_port_var, 0, 4),
        ]:
            Label(connection_frame, text=label).grid(row=row, column=column, sticky="w", padx=6, pady=6)
            Entry(connection_frame, textvariable=variable).grid(row=row, column=column + 1, sticky="ew", padx=6, pady=6)
        Label(connection_frame, text="SSH key").grid(row=1, column=0, sticky="w", padx=6, pady=6)
        Entry(connection_frame, textvariable=self.cluster_ssh_key_var).grid(row=1, column=1, columnspan=3, sticky="ew", padx=6, pady=6)
        Button(connection_frame, text="Choose key", command=self.choose_cluster_ssh_key).grid(row=1, column=4, padx=3, pady=6)
        connection_frame.columnconfigure(1, weight=1)
        connection_frame.columnconfigure(3, weight=1)

        auth_frame = LabelFrame(body, text="Authentication")
        auth_frame.pack(fill="x", pady=(0, 8))
        Label(auth_frame, text="Method").grid(row=0, column=0, sticky="w", padx=6, pady=6)
        ttk.Combobox(auth_frame, textvariable=self.auth_method_var, values=["ssh_key", "password"], state="readonly", width=14).grid(row=0, column=1, sticky="w", padx=6, pady=6)
        Label(auth_frame, text="Password").grid(row=0, column=2, sticky="w", padx=6, pady=6)
        Entry(auth_frame, textvariable=self.password_var, show="*").grid(row=0, column=3, sticky="ew", padx=6, pady=6)
        Checkbutton(auth_frame, text="Remember securely", variable=self.remember_password_var).grid(row=0, column=4, sticky="w", padx=6, pady=6)
        Button(auth_frame, text="Load saved", command=self.load_saved_password).grid(row=0, column=5, padx=3, pady=6)
        Button(auth_frame, text="Forget saved", command=self.forget_saved_password).grid(row=0, column=6, padx=3, pady=6)
        auth_frame.columnconfigure(3, weight=1)

        paths_frame = LabelFrame(body, text="Workflow paths and local tools")
        paths_frame.pack(fill="x", pady=(0, 8))
        for row, (label, variable) in enumerate([
            ("Remote upload parent", self.cluster_remote_upload_parent_var),
            ("Remote workflow dir fallback", self.cluster_remote_workflow_dir_var),
            ("Local download dir fallback", self.cluster_local_download_dir_var),
            ("PyMOL executable", self.cluster_pymol_var),
            ("bjobs user", self.cluster_bjobs_user_var),
        ]):
            Label(paths_frame, text=label).grid(row=row, column=0, sticky="w", padx=6, pady=5)
            Entry(paths_frame, textvariable=variable).grid(row=row, column=1, sticky="ew", padx=6, pady=5)
            if label == "PyMOL executable":
                Button(paths_frame, text="Choose", command=self.choose_pymol_executable).grid(row=row, column=2, padx=3, pady=5)
        paths_frame.columnconfigure(1, weight=1)

        scan_frame = LabelFrame(body, text="Scan and SSH options")
        scan_frame.pack(fill="x", pady=(0, 8))
        Label(scan_frame, text="Result scan max depth").grid(row=0, column=0, sticky="w", padx=6, pady=6)
        Entry(scan_frame, textvariable=self.cluster_scan_depth_var, width=8).grid(row=0, column=1, sticky="w", padx=6, pady=6)
        Label(scan_frame, text="Job log scan max depth").grid(row=0, column=2, sticky="w", padx=6, pady=6)
        Entry(scan_frame, textvariable=self.cluster_log_depth_var, width=8).grid(row=0, column=3, sticky="w", padx=6, pady=6)
        Label(scan_frame, text="Result file patterns").grid(row=1, column=0, sticky="w", padx=6, pady=6)
        Entry(scan_frame, textvariable=self.cluster_result_patterns_var).grid(row=1, column=1, columnspan=3, sticky="ew", padx=6, pady=6)
        Label(scan_frame, text="SSH extra args").grid(row=2, column=0, sticky="w", padx=6, pady=6)
        Entry(scan_frame, textvariable=self.cluster_ssh_extra_args_var).grid(row=2, column=1, columnspan=3, sticky="ew", padx=6, pady=6)
        scan_frame.columnconfigure(1, weight=1)
        scan_frame.columnconfigure(3, weight=1)

        hint = LabelFrame(body, text="Notes")
        hint.pack(fill="both", expand=True, pady=(0, 8))
        text = Text(hint, height=7, wrap="word")
        text.pack(fill="both", expand=True, padx=6, pady=6)
        text.insert(
            END,
            "Passwords are never written to cluster_profile.local.json. If Remember securely is enabled, the password is stored separately in the local encrypted credential store.\n"
            "Remote upload parent is the recommended setting. The concrete remote_workflow_dir is derived from the active target config before upload/submit.\n"
            "Result file patterns accept comma or whitespace separated patterns, for example: *.csv, *.pdb, *.out, *.err, *.json, *.txt.\n",
        )
        text.configure(state="disabled")

        self._build_status_bar(body, self.cluster_settings_status_var, "Settings").pack(fill="x", pady=(0, 8))
        self._load_cluster_settings_vars()

    def _build_cluster_dashboard(self, parent: Frame):
        body = Frame(parent)
        body.pack(fill="both", expand=True, padx=12, pady=12)

        self.stage_var = StringVar(value="RFDiffusion")
        self.auto_refresh_var = BooleanVar(value=False)
        self.job_entries = []
        self.job_entries_by_iid = {}
        self.job_filter_var = StringVar()
        self.job_status_filter_var = StringVar(value="All")
        self.job_summary_var = StringVar(value="Jobs not loaded")
        self.job_selection_var = StringVar(value="Selected: 0")
        self.job_sort_column = "jobid"
        self.job_sort_reverse = False
        self.queue_name_var = StringVar()
        self.queue_status_var = StringVar(value="Queue status not queried")

        settings_link = LabelFrame(body, text="Cluster settings")
        settings_link.pack(fill="x", pady=(0, 8))
        Label(settings_link, text="Active profile").grid(row=0, column=0, sticky="w", padx=6, pady=6)
        Label(settings_link, textvariable=self.cluster_profile_path_var, anchor="w").grid(row=0, column=1, sticky="ew", padx=6, pady=6)
        Button(settings_link, text="Open settings", command=self.show_cluster_settings).grid(row=0, column=2, padx=3, pady=6)
        settings_link.columnconfigure(1, weight=1)

        summary_frame = LabelFrame(body, text="Current target parameters")
        summary_frame.pack(fill="x", pady=(0, 8))
        self.parameter_summary = Text(summary_frame, height=8, wrap="word")
        self.parameter_summary.pack(fill="x", padx=6, pady=6)

        actions = LabelFrame(body, text="Actions")
        actions.pack(fill="x", pady=(0, 8))
        Button(actions, text="Test connection", command=self.test_cluster_connection).grid(row=0, column=0, padx=4, pady=6, sticky="ew")
        Button(actions, text="Upload workflow", command=self.upload_workflow).grid(row=0, column=1, padx=4, pady=6, sticky="ew")
        self.stage_box = ttk.Combobox(actions, textvariable=self.stage_var, values=list(self.cluster_profile.get("stage_scripts", {}).keys()), state="readonly", width=18)
        self.stage_box.grid(row=0, column=2, padx=4, pady=6, sticky="ew")
        Button(actions, text="Sync + submit stage", command=self.submit_selected_stage).grid(row=0, column=3, padx=4, pady=6, sticky="ew")
        self._build_status_bar(actions, self.status_var, "Cluster").grid(row=1, column=0, columnspan=4, sticky="ew", padx=6, pady=(0, 6))
        for column in range(4):
            actions.columnconfigure(column, weight=1)

        queue_frame = LabelFrame(body, text="Queue status")
        queue_frame.pack(fill="x", pady=(0, 8))
        Button(queue_frame, text="queueinfo", command=lambda: self.query_queue_info("all")).grid(row=0, column=0, padx=4, pady=6, sticky="ew")
        Button(queue_frame, text="queueinfo -gpu", command=lambda: self.query_queue_info("gpu")).grid(row=0, column=1, padx=4, pady=6, sticky="ew")
        Label(queue_frame, text="Queue name").grid(row=0, column=2, sticky="e", padx=(10, 4), pady=6)
        Entry(queue_frame, textvariable=self.queue_name_var, width=24).grid(row=0, column=3, padx=4, pady=6, sticky="ew")
        Button(queue_frame, text="queueinfo -l", command=lambda: self.query_queue_info("detail")).grid(row=0, column=4, padx=4, pady=6, sticky="ew")
        self._build_status_bar(queue_frame, self.queue_status_var, "Queue").grid(row=1, column=0, columnspan=5, sticky="ew", padx=6, pady=(0, 6))
        for column in (0, 1, 4):
            queue_frame.columnconfigure(column, weight=1)
        queue_frame.columnconfigure(3, weight=2)

        jobs_frame = LabelFrame(body, text="Jobs")
        jobs_frame.pack(fill="both", expand=True)
        jobs_toolbar = Frame(jobs_frame)
        jobs_toolbar.pack(fill="x", padx=6, pady=6)
        Button(jobs_toolbar, text="Refresh jobs", command=self.refresh_jobs).pack(side="left", padx=(0, 4))
        Button(jobs_toolbar, text="View bjobs -l", command=self.view_selected_job_details).pack(side="left", padx=4)
        Button(jobs_toolbar, text="Kill selected", command=self.kill_selected_jobs).pack(side="left", padx=4)
        Button(jobs_toolbar, text="Kill all visible", command=self.kill_all_visible_jobs).pack(side="left", padx=4)
        Checkbutton(jobs_toolbar, text="Auto refresh every 10s", variable=self.auto_refresh_var, command=self._schedule_auto_refresh).pack(side="left", padx=12)
        Label(jobs_toolbar, text="Search").pack(side="left", padx=(10, 4))
        Entry(jobs_toolbar, textvariable=self.job_filter_var, width=24).pack(side="left", padx=(0, 8))
        Label(jobs_toolbar, text="Status").pack(side="left", padx=(0, 4))
        self.job_status_box = ttk.Combobox(jobs_toolbar, textvariable=self.job_status_filter_var, values=["All"], state="readonly", width=10)
        self.job_status_box.pack(side="left")
        Label(jobs_toolbar, textvariable=self.job_selection_var, fg=self.colors["muted"]).pack(side="right", padx=6)
        self.job_filter_var.trace_add("write", lambda *_: self._apply_job_filter())
        self.job_status_filter_var.trace_add("write", lambda *_: self._apply_job_filter())

        self._build_status_bar(jobs_frame, self.job_summary_var, "Jobs").pack(fill="x", padx=6, pady=(0, 6))

        jobs_table = Frame(jobs_frame)
        jobs_table.pack(fill="both", expand=True, padx=6, pady=(0, 6))
        jobs_table.rowconfigure(0, weight=1)
        jobs_table.columnconfigure(0, weight=1)

        columns = ("jobid", "stat", "queue", "job_name", "exec_host", "submit_time", "user")
        labels = {
            "jobid": "Job ID",
            "stat": "Status",
            "queue": "Queue",
            "job_name": "Job name",
            "exec_host": "Exec host",
            "submit_time": "Submit time",
            "user": "User",
        }
        widths = {
            "jobid": 95,
            "stat": 78,
            "queue": 120,
            "job_name": 320,
            "exec_host": 220,
            "submit_time": 150,
            "user": 110,
        }
        self.jobs_tree = ttk.Treeview(jobs_table, columns=columns, show="headings", selectmode="extended", height=12)
        jobs_y_scroll = ttk.Scrollbar(jobs_table, orient="vertical", command=self.jobs_tree.yview)
        jobs_x_scroll = ttk.Scrollbar(jobs_table, orient="horizontal", command=self.jobs_tree.xview)
        self.jobs_tree.configure(yscrollcommand=jobs_y_scroll.set, xscrollcommand=jobs_x_scroll.set)
        for column in columns:
            self.jobs_tree.heading(column, text=labels[column], command=lambda selected_column=column: self._sort_jobs(selected_column))
            self.jobs_tree.column(column, width=widths[column], minwidth=70, anchor="w", stretch=False)
        self.jobs_tree.tag_configure("job_run", background="#e7f6e7")
        self.jobs_tree.tag_configure("job_pend", background="#fff4d6")
        self.jobs_tree.tag_configure("job_exit", background="#fde2e1")
        self.jobs_tree.tag_configure("job_done", background="#eef2f7")
        self.jobs_tree.grid(row=0, column=0, sticky="nsew")
        jobs_y_scroll.grid(row=0, column=1, sticky="ns")
        jobs_x_scroll.grid(row=1, column=0, sticky="ew")
        self.jobs_tree.bind("<<TreeviewSelect>>", self._update_job_selection_detail)

        detail_frame = LabelFrame(jobs_frame, text="Selected job summary")
        detail_frame.pack(fill="x", padx=6, pady=(0, 6))
        self.job_detail_text = Text(detail_frame, height=4, wrap="word")
        self.job_detail_text.pack(fill="x", padx=6, pady=6)
        self.job_detail_text.insert(END, "Select a job to inspect its parsed bjobs row. Use View bjobs -l for full pending reason and LSF details.")

    def _build_results_browser(self, parent: Frame):
        body = Frame(parent)
        body.pack(fill="both", expand=True, padx=12, pady=12)

        controls = LabelFrame(body, text="Remote result scan")
        controls.pack(fill="x", pady=(0, 8))
        self.result_search_var = StringVar()
        self.result_target_filter_var = StringVar(value="All")
        self.result_pilot_filter_var = StringVar(value="All")
        self.result_stage_filter_var = StringVar(value="All")
        self.result_shard_filter_var = StringVar(value="All")
        self.result_scratch_root_var = StringVar()
        self.result_scan_date_var = StringVar()
        self.result_scan_user_var = StringVar()
        self.result_scan_status_var = StringVar(value="Not scanned")
        self.result_status_display_var = StringVar(value="Not scanned | Selected: 0 files; Visible: 0; Scanned: 0")
        self.result_scan_status_var.trace_add("write", lambda *_: self._update_result_status_display())
        self.result_filter_boxes = {}
        Label(controls, text="Search").grid(row=0, column=0, sticky="w", padx=6, pady=6)
        Entry(controls, textvariable=self.result_search_var).grid(row=0, column=1, columnspan=2, sticky="ew", padx=6, pady=6)
        self.result_search_var.trace_add("write", lambda *_: self._apply_result_filter())
        Button(controls, text="Scan scratch", command=self.scan_results).grid(row=0, column=3, padx=3, pady=6)
        Button(controls, text="Download selected", command=self.download_selected_results).grid(row=0, column=4, padx=3, pady=6)
        Button(controls, text="Download all filtered", command=self.download_filtered_results).grid(row=0, column=5, padx=3, pady=6)
        Button(controls, text="Preview CSV", command=self.preview_selected_csv).grid(row=0, column=6, padx=3, pady=6)
        Button(controls, text="Open PDB in PyMOL", command=self.open_selected_pdb_in_pymol).grid(row=0, column=7, padx=3, pady=6)

        filters = [
            ("Target", "target", self.result_target_filter_var),
            ("Pilot", "pilot", self.result_pilot_filter_var),
            ("Stage", "stage", self.result_stage_filter_var),
            ("Shard", "shard", self.result_shard_filter_var),
        ]
        for offset, (label, key, variable) in enumerate(filters):
            column = offset * 2
            Label(controls, text=label).grid(row=1, column=column, sticky="w", padx=6, pady=(0, 6))
            box = ttk.Combobox(controls, textvariable=variable, values=["All"], state="readonly", width=16)
            box.grid(row=1, column=column + 1, sticky="ew", padx=6, pady=(0, 6))
            box.bind("<<ComboboxSelected>>", lambda _event: self._apply_result_filter())
            self.result_filter_boxes[key] = box
        scan_fields = [
            ("Scratch root", self.result_scratch_root_var, 18),
            ("Date", self.result_scan_date_var, 14),
            ("User folder", self.result_scan_user_var, 18),
        ]
        start_column = 0
        for label, variable, width in scan_fields:
            Label(controls, text=label).grid(row=2, column=start_column, sticky="w", padx=6, pady=(0, 6))
            Entry(controls, textvariable=variable, width=width).grid(row=2, column=start_column + 1, sticky="ew", padx=6, pady=(0, 6))
            start_column += 2
        Button(controls, text="Use current target", command=self.use_current_target_result_scan_settings).grid(row=2, column=6, padx=3, pady=(0, 6))
        self._build_status_bar(controls, self.result_status_display_var, "Results").grid(row=4, column=0, columnspan=8, sticky="ew", padx=6, pady=(0, 6))
        Button(controls, text="Scan logs/errors", command=self.scan_job_logs).grid(row=3, column=0, columnspan=2, sticky="ew", padx=6, pady=(0, 6))
        Button(controls, text="Preview log text", command=self.preview_selected_log_text).grid(row=3, column=2, columnspan=2, sticky="ew", padx=6, pady=(0, 6))
        controls.columnconfigure(1, weight=1)
        controls.columnconfigure(2, weight=1)
        controls.columnconfigure(5, weight=1)

        columns = ("target", "pilot", "stage", "shard", "type", "size", "mtime", "name", "path")
        results_table_frame = Frame(body)
        results_table_frame.pack(fill="both", expand=True)
        self.results_tree = ttk.Treeview(results_table_frame, columns=columns, show="headings", selectmode="extended")
        results_y_scroll = ttk.Scrollbar(results_table_frame, orient="vertical", command=self.results_tree.yview)
        results_x_scroll = ttk.Scrollbar(results_table_frame, orient="horizontal", command=self.results_tree.xview)
        self.results_tree.configure(yscrollcommand=results_y_scroll.set, xscrollcommand=results_x_scroll.set)
        self.results_tree.bind("<<TreeviewSelect>>", lambda _event: self._update_result_status_display())
        for column in columns:
            self.results_tree.heading(column, text=column, command=lambda col=column: self._sort_results(col))
            width = 90
            if column == "path":
                width = 760
            elif column == "name":
                width = 360
            elif column == "stage":
                width = 170
            elif column == "mtime":
                width = 130
            self.results_tree.column(column, width=width, minwidth=60, anchor="w", stretch=False)
        self.results_tree.grid(row=0, column=0, sticky="nsew")
        results_y_scroll.grid(row=0, column=1, sticky="ns")
        results_x_scroll.grid(row=1, column=0, sticky="ew")
        results_table_frame.rowconfigure(0, weight=1)
        results_table_frame.columnconfigure(0, weight=1)

        hint = Label(
            body,
            text="扫描目录为 <scratch root>/<date>/<user folder>；结果按 target / pilot / stage / shard 归类。下载文件会进入当前 workflow 的 retrieved_results 文件夹。",
        )
        hint.pack(fill="x", pady=(6, 0))
        self._load_result_scan_vars()

    def _build_data_processing(self, parent: Frame):
        body = Frame(parent)
        body.pack(fill="both", expand=True, padx=12, pady=12)

        self.data_csv_entries = []
        self.data_csv_all_entries = []
        self.data_csv_search_var = StringVar()
        self.data_target_filter_var = StringVar(value="All")
        self.data_pilot_filter_var = StringVar(value="All")
        self.data_stage_filter_var = StringVar(value="All")
        self.data_shard_filter_var = StringVar(value="All")
        self.data_filter_boxes = {}
        self.data_loaded_csv_path = None
        self.data_loaded_headers = []
        self.data_metric_var = StringVar()
        self.data_metric_direction_var = StringVar(value="higher")
        self.data_percentile_var = StringVar(value="90")
        self.data_threshold_var = StringVar()
        self.data_bins_var = StringVar(value="40")
        self.data_hist_min_var = StringVar()
        self.data_hist_max_var = StringVar()
        self.data_robust_hist_var = BooleanVar(value=True)
        self.data_status_var = StringVar(value="No merged CSV loaded. Steps: scan -> load for analysis -> plot histogram + scatter.")
        self.hit_afcyc_csv_var = StringVar()
        self.hit_pyro_csv_var = StringVar()
        self.hit_status_var = StringVar(value="No hit screening run")
        self.hit_rows = []
        self.hit_afcyc_ipae_col_var = StringVar(value=HIT_SCREEN_DEFAULTS["afcyc_ipae_column"])
        self.hit_afcyc_ipae_max_var = StringVar(value=HIT_SCREEN_DEFAULTS["afcyc_ipae_max"])
        self.hit_afcyc_rmsd_col_var = StringVar(value=HIT_SCREEN_DEFAULTS["afcyc_rmsd_column"])
        self.hit_afcyc_rmsd_max_var = StringVar(value=HIT_SCREEN_DEFAULTS["afcyc_rmsd_max"])
        self.hit_afcyc_plddt_col_var = StringVar(value=HIT_SCREEN_DEFAULTS["afcyc_plddt_column"])
        self.hit_afcyc_plddt_min_var = StringVar(value=HIT_SCREEN_DEFAULTS["afcyc_plddt_min"])
        self.hit_pyro_dg_col_var = StringVar(value=HIT_SCREEN_DEFAULTS["pyro_dg_column"])
        self.hit_pyro_dg_max_var = StringVar(value=HIT_SCREEN_DEFAULTS["pyro_dg_max"])
        self.hit_pyro_sap_col_var = StringVar(value=HIT_SCREEN_DEFAULTS["pyro_sap_column"])
        self.hit_pyro_sap_max_var = StringVar(value=HIT_SCREEN_DEFAULTS["pyro_sap_max"])
        self.hit_pyro_cms_col_var = StringVar(value=HIT_SCREEN_DEFAULTS["pyro_cms_column"])
        self.hit_pyro_cms_min_var = StringVar(value=HIT_SCREEN_DEFAULTS["pyro_cms_min"])
        self.hit_pyro_sasa_col_var = StringVar(value=HIT_SCREEN_DEFAULTS["pyro_sasa_column"])
        self.hit_pyro_sasa_min_var = StringVar(value=HIT_SCREEN_DEFAULTS["pyro_sasa_min"])
        self.hit_pyro_packstat_col_var = StringVar(value=HIT_SCREEN_DEFAULTS["pyro_packstat_column"])
        self.hit_pyro_packstat_min_var = StringVar(value=HIT_SCREEN_DEFAULTS["pyro_packstat_min"])
        self.hit_limit_var = StringVar(value="200")

        scan_frame = LabelFrame(body, text="Merged CSV scan")
        scan_frame.pack(fill="x", pady=(0, 8))
        Label(scan_frame, text="Scan root").grid(row=0, column=0, sticky="w", padx=6, pady=6)
        Label(scan_frame, textvariable=self.result_scan_status_var).grid(row=0, column=1, sticky="w", padx=6, pady=6)
        Button(scan_frame, text="Use current target", command=self.use_current_target_result_scan_settings).grid(row=0, column=2, padx=3, pady=6)
        Button(scan_frame, text="Scan merged CSV", command=self.scan_data_merged_csvs).grid(row=0, column=3, padx=3, pady=6)
        Button(scan_frame, text="Preview selected CSV", command=self.preview_data_csv).grid(row=0, column=4, padx=3, pady=6)
        Button(scan_frame, text="Load for analysis", command=self.load_data_csv_for_analysis).grid(row=0, column=5, padx=3, pady=6)
        Label(scan_frame, text="Search").grid(row=1, column=0, sticky="w", padx=6, pady=(0, 6))
        Entry(scan_frame, textvariable=self.data_csv_search_var).grid(row=1, column=1, columnspan=2, sticky="ew", padx=6, pady=(0, 6))
        self.data_csv_search_var.trace_add("write", lambda *_: self._apply_data_csv_filter())
        data_filters = [
            ("Target", "target", self.data_target_filter_var),
            ("Pilot", "pilot", self.data_pilot_filter_var),
            ("Stage", "stage", self.data_stage_filter_var),
            ("Shard", "shard", self.data_shard_filter_var),
        ]
        for offset, (label, key, variable) in enumerate(data_filters):
            column = offset * 2
            Label(scan_frame, text=label).grid(row=2, column=column, sticky="w", padx=6, pady=(0, 6))
            box = ttk.Combobox(scan_frame, textvariable=variable, values=["All"], state="readonly", width=16)
            box.grid(row=2, column=column + 1, sticky="ew", padx=6, pady=(0, 6))
            box.bind("<<ComboboxSelected>>", lambda _event: self._apply_data_csv_filter())
            self.data_filter_boxes[key] = box
        self._build_status_bar(scan_frame, self.data_status_var, "Data").grid(row=3, column=0, columnspan=8, sticky="ew", padx=6, pady=(0, 6))
        scan_frame.columnconfigure(1, weight=1)

        csv_frame = Frame(body)
        csv_frame.pack(fill="both", expand=False, pady=(0, 8))
        csv_columns = ("stage", "size", "mtime", "name", "path")
        self.data_csv_tree = ttk.Treeview(csv_frame, columns=csv_columns, show="headings", selectmode="browse", height=6)
        csv_y_scroll = ttk.Scrollbar(csv_frame, orient="vertical", command=self.data_csv_tree.yview)
        csv_x_scroll = ttk.Scrollbar(csv_frame, orient="horizontal", command=self.data_csv_tree.xview)
        self.data_csv_tree.configure(yscrollcommand=csv_y_scroll.set, xscrollcommand=csv_x_scroll.set)
        for column in csv_columns:
            self.data_csv_tree.heading(column, text=column)
            width = 100
            if column == "path":
                width = 720
            elif column == "name":
                width = 260
            elif column == "mtime":
                width = 130
            self.data_csv_tree.column(column, width=width, minwidth=70, anchor="w", stretch=False)
        self.data_csv_tree.grid(row=0, column=0, sticky="nsew")
        csv_y_scroll.grid(row=0, column=1, sticky="ns")
        csv_x_scroll.grid(row=1, column=0, sticky="ew")
        csv_frame.rowconfigure(0, weight=1)
        csv_frame.columnconfigure(0, weight=1)

        analysis_frame = LabelFrame(body, text="PyRosetta metric analysis")
        analysis_frame.pack(fill="x", pady=(0, 8))
        Label(analysis_frame, text="Metric").grid(row=0, column=0, sticky="w", padx=6, pady=6)
        self.data_metric_box = ttk.Combobox(analysis_frame, textvariable=self.data_metric_var, values=[], state="readonly", width=32)
        self.data_metric_box.grid(row=0, column=1, sticky="ew", padx=6, pady=6)
        Button(analysis_frame, text="AfCyc defaults", command=self.use_afcyc_metric_defaults).grid(row=0, column=2, sticky="ew", padx=3, pady=6)
        Button(analysis_frame, text="PyRosetta defaults", command=self.use_pyrosetta_metric_defaults).grid(row=0, column=3, sticky="ew", padx=3, pady=6)
        Label(analysis_frame, text="Direction").grid(row=0, column=4, sticky="w", padx=6, pady=6)
        ttk.Combobox(analysis_frame, textvariable=self.data_metric_direction_var, values=["higher", "lower"], state="readonly", width=9).grid(row=0, column=5, sticky="w", padx=6, pady=6)
        Label(analysis_frame, text="Percentile").grid(row=0, column=6, sticky="w", padx=6, pady=6)
        Entry(analysis_frame, textvariable=self.data_percentile_var, width=8).grid(row=0, column=7, sticky="w", padx=6, pady=6)
        Label(analysis_frame, text="Custom threshold").grid(row=1, column=0, sticky="w", padx=6, pady=(0, 6))
        Entry(analysis_frame, textvariable=self.data_threshold_var, width=12).grid(row=1, column=1, sticky="w", padx=6, pady=(0, 6))
        Label(analysis_frame, text="Bins").grid(row=1, column=2, sticky="w", padx=6, pady=(0, 6))
        Entry(analysis_frame, textvariable=self.data_bins_var, width=8).grid(row=1, column=3, sticky="w", padx=6, pady=(0, 6))
        Label(analysis_frame, text="Hist min").grid(row=1, column=4, sticky="w", padx=6, pady=(0, 6))
        Entry(analysis_frame, textvariable=self.data_hist_min_var, width=12).grid(row=1, column=5, sticky="w", padx=6, pady=(0, 6))
        Label(analysis_frame, text="Hist max").grid(row=1, column=6, sticky="w", padx=6, pady=(0, 6))
        Entry(analysis_frame, textvariable=self.data_hist_max_var, width=12).grid(row=1, column=7, sticky="w", padx=6, pady=(0, 6))
        Checkbutton(analysis_frame, text="Auto robust range 1-99% when min/max blank", variable=self.data_robust_hist_var).grid(row=2, column=0, columnspan=8, sticky="w", padx=6, pady=(0, 6))
        Button(
            analysis_frame,
            text="Plot histogram + scatter",
            command=self.analyze_data_metric,
            bg=self.colors["primary"],
            fg="#ffffff",
            activebackground=self.colors["primary_dark"],
        ).grid(row=3, column=0, columnspan=2, sticky="ew", padx=6, pady=(0, 6))
        self._build_status_bar(analysis_frame, self.data_status_var, "Data").grid(row=3, column=2, columnspan=6, sticky="ew", padx=6, pady=(0, 6))
        analysis_frame.columnconfigure(1, weight=1)

        plot_frame = LabelFrame(body, text="Metric plots")
        plot_frame.pack(fill="both", expand=True, pady=(0, 8))
        plot_toolbar = Frame(plot_frame)
        plot_toolbar.pack(fill="x", padx=6, pady=(6, 0))
        Button(
            plot_toolbar,
            text="Draw / refresh histogram + scatter",
            command=self.analyze_data_metric,
            bg=self.colors["primary"],
            fg="#ffffff",
            activebackground=self.colors["primary_dark"],
        ).pack(side="left", padx=(0, 8))
        Label(
            plot_toolbar,
            text="Load a merged CSV, choose a metric, then click this button to render both plots.",
            fg=self.colors["muted"],
        ).pack(side="left")
        hist_frame = LabelFrame(plot_frame, text="Histogram")
        hist_frame.pack(fill="both", expand=True, padx=6, pady=(6, 3))
        self.data_hist_canvas = Canvas(hist_frame, height=260, bg="white")
        self.data_hist_canvas.pack(fill="both", expand=True, padx=6, pady=6)
        self.data_hist_canvas.create_text(240, 130, text="Histogram will appear after plotting.", fill="#64748b")
        scatter_frame = LabelFrame(plot_frame, text="Scatter by CSV order")
        scatter_frame.pack(fill="both", expand=True, padx=6, pady=(3, 6))
        self.data_scatter_canvas = Canvas(scatter_frame, height=260, bg="white")
        self.data_scatter_canvas.pack(fill="both", expand=True, padx=6, pady=6)
        self.data_scatter_canvas.create_text(240, 130, text="Scatter plot will appear after plotting.", fill="#64748b")

        candidates_frame = LabelFrame(body, text="Candidate rows within display range above percentile/custom threshold")
        candidates_frame.pack(fill="both", expand=True)
        candidate_columns = ("rank", "row", "value", "candidate", "summary")
        self.data_candidate_tree = ttk.Treeview(candidates_frame, columns=candidate_columns, show="headings", height=8)
        candidate_y_scroll = ttk.Scrollbar(candidates_frame, orient="vertical", command=self.data_candidate_tree.yview)
        candidate_x_scroll = ttk.Scrollbar(candidates_frame, orient="horizontal", command=self.data_candidate_tree.xview)
        self.data_candidate_tree.configure(yscrollcommand=candidate_y_scroll.set, xscrollcommand=candidate_x_scroll.set)
        for column in candidate_columns:
            self.data_candidate_tree.heading(column, text=column)
            width = 90
            if column == "candidate":
                width = 260
            elif column == "summary":
                width = 720
            self.data_candidate_tree.column(column, width=width, minwidth=60, anchor="w", stretch=False)
        self.data_candidate_tree.grid(row=0, column=0, sticky="nsew")
        candidate_y_scroll.grid(row=0, column=1, sticky="ns")
        candidate_x_scroll.grid(row=1, column=0, sticky="ew")
        candidates_frame.rowconfigure(0, weight=1)
        candidates_frame.columnconfigure(0, weight=1)

        self._build_hit_screening(body)

    def _build_hit_screening(self, parent: Frame):
        hit_frame = LabelFrame(parent, text="Hit screening from AfCycDesign + PyRosetta merged CSV")
        hit_frame.pack(fill="both", expand=True, pady=(0, 8))

        source_frame = LabelFrame(hit_frame, text="CSV inputs")
        source_frame.pack(fill="x", padx=6, pady=6)
        Label(source_frame, text="AfCyc CSV").grid(row=0, column=0, sticky="w", padx=6, pady=6)
        Entry(source_frame, textvariable=self.hit_afcyc_csv_var).grid(row=0, column=1, sticky="ew", padx=6, pady=6)
        Button(source_frame, text="Use selected scanned CSV", command=lambda: self.set_hit_csv_from_selected("afcyc")).grid(row=0, column=2, padx=3, pady=6)
        Button(source_frame, text="Browse", command=lambda: self.browse_hit_csv("afcyc")).grid(row=0, column=3, padx=3, pady=6)
        Label(source_frame, text="PyRosetta CSV").grid(row=1, column=0, sticky="w", padx=6, pady=(0, 6))
        Entry(source_frame, textvariable=self.hit_pyro_csv_var).grid(row=1, column=1, sticky="ew", padx=6, pady=(0, 6))
        Button(source_frame, text="Use selected scanned CSV", command=lambda: self.set_hit_csv_from_selected("pyro")).grid(row=1, column=2, padx=3, pady=(0, 6))
        Button(source_frame, text="Browse", command=lambda: self.browse_hit_csv("pyro")).grid(row=1, column=3, padx=3, pady=(0, 6))
        source_frame.columnconfigure(1, weight=1)

        threshold_frame = LabelFrame(hit_frame, text="RFpeptides-inspired thresholds")
        threshold_frame.pack(fill="x", padx=6, pady=(0, 6))
        hit_threshold_rows = [
            ("AfCyc iPAE col", self.hit_afcyc_ipae_col_var, "<=", self.hit_afcyc_ipae_max_var),
            ("AfCyc RMSD col", self.hit_afcyc_rmsd_col_var, "<=", self.hit_afcyc_rmsd_max_var),
            ("binder pLDDT col", self.hit_afcyc_plddt_col_var, ">=", self.hit_afcyc_plddt_min_var),
            ("interface dG col", self.hit_pyro_dg_col_var, "<=", self.hit_pyro_dg_max_var),
            ("SAP col", self.hit_pyro_sap_col_var, "<=", self.hit_pyro_sap_max_var),
            ("CMS col", self.hit_pyro_cms_col_var, ">=", self.hit_pyro_cms_min_var),
            ("interface SASA col", self.hit_pyro_sasa_col_var, ">=", self.hit_pyro_sasa_min_var),
            ("packstat col", self.hit_pyro_packstat_col_var, ">=", self.hit_pyro_packstat_min_var),
        ]
        for row, (label, column_var, operator, threshold_var) in enumerate(hit_threshold_rows):
            Label(threshold_frame, text=label).grid(row=row, column=0, sticky="w", padx=6, pady=4)
            Entry(threshold_frame, textvariable=column_var, width=26).grid(row=row, column=1, sticky="ew", padx=6, pady=4)
            Label(threshold_frame, text=operator).grid(row=row, column=2, sticky="w")
            Entry(threshold_frame, textvariable=threshold_var, width=10).grid(row=row, column=3, sticky="w", padx=6, pady=4)

        limit_row = len(hit_threshold_rows)
        Label(threshold_frame, text="Max hits").grid(row=limit_row, column=0, sticky="w", padx=6, pady=(2, 6))
        Entry(threshold_frame, textvariable=self.hit_limit_var, width=10).grid(row=limit_row, column=1, sticky="w", padx=6, pady=(2, 6))
        Button(threshold_frame, text="Reset paper defaults", command=self.reset_hit_thresholds).grid(row=limit_row, column=2, columnspan=2, sticky="ew", padx=6, pady=(2, 6))
        Label(
            threshold_frame,
            text="Defaults follow RFpeptides-style filters: low iPAE/RMSD/ddG/SAP and high CMS; pLDDT is optional if blank.",
            fg=self.colors["muted"],
            wraplength=560,
            justify="left",
        ).grid(row=limit_row + 1, column=0, columnspan=4, sticky="w", padx=6, pady=(0, 6))
        threshold_frame.columnconfigure(1, weight=1)

        actions = Frame(hit_frame)
        actions.pack(fill="x", padx=6, pady=(0, 6))
        hit_action_buttons = [
            ("Run hit screening", self.run_hit_screening, True),
            ("Save hit CSV", self.save_hit_screening_csv, False),
            ("Download selected hit structures", lambda: self.download_hit_structures(selected_only=True), False),
            ("Download all hit structures", lambda: self.download_hit_structures(selected_only=False), False),
        ]
        for index, (text, command, primary) in enumerate(hit_action_buttons):
            options = {
                "text": text,
                "command": command,
            }
            if primary:
                options.update({"bg": self.colors["primary"], "fg": "#ffffff", "activebackground": self.colors["primary_dark"]})
            Button(actions, **options).grid(row=index // 2, column=index % 2, sticky="ew", padx=3, pady=3)
        actions.columnconfigure(0, weight=1)
        actions.columnconfigure(1, weight=1)
        self._build_status_bar(hit_frame, self.hit_status_var, "Hits").pack(fill="x", padx=6, pady=(0, 6))

        table_frame = Frame(hit_frame)
        table_frame.pack(fill="both", expand=True, padx=6, pady=(0, 6))
        table_frame.rowconfigure(0, weight=1)
        table_frame.columnconfigure(0, weight=1)
        hit_columns = ("rank", "description", "i_pae", "rmsd", "plddt", "interface_dG", "sap", "cms", "binder_seq", "afcyc_pdb", "mpnn_pdb")
        self.hit_tree = ttk.Treeview(table_frame, columns=hit_columns, show="headings", selectmode="extended", height=8)
        hit_y_scroll = ttk.Scrollbar(table_frame, orient="vertical", command=self.hit_tree.yview)
        hit_x_scroll = ttk.Scrollbar(table_frame, orient="horizontal", command=self.hit_tree.xview)
        self.hit_tree.configure(yscrollcommand=hit_y_scroll.set, xscrollcommand=hit_x_scroll.set)
        widths = {
            "rank": 60,
            "description": 300,
            "i_pae": 90,
            "rmsd": 90,
            "plddt": 90,
            "interface_dG": 110,
            "sap": 90,
            "cms": 90,
            "binder_seq": 160,
            "afcyc_pdb": 520,
            "mpnn_pdb": 520,
        }
        for column in hit_columns:
            self.hit_tree.heading(column, text=column)
            self.hit_tree.column(column, width=widths.get(column, 100), minwidth=60, anchor="w", stretch=False)
        self.hit_tree.grid(row=0, column=0, sticky="nsew")
        hit_y_scroll.grid(row=0, column=1, sticky="ns")
        hit_x_scroll.grid(row=1, column=0, sticky="ew")

    def _build_scratch_safety(self, parent: Frame):
        body = Frame(parent)
        body.pack(fill="both", expand=True, padx=12, pady=12)

        info = LabelFrame(body, text="Move scratch user directory")
        info.pack(fill="x", pady=(0, 8))
        self.scratch_root_var = StringVar()
        self.scratch_source_date_var = StringVar()
        self.scratch_target_date_var = StringVar()
        self.scratch_user_dir_var = StringVar()
        self.scratch_status_var = StringVar(value="Scratch safety idle.")
        fields = [
            ("Scratch root", self.scratch_root_var),
            ("Current date folder", self.scratch_source_date_var),
            ("Target date folder", self.scratch_target_date_var),
            ("User folder", self.scratch_user_dir_var),
        ]
        for row, (label, variable) in enumerate(fields):
            Label(info, text=label).grid(row=row, column=0, sticky="w", padx=6, pady=5)
            Entry(info, textvariable=variable).grid(row=row, column=1, sticky="ew", padx=6, pady=5)
        info.columnconfigure(1, weight=1)

        actions = Frame(body)
        actions.pack(fill="x", pady=(0, 8))
        Button(actions, text="Use current target settings", command=self.use_current_target_scratch_settings).pack(side="left", padx=(0, 6))
        Button(actions, text="Preview move", command=self.preview_scratch_move).pack(side="left", padx=(0, 6))
        Button(actions, text="Move now", command=self.execute_scratch_move).pack(side="left", padx=6)
        Label(actions, text="target_date=auto uses the cluster date from `date +%F`.").pack(side="left", padx=12)
        self._build_status_bar(body, self.scratch_status_var, "Scratch").pack(fill="x", pady=(0, 8))

        explanation = LabelFrame(body, text="Safety rules")
        explanation.pack(fill="both", expand=True)
        text = Text(explanation, height=12, wrap="word")
        text.pack(fill="both", expand=True, padx=6, pady=6)
        text.insert(
            END,
            "This feature is isolated from job submission and result download.\n\n"
            "It moves exactly one directory:\n"
            "  <scratch_root>/<current date folder>/<user folder>\n"
            "to:\n"
            "  <scratch_root>/<target date folder>/<user folder>\n\n"
            "It refuses to run if the source directory does not exist.\n"
            "It refuses to run if the target user directory already exists, so it will not merge or overwrite data.\n"
            "Use Preview move first and check Raw SSH/LSF output before Move now.\n",
        )
        text.configure(state="disabled")
        self._load_scratch_migration_vars()

    def _load_initial_cluster_profile(self) -> dict:
        if cluster_ops.DEFAULT_PROFILE.exists():
            self.cluster_profile_path = cluster_ops.DEFAULT_PROFILE
            return cluster_ops.load_profile(cluster_ops.DEFAULT_PROFILE)
        self.cluster_profile_path = cluster_ops.DEFAULT_EXAMPLE_PROFILE
        return cluster_ops.default_profile()

    def show_cluster_settings(self):
        self._show_page("Cluster Settings")

    def _load_cluster_settings_vars(self):
        if not hasattr(self, "cluster_host_var"):
            return
        profile = self.cluster_profile
        self.cluster_profile_path_var.set(str(self.cluster_profile_path))
        self.cluster_host_var.set(str(profile.get("host") or ""))
        self.cluster_user_var.set(str(profile.get("user") or ""))
        self.cluster_port_var.set(str(profile.get("port") or "22"))
        self.auth_method_var.set(str(profile.get("auth_method") or "ssh_key"))
        self.password_var.set("")
        self.remember_password_var.set(bool(profile.get("remember_password")))
        self.cluster_ssh_key_var.set(str(profile.get("ssh_key") or ""))
        self.cluster_remote_upload_parent_var.set(str(profile.get("remote_upload_parent") or ""))
        self.cluster_remote_workflow_dir_var.set(str(profile.get("remote_workflow_dir") or ""))
        self.cluster_local_download_dir_var.set(str(profile.get("local_download_dir") or ""))
        self.cluster_pymol_var.set(str(profile.get("pymol_executable") or ""))
        self.cluster_bjobs_user_var.set(str(profile.get("bjobs_user") or "$USER"))
        self.cluster_scan_depth_var.set(str(profile.get("scratch_scan_max_depth") or 8))
        self.cluster_log_depth_var.set(str(profile.get("job_log_scan_max_depth") or profile.get("scratch_scan_max_depth") or 8))
        self.cluster_result_patterns_var.set(", ".join(str(item) for item in profile.get("result_file_patterns", [])))
        self.cluster_ssh_extra_args_var.set(" ".join(shlex.quote(str(item)) for item in profile.get("ssh_extra_args", [])))
        if self.remember_password_var.get():
            self.load_saved_password(silent=True)
        self.cluster_settings_status_var.set(f"Loaded cluster profile: {self.cluster_profile_path}")

    def choose_cluster_ssh_key(self):
        path = filedialog.askopenfilename(
            title="Choose SSH private key",
            initialdir=str(Path.home() / ".ssh"),
            filetypes=[("SSH key files", "id_* *.pem *.key"), ("All files", "*.*")],
        )
        if path:
            self.cluster_ssh_key_var.set(path)

    def choose_pymol_executable(self):
        path = filedialog.askopenfilename(
            title="Choose PyMOL executable",
            initialdir=str(Path.home()),
            filetypes=[("Executable files", "*.exe *.bat *.cmd"), ("All files", "*.*")],
        )
        if path:
            self.cluster_pymol_var.set(path)

    def _load_scratch_migration_vars(self):
        if not hasattr(self, "scratch_root_var"):
            return
        migration = self.cluster_profile.get("scratch_migration", {})
        self.scratch_root_var.set(str(migration.get("scratch_root") or "/scratch"))
        self.scratch_source_date_var.set(str(migration.get("source_date") or ""))
        self.scratch_target_date_var.set(str(migration.get("target_date") or "auto"))
        self.scratch_user_dir_var.set(str(migration.get("user_dir") or self.cluster_profile.get("user") or ""))
        if hasattr(self, "auth_method_var"):
            self.auth_method_var.set(str(self.cluster_profile.get("auth_method") or "ssh_key"))
        if hasattr(self, "remember_password_var"):
            self.remember_password_var.set(bool(self.cluster_profile.get("remember_password")))
            if self.remember_password_var.get():
                self.load_saved_password(silent=True)

    def _load_result_scan_vars(self):
        if not hasattr(self, "result_scan_date_var"):
            return
        data = deepcopy(self.config_data)
        roots = self.cluster_profile.get("scratch_scan_roots") or []
        parsed_root = self._parse_scratch_user_root(str(roots[0])) if roots else {}
        self.result_scratch_root_var.set(parsed_root.get("scratch_root") or "/scratch")
        self.result_scan_date_var.set(parsed_root.get("date") or str(data.get("scratch_date") or ""))
        self.result_scan_user_var.set(parsed_root.get("user") or str(data.get("cluster_user") or self.cluster_profile.get("user") or ""))
        root = self._remote_result_scan_root()
        self.result_scan_status_var.set(f"Scan root: {root}" if root else "Set scan date/user")

    def _parse_scratch_user_root(self, root: str) -> dict:
        parts = [part for part in root.replace("\\", "/").split("/") if part]
        if len(parts) >= 3 and parts[0] == "scratch":
            return {"scratch_root": "/scratch", "date": parts[1], "user": parts[2]}
        return {}

    def _remote_result_scan_root(self) -> str:
        scratch_root = self.result_scratch_root_var.get().strip().rstrip("/") or "/scratch"
        date = self.result_scan_date_var.get().strip()
        user = self.result_scan_user_var.get().strip()
        if not date or not user:
            return ""
        return f"{scratch_root}/{date}/{user}"

    def use_current_target_result_scan_settings(self):
        data = self._collect_config()
        self.result_scratch_root_var.set("/scratch")
        self.result_scan_date_var.set(str(data.get("scratch_date") or ""))
        self.result_scan_user_var.set(str(data.get("cluster_user") or self.cluster_profile.get("user") or ""))
        root = self._remote_result_scan_root()
        self.result_scan_status_var.set(f"Scan root: {root}" if root else "Set scan date/user")

    def _build_group(self, parent: Frame, group_name: str, fields):
        body = Frame(parent)
        body.pack(fill="both", expand=True, padx=12, pady=12)

        for row, (key, label, value_type) in enumerate(fields):
            Label(body, text=label).grid(row=row, column=0, sticky="w", padx=(0, 10), pady=5)
            if value_type is bool:
                var = BooleanVar()
                widget = Checkbutton(body, variable=var)
                widget.grid(row=row, column=1, sticky="w", pady=5)
            else:
                var = StringVar()
                widget = Entry(body, textvariable=var)
                widget.grid(row=row, column=1, sticky="ew", pady=5)
            self.variables[key] = (var, value_type)
        body.columnconfigure(1, weight=1)

        if group_name == "Project":
            self.project_sync_status_var = StringVar(value="Project-derived fields have not been synchronized.")
            sync_frame = LabelFrame(body, text="Project configuration helper")
            sync_frame.grid(row=len(fields), column=0, columnspan=2, sticky="ew", pady=(14, 0))
            Label(
                sync_frame,
                text=(
                    "After entering target, pilot, cluster user/date and project directories, synchronize derived paths "
                    "to RFDiffusion, AfCycDesign, PyRosetta, Results Browser, Scratch Safety and Dashboard."
                ),
                wraplength=920,
                justify="left",
                fg=self.colors["muted"],
            ).grid(row=0, column=0, columnspan=2, sticky="w", padx=8, pady=(8, 4))
            Button(
                sync_frame,
                text="Apply project info to all pages",
                command=self.apply_project_info_to_all_pages,
                bg=self.colors["primary"],
                fg="#ffffff",
                activebackground=self.colors["primary_dark"],
            ).grid(row=1, column=0, sticky="ew", padx=8, pady=(4, 8))
            self._build_status_bar(sync_frame, self.project_sync_status_var, "Project").grid(
                row=1,
                column=1,
                sticky="ew",
                padx=(0, 8),
                pady=(4, 8),
            )
            sync_frame.columnconfigure(0, weight=1)
            sync_frame.columnconfigure(1, weight=3)

    def apply_project_info_to_all_pages(self):
        try:
            target = self._project_field("target")
            if not target:
                raise ValueError("Target name is required")
            pilot = self._project_field("pilot") or "pilot0"
            project_dir = self._project_field("project_dir_name") or f"{target}_test"
            home_project_dir = self._project_field("home_project_dir") or f"$HOME/{target}"
            cluster_user = self._project_field("cluster_user")
            scratch_date = self._project_field("scratch_date")

            self._set_config_variable("pilot", pilot)
            self._set_config_variable("project_dir_name", project_dir)
            self._set_config_variable("home_project_dir", home_project_dir)

            changed = [
                f"pilot={pilot}",
                f"project_dir_name={project_dir}",
                f"home_project_dir={home_project_dir}",
            ]
            self._replace_demo_bundle_sources(target, changed)

            preview_profile = deepcopy(self.cluster_profile)
            if hasattr(self, "cluster_remote_upload_parent_var"):
                preview_profile["remote_upload_parent"] = self.cluster_remote_upload_parent_var.get().strip()
            if hasattr(self, "cluster_remote_workflow_dir_var"):
                preview_profile["remote_workflow_dir"] = self.cluster_remote_workflow_dir_var.get().strip()
            base_data = self._collect_config()
            remote_workflow_dir = self._remote_workflow_dir_for_config(preview_profile, base_data).rstrip("/")
            if hasattr(self, "cluster_remote_workflow_dir_var"):
                self.cluster_remote_workflow_dir_var.set(remote_workflow_dir)

            local_pdb = self.top_local_pdb_var.get().strip() if hasattr(self, "top_local_pdb_var") else ""
            current_input_pdb = self._project_field("input_pdb")
            input_name = self._path_basename(local_pdb or current_input_pdb)
            changed.append(f"remote_workflow_dir={remote_workflow_dir}")
            if input_name:
                remote_input_pdb = f"{remote_workflow_dir}/inputs/{input_name}"
                self._set_config_variable("input_pdb", remote_input_pdb)
                changed.append(f"input_pdb={remote_input_pdb}")

            bundle_names = self._workflow_bundle_names(base_data)
            derived_paths = {
                "afcyc.afcyc_script": f"{remote_workflow_dir}/scripts/{bundle_names['afcyc_script']}",
                "afcyc.rmsd_script": f"{remote_workflow_dir}/scripts/{bundle_names['afcyc_rmsd_script']}",
                "afcyc.merge_script": f"{remote_workflow_dir}/scripts/{bundle_names['afcyc_merge_script']}",
                "pyrosetta.script": f"{remote_workflow_dir}/scripts/{bundle_names['pyro_script']}",
                "pyrosetta.merge_script": f"{remote_workflow_dir}/scripts/{bundle_names['pyro_merge_script']}",
            }
            for key, value in derived_paths.items():
                self._set_config_variable(key, value)
                changed.append(f"{key}={value}")

            if hasattr(self, "result_scratch_root_var"):
                self.result_scratch_root_var.set("/scratch")
                self.result_scan_date_var.set(scratch_date)
                self.result_scan_user_var.set(cluster_user)
                root = self._remote_result_scan_root()
                self.result_scan_status_var.set(f"Scan root: {root}" if root else "Set scan date/user")
            if hasattr(self, "scratch_source_date_var"):
                self.scratch_source_date_var.set(scratch_date)
                self.scratch_user_dir_var.set(cluster_user)

            self._refresh_stage_selector()
            self._update_current_project_banner()
            self._update_parameter_summary()
            self.project_sync_status_var.set(
                f"Synchronization complete for {target}/{pilot}: stage paths, PDB path, scan and scratch context updated. Save config when ready."
            )
            self._log("Applied project information to dependent pages:")
            for item in changed:
                self._log(f"  {item}")
        except Exception as exc:
            self._show_error("Failed to synchronize project configuration", exc)

    def _project_field(self, key: str) -> str:
        variable_entry = self.variables.get(key)
        if not variable_entry:
            return ""
        return str(variable_entry[0].get() or "").strip()

    def _set_config_variable(self, key: str, value):
        variable_entry = self.variables.get(key)
        if not variable_entry:
            return
        variable, value_type = variable_entry
        if value_type is bool:
            variable.set(bool(value))
        else:
            variable.set(str(value))

    def _path_basename(self, path_text: str) -> str:
        normalized = str(path_text or "").strip().replace("\\", "/").rstrip("/")
        return normalized.rsplit("/", 1)[-1] if normalized else ""

    def _workflow_bundle_names(self, data: dict) -> dict[str, str]:
        afcyc = data.get("afcyc", {})
        pyrosetta = data.get("pyrosetta", {})
        return {
            "afcyc_script": self._path_basename(afcyc.get("local_afcyc_script") or "afcyc_predict_batch.py"),
            "afcyc_rmsd_script": self._path_basename(afcyc.get("local_rmsd_script") or "rmsd_from_afcyc.py"),
            "afcyc_merge_script": self._path_basename(afcyc.get("local_merge_script") or "merge_afcyc_csvs.py"),
            "pyro_script": self._path_basename(pyrosetta.get("local_script") or "PyRosetta_fullScoring_v4_debug.py"),
            "pyro_merge_script": self._path_basename(pyrosetta.get("local_merge_script") or "merge_pyrosetta_csvs.py"),
        }

    def _replace_demo_bundle_sources(self, target: str, changed: list[str]):
        if target.upper().startswith("DUMMY"):
            return
        candidates = {
            "afcyc.local_afcyc_script": "../Test3_PGLYRP1/3.AfCycDesign/v3/afcyc_predict_batch.py",
            "afcyc.local_rmsd_script": "../Test3_PGLYRP1/3.AfCycDesign/v3/rmsd_from_afcyc.py",
            "afcyc.local_merge_script": "../Test3_PGLYRP1/3.AfCycDesign/v3/merge_afcyc_csvs.py",
            "pyrosetta.local_script": "../Test3_PGLYRP1/4.PyRosetta/PyRosetta_fullScoring_v4_debug.py",
            "pyrosetta.local_merge_script": "../Exercise_Phase2_design/4.PyRosetta/merge_pyrosetta_csvs.py",
        }
        for key, candidate in candidates.items():
            current = self._project_field(key)
            if current and not self._is_demo_bundle_path(current):
                continue
            candidate_path = generate_workflow.resolve_local_path(candidate)
            if not candidate_path.is_file():
                continue
            self._set_config_variable(key, candidate)
            changed.append(f"{key}={candidate}")

    def _is_demo_bundle_path(self, path_text: str) -> bool:
        normalized = str(path_text or "").replace("\\", "/").lower()
        if "/examples/scripts/" in f"/{normalized.lstrip('/')}":
            return True
        try:
            path = generate_workflow.resolve_local_path(path_text)
            if path.is_file():
                text = path.read_text(encoding="utf-8", errors="replace")[:2000].lower()
                return "dummy " in text and "entrypoint" in text
        except Exception:
            pass
        return False

    def _load_values(self, config: dict):
        for key, (var, value_type) in self.variables.items():
            value = get_nested(config, key)
            if (value is None or (value == "" and value_type is int)) and key in FIELD_DEFAULTS:
                value = FIELD_DEFAULTS[key]
            if value_type is bool:
                var.set(bool(value))
            else:
                var.set("" if value is None else str(value))
        if hasattr(self, "top_local_pdb_var"):
            local_input = str(config.get("local_input_pdb") or "")
            self.top_local_pdb_var.set(local_input)
            if hasattr(self, "pdb_raw_path_var") and not self.pdb_raw_path_var.get().strip():
                self.pdb_raw_path_var.set(local_input)
                if local_input:
                    try:
                        self.pdb_clean_output_var.set(str(self._default_clean_pdb_output_path(self._resolve_local_config_path(local_input))))
                    except Exception:
                        pass
        if hasattr(self, "parameter_summary"):
            self._update_parameter_summary()
        self._update_current_project_banner()
        self._load_result_scan_vars()

    def _collect_config(self) -> dict:
        data = deepcopy(self.config_data)
        for key, (var, value_type) in self.variables.items():
            raw_value = var.get()
            if value_type is bool:
                value = bool(raw_value)
            elif value_type is int:
                raw_text = str(raw_value).strip()
                try:
                    value = int(raw_text)
                except ValueError as exc:
                    raise ValueError(f"{key} must be an integer") from exc
                if value < 0:
                    raise ValueError(f"{key} must be non-negative")
            else:
                value = str(raw_value).strip()
            set_nested(data, key, value)
        if hasattr(self, "top_local_pdb_var"):
            data["local_input_pdb"] = self.top_local_pdb_var.get().strip()
        return data

    def open_config(self):
        path = filedialog.askopenfilename(
            title="Open workflow config",
            initialdir=str(SCRIPT_DIR / "configs"),
            filetypes=[("JSON config", "*.json"), ("All files", "*.*")],
        )
        if not path:
            return
        try:
            self.config_path = Path(path)
            self.config_data = generate_workflow.load_config(self.config_path)
            self.config_path_var.set(str(self.config_path))
            self._load_values(self.config_data)
            self._update_parameter_summary()
            self._load_result_scan_vars()
            self._refresh_stage_selector()
            self._log(f"Loaded config: {self.config_path}")
        except Exception as exc:
            self._show_error("Failed to open config", exc)

    def save_config(self):
        try:
            data = self._collect_config()
            path = Path(self.config_path_var.get()).expanduser()
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_text(json.dumps(data, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
            self.config_path = path
            self.config_data = data
            self._update_parameter_summary()
            self._update_current_project_banner()
            self._log(f"Saved config: {path}")
        except Exception as exc:
            self._show_error("Failed to save config", exc)

    def save_config_as(self):
        path = filedialog.asksaveasfilename(
            title="Save workflow config",
            initialdir=str(SCRIPT_DIR / "configs"),
            initialfile=f"{self.variables['target'][0].get() or 'NEW_TARGET'}.json",
            defaultextension=".json",
            filetypes=[("JSON config", "*.json"), ("All files", "*.*")],
        )
        if not path:
            return
        self.config_path_var.set(path)
        self.save_config()

    def choose_output_dir(self):
        path = filedialog.askdirectory(
            title="Choose workflow output directory",
            initialdir=str(Path(self.output_dir_var.get()).expanduser()),
        )
        if path:
            self.output_dir_var.set(path)

    def choose_local_input_pdb(self):
        path = filedialog.askopenfilename(
            title="Choose local input PDB",
            initialdir=str(SCRIPT_DIR),
            filetypes=[("PDB files", "*.pdb *.ent"), ("All files", "*.*")],
        )
        if not path:
            return
        source_path = Path(path)
        self._set_workflow_local_pdb(source_path)
        self._set_preprocess_source(source_path, register_original=True)
        self._log(f"Selected local PDB: {path}")

    def preview_top_local_pdb(self):
        path_text = self.top_local_pdb_var.get().strip()
        if not path_text:
            self.choose_local_input_pdb()
            path_text = self.top_local_pdb_var.get().strip()
            if not path_text:
                return
        self.pdb_raw_path_var.set(path_text)
        self._show_page("PDB Preprocess")
        self.preview_preprocess_pdb()

    def choose_preprocess_pdb(self):
        path = filedialog.askopenfilename(
            title="Choose raw PDB for preprocessing",
            initialdir=str(SCRIPT_DIR),
            filetypes=[("PDB files", "*.pdb *.ent"), ("All files", "*.*")],
        )
        if not path:
            return
        self._set_preprocess_source(Path(path), register_original=True)
        self.preview_preprocess_pdb()

    def fetch_rcsb_pdb(self):
        pdb_id = self.pdb_rcsb_id_var.get().strip()
        try:
            normalized_id = pdb_preprocess.normalize_pdb_id(pdb_id)
        except Exception as exc:
            self._show_error("Invalid PDB ID", exc)
            return
        destination = SCRIPT_DIR / "inputs_preprocessed" / "rcsb"
        self.pdb_preprocess_status_var.set(f"Downloading {normalized_id} from RCSB PDB ...")

        def on_success(result):
            path = Path(result["path"])
            self._set_preprocess_source(path, register_original=True)
            self.pdb_preprocess_status_var.set(f"Downloaded {normalized_id}: {path}")
            self._log(f"Fetched RCSB PDB {normalized_id} from {result['url']} -> {path}")
            self.preview_preprocess_pdb()

        self._run_background(
            f"Fetch RCSB PDB {normalized_id}",
            lambda: pdb_preprocess.fetch_rcsb_pdb(normalized_id, destination),
            on_success,
            lambda exc: self.pdb_preprocess_status_var.set(f"RCSB fetch failed: {exc}"),
        )

    def _set_preprocess_source(self, path: Path, register_original: bool):
        source_path = Path(path).resolve()
        if register_original:
            original_path = pdb_preprocess.snapshot_original_pdb(
                source_path,
                SCRIPT_DIR / "inputs_preprocessed" / "originals",
            )
            self.pdb_original_path_var.set(str(original_path))
        self.pdb_raw_path_var.set(str(source_path))
        self.pdb_clean_output_var.set(str(self._default_clean_pdb_output_path(source_path)))

    def restore_original_pdb(self):
        try:
            original_text = self.pdb_original_path_var.get().strip()
            if not original_text:
                raise ValueError("No original PDB snapshot is available")
            original_path = self._resolve_local_config_path(original_text)
            if not original_path.is_file():
                raise FileNotFoundError(f"Original PDB snapshot not found: {original_path}")
            self.pdb_raw_path_var.set(str(original_path))
            self.pdb_clean_output_var.set(str(self._default_clean_pdb_output_path(original_path)))
            self._set_workflow_local_pdb(original_path)
            self.preview_preprocess_pdb()
            self.pdb_preprocess_status_var.set(f"Restored original PDB: {original_path.name}")
            self._log(f"Restored original PDB snapshot: {original_path}")
        except Exception as exc:
            self._show_error("Failed to restore original PDB", exc)

    def choose_clean_pdb_output(self):
        current = self.pdb_clean_output_var.get().strip()
        initialdir = str(Path(current).expanduser().parent) if current else str(SCRIPT_DIR / "inputs_preprocessed")
        path = filedialog.asksaveasfilename(
            title="Save cleaned PDB as",
            initialdir=initialdir,
            initialfile=Path(current).name if current else "cleaned_input.pdb",
            defaultextension=".pdb",
            filetypes=[("PDB files", "*.pdb"), ("All files", "*.*")],
        )
        if path:
            self.pdb_clean_output_var.set(path)

    def preview_preprocess_pdb(self):
        try:
            path = self._preprocess_input_path()
            if not self.pdb_original_path_var.get().strip():
                original_path = pdb_preprocess.snapshot_original_pdb(
                    path,
                    SCRIPT_DIR / "inputs_preprocessed" / "originals",
                )
                self.pdb_original_path_var.set(str(original_path))
            if not self.pdb_clean_output_var.get().strip():
                self.pdb_clean_output_var.set(str(self._default_clean_pdb_output_path(path)))
            summary = pdb_preprocess.parse_pdb(path)
            self._populate_pdb_preview(summary)
            self.pdb_preprocess_status_var.set(
                f"Previewed {path.name}: {len(summary['protein_residues'])} protein residues, "
                f"{summary['nonprotein_atoms']} non-protein atoms"
            )
            self._log(f"Previewed PDB: {path}")
        except Exception as exc:
            self._show_error("Failed to preview PDB", exc)

    def use_preprocess_pdb_as_input(self):
        try:
            path = self._preprocess_input_path()
            self._set_workflow_local_pdb(path)
            self.pdb_preprocess_status_var.set(f"Workflow local input set to: {path}")
            self._log(f"Workflow local input PDB set to: {path}")
        except Exception as exc:
            self._show_error("Failed to set workflow input PDB", exc)

    def clean_preprocess_pdb(self):
        try:
            input_path = self._preprocess_input_path()
            output_text = self.pdb_clean_output_var.get().strip()
            output_path = self._resolve_local_config_path(output_text) if output_text else self._default_clean_pdb_output_path(input_path)
            if input_path == output_path:
                raise ValueError("Cleaned output must be different from the current input PDB")
            selected_chains = self._selected_preprocess_chains()
            if self.pdb_chain_ids and not selected_chains:
                raise ValueError("Select at least one chain to keep")
            if self.pdb_renumber_residues_var.get():
                if not messagebox.askyesno(
                    "Confirm residue renumbering",
                    "Continuous residue renumbering changes residue IDs. RFdiffusion contigs must match the cleaned PDB.\n\nContinue?",
                ):
                    return
            result = pdb_preprocess.clean_pdb(
                input_path=input_path,
                output_path=output_path,
                keep_protein_only=bool(self.pdb_keep_protein_only_var.get()),
                renumber_atoms=bool(self.pdb_renumber_atoms_var.get()),
                renumber_residues=bool(self.pdb_renumber_residues_var.get()),
                keep_chains=selected_chains if self.pdb_chain_ids else None,
                renumber_chains=bool(self.pdb_renumber_chains_var.get()),
            )
            self.pdb_clean_output_var.set(str(output_path))
            self.pdb_raw_path_var.set(str(output_path))
            summary = pdb_preprocess.parse_pdb(output_path)
            self._populate_pdb_preview(summary)
            if self.pdb_set_as_input_var.get():
                self._set_workflow_local_pdb(output_path)
            self.pdb_preprocess_status_var.set(
                f"Cleaned PDB: kept {result['kept_atoms']} atoms in {', '.join(result['kept_chains']) or 'all chains'}, "
                f"chain map {self._format_chain_map(result['chain_map'])}, removed {result['removed_atoms']} atoms -> {output_path.name}"
            )
            self._log(f"Cleaned PDB written: {output_path}")
        except Exception as exc:
            self._show_error("Failed to clean PDB", exc)

    def _preprocess_input_path(self) -> Path:
        path_text = self.pdb_raw_path_var.get().strip()
        if not path_text:
            raise ValueError("Select a raw PDB file first")
        path = self._resolve_local_config_path(path_text)
        if not path.is_file():
            raise FileNotFoundError(f"PDB file not found: {path}")
        return path

    def _resolve_local_config_path(self, path_text: str) -> Path:
        path = Path(path_text).expanduser()
        if not path.is_absolute():
            path = SCRIPT_DIR / path
        return path.resolve()

    def _default_clean_pdb_output_path(self, input_path: Path) -> Path:
        return SCRIPT_DIR / "inputs_preprocessed" / f"{Path(input_path).stem}_protein_clean.pdb"

    def _config_path_for_local_file(self, path: Path) -> str:
        resolved = Path(path).resolve()
        try:
            return resolved.relative_to(SCRIPT_DIR).as_posix()
        except ValueError:
            return str(resolved)

    def _set_workflow_local_pdb(self, path: Path):
        path = Path(path).resolve()
        self.top_local_pdb_var.set(self._config_path_for_local_file(path))
        self._update_parameter_summary()

    def _populate_pdb_preview(self, summary: dict):
        self.pdb_summary_text.delete("1.0", END)
        self.pdb_summary_text.insert(END, pdb_preprocess.summary_text(summary))
        self.pdb_component_text.delete("1.0", END)
        self.pdb_component_text.insert(END, pdb_preprocess.component_summary_text(summary))
        previous_selected = {
            self.pdb_chain_ids[index]
            for index in self.pdb_chain_listbox.curselection()
            if index < len(self.pdb_chain_ids)
        }
        self.pdb_chain_listbox.delete(0, END)
        self.pdb_chain_info = summary.get("chains", {})
        self.pdb_chain_ids = sorted(self.pdb_chain_info)
        for index, chain in enumerate(self.pdb_chain_ids):
            info = self.pdb_chain_info[chain]
            references = ", ".join(info.get("dbrefs", [])) or "no database reference"
            label = (
                f"Chain {chain} | {info.get('description') or 'Unspecified molecule'} | "
                f"protein={info.get('protein_residues', 0)}, other={info.get('nonprotein_residues', 0)}, atoms={info.get('atoms', 0)} | "
                f"{references}"
            )
            self.pdb_chain_listbox.insert(END, label)
            if not previous_selected or chain in previous_selected:
                self.pdb_chain_listbox.selection_set(index)
        for item in self.pdb_residue_tree.get_children():
            self.pdb_residue_tree.delete(item)
        for index, residue in enumerate(summary["protein_residues"]):
            self.pdb_residue_tree.insert(
                "",
                END,
                iid=str(index),
                values=(
                    residue.chain,
                    residue.resseq,
                    residue.icode,
                    residue.resname,
                    residue.aa,
                    residue.atom_count,
                ),
            )

    def _selected_preprocess_chains(self) -> set[str]:
        return {
            self.pdb_chain_ids[index]
            for index in self.pdb_chain_listbox.curselection()
            if index < len(self.pdb_chain_ids)
        }

    def select_all_pdb_chains(self):
        self.pdb_chain_listbox.selection_set(0, END)

    def select_protein_pdb_chains(self):
        self.pdb_chain_listbox.selection_clear(0, END)
        for index, chain in enumerate(self.pdb_chain_ids):
            if self.pdb_chain_info.get(chain, {}).get("protein_residues", 0):
                self.pdb_chain_listbox.selection_set(index)

    def _format_chain_map(self, chain_map: dict[str, str]) -> str:
        changes = [f"{old}->{new}" for old, new in chain_map.items() if old != new]
        return ", ".join(changes) if changes else "unchanged"

    def generate(self):
        try:
            data = self._collect_config()
            output_dir = Path(self.output_dir_var.get()).expanduser()
            written = generate_workflow.generate(data, output_dir, force=True)
            self.config_data = data
            self.output_dir = output_dir
            self._update_parameter_summary()
            self._update_current_project_banner()
            self._log("Generated workflow files:")
            for path in written:
                self._log(f"  {path}")
            messagebox.showinfo("Workflow generated", f"Generated {len(written)} submit scripts.")
        except Exception as exc:
            self._show_error("Failed to generate workflow", exc)

    def _apply_cluster_settings_to_profile(self, profile: dict) -> dict:
        if not hasattr(self, "cluster_host_var"):
            return profile
        profile["host"] = self.cluster_host_var.get().strip()
        profile["user"] = self.cluster_user_var.get().strip()
        port_text = self.cluster_port_var.get().strip()
        if port_text:
            try:
                profile["port"] = int(port_text)
            except ValueError as exc:
                raise ValueError("Cluster port must be an integer") from exc
        else:
            profile.pop("port", None)
        profile["auth_method"] = self.auth_method_var.get().strip() or "ssh_key"
        profile["remember_password"] = bool(self.remember_password_var.get())
        for key, variable in [
            ("ssh_key", self.cluster_ssh_key_var),
            ("remote_upload_parent", self.cluster_remote_upload_parent_var),
            ("remote_workflow_dir", self.cluster_remote_workflow_dir_var),
            ("local_download_dir", self.cluster_local_download_dir_var),
            ("pymol_executable", self.cluster_pymol_var),
            ("bjobs_user", self.cluster_bjobs_user_var),
        ]:
            value = variable.get().strip()
            if value:
                profile[key] = value
            else:
                profile.pop(key, None)
        for key, variable in [
            ("scratch_scan_max_depth", self.cluster_scan_depth_var),
            ("job_log_scan_max_depth", self.cluster_log_depth_var),
        ]:
            value = variable.get().strip()
            if value:
                try:
                    profile[key] = int(value)
                except ValueError as exc:
                    raise ValueError(f"{key} must be an integer") from exc
            else:
                profile.pop(key, None)
        patterns = self._split_profile_list(self.cluster_result_patterns_var.get())
        if patterns:
            profile["result_file_patterns"] = patterns
        else:
            profile.pop("result_file_patterns", None)
        extra_args = shlex.split(self.cluster_ssh_extra_args_var.get().strip()) if self.cluster_ssh_extra_args_var.get().strip() else []
        profile["ssh_extra_args"] = extra_args
        password = self.password_var.get()
        if profile["auth_method"] == "password" and password:
            profile["password"] = password
        else:
            profile.pop("password", None)
        return profile

    def _split_profile_list(self, text: str) -> list[str]:
        return [item.strip() for item in text.replace(",", " ").split() if item.strip()]

    def _collect_credential_profile(self) -> dict:
        profile = deepcopy(self.cluster_profile)
        return self._apply_cluster_settings_to_profile(profile)

    def _collect_cluster_profile(self) -> dict:
        profile = deepcopy(self.cluster_profile)
        profile = self._apply_cluster_settings_to_profile(profile)
        data = self._collect_config()
        profile["stage_scripts"] = self._stage_scripts_for_config(data)
        profile["remote_workflow_dir"] = self._remote_workflow_dir_for_config(profile, data)
        profile.setdefault("local_download_dir", f"downloads/{data.get('target', 'target')}")
        if hasattr(self, "scratch_root_var"):
            profile["scratch_migration"] = {
                "scratch_root": self.scratch_root_var.get().strip() or "/scratch",
                "source_date": self.scratch_source_date_var.get().strip(),
                "target_date": self.scratch_target_date_var.get().strip() or "auto",
                "user_dir": self.scratch_user_dir_var.get().strip(),
            }
        if hasattr(self, "result_scan_date_var"):
            scan_root = self._remote_result_scan_root()
            if scan_root:
                profile["scratch_scan_roots"] = [scan_root]
                profile["job_log_scan_roots"] = [scan_root]
        return profile

    def _stage_scripts_for_config(self, data: dict) -> dict:
        target = str(data.get("target") or "TARGET")
        pilot = str(data.get("pilot") or "pilot0")
        relax_cycles = data.get("proteinmpnn", {}).get("relax_cycles") or 4
        return {
            "RFDiffusion": f"1.RFDiffusion/submit_{target}_rfdiffusion_{pilot}.sh",
            "ProteinMPNN": f"2.ProteinMPNN/submit_{target}_mpnn_relax{relax_cycles}_shards.sh",
            "AfCycDesign": "3.AfCycDesign/submit_afcyc_shards_integrated.sh",
            "PyRosetta": "4.PyRosetta/submit_pyrosetta_shards.sh",
        }

    def _remote_workflow_dir_for_config(self, profile: dict, data: dict) -> str:
        workflow_name = f"{data.get('target')}_workflow"
        parent = str(profile.get("remote_upload_parent") or "").strip().rstrip("/")
        if parent:
            return f"{parent}/{workflow_name}"
        remote_workflow_dir = str(profile.get("remote_workflow_dir") or "").strip().rstrip("/")
        if remote_workflow_dir:
            if "/" not in remote_workflow_dir:
                return workflow_name
            parent = remote_workflow_dir.rsplit("/", 1)[0]
            return f"{parent}/{workflow_name}"
        return workflow_name

    def open_cluster_profile(self):
        path = filedialog.askopenfilename(
            title="Open cluster profile",
            initialdir=str(SCRIPT_DIR),
            filetypes=[("JSON profile", "*.json"), ("All files", "*.*")],
        )
        if not path:
            return
        try:
            self.cluster_profile_path = Path(path)
            self.cluster_profile = cluster_ops.load_profile(self.cluster_profile_path)
            self.cluster_profile_path_var.set(str(self.cluster_profile_path))
            if hasattr(self, "password_var"):
                self.password_var.set("")
            self._load_cluster_settings_vars()
            self._refresh_stage_selector()
            self._load_scratch_migration_vars()
            self._load_result_scan_vars()
            if hasattr(self, "cluster_settings_status_var"):
                self.cluster_settings_status_var.set(f"Loaded cluster profile: {self.cluster_profile_path}")
            self._log(f"Loaded cluster profile: {self.cluster_profile_path}")
        except Exception as exc:
            self._show_error("Failed to open cluster profile", exc)

    def save_cluster_profile(self):
        try:
            path = Path(self.cluster_profile_path_var.get()).expanduser()
            self.cluster_profile = self._collect_cluster_profile()
            profile_to_save = deepcopy(self.cluster_profile)
            profile_to_save.pop("password", None)
            cluster_ops.save_profile(path, profile_to_save)
            self.cluster_profile_path = path
            self._save_remembered_password_if_requested()
            self._load_cluster_settings_vars()
            if hasattr(self, "cluster_settings_status_var"):
                self.cluster_settings_status_var.set(f"Saved cluster profile: {path}")
            self._log(f"Saved cluster profile: {path}")
        except Exception as exc:
            self._show_error("Failed to save cluster profile", exc)

    def save_cluster_profile_as(self):
        path = filedialog.asksaveasfilename(
            title="Save cluster profile",
            initialdir=str(SCRIPT_DIR),
            initialfile="cluster_profile.local.json",
            defaultextension=".json",
            filetypes=[("JSON profile", "*.json"), ("All files", "*.*")],
        )
        if not path:
            return
        self.cluster_profile_path_var.set(path)
        self.save_cluster_profile()

    def _refresh_stage_selector(self):
        try:
            stage_names = list(self._collect_cluster_profile().get("stage_scripts", {}).keys())
        except Exception:
            stage_names = list(self.cluster_profile.get("stage_scripts", {}).keys())
        if hasattr(self, "stage_box"):
            self.stage_box.configure(values=stage_names)
        if stage_names and self.stage_var.get() not in stage_names:
            self.stage_var.set(stage_names[0])

    def _current_local_workflow_dir(self) -> Path:
        data = self._collect_config()
        return Path(self.output_dir_var.get()).expanduser() / f"{data['target']}_workflow"

    def _update_current_project_banner(self):
        if not hasattr(self, "current_project_var"):
            return
        try:
            def field_value(key: str) -> str:
                item = self.variables.get(key)
                return item[0].get().strip() if item else ""

            target = field_value("target") or "未设置"
            pilot = field_value("pilot") or "未设置"
            scratch_date = field_value("scratch_date") or "未设置"
            project_dir = field_value("project_dir_name") or "未设置"
            config_name = Path(self.config_path_var.get()).name if hasattr(self, "config_path_var") else Path(self.config_path).name
            self.current_project_var.set(
                f"Target={target}  |  Pilot={pilot}  |  Config={config_name}  |  Scratch={scratch_date}/{project_dir}"
            )
        except Exception:
            self.current_project_var.set("Current project: unable to read config fields")

    def _update_parameter_summary(self):
        try:
            data = self._collect_config()
            lines = [
                f"Target: {data.get('target')}    Pilot: {data.get('pilot')}",
                f"Remarks: {data.get('remarks') or '-'}",
                f"Input PDB: {data.get('input_pdb')}",
                f"Local input PDB: {data.get('local_input_pdb')}",
                f"Chains: target={data.get('target_chain')} binder={data.get('binder_chain')}",
                f"Contigs: {data.get('contigs')}",
                f"Hotspots: {data['rfdiffusion'].get('hotspot_res') or 'none'}",
                f"Scratch: /scratch/{data.get('scratch_date')}/{data.get('cluster_user')}/{data.get('project_dir_name')}",
                f"RFDiffusion: {data['rfdiffusion'].get('n_shards')} shards × {data['rfdiffusion'].get('designs_per_shard')} designs, GPU queue={data['rfdiffusion'].get('queue')}, GPU cores={data['rfdiffusion'].get('gpu_ncpu')}",
                f"ProteinMPNN: {data['proteinmpnn'].get('n_shards')} shards, seqs={data['proteinmpnn'].get('seqs_per_struct')}, relax={data['proteinmpnn'].get('relax_cycles')}, LSF -n={data['proteinmpnn'].get('ncpu')}, LSF -R={data['proteinmpnn'].get('resource_req') or 'omitted'}, queue={data['proteinmpnn'].get('queue')}",
                f"AfCycDesign: {data['afcyc'].get('n_shards')} shards, GPU queue={data['afcyc'].get('gpu_queue')}, GPU cores={data['afcyc'].get('gpu_ncpu')}, CPU queue={data['afcyc'].get('cpu_queue')}",
                f"PyRosetta: {data['pyrosetta'].get('n_shards')} shards, LSF -n={data['pyrosetta'].get('ncpu')}, LSF -R={data['pyrosetta'].get('resource_req') or 'omitted'}, queue={data['pyrosetta'].get('queue')}",
            ]
            self.parameter_summary.delete("1.0", END)
            self.parameter_summary.insert(END, "\n".join(lines))
        except Exception:
            pass

    def _run_background(self, description: str, worker, on_success=None, on_error=None):
        self._set_status(f"Running: {description}")
        self._log(f"Started: {description}")

        def run():
            try:
                result = worker()
                self.root.after(
                    0,
                    lambda desc=description, res=result, callback=on_success: self._finish_background(desc, res, callback),
                )
            except Exception as exc:
                self.root.after(
                    0,
                    lambda desc=description, error=exc, callback=on_error: self._fail_background(desc, error, callback),
                )

        threading.Thread(target=run, daemon=True).start()

    def _finish_background(self, description: str, result, on_success=None):
        if isinstance(result, cluster_ops.CommandResult):
            self._record_command_result(description, result)
            self._log(f"Finished: {description} (exit {result.returncode})")
            if result.stdout.strip():
                self._log(result.stdout.strip())
            if result.stderr.strip():
                self._log(result.stderr.strip())
        elif isinstance(result, tuple) and result and isinstance(result[0], cluster_ops.CommandResult):
            self._record_command_result(description, result[0])
            self._log(f"Finished: {description} (exit {result[0].returncode})")
        elif isinstance(result, list) and result and all(isinstance(item, cluster_ops.CommandResult) for item in result):
            for index, item in enumerate(result, 1):
                self._record_command_result(f"{description} #{index}", item)
            self._log(f"Finished: {description} ({len(result)} commands)")
        else:
            self._log(f"Finished: {description}")
        self._set_status(f"Finished: {description}")
        if on_success:
            on_success(result)

    def _fail_background(self, description: str, error: Exception, on_error=None):
        self._show_error(description, error)
        if on_error:
            on_error(error)

    def test_cluster_connection(self):
        def on_success(result):
            self._save_remembered_password_if_requested()
            if result.returncode == 0:
                summary = result.stdout.strip() or "Connection succeeded."
                messagebox.showinfo("Connection succeeded", summary[:1500])
                self._set_status("Connection succeeded")
            else:
                summary = (result.stderr.strip() or result.stdout.strip() or "Connection failed.")
                messagebox.showerror("Connection failed", summary[:1500])
                self._set_status("Connection failed")

        self._run_background("Test cluster connection", lambda: cluster_ops.test_connection(self._collect_cluster_profile()), on_success)

    def upload_workflow(self):
        data = self._collect_config()
        output_dir = Path(self.output_dir_var.get()).expanduser()
        profile = self._collect_cluster_profile()

        def worker():
            generate_workflow.generate(data, output_dir, force=True)
            local_dir = output_dir / f"{data['target']}_workflow"
            return cluster_ops.upload_workflow(profile, local_dir)

        self.config_data = data
        self.output_dir = output_dir
        self._update_parameter_summary()
        self._run_background("Generate and upload workflow", worker)

    def submit_selected_stage(self):
        stage_name = self.stage_var.get()
        data = self._collect_config()
        output_dir = Path(self.output_dir_var.get()).expanduser()
        profile = self._collect_cluster_profile()
        script = profile.get("stage_scripts", {}).get(stage_name, "")
        local_dir = output_dir / f"{data['target']}_workflow"
        self.config_data = data
        self.output_dir = output_dir
        self._update_parameter_summary()
        self._update_current_project_banner()
        self._log(
            f"Sync and submit current config: target={data.get('target')} pilot={data.get('pilot')} "
            f"-> {profile.get('remote_workflow_dir')} :: {stage_name} -> {script}"
        )

        def worker():
            generate_workflow.generate(data, output_dir, force=True)
            upload_result = cluster_ops.upload_workflow(profile, local_dir)
            if upload_result.returncode != 0:
                return [upload_result]
            submit_result = cluster_ops.submit_stage(profile, stage_name)
            return [upload_result, submit_result]

        def on_success(results):
            if isinstance(results, list) and results:
                if len(results) == 1 and results[0].returncode != 0:
                    self._set_status(f"Upload failed; {stage_name} was not submitted")
                    return
                if results[-1].returncode != 0:
                    self._set_status(f"Submit failed: {stage_name}")
                    return
            self.refresh_jobs()

        self._run_background(f"Sync and submit {stage_name}", worker, on_success)

    def refresh_jobs(self):
        if hasattr(self, "job_summary_var"):
            self.job_summary_var.set("Running: refreshing jobs from cluster...")

        def worker():
            return cluster_ops.list_jobs(self._collect_cluster_profile())

        def on_success(result):
            command_result, jobs = result
            self._log(f"Refresh jobs exit {command_result.returncode}; {len(jobs)} jobs visible.")
            self._populate_jobs(jobs)
            if not jobs:
                detail = (command_result.stdout.strip() or command_result.stderr.strip())
                if detail:
                    self._log("No jobs were parsed from bjobs output; check Raw SSH/LSF output.")
                    self.job_summary_var.set("Jobs loaded: no parsed jobs; check Raw SSH/LSF output for bjobs text.")
                else:
                    self._log("bjobs returned no visible jobs for the configured user.")
                    self.job_summary_var.set("Jobs loaded: no visible jobs for the configured user.")

        self._run_background("Refresh jobs", worker, on_success)

    def view_selected_job_details(self):
        job_ids = self._selected_job_ids()
        if not job_ids:
            messagebox.showinfo("No jobs selected", "Select one or more jobs first.")
            return

        def on_success(result):
            content = result.stdout.strip() or result.stderr.strip() or "No bjobs -l output returned."
            self._show_text_window("bjobs -l details", content)

        self._run_background(
            "View bjobs -l",
            lambda: cluster_ops.describe_jobs(self._collect_cluster_profile(), job_ids),
            on_success,
        )

    def query_queue_info(self, mode: str):
        queue_name = self.queue_name_var.get().strip()
        if mode == "detail" and not queue_name:
            messagebox.showinfo("Queue name required", "Enter a queue name before running queueinfo -l.")
            return
        label = "queueinfo"
        if mode == "gpu":
            label = "queueinfo -gpu"
        elif mode == "detail":
            label = f"queueinfo -l {queue_name}"
        self.queue_status_var.set(f"Running: {label}")

        def on_success(result):
            content = result.stdout.strip() or result.stderr.strip() or f"{label} returned no output."
            if result.returncode == 0:
                self.queue_status_var.set(f"Finished: {label}")
            else:
                self.queue_status_var.set(f"Failed: {label} exit {result.returncode}")
            self._show_text_window(label, content)

        self._run_background(
            label,
            lambda: cluster_ops.queue_info(self._collect_cluster_profile(), mode=mode, queue_name=queue_name),
            on_success,
            lambda exc: self.queue_status_var.set(f"Failed: {label}: {exc}"),
        )

    def scan_results(self):
        scan_root = self._remote_result_scan_root()
        if not scan_root:
            messagebox.showinfo("Scan root required", "Set scratch date and user folder first.")
            return
        self.result_scan_status_var.set(f"Scanning results: {scan_root}")

        def worker():
            return cluster_ops.scan_results(self._collect_cluster_profile())

        def on_success(result):
            command_result, entries = result
            self._log(f"Result scan exit {command_result.returncode}; {len(entries)} files found.")
            if entries:
                self.result_scan_status_var.set(f"Scan complete: {len(entries)} files")
            else:
                self.result_scan_status_var.set("Scan complete: 0 files")
                self._log(f"No result files found under {scan_root}.")
            self.result_entries = entries
            self._refresh_result_filter_values()
            self._apply_result_filter()

        self._run_background("Scan scratch results", worker, on_success)

    def scan_job_logs(self):
        scan_root = self._remote_result_scan_root()
        if not scan_root:
            messagebox.showinfo("Scan root required", "Set scratch date and user folder first.")
            return
        self.result_scan_status_var.set(f"Scanning logs: {scan_root}")

        def worker():
            return cluster_ops.scan_job_logs(self._collect_cluster_profile())

        def on_success(result):
            command_result, entries = result
            self._log(f"Job log scan exit {command_result.returncode}; {len(entries)} .out/.err files found.")
            if entries:
                self.result_scan_status_var.set(f"Log scan complete: {len(entries)} files")
            else:
                self.result_scan_status_var.set("Log scan complete: 0 files")
                self._log(f"No .out/.err files found under {scan_root}.")
            self.result_entries = entries
            self._refresh_result_filter_values()
            self._apply_result_filter()

        self._run_background("Scan job logs/errors", worker, on_success)

    def scan_data_merged_csvs(self):
        scan_root = self._remote_result_scan_root()
        if not scan_root:
            messagebox.showinfo("Scan root required", "Set scratch date and user folder first.")
            return
        self.data_status_var.set(f"Scanning merged CSVs: {scan_root}")

        def worker():
            return cluster_ops.scan_merged_csvs(self._collect_cluster_profile())

        def on_success(result):
            command_result, entries = result
            self.data_csv_all_entries = entries
            self._refresh_data_csv_filter_values()
            self._apply_data_csv_filter()
            self.data_status_var.set(f"Found {len(entries)} merged CSV files")
            self._log(f"Merged CSV scan exit {command_result.returncode}; {len(entries)} files found.")

        self._run_background("Scan merged CSV files", worker, on_success)

    def _refresh_data_csv_filter_values(self):
        if not hasattr(self, "data_filter_boxes"):
            return
        for key, box in self.data_filter_boxes.items():
            values = sorted({str(entry.get(key) or "unknown") for entry in self.data_csv_all_entries})
            selected = getattr(self, f"data_{key}_filter_var").get()
            box.configure(values=["All"] + values)
            if selected not in ["All"] + values:
                getattr(self, f"data_{key}_filter_var").set("All")

    def _entry_matches_data_filters(self, entry: dict) -> bool:
        filters = {
            "target": self.data_target_filter_var.get(),
            "pilot": self.data_pilot_filter_var.get(),
            "stage": self.data_stage_filter_var.get(),
            "shard": self.data_shard_filter_var.get(),
        }
        for key, selected in filters.items():
            if selected != "All" and str(entry.get(key) or "unknown") != selected:
                return False
        query = self.data_csv_search_var.get().strip().lower()
        if query:
            haystack = " ".join(str(entry.get(key, "")) for key in ("target", "pilot", "stage", "shard", "name", "path")).lower()
            if query not in haystack:
                return False
        return True

    def _apply_data_csv_filter(self):
        if not hasattr(self, "data_csv_tree"):
            return
        self.data_csv_entries = [
            entry
            for entry in self.data_csv_all_entries
            if self._entry_matches_data_filters(entry)
        ]
        self._populate_data_csv_tree(self.data_csv_entries)
        if hasattr(self, "data_status_var") and self.data_csv_all_entries:
            self.data_status_var.set(f"Merged CSV visible: {len(self.data_csv_entries)} / scanned: {len(self.data_csv_all_entries)}")

    def _populate_data_csv_tree(self, entries: list[dict]):
        for item in self.data_csv_tree.get_children():
            self.data_csv_tree.delete(item)
        for index, entry in enumerate(entries):
            self.data_csv_tree.insert(
                "",
                END,
                iid=str(index),
                values=(
                    entry.get("stage", ""),
                    entry.get("size", 0),
                    entry.get("mtime", ""),
                    entry.get("name", ""),
                    entry.get("path", ""),
                ),
            )

    def _selected_data_csv_entry(self) -> dict | None:
        selected = list(self.data_csv_tree.selection())
        if len(selected) != 1:
            messagebox.showinfo("Select one CSV", "Select exactly one merged CSV first.")
            return None
        return self.data_csv_entries[int(selected[0])]

    def preview_data_csv(self):
        entry = self._selected_data_csv_entry()
        if not entry:
            return
        profile = self._collect_cluster_profile()
        download_root = self._result_download_root([entry])
        local_path = cluster_ops.local_download_path_for_entry(profile, entry, download_root, layout="stage")
        if local_path.exists():
            self.data_status_var.set(f"Using local merged CSV: {local_path}")
            self._show_csv_file(local_path)
            return

        def on_success(results):
            failures = sum(1 for result in results if result.returncode != 0)
            if local_path.exists():
                if failures:
                    self.data_status_var.set(f"Download finished with {failures} failed command(s); opened CSV: {local_path}")
                else:
                    self.data_status_var.set(f"Download finished: merged CSV preview ready -> {local_path}")
                self._show_csv_file(local_path)
            else:
                self.data_status_var.set(f"Download finished but merged CSV was not found: {local_path}")
                messagebox.showerror("CSV not downloaded", f"Expected file was not found:\n{local_path}")

        self.data_status_var.set(f"Downloading merged CSV for preview to {local_path.parent} ...")
        self._run_background(
            "Download merged CSV for preview",
            lambda: cluster_ops.download_files(profile, [entry], local_dir=download_root, layout="stage"),
            on_success,
            lambda exc: self.data_status_var.set(f"Download failed: {exc}"),
        )

    def load_data_csv_for_analysis(self):
        entry = self._selected_data_csv_entry()
        if not entry:
            return
        profile = self._collect_cluster_profile()
        download_root = self._result_download_root([entry])
        local_path = cluster_ops.local_download_path_for_entry(profile, entry, download_root, layout="stage")

        def on_success(results):
            failures = sum(1 for result in results if result.returncode != 0)
            if local_path.exists():
                if failures:
                    self.data_status_var.set(f"Download finished with {failures} failed command(s); loading CSV: {local_path}")
                else:
                    self.data_status_var.set(f"Download finished: loading merged CSV -> {local_path}")
                self._load_local_data_csv(local_path)
            else:
                self.data_status_var.set(f"Download finished but merged CSV was not found: {local_path}")
                messagebox.showerror("CSV not downloaded", f"Expected file was not found:\n{local_path}")

        refresh_note = "refreshing existing local copy" if local_path.exists() else "downloading"
        self.data_status_var.set(f"{refresh_note.capitalize()} merged CSV for analysis to {local_path.parent} ...")
        self._run_background(
            "Refresh merged CSV for analysis",
            lambda: cluster_ops.download_files(profile, [entry], local_dir=download_root, layout="stage"),
            on_success,
            lambda exc: self.data_status_var.set(f"Download failed: {exc}"),
        )

    def _load_local_data_csv(self, path: Path):
        try:
            numeric_headers, sampled_rows = self._detect_numeric_csv_headers(path)
        except Exception as exc:
            self._show_error("Failed to inspect CSV metrics", exc)
            return
        if not numeric_headers:
            messagebox.showinfo("No numeric metrics", f"No numeric columns were detected in:\n{path}")
            return
        self.data_loaded_csv_path = path
        self.data_loaded_headers = numeric_headers
        self.data_metric_box.configure(values=numeric_headers)
        if self.data_metric_var.get() not in numeric_headers:
            self.data_metric_var.set(numeric_headers[0])
        total_rows = self._count_csv_data_rows(path)
        self.data_status_var.set(
            f"Loaded {path.name}; rows={total_rows}, {len(numeric_headers)} numeric columns from {sampled_rows} sampled rows. "
            "Select a metric, then click Plot histogram + scatter."
        )

    def _count_csv_data_rows(self, path: Path) -> int:
        with path.open("r", encoding="utf-8-sig", errors="replace", newline="") as handle:
            reader = csv.reader(handle)
            try:
                next(reader)
            except StopIteration:
                return 0
            return sum(1 for _row in reader)

    def use_afcyc_metric_defaults(self):
        self._apply_metric_preset(AFCYC_SCORE_PRESETS, "AfCycDesign")

    def use_pyrosetta_metric_defaults(self):
        self._apply_metric_preset(PYROSETTA_SCORE_PRESETS, "PyRosetta")

    def _apply_metric_preset(self, presets: list[tuple[str, str, str]], label: str):
        if not self.data_loaded_headers:
            self.data_status_var.set(f"Load a merged CSV before applying {label} metric defaults.")
            return
        for metric, direction, percentile in presets:
            if metric in self.data_loaded_headers:
                self.data_metric_var.set(metric)
                self.data_metric_direction_var.set(direction)
                self.data_percentile_var.set(percentile)
                self.data_status_var.set(f"{label} default metric selected: {metric} ({direction} is better)")
                return
        self.data_status_var.set(f"No {label} default metric was found in the loaded CSV.")

    def _detect_numeric_csv_headers(self, path: Path, sample_limit: int = 5000) -> tuple[list[str], int]:
        with path.open("r", encoding="utf-8-sig", errors="replace", newline="") as handle:
            reader = csv.DictReader(handle)
            if not reader.fieldnames:
                return [], 0
            numeric_counts = {name: 0 for name in reader.fieldnames if name}
            sampled_rows = 0
            for row in reader:
                sampled_rows += 1
                for name in numeric_counts:
                    value = self._to_float(row.get(name, ""))
                    if value is not None:
                        numeric_counts[name] += 1
                if sampled_rows >= sample_limit:
                    break
        numeric_headers = [name for name, count in numeric_counts.items() if count > 0]
        return numeric_headers, sampled_rows

    def analyze_data_metric(self):
        if not self.data_loaded_csv_path:
            messagebox.showinfo("No CSV loaded", "Load a merged CSV for analysis first.")
            return
        metric = self.data_metric_var.get().strip()
        if not metric:
            messagebox.showinfo("No metric selected", "Select a numeric metric first.")
            return
        try:
            higher_is_better = self.data_metric_direction_var.get() != "lower"
            percentile = max(0.0, min(100.0, float(self.data_percentile_var.get().strip() or "90")))
            bins = max(5, min(200, int(self.data_bins_var.get().strip() or "40")))
            custom_threshold = self._to_float(self.data_threshold_var.get().strip())
            manual_min = self._optional_float_entry(self.data_hist_min_var.get(), "Histogram min")
            manual_max = self._optional_float_entry(self.data_hist_max_var.get(), "Histogram max")
            values = self._read_metric_values(self.data_loaded_csv_path, metric)
            if not values:
                messagebox.showinfo("No numeric values", f"No numeric values found for metric:\n{metric}")
                return
            low, high, outside_count, range_note = self._metric_display_range(
                values,
                robust=self.data_robust_hist_var.get(),
                manual_min=manual_min,
                manual_max=manual_max,
            )
            display_values = [value for value in values if low <= value <= high]
            if not display_values:
                messagebox.showinfo("No displayed values", f"No values are inside the selected display range:\n{low:g} to {high:g}")
                return
            percentile_value = self._percentile(display_values, percentile)
            candidate_threshold = custom_threshold if custom_threshold is not None else percentile_value
            candidates, total_passing = self._read_candidate_rows(self.data_loaded_csv_path, metric, candidate_threshold, low, high, higher_is_better=higher_is_better)
            self._draw_metric_histogram(values, bins, percentile_value, custom_threshold, percentile, low, high, outside_count, range_note)
            self._draw_metric_scatter(values, percentile_value, custom_threshold, percentile, low, high, outside_count, range_note)
            self._populate_candidate_rows(candidates)
            threshold_label = "custom threshold" if custom_threshold is not None else f"P{percentile:g}"
            operator = ">=" if higher_is_better else "<="
            self.data_status_var.set(
                f"{metric}: n={len(values)}, displayed={len(display_values)}, P{percentile:g} within range={percentile_value:.4g}, "
                f"{threshold_label}{operator}{candidate_threshold:.4g}, direction={self.data_metric_direction_var.get()}, range=[{low:.4g}, {high:.4g}], "
                f"outside={outside_count}, candidates={total_passing} (table shows top {len(candidates)}, capped at 500)"
            )
        except Exception as exc:
            self._show_error("Failed to analyze metric", exc)

    def browse_hit_csv(self, kind: str):
        initial_dir = self._default_hit_csv_browse_dir(kind)
        path = filedialog.askopenfilename(
            title="Choose merged CSV",
            initialdir=str(initial_dir),
            filetypes=[("CSV files", "*.csv"), ("All files", "*.*")],
        )
        if not path:
            return
        if kind == "afcyc":
            self.hit_afcyc_csv_var.set(path)
        else:
            self.hit_pyro_csv_var.set(path)

    def _default_hit_csv_browse_dir(self, kind: str) -> Path:
        stage = "AfCycDesign" if kind == "afcyc" else "PyRosetta"
        target = self._project_field("target") or self._download_target_for_entries(None)
        output_dir = Path(self.output_dir_var.get()).expanduser()
        browse_dir = output_dir / f"{cluster_ops.safe_local_name(target)}_workflow" / "retrieved_results" / stage
        browse_dir.mkdir(parents=True, exist_ok=True)
        return browse_dir

    def set_hit_csv_from_selected(self, kind: str):
        entry = self._selected_data_csv_entry()
        if not entry:
            return
        profile = self._collect_cluster_profile()
        download_root = self._result_download_root([entry])
        local_path = cluster_ops.local_download_path_for_entry(profile, entry, download_root, layout="stage")
        label = "AfCycDesign" if kind == "afcyc" else "PyRosetta"
        if local_path.exists():
            self._set_hit_csv_path(kind, local_path)
            self.hit_status_var.set(f"{label} CSV selected: {local_path}")
            return

        def on_success(results):
            failures = sum(1 for result in results if result.returncode != 0)
            if local_path.exists():
                self._set_hit_csv_path(kind, local_path)
                if failures:
                    self.hit_status_var.set(f"{label} CSV downloaded with {failures} failed command(s): {local_path}")
                else:
                    self.hit_status_var.set(f"{label} CSV downloaded and selected: {local_path}")
            else:
                self.hit_status_var.set(f"{label} CSV download finished but file was not found: {local_path}")

        self.hit_status_var.set(f"Downloading selected {label} merged CSV to {local_path.parent} ...")
        self._run_background(
            f"Download {label} hit-screening CSV",
            lambda: cluster_ops.download_files(profile, [entry], local_dir=download_root, layout="stage"),
            on_success,
            lambda exc: self.hit_status_var.set(f"Download failed: {exc}"),
        )

    def _set_hit_csv_path(self, kind: str, path: Path):
        if kind == "afcyc":
            self.hit_afcyc_csv_var.set(str(path))
        else:
            self.hit_pyro_csv_var.set(str(path))

    def reset_hit_thresholds(self):
        self.hit_afcyc_ipae_col_var.set(HIT_SCREEN_DEFAULTS["afcyc_ipae_column"])
        self.hit_afcyc_ipae_max_var.set(HIT_SCREEN_DEFAULTS["afcyc_ipae_max"])
        self.hit_afcyc_rmsd_col_var.set(HIT_SCREEN_DEFAULTS["afcyc_rmsd_column"])
        self.hit_afcyc_rmsd_max_var.set(HIT_SCREEN_DEFAULTS["afcyc_rmsd_max"])
        self.hit_afcyc_plddt_col_var.set(HIT_SCREEN_DEFAULTS["afcyc_plddt_column"])
        self.hit_afcyc_plddt_min_var.set(HIT_SCREEN_DEFAULTS["afcyc_plddt_min"])
        self.hit_pyro_dg_col_var.set(HIT_SCREEN_DEFAULTS["pyro_dg_column"])
        self.hit_pyro_dg_max_var.set(HIT_SCREEN_DEFAULTS["pyro_dg_max"])
        self.hit_pyro_sap_col_var.set(HIT_SCREEN_DEFAULTS["pyro_sap_column"])
        self.hit_pyro_sap_max_var.set(HIT_SCREEN_DEFAULTS["pyro_sap_max"])
        self.hit_pyro_cms_col_var.set(HIT_SCREEN_DEFAULTS["pyro_cms_column"])
        self.hit_pyro_cms_min_var.set(HIT_SCREEN_DEFAULTS["pyro_cms_min"])
        self.hit_pyro_sasa_col_var.set(HIT_SCREEN_DEFAULTS["pyro_sasa_column"])
        self.hit_pyro_sasa_min_var.set(HIT_SCREEN_DEFAULTS["pyro_sasa_min"])
        self.hit_pyro_packstat_col_var.set(HIT_SCREEN_DEFAULTS["pyro_packstat_column"])
        self.hit_pyro_packstat_min_var.set(HIT_SCREEN_DEFAULTS["pyro_packstat_min"])
        self.hit_status_var.set("Hit thresholds reset to RFpeptides-inspired defaults.")

    def run_hit_screening(self):
        try:
            afcyc_path = self._required_csv_path(self.hit_afcyc_csv_var.get(), "AfCycDesign CSV")
            pyro_path = self._required_csv_path(self.hit_pyro_csv_var.get(), "PyRosetta CSV")
            thresholds = self._collect_hit_thresholds()
            max_hits = max(1, min(5000, int(str(self.hit_limit_var.get()).strip() or "200")))
            afcyc_rows = self._load_csv_rows(afcyc_path)
            pyro_rows = self._load_csv_rows(pyro_path)
            hits, stats = self._screen_hit_rows(afcyc_rows, pyro_rows, thresholds, max_hits)
            self.hit_rows = hits
            self._populate_hit_rows(hits)
            self.hit_status_var.set(
                f"Hit screening complete: {len(hits)} shown / {stats['passing']} passing; "
                f"matched={stats['matched']}, AfCyc rows={len(afcyc_rows)}, PyRosetta rows={len(pyro_rows)}"
            )
        except Exception as exc:
            self._show_error("Failed to run hit screening", exc)

    def _required_csv_path(self, value: str, label: str) -> Path:
        text = str(value or "").strip()
        if not text:
            raise ValueError(f"{label} is required")
        path = Path(text).expanduser()
        if not path.exists():
            raise ValueError(f"{label} does not exist: {path}")
        return path

    def _load_csv_rows(self, path: Path) -> list[dict]:
        with path.open("r", encoding="utf-8-sig", errors="replace", newline="") as handle:
            return list(csv.DictReader(handle))

    def _collect_hit_thresholds(self) -> dict:
        return {
            "afcyc_ipae": (self.hit_afcyc_ipae_col_var.get().strip(), self._optional_float_entry(self.hit_afcyc_ipae_max_var.get(), "AfCyc iPAE max"), "max"),
            "afcyc_rmsd": (self.hit_afcyc_rmsd_col_var.get().strip(), self._optional_float_entry(self.hit_afcyc_rmsd_max_var.get(), "AfCyc RMSD max"), "max"),
            "afcyc_plddt": (self.hit_afcyc_plddt_col_var.get().strip(), self._optional_float_entry(self.hit_afcyc_plddt_min_var.get(), "binder pLDDT min"), "min"),
            "pyro_dg": (self.hit_pyro_dg_col_var.get().strip(), self._optional_float_entry(self.hit_pyro_dg_max_var.get(), "interface dG max"), "max"),
            "pyro_sap": (self.hit_pyro_sap_col_var.get().strip(), self._optional_float_entry(self.hit_pyro_sap_max_var.get(), "SAP max"), "max"),
            "pyro_cms": (self.hit_pyro_cms_col_var.get().strip(), self._optional_float_entry(self.hit_pyro_cms_min_var.get(), "CMS min"), "min"),
            "pyro_sasa": (self.hit_pyro_sasa_col_var.get().strip(), self._optional_float_entry(self.hit_pyro_sasa_min_var.get(), "interface SASA min"), "min"),
            "pyro_packstat": (self.hit_pyro_packstat_col_var.get().strip(), self._optional_float_entry(self.hit_pyro_packstat_min_var.get(), "packstat min"), "min"),
        }

    def _screen_hit_rows(self, afcyc_rows: list[dict], pyro_rows: list[dict], thresholds: dict, max_hits: int) -> tuple[list[dict], dict]:
        pyro_by_key = {}
        for row in pyro_rows:
            key = self._hit_key(row)
            if key:
                pyro_by_key.setdefault(key, row)

        passing = []
        matched = 0
        for afcyc_row in afcyc_rows:
            key = self._hit_key(afcyc_row)
            if not key:
                continue
            pyro_row = pyro_by_key.get(key)
            if not pyro_row:
                continue
            matched += 1
            values = {
                "afcyc_ipae": self._hit_threshold_value(afcyc_row, thresholds["afcyc_ipae"]),
                "afcyc_rmsd": self._hit_threshold_value(afcyc_row, thresholds["afcyc_rmsd"]),
                "afcyc_plddt": self._hit_threshold_value(afcyc_row, thresholds["afcyc_plddt"]),
                "pyro_dg": self._hit_threshold_value(pyro_row, thresholds["pyro_dg"]),
                "pyro_sap": self._hit_threshold_value(pyro_row, thresholds["pyro_sap"]),
                "pyro_cms": self._hit_threshold_value(pyro_row, thresholds["pyro_cms"]),
                "pyro_sasa": self._hit_threshold_value(pyro_row, thresholds["pyro_sasa"]),
                "pyro_packstat": self._hit_threshold_value(pyro_row, thresholds["pyro_packstat"]),
            }
            if not self._hit_passes_thresholds(values, thresholds):
                continue
            passing.append({
                "description": key,
                "binder_seq": afcyc_row.get("binder_seq") or pyro_row.get("binder_seq") or "",
                "i_pae": values["afcyc_ipae"],
                "rmsd": values["afcyc_rmsd"],
                "plddt": values["afcyc_plddt"],
                "interface_dG": values["pyro_dg"],
                "sap": values["pyro_sap"],
                "sap_bound": values["pyro_sap"],
                "cms": values["pyro_cms"],
                "interface_delta_sasa": values["pyro_sasa"],
                "interface_packstat": values["pyro_packstat"],
                "afcyc_pdb": self._first_nonempty(afcyc_row, ["pred_pdb", "pred_pdb_rmsd"]),
                "mpnn_pdb": self._first_nonempty(afcyc_row, ["input_pdb", "input_pdb_rmsd"]) or self._first_nonempty(pyro_row, ["pdb_path", "input_pdb"]),
                "afcyc_row": afcyc_row,
                "pyro_row": pyro_row,
            })

        passing.sort(key=self._hit_sort_key)
        return passing[:max_hits], {"matched": matched, "passing": len(passing)}

    def _hit_key(self, row: dict) -> str:
        for key in ("description", "name", "tag", "model", "design"):
            value = str(row.get(key, "")).strip()
            if value:
                return value
        for key in ("pred_pdb", "input_pdb", "pdb_path"):
            value = str(row.get(key, "")).strip()
            if value:
                return Path(value).stem.replace("_afcyc", "")
        return ""

    def _hit_threshold_value(self, row: dict, threshold_spec: tuple[str, float | None, str]) -> float | None:
        column, threshold, _mode = threshold_spec
        if threshold is None:
            return self._to_float(row.get(column, "")) if column else None
        if not column:
            raise ValueError("A threshold was set but its CSV column is blank")
        if column not in row:
            raise ValueError(f"Column not found in CSV: {column}")
        return self._to_float(row.get(column, ""))

    def _hit_passes_thresholds(self, values: dict, thresholds: dict) -> bool:
        for key, (_column, threshold, mode) in thresholds.items():
            if threshold is None:
                continue
            value = values.get(key)
            if value is None:
                return False
            if mode == "max" and value > threshold:
                return False
            if mode == "min" and value < threshold:
                return False
        return True

    def _hit_sort_key(self, hit: dict):
        return (
            self._sort_float(hit.get("interface_dG"), default=math.inf),
            self._sort_float(hit.get("i_pae"), default=math.inf),
            self._sort_float(hit.get("rmsd"), default=math.inf),
            self._sort_float(hit.get("sap"), default=math.inf),
            -self._sort_float(hit.get("cms"), default=-math.inf),
        )

    def _sort_float(self, value, default: float) -> float:
        return value if isinstance(value, (int, float)) and math.isfinite(value) else default

    def _first_nonempty(self, row: dict, keys: list[str]) -> str:
        for key in keys:
            value = str(row.get(key, "")).strip()
            if value:
                return value
        return ""

    def _populate_hit_rows(self, hits: list[dict]):
        for item in self.hit_tree.get_children():
            self.hit_tree.delete(item)
        for rank, hit in enumerate(hits, start=1):
            self.hit_tree.insert(
                "",
                END,
                iid=str(rank - 1),
                values=(
                    rank,
                    hit.get("description", ""),
                    self._format_optional_float(hit.get("i_pae")),
                    self._format_optional_float(hit.get("rmsd")),
                    self._format_optional_float(hit.get("plddt")),
                    self._format_optional_float(hit.get("interface_dG")),
                    self._format_optional_float(hit.get("sap")),
                    self._format_optional_float(hit.get("cms")),
                    hit.get("binder_seq", ""),
                    hit.get("afcyc_pdb", ""),
                    hit.get("mpnn_pdb", ""),
                ),
            )

    def _format_optional_float(self, value) -> str:
        if isinstance(value, (int, float)) and math.isfinite(value):
            return f"{value:.6g}"
        return ""

    def save_hit_screening_csv(self):
        try:
            if not self.hit_rows:
                messagebox.showinfo("No hits", "Run hit screening first.")
                return
            output_path = self._timestamped_hit_csv_path()
            output_path.parent.mkdir(parents=True, exist_ok=True)
            columns = self._hit_csv_columns()
            with output_path.open("w", encoding="utf-8", newline="") as handle:
                writer = csv.DictWriter(handle, fieldnames=columns)
                writer.writeheader()
                for rank, hit in enumerate(self.hit_rows, start=1):
                    writer.writerow(self._hit_csv_row(rank, hit, columns))
            self.hit_status_var.set(f"Hit CSV saved: {output_path}")
            self._log(f"Hit screening CSV saved: {output_path}")
            messagebox.showinfo("Hit CSV saved", f"Saved {len(self.hit_rows)} hit row(s):\n{output_path}")
        except Exception as exc:
            self.hit_status_var.set(f"Hit CSV save failed: {exc}")
            self._show_error("Failed to save hit CSV", exc)

    def _hit_csv_columns(self) -> list[str]:
        leading = ["rank", "description", "binder_seq"] + HIT_FILTER_OUTPUT_COLUMNS
        remaining_pyro = [column for column in PYROSETTA_HIT_OUTPUT_COLUMNS if column not in HIT_FILTER_OUTPUT_COLUMNS]
        trailing = ["afcyc_pdb", "mpnn_pdb"]
        return leading + remaining_pyro + trailing

    def _hit_csv_row(self, rank: int, hit: dict, columns: list[str]) -> dict:
        pyro_row = hit.get("pyro_row") or {}
        row = {
            "rank": rank,
            "description": hit.get("description", ""),
            "binder_seq": hit.get("binder_seq", ""),
            "i_pae": self._format_optional_float(hit.get("i_pae")),
            "rmsd": self._format_optional_float(hit.get("rmsd")),
            "plddt": self._format_optional_float(hit.get("plddt")),
            "interface_dG": self._format_optional_float(hit.get("interface_dG")),
            "sap_bound": self._format_optional_float(hit.get("sap_bound", hit.get("sap"))),
            "cms": self._format_optional_float(hit.get("cms")),
            "interface_delta_sasa": self._format_optional_float(hit.get("interface_delta_sasa")),
            "interface_packstat": self._format_optional_float(hit.get("interface_packstat")),
            "afcyc_pdb": hit.get("afcyc_pdb", ""),
            "mpnn_pdb": hit.get("mpnn_pdb", ""),
        }
        for column in PYROSETTA_HIT_OUTPUT_COLUMNS:
            if column not in row:
                row[column] = pyro_row.get(column, "")
        return {column: row.get(column, "") for column in columns}

    def download_hit_structures(self, selected_only: bool):
        if not self.hit_rows:
            messagebox.showinfo("No hits", "Run hit screening first.")
            return
        if selected_only:
            selected = list(self.hit_tree.selection())
            if not selected:
                messagebox.showinfo("No hits selected", "Select one or more hit rows first.")
                return
            hits = [self.hit_rows[int(item)] for item in selected]
        else:
            hits = self.hit_rows
        entries = self._hit_download_entries(hits)
        if not entries:
            messagebox.showinfo("No structures", "No downloadable AfCyc or MPNN PDB paths were found in the hit rows.")
            return
        if not messagebox.askyesno("Confirm hit download", f"Download {len(entries)} structure file(s) for {len(hits)} hit(s)?"):
            return
        profile = self._collect_cluster_profile()
        download_root = self._timestamped_hit_download_root(selected_only)

        def on_success(results):
            failures = sum(1 for result in results if result.returncode != 0)
            self._log_download_results(results, download_root)
            if failures:
                self.hit_status_var.set(f"Hit structure download finished with {failures} failed command(s): {download_root}")
            else:
                self.hit_status_var.set(f"Hit structure download finished: {len(entries)} file(s) -> {download_root}")

        self.hit_status_var.set(f"Downloading {len(entries)} hit structure file(s) to {download_root} ...")
        self._run_background(
            "Download hit structures",
            lambda: cluster_ops.download_files(profile, entries, local_dir=download_root, layout="stage"),
            on_success,
            lambda exc: self.hit_status_var.set(f"Download failed: {exc}"),
        )

    def _hit_download_entries(self, hits: list[dict]) -> list[dict]:
        entries = []
        seen = set()
        for hit in hits:
            for stage, key in (("AfCycDesign_predicted", "afcyc_pdb"), ("MPNN_input", "mpnn_pdb")):
                path = str(hit.get(key, "")).strip()
                if not path or path in seen:
                    continue
                seen.add(path)
                info = cluster_ops.classify_result_path(path)
                entries.append({
                    "path": path,
                    "name": Path(path).name,
                    "size": 0,
                    "mtime": "",
                    "target": info.get("target") or self._download_target_for_entries(None),
                    "pilot": info.get("pilot") or self._current_pilot(),
                    "stage": stage,
                    "shard": info.get("shard") or "all",
                    "type": "pdb",
                })
        return entries

    def _hit_output_root(self) -> Path:
        target = self._download_target_for_entries(None)
        output_dir = Path(self.output_dir_var.get()).expanduser()
        return output_dir / f"{cluster_ops.safe_local_name(target)}_workflow" / "retrieved_hits"

    def _timestamped_hit_download_root(self, selected_only: bool) -> Path:
        suffix = "selected" if selected_only else "all"
        return self._timestamped_hit_run_root(suffix)

    def _timestamped_hit_csv_path(self) -> Path:
        return self._timestamped_hit_run_root("csv") / "hit_screening_filtered.csv"

    def _timestamped_hit_run_root(self, suffix: str) -> Path:
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        return self._hit_output_root() / f"{timestamp}_{suffix}"

    def _to_float(self, value):
        try:
            if value is None:
                return None
            text = str(value).strip()
            if not text or text.lower() in {"nan", "none", "null"}:
                return None
            number = float(text)
            if not math.isfinite(number):
                return None
            return number
        except ValueError:
            return None

    def _optional_float_entry(self, value, label: str) -> float | None:
        text = str(value or "").strip()
        if not text:
            return None
        number = self._to_float(text)
        if number is None:
            raise ValueError(f"{label} must be a finite number")
        return number

    def _read_metric_values(self, path: Path, metric: str) -> list[float]:
        values = []
        with path.open("r", encoding="utf-8-sig", errors="replace", newline="") as handle:
            reader = csv.DictReader(handle)
            if metric not in (reader.fieldnames or []):
                raise ValueError(f"Metric not found in CSV: {metric}")
            for row in reader:
                value = self._to_float(row.get(metric, ""))
                if value is not None:
                    values.append(value)
        return values

    def _percentile(self, values: list[float], percentile: float) -> float:
        sorted_values = sorted(values)
        if len(sorted_values) == 1:
            return sorted_values[0]
        position = (len(sorted_values) - 1) * percentile / 100.0
        lower = math.floor(position)
        upper = math.ceil(position)
        if lower == upper:
            return sorted_values[int(position)]
        fraction = position - lower
        return sorted_values[lower] * (1 - fraction) + sorted_values[upper] * fraction

    def _read_candidate_rows(
        self,
        path: Path,
        metric: str,
        threshold: float,
        low: float | None = None,
        high: float | None = None,
        limit: int = 500,
        higher_is_better: bool = True,
    ) -> tuple[list[tuple], int]:
        heap = []
        total = 0
        with path.open("r", encoding="utf-8-sig", errors="replace", newline="") as handle:
            reader = csv.DictReader(handle)
            for data_row_index, row in enumerate(reader, start=2):
                value = self._to_float(row.get(metric, ""))
                if value is None:
                    continue
                if low is not None and value < low:
                    continue
                if high is not None and value > high:
                    continue
                if higher_is_better and value < threshold:
                    continue
                if not higher_is_better and value > threshold:
                    continue
                total += 1
                candidate_id = self._candidate_identifier(row)
                summary = self._candidate_summary(row, metric)
                score = value if higher_is_better else -value
                item = (score, value, data_row_index, candidate_id, summary)
                if len(heap) < limit:
                    heapq.heappush(heap, item)
                elif score > heap[0][0]:
                    heapq.heapreplace(heap, item)
        candidates = [
            (value, row_index, candidate_id, summary)
            for _score, value, row_index, candidate_id, summary in sorted(heap, key=lambda item: item[0], reverse=True)
        ]
        return candidates, total

    def _candidate_identifier(self, row: dict) -> str:
        preferred = [
            "description",
            "tag",
            "name",
            "pdb",
            "pdb_name",
            "filename",
            "file",
            "path",
            "input_pdb",
            "model",
            "design",
        ]
        lower_map = {str(key).lower(): key for key in row}
        for key in preferred:
            original = lower_map.get(key)
            if original and str(row.get(original, "")).strip():
                return str(row[original]).strip()
        for key, value in row.items():
            if str(value).strip():
                return f"{key}={value}"
        return "unknown"

    def _candidate_summary(self, row: dict, metric: str, max_fields: int = 8) -> str:
        parts = []
        for key, value in row.items():
            if key == metric or not str(value).strip():
                continue
            parts.append(f"{key}={value}")
            if len(parts) >= max_fields:
                break
        return "; ".join(parts)

    def _populate_candidate_rows(self, candidates: list[tuple]):
        for item in self.data_candidate_tree.get_children():
            self.data_candidate_tree.delete(item)
        for rank, (value, row_index, candidate_id, summary) in enumerate(candidates, start=1):
            self.data_candidate_tree.insert(
                "",
                END,
                values=(rank, row_index, f"{value:.6g}", candidate_id, summary),
            )

    def _metric_display_range(
        self,
        values: list[float],
        robust: bool,
        manual_min: float | None,
        manual_max: float | None,
    ) -> tuple[float, float, int, str]:
        sorted_values = sorted(values)
        low = sorted_values[0]
        high = sorted_values[-1]
        range_note = "full range"
        if robust and len(sorted_values) >= 10:
            low = self._percentile(sorted_values, 1)
            high = self._percentile(sorted_values, 99)
            range_note = "auto 1-99% range"
        if manual_min is not None:
            low = manual_min
            range_note = "manual range" if manual_max is not None else "manual min + auto max"
        if manual_max is not None:
            high = manual_max
            range_note = "manual range" if manual_min is not None else "auto min + manual max"
        if high <= low:
            raise ValueError(f"Histogram max must be greater than min: {high:g} <= {low:g}")
        outside_count = sum(1 for value in values if value < low or value > high)
        return low, high, outside_count, range_note

    def _draw_metric_histogram(
        self,
        values: list[float],
        bins: int,
        percentile_value: float,
        custom_threshold: float | None,
        percentile: float,
        low: float,
        high: float,
        outside_count: int,
        range_note: str,
    ):
        canvas = self.data_hist_canvas
        canvas.delete("all")
        canvas.update_idletasks()
        width, height = self._plot_canvas_size(canvas)
        left, right, top, bottom = 70, 30, 30, 50
        plot_width = width - left - right
        plot_height = height - top - bottom

        display_values = [value for value in values if low <= value <= high]
        if not display_values:
            canvas.create_text(width / 2, height / 2, text="No values inside histogram display range", fill="black")
            return

        counts = [0] * bins
        span = high - low
        for value in display_values:
            index = int((value - low) / span * bins)
            index = max(0, min(bins - 1, index))
            counts[index] += 1
        max_count = max(counts) or 1

        canvas.create_line(left, top + plot_height, left + plot_width, top + plot_height, fill="black")
        canvas.create_line(left, top, left, top + plot_height, fill="black")
        bar_width = plot_width / bins
        for index, count in enumerate(counts):
            x0 = left + index * bar_width
            x1 = left + (index + 1) * bar_width - 1
            y0 = top + plot_height - (count / max_count * plot_height)
            y1 = top + plot_height
            canvas.create_rectangle(x0, y0, x1, y1, fill="#7aa6c2", outline="")

        def draw_line(value: float, color: str, label: str):
            clamped = max(low, min(high, value))
            x = left + (clamped - low) / (high - low) * plot_width
            off_scale = value < low or value > high
            line_options = {"fill": color, "width": 2}
            if off_scale:
                line_options["dash"] = (4, 3)
            canvas.create_line(x, top, x, top + plot_height, **line_options)
            suffix = " off scale" if off_scale else ""
            text_x = min(x + 4, left + plot_width - 120)
            canvas.create_text(text_x, top + 12, text=f"{label}{suffix}", fill=color, anchor="w")

        draw_line(percentile_value, "red", f"range P{percentile:g}={percentile_value:.4g}")
        if custom_threshold is not None:
            draw_line(custom_threshold, "blue", f"threshold={custom_threshold:.4g}")

        canvas.create_text(left, height - 22, text=f"{low:.4g}", anchor="w")
        canvas.create_text(left + plot_width, height - 22, text=f"{high:.4g}", anchor="e")
        note = f"n={len(values)}, displayed={len(display_values)}, outside={outside_count}, {range_note}"
        canvas.create_text(left, 14, text=note, anchor="w", fill="black")

    def _draw_metric_scatter(
        self,
        values: list[float],
        percentile_value: float,
        custom_threshold: float | None,
        percentile: float,
        low: float,
        high: float,
        outside_count: int,
        range_note: str,
    ):
        canvas = self.data_scatter_canvas
        canvas.delete("all")
        canvas.update_idletasks()
        width, height = self._plot_canvas_size(canvas)
        left, right, top, bottom = 70, 30, 30, 50
        plot_width = width - left - right
        plot_height = height - top - bottom

        indexed_values = [(index, value) for index, value in enumerate(values) if low <= value <= high]
        if not indexed_values:
            canvas.create_text(width / 2, height / 2, text="No values inside scatter display range", fill="black")
            return

        max_points = 8000
        step = max(1, math.ceil(len(indexed_values) / max_points))
        sampled = indexed_values[::step]
        denominator = max(1, len(values) - 1)

        canvas.create_line(left, top + plot_height, left + plot_width, top + plot_height, fill="black")
        canvas.create_line(left, top, left, top + plot_height, fill="black")

        for index, value in sampled:
            x = left + index / denominator * plot_width
            y = top + plot_height - (value - low) / (high - low) * plot_height
            canvas.create_rectangle(x, y, x + 1, y + 1, fill="#3f6f9f", outline="")

        def draw_horizontal(value: float, color: str, label: str, y_offset: int):
            clamped = max(low, min(high, value))
            y = top + plot_height - (clamped - low) / (high - low) * plot_height
            off_scale = value < low or value > high
            line_options = {"fill": color, "width": 2}
            if off_scale:
                line_options["dash"] = (4, 3)
            canvas.create_line(left, y, left + plot_width, y, **line_options)
            suffix = " off scale" if off_scale else ""
            canvas.create_text(left + 4, y + y_offset, text=f"{label}{suffix}", fill=color, anchor="w")

        draw_horizontal(percentile_value, "red", f"range P{percentile:g}={percentile_value:.4g}", -10)
        if custom_threshold is not None:
            draw_horizontal(custom_threshold, "blue", f"threshold={custom_threshold:.4g}", 12)

        canvas.create_text(left, height - 22, text="first", anchor="w")
        canvas.create_text(left + plot_width, height - 22, text="last", anchor="e")
        canvas.create_text(left - 6, top, text=f"{high:.4g}", anchor="e")
        canvas.create_text(left - 6, top + plot_height, text=f"{low:.4g}", anchor="e")
        note = f"plotted={len(sampled)}/{len(indexed_values)}, outside={outside_count}, {range_note}"
        if step > 1:
            note += f", sampled every {step}"
        canvas.create_text(left, 14, text=note, anchor="w", fill="black")

    def _plot_canvas_size(self, canvas: Canvas) -> tuple[int, int]:
        width = canvas.winfo_width()
        height = canvas.winfo_height()
        if width <= 2:
            width = canvas.winfo_reqwidth() or 460
        if height <= 2:
            height = canvas.winfo_reqheight() or 240
        return max(width, 260), max(height, 200)

    def _populate_jobs(self, jobs: list[dict]):
        self.job_entries = list(jobs)
        self._refresh_job_status_filter_values()
        self._apply_job_filter()

    def _refresh_job_status_filter_values(self):
        statuses = sorted({str(job.get("stat") or "").strip() for job in self.job_entries if str(job.get("stat") or "").strip()})
        values = ["All"] + statuses
        if hasattr(self, "job_status_box"):
            self.job_status_box.configure(values=values)
        if hasattr(self, "job_status_filter_var") and self.job_status_filter_var.get() not in values:
            self.job_status_filter_var.set("All")

    def _apply_job_filter(self):
        if not hasattr(self, "jobs_tree"):
            return
        for item in self.jobs_tree.get_children():
            self.jobs_tree.delete(item)
        self.job_entries_by_iid = {}
        query = self.job_filter_var.get().strip().lower() if hasattr(self, "job_filter_var") else ""
        status_filter = self.job_status_filter_var.get() if hasattr(self, "job_status_filter_var") else "All"
        visible_jobs = []
        for index, job in enumerate(self.job_entries):
            if status_filter != "All" and str(job.get("stat") or "") != status_filter:
                continue
            haystack = " ".join(str(job.get(key, "")) for key in ("jobid", "stat", "queue", "job_name", "exec_host", "submit_time", "user", "raw")).lower()
            if query and query not in haystack:
                continue
            visible_jobs.append(job)
            iid = self._job_tree_iid(job, index)
            self.job_entries_by_iid[iid] = job
            status_tag = self._job_status_tag(str(job.get("stat") or ""))
            self.jobs_tree.insert(
                "",
                END,
                iid=iid,
                values=(
                    job.get("jobid", ""),
                    job.get("stat", ""),
                    job.get("queue", ""),
                    job.get("job_name", ""),
                    job.get("exec_host", ""),
                    job.get("submit_time", ""),
                    job.get("user", ""),
                ),
                tags=(status_tag,) if status_tag else (),
            )
        self._update_job_summary(visible_jobs)
        self._update_job_selection_detail()

    def _job_tree_iid(self, job: dict, index: int) -> str:
        base = str(job.get("jobid") or f"job_{index}")
        if base not in self.job_entries_by_iid:
            return base
        return f"{base}#{index}"

    def _job_status_tag(self, status: str) -> str:
        normalized = status.upper()
        if normalized == "RUN":
            return "job_run"
        if normalized == "PEND":
            return "job_pend"
        if normalized in {"EXIT", "ZOMBI"}:
            return "job_exit"
        if normalized in {"DONE", "PSUSP", "USUSP", "SSUSP"}:
            return "job_done"
        return ""

    def _sort_jobs(self, column: str):
        if self.job_sort_column == column:
            self.job_sort_reverse = not self.job_sort_reverse
        else:
            self.job_sort_column = column
            self.job_sort_reverse = False
        self.job_entries.sort(key=lambda job: self._job_sort_value(job, column), reverse=self.job_sort_reverse)
        self._apply_job_filter()

    def _job_sort_value(self, job: dict, column: str):
        value = str(job.get(column) or "")
        if column == "jobid":
            number = value.split("[", 1)[0].split("#", 1)[0]
            if number.isdigit():
                return int(number)
        return value.lower()

    def _update_job_summary(self, visible_jobs: list[dict]):
        if not hasattr(self, "job_summary_var"):
            return
        counts = {}
        for job in self.job_entries:
            status = str(job.get("stat") or "UNKNOWN")
            counts[status] = counts.get(status, 0) + 1
        count_text = ", ".join(f"{key}={counts[key]}" for key in sorted(counts)) if counts else "no jobs"
        timestamp = datetime.now().strftime("%H:%M:%S")
        self.job_summary_var.set(f"Jobs loaded: total={len(self.job_entries)}, visible={len(visible_jobs)} | {count_text} | refreshed {timestamp}")

    def _update_job_selection_detail(self, _event=None):
        if not hasattr(self, "jobs_tree"):
            return
        selected_iids = list(self.jobs_tree.selection())
        if hasattr(self, "job_selection_var"):
            self.job_selection_var.set(f"Selected: {len(selected_iids)}")
        if not hasattr(self, "job_detail_text"):
            return
        self.job_detail_text.delete("1.0", END)
        if not selected_iids:
            self.job_detail_text.insert(END, "Select a job to inspect its parsed bjobs row. Use View bjobs -l for full pending reason and LSF details.")
            return
        if len(selected_iids) > 1:
            job_ids = [str(self.job_entries_by_iid.get(iid, {}).get("jobid", iid)) for iid in selected_iids]
            self.job_detail_text.insert(END, f"{len(selected_iids)} jobs selected:\n" + ", ".join(job_ids[:40]))
            if len(job_ids) > 40:
                self.job_detail_text.insert(END, f"\n... {len(job_ids) - 40} more")
            return
        job = self.job_entries_by_iid.get(selected_iids[0], {})
        lines = [
            f"Job ID: {job.get('jobid', '')}    Status: {job.get('stat', '')}    Queue: {job.get('queue', '')}",
            f"Name: {job.get('job_name', '')}",
            f"User: {job.get('user', '')}    From: {job.get('from_host', '')}    Exec: {job.get('exec_host', '')}    Submit: {job.get('submit_time', '')}",
            f"Raw: {job.get('raw', '')}",
        ]
        self.job_detail_text.insert(END, "\n".join(lines))

    def _selected_job_ids(self) -> list[str]:
        selected = []
        for iid in self.jobs_tree.selection():
            job = self.job_entries_by_iid.get(iid)
            selected.append(str(job.get("jobid") if job else iid).split("#", 1)[0])
        return selected

    def kill_selected_jobs(self):
        job_ids = self._selected_job_ids()
        if not job_ids:
            messagebox.showinfo("No jobs selected", "Select one or more jobs first.")
            return
        if not messagebox.askyesno("Confirm bkill", f"Kill selected jobs?\n\n{', '.join(job_ids)}"):
            return
        self._run_background("Kill selected jobs", lambda: cluster_ops.kill_jobs(self._collect_cluster_profile(), job_ids), lambda _: self.refresh_jobs())

    def kill_all_visible_jobs(self):
        job_ids = []
        for iid in self.jobs_tree.get_children():
            job = self.job_entries_by_iid.get(iid)
            job_ids.append(str(job.get("jobid") if job else iid).split("#", 1)[0])
        if not job_ids:
            messagebox.showinfo("No jobs visible", "Refresh jobs first.")
            return
        if not messagebox.askyesno("Confirm bkill", f"Kill all visible jobs?\n\n{', '.join(job_ids)}"):
            return
        self._run_background("Kill all visible jobs", lambda: cluster_ops.kill_jobs(self._collect_cluster_profile(), job_ids), lambda _: self.refresh_jobs())

    def download_selected_results(self):
        selected = list(self.results_tree.selection())
        if not selected:
            messagebox.showinfo("No result files selected", "Select one or more result files first.")
            return
        entries = [self.result_entries[int(item)] for item in selected]
        download_root = self._result_download_summary(entries)

        def on_success(results):
            self._log_download_results(results, download_root)
            failures = sum(1 for result in results if result.returncode != 0)
            if failures:
                self.result_scan_status_var.set(f"Download finished with {failures} failed command(s); local: {download_root}")
            else:
                self.result_scan_status_var.set(f"Download finished: {len(entries)} selected file(s) -> {download_root}")

        self.result_scan_status_var.set(f"Downloading {len(entries)} selected file(s) to {download_root} ...")
        self._run_background(
            "Download selected result files",
            lambda: self._download_result_entries(self._collect_cluster_profile(), entries),
            on_success,
            lambda exc: self.result_scan_status_var.set(f"Download failed: {exc}"),
        )

    def download_filtered_results(self):
        entries = [self.result_entries[int(item)] for item in self.results_tree.get_children()]
        if not entries:
            messagebox.showinfo("No result files visible", "Scan or adjust the search filter first.")
            return
        if not messagebox.askyesno("Confirm download", f"Download all {len(entries)} filtered files?"):
            return
        download_root = self._result_download_summary(entries)

        def on_success(results):
            self._log_download_results(results, download_root)
            failures = sum(1 for result in results if result.returncode != 0)
            if failures:
                self.result_scan_status_var.set(f"Download finished with {failures} failed command(s); local: {download_root}")
            else:
                self.result_scan_status_var.set(f"Download finished: {len(entries)} filtered file(s) -> {download_root}")

        self.result_scan_status_var.set(f"Downloading {len(entries)} filtered file(s) to {download_root} ...")
        self._run_background(
            "Download filtered result files",
            lambda: self._download_result_entries(self._collect_cluster_profile(), entries),
            on_success,
            lambda exc: self.result_scan_status_var.set(f"Download failed: {exc}"),
        )

    def preview_selected_csv(self):
        entry = self._single_selected_result("csv")
        if not entry:
            return
        profile = self._collect_cluster_profile()
        download_root = self._result_download_root([entry])
        local_path = cluster_ops.local_download_path_for_entry(profile, entry, download_root, layout="stage")
        if local_path.exists():
            self.result_scan_status_var.set(f"Using local CSV: {local_path}")
            self._show_csv_file(local_path)
            return

        def on_success(results):
            failures = sum(1 for result in results if result.returncode != 0)
            if local_path.exists():
                if failures:
                    self.result_scan_status_var.set(f"Download finished with {failures} failed command(s); opened local CSV: {local_path}")
                else:
                    self.result_scan_status_var.set(f"Download finished: CSV preview ready -> {local_path}")
                self._show_csv_file(local_path)
            else:
                self.result_scan_status_var.set(f"Download finished but CSV was not found: {local_path}")
                messagebox.showerror("CSV not downloaded", f"Expected file was not found:\n{local_path}")

        self.result_scan_status_var.set(f"Downloading CSV for preview to {local_path.parent} ...")
        self._run_background(
            "Download CSV for preview",
            lambda: cluster_ops.download_files(profile, [entry], local_dir=download_root, layout="stage"),
            on_success,
            lambda exc: self.result_scan_status_var.set(f"Download failed: {exc}"),
        )

    def open_selected_pdb_in_pymol(self):
        entry = self._single_selected_result("pdb")
        if not entry:
            return
        profile = self._collect_cluster_profile()
        download_root = self._result_download_root([entry])
        local_path = cluster_ops.local_download_path_for_entry(profile, entry, download_root, layout="stage")
        if local_path.exists():
            self.result_scan_status_var.set(f"Using local PDB: {local_path}")
            self._open_pdb_file_in_pymol(local_path)
            return

        def on_success(results):
            failures = sum(1 for result in results if result.returncode != 0)
            if local_path.exists():
                if failures:
                    self.result_scan_status_var.set(f"Download finished with {failures} failed command(s); opened local PDB: {local_path}")
                else:
                    self.result_scan_status_var.set(f"Download finished: PDB ready for PyMOL -> {local_path}")
                self._open_pdb_file_in_pymol(local_path)
            else:
                self.result_scan_status_var.set(f"Download finished but PDB was not found: {local_path}")
                messagebox.showerror("PDB not downloaded", f"Expected file was not found:\n{local_path}")

        self.result_scan_status_var.set(f"Downloading PDB for PyMOL to {local_path.parent} ...")
        self._run_background(
            "Download PDB for PyMOL",
            lambda: cluster_ops.download_files(profile, [entry], local_dir=download_root, layout="stage"),
            on_success,
            lambda exc: self.result_scan_status_var.set(f"Download failed: {exc}"),
        )

    def preview_selected_log_text(self):
        entry = self._single_selected_result(("out", "err", "txt", "log"))
        if not entry:
            return
        profile = self._collect_cluster_profile()
        download_root = self._result_download_root([entry])
        local_path = cluster_ops.local_download_path_for_entry(profile, entry, download_root, layout="stage")
        if local_path.exists():
            self.result_scan_status_var.set(f"Using local log/text: {local_path}")
            self._show_text_file(local_path)
            return

        def on_success(results):
            failures = sum(1 for result in results if result.returncode != 0)
            if local_path.exists():
                if failures:
                    self.result_scan_status_var.set(f"Download finished with {failures} failed command(s); opened log/text: {local_path}")
                else:
                    self.result_scan_status_var.set(f"Download finished: log/text preview ready -> {local_path}")
                self._show_text_file(local_path)
            else:
                self.result_scan_status_var.set(f"Download finished but log/text file was not found: {local_path}")
                messagebox.showerror("Log not downloaded", f"Expected file was not found:\n{local_path}")

        self.result_scan_status_var.set(f"Downloading log/text for preview to {local_path.parent} ...")
        self._run_background(
            "Download log for preview",
            lambda: cluster_ops.download_files(profile, [entry], local_dir=download_root, layout="stage"),
            on_success,
            lambda exc: self.result_scan_status_var.set(f"Download failed: {exc}"),
        )

    def preview_scratch_move(self):
        self.scratch_status_var.set("Previewing scratch move on cluster ...")

        def on_success(result):
            output = (result.stdout or "") + "\n" + (result.stderr or "")
            if result.returncode != 0:
                self.scratch_status_var.set(f"Preview failed with exit {result.returncode}; check Raw SSH/LSF output.")
            elif "SOURCE_EXISTS=no" in output:
                self.scratch_status_var.set("Preview complete: source directory does not exist; move would fail.")
            elif "TARGET_EXISTS=yes" in output:
                self.scratch_status_var.set("Preview complete: target already exists; move will refuse to overwrite.")
            else:
                self.scratch_status_var.set("Preview complete: source exists and target is available.")

        self._run_background(
            "Preview scratch directory move",
            lambda: cluster_ops.preview_scratch_migration(self._collect_cluster_profile()),
            on_success,
            lambda exc: self.scratch_status_var.set(f"Preview failed: {exc}"),
        )

    def use_current_target_scratch_settings(self):
        try:
            data = self._collect_config()
            self.scratch_source_date_var.set(str(data.get("scratch_date") or ""))
            self.scratch_user_dir_var.set(str(data.get("cluster_user") or self.cluster_profile.get("user") or ""))
            if not self.scratch_root_var.get().strip():
                self.scratch_root_var.set("/scratch")
            if not self.scratch_target_date_var.get().strip():
                self.scratch_target_date_var.set("auto")
            self._log("Scratch migration fields updated from current target settings.")
            self.scratch_status_var.set("Scratch fields updated from current target settings.")
        except Exception as exc:
            self._show_error("Failed to use current target settings", exc)

    def _save_remembered_password_if_requested(self):
        try:
            profile = self._collect_credential_profile()
            if profile.get("auth_method") != "password" or not profile.get("remember_password"):
                return
            password = self.password_var.get()
            if not password:
                return
            key = credential_store.save_password(profile, password)
            self._log(f"Saved encrypted password for {key}.")
        except Exception as exc:
            self._show_error("Failed to remember password", exc)

    def load_saved_password(self, silent: bool = False):
        try:
            profile = self._collect_credential_profile()
            password = credential_store.load_password(profile)
            if password:
                self.password_var.set(password)
                if not silent:
                    self._log("Loaded saved encrypted password.")
            elif not silent:
                messagebox.showinfo("No saved password", "No saved password found for this host/user/port.")
        except Exception as exc:
            if not silent:
                self._show_error("Failed to load saved password", exc)

    def forget_saved_password(self):
        try:
            profile = self._collect_credential_profile()
            existed = credential_store.delete_password(profile)
            self.password_var.set("")
            if existed:
                self._log("Deleted saved encrypted password.")
            else:
                messagebox.showinfo("No saved password", "No saved password existed for this host/user/port.")
        except Exception as exc:
            self._show_error("Failed to forget saved password", exc)

    def execute_scratch_move(self):
        try:
            profile = self._collect_cluster_profile()
            paths = cluster_ops.scratch_migration_paths(profile)
        except Exception as exc:
            self._show_error("Invalid scratch migration settings", exc)
            return
        message = (
            "Move scratch user directory now?\n\n"
            f"From: {paths['source']}\n"
            f"To:   {paths['target']}\n\n"
            "This refuses to overwrite an existing target directory, but it still moves data on the cluster."
        )
        if not messagebox.askyesno("Confirm scratch move", message):
            return
        self.scratch_status_var.set(f"Moving scratch directory: {paths['source']} -> {paths['target']} ...")

        def on_success(result):
            output = (result.stdout or "") + "\n" + (result.stderr or "")
            if result.returncode == 0 and "MOVE_DONE" in output:
                self.scratch_status_var.set(f"Move finished: {paths['target']}")
            elif result.returncode == 0:
                self.scratch_status_var.set("Move command finished; check Raw SSH/LSF output for details.")
            else:
                self.scratch_status_var.set(f"Move failed with exit {result.returncode}; check Raw SSH/LSF output.")

        self._run_background(
            "Move scratch user directory",
            lambda: cluster_ops.move_scratch_user_dir(profile),
            on_success,
            lambda exc: self.scratch_status_var.set(f"Move failed: {exc}"),
        )

    def _log_download_results(self, results, download_root: Path | str | None = None):
        for result in results:
            self._log(f"Download command exit {result.returncode}")
            if result.stdout.strip():
                self._log(result.stdout.strip())
            if result.stderr.strip():
                self._log(result.stderr.strip())
        self._log(f"Downloaded result files are under: {download_root or self._result_download_root()}")

    def _update_result_status_display(self):
        if not hasattr(self, "result_status_display_var"):
            return
        base_status = self.result_scan_status_var.get() if hasattr(self, "result_scan_status_var") else ""
        selected_count = len(self.results_tree.selection()) if hasattr(self, "results_tree") else 0
        visible_count = len(self.results_tree.get_children()) if hasattr(self, "results_tree") else 0
        scanned_count = len(self.result_entries) if hasattr(self, "result_entries") else 0
        count_text = f"Selected: {selected_count} file(s); Visible: {visible_count}; Scanned: {scanned_count}"
        self.result_status_display_var.set(f"{base_status} | {count_text}" if base_status else count_text)

    def _download_result_entries(self, profile: dict, entries: list[dict]):
        results = []
        for target, target_entries in self._download_entry_groups(entries).items():
            results.extend(
                cluster_ops.download_files(
                    profile,
                    target_entries,
                    local_dir=self._result_download_root_for_target(target),
                    layout="stage",
                )
            )
        return results

    def _result_download_root(self, entries: list[dict] | None = None) -> Path:
        return self._result_download_root_for_target(self._download_target_for_entries(entries))

    def _result_download_root_for_target(self, target: str) -> Path:
        output_dir = Path(self.output_dir_var.get()).expanduser()
        return output_dir / f"{cluster_ops.safe_local_name(target)}_workflow" / "retrieved_results"

    def _result_download_summary(self, entries: list[dict]) -> str:
        roots = [self._result_download_root_for_target(target) for target in self._download_entry_groups(entries)]
        if len(roots) == 1:
            return str(roots[0])
        return f"{len(roots)} target workflow folders: " + "; ".join(str(root) for root in roots)

    def _download_entry_groups(self, entries: list[dict]) -> dict[str, list[dict]]:
        groups: dict[str, list[dict]] = {}
        for entry in entries:
            target = self._download_target_for_entry(entry)
            groups.setdefault(target, []).append(entry)
        return groups

    def _download_target_for_entry(self, entry: dict) -> str:
        target = str(entry.get("target") or "").strip()
        if target and target.lower() != "unknown":
            return target
        return self._download_target_for_entries(None)

    def _download_target_for_entries(self, entries: list[dict] | None = None) -> str:
        targets = sorted(
            {
                str(entry.get("target") or "").strip()
                for entry in (entries or [])
                if str(entry.get("target") or "").strip().lower() != "unknown"
            }
        )
        if len(targets) == 1:
            return targets[0]
        try:
            current_target = str(self._collect_config().get("target") or "").strip()
        except Exception:
            current_target = ""
        if current_target:
            return current_target
        return targets[0] if targets else "target"

    def _current_pilot(self) -> str:
        try:
            pilot = str(self._collect_config().get("pilot") or "").strip()
        except Exception:
            pilot = ""
        return pilot or "pilot"

    def _refresh_result_filter_values(self):
        if not hasattr(self, "result_filter_boxes"):
            return
        for key, box in self.result_filter_boxes.items():
            values = sorted({str(entry.get(key, "") or "unknown") for entry in self.result_entries})
            current = box.get() or "All"
            options = ["All"] + values
            box.configure(values=options)
            if current not in options:
                box.set("All")

    def _entry_matches_result_filters(self, entry: dict) -> bool:
        filters = {
            "target": self.result_target_filter_var.get(),
            "pilot": self.result_pilot_filter_var.get(),
            "stage": self.result_stage_filter_var.get(),
            "shard": self.result_shard_filter_var.get(),
        }
        for key, selected in filters.items():
            if selected and selected != "All" and str(entry.get(key, "")) != selected:
                return False
        return True

    def _single_selected_result(self, file_type: str | tuple[str, ...]) -> dict | None:
        allowed_types = (file_type,) if isinstance(file_type, str) else file_type
        selected = list(self.results_tree.selection())
        if len(selected) != 1:
            messagebox.showinfo("Select one file", f"Select exactly one {', '.join('.' + item for item in allowed_types)} file first.")
            return None
        entry = self.result_entries[int(selected[0])]
        if str(entry.get("type", "")).lower() not in {item.lower() for item in allowed_types}:
            messagebox.showinfo("Wrong file type", f"Selected file is not one of: {', '.join('.' + item for item in allowed_types)}.")
            return None
        return entry

    def _show_text_file(self, path: Path):
        try:
            content = path.read_text(encoding="utf-8", errors="replace")
        except Exception as exc:
            self._show_error("Failed to open text file", exc)
            return
        if len(content) > 2_000_000:
            content = content[:2_000_000] + "\n\n[Preview truncated at 2 MB.]"
        self._show_text_window(path.name, content, subtitle=str(path))

    def _show_text_window(self, title: str, content: str, subtitle: str = ""):
        window = Toplevel(self.root)
        window.title(title)
        window.geometry("1000x650")
        if subtitle:
            Label(window, text=subtitle).pack(fill="x", padx=8, pady=6)
        frame = Frame(window)
        frame.pack(fill="both", expand=True, padx=8, pady=8)
        text = Text(frame, wrap="none")
        y_scroll = ttk.Scrollbar(frame, orient="vertical", command=text.yview)
        x_scroll = ttk.Scrollbar(frame, orient="horizontal", command=text.xview)
        text.configure(yscrollcommand=y_scroll.set, xscrollcommand=x_scroll.set)
        text.grid(row=0, column=0, sticky="nsew")
        y_scroll.grid(row=0, column=1, sticky="ns")
        x_scroll.grid(row=1, column=0, sticky="ew")
        frame.rowconfigure(0, weight=1)
        frame.columnconfigure(0, weight=1)
        text.insert(END, content)
        text.configure(state="disabled")

    def _read_csv_page(self, path: Path, start_index: int, page_size: int):
        with path.open("r", encoding="utf-8-sig", errors="replace", newline="") as handle:
            reader = csv.reader(handle)
            try:
                headers = list(next(reader))
            except StopIteration:
                return [], [], False

            skipped = 0
            while skipped < start_index:
                try:
                    next(reader)
                except StopIteration:
                    break
                skipped += 1

            rows = []
            for _ in range(page_size + 1):
                try:
                    rows.append(next(reader))
                except StopIteration:
                    break
            has_more = len(rows) > page_size
            if has_more:
                rows = rows[:page_size]

        max_columns = max([len(headers)] + [len(row) for row in rows], default=0)
        if max_columns and len(headers) < max_columns:
            headers.extend(f"col_{index + 1}" for index in range(len(headers), max_columns))
        return headers, rows, has_more

    def _show_csv_file(self, path: Path):
        window = Toplevel(self.root)
        window.title(f"CSV Viewer - {path.name}")
        window.geometry("1280x760")

        page_size_var = StringVar(value="200")
        column_search_var = StringVar()
        status_var = StringVar(value=f"Opening {path.name} ...")
        column_matches = []
        state = {"start": 0, "has_more": False, "headers": []}

        header_frame = Frame(window)
        header_frame.pack(fill="x", padx=8, pady=6)
        Label(header_frame, text=str(path)).grid(row=0, column=0, columnspan=8, sticky="w")
        Label(header_frame, text="Rows/page").grid(row=1, column=0, sticky="w", pady=(6, 0))
        Entry(header_frame, textvariable=page_size_var, width=8).grid(row=1, column=1, sticky="w", padx=4, pady=(6, 0))
        prev_button = Button(header_frame, text="Prev page")
        next_button = Button(header_frame, text="Next page")
        reload_button = Button(header_frame, text="Reload page")
        copy_columns_button = Button(header_frame, text="Copy column names")
        prev_button.grid(row=1, column=2, padx=4, pady=(6, 0))
        next_button.grid(row=1, column=3, padx=4, pady=(6, 0))
        reload_button.grid(row=1, column=4, padx=4, pady=(6, 0))
        copy_columns_button.grid(row=1, column=5, padx=4, pady=(6, 0))
        Label(header_frame, textvariable=status_var).grid(row=1, column=6, sticky="w", padx=10, pady=(6, 0))
        header_frame.columnconfigure(6, weight=1)

        body = Frame(window)
        body.pack(fill="both", expand=True, padx=8, pady=(0, 8))

        column_frame = LabelFrame(body, text="Columns")
        column_frame.grid(row=0, column=0, sticky="ns", padx=(0, 8))
        Label(column_frame, text="Search").pack(fill="x", padx=6, pady=(6, 0))
        Entry(column_frame, textvariable=column_search_var).pack(fill="x", padx=6, pady=4)
        column_list_frame = Frame(column_frame)
        column_list_frame.pack(fill="both", expand=True, padx=6, pady=(0, 6))
        column_list = Listbox(column_list_frame, width=34, exportselection=False)
        column_list_scroll = ttk.Scrollbar(column_list_frame, orient="vertical", command=column_list.yview)
        column_list.configure(yscrollcommand=column_list_scroll.set)
        column_list.grid(row=0, column=0, sticky="nsew")
        column_list_scroll.grid(row=0, column=1, sticky="ns")
        column_list_frame.rowconfigure(0, weight=1)
        column_list_frame.columnconfigure(0, weight=1)

        table_frame = Frame(body)
        table_frame.grid(row=0, column=1, sticky="nsew")
        body.rowconfigure(0, weight=1)
        body.columnconfigure(1, weight=1)

        tree = ttk.Treeview(table_frame, show="headings")
        y_scroll = ttk.Scrollbar(table_frame, orient="vertical", command=tree.yview)
        x_scroll = ttk.Scrollbar(table_frame, orient="horizontal", command=tree.xview)
        tree.configure(yscrollcommand=y_scroll.set, xscrollcommand=x_scroll.set)
        tree.grid(row=0, column=0, sticky="nsew")
        y_scroll.grid(row=0, column=1, sticky="ns")
        x_scroll.grid(row=1, column=0, sticky="ew")
        table_frame.rowconfigure(0, weight=1)
        table_frame.columnconfigure(0, weight=1)

        def page_size() -> int:
            try:
                return max(1, min(5000, int(page_size_var.get().strip())))
            except ValueError:
                page_size_var.set("200")
                return 200

        def refresh_column_list():
            nonlocal column_matches
            query = column_search_var.get().strip().lower()
            column_list.delete(0, END)
            column_matches = []
            for index, header in enumerate(state["headers"]):
                label = f"{index + 1}: {header or f'col_{index + 1}'}"
                if query and query not in label.lower():
                    continue
                column_matches.append(index)
                column_list.insert(END, label)

        def focus_selected_column(_event=None):
            selection = column_list.curselection()
            if not selection or not state["headers"]:
                return
            column_index = column_matches[selection[0]]
            fraction = column_index / max(1, len(state["headers"]) - 1)
            tree.xview_moveto(fraction)
            status_var.set(f"Focused column {column_index + 1}: {state['headers'][column_index]}")

        def copy_column_names():
            names = [f"{index + 1}\t{name}" for index, name in enumerate(state["headers"])]
            window.clipboard_clear()
            window.clipboard_append("\n".join(names))
            status_var.set(f"Copied {len(names)} column names")

        def load_page(start_index: int):
            start_index = max(0, start_index)
            size = page_size()
            status_var.set(f"Loading rows {start_index + 1}-{start_index + size} ...")
            window.update_idletasks()
            try:
                headers, rows, has_more = self._read_csv_page(path, start_index, size)
            except Exception as exc:
                self._show_error("Failed to open CSV", exc)
                return
            if not headers and not rows:
                messagebox.showinfo("Empty CSV", f"No rows found in:\n{path}")
                window.destroy()
                return

            state["start"] = start_index
            state["has_more"] = has_more
            state["headers"] = headers
            column_ids = [f"col_{index}" for index, _header in enumerate(headers)]

            for item in tree.get_children():
                tree.delete(item)
            tree.configure(columns=column_ids)
            for column_id, header in zip(column_ids, headers):
                tree.heading(column_id, text=header or column_id)
                tree.column(column_id, width=max(120, min(260, len(str(header)) * 9 + 40)), minwidth=80, anchor="w", stretch=False)
            for row in rows:
                padded = row + [""] * max(0, len(headers) - len(row))
                tree.insert("", END, values=padded[:len(headers)])

            refresh_column_list()
            end_row = start_index + len(rows)
            more_text = "more rows available" if has_more else "end of file"
            status_var.set(f"Showing rows {start_index + 1}-{end_row}; {len(headers)} columns; {more_text}")
            prev_button.configure(state="normal" if start_index > 0 else "disabled")
            next_button.configure(state="normal" if has_more else "disabled")

        prev_button.configure(command=lambda: load_page(max(0, state["start"] - page_size())))
        next_button.configure(command=lambda: load_page(state["start"] + page_size()))
        reload_button.configure(command=lambda: load_page(state["start"]))
        copy_columns_button.configure(command=copy_column_names)
        column_search_var.trace_add("write", lambda *_: refresh_column_list())
        column_list.bind("<<ListboxSelect>>", focus_selected_column)

        load_page(0)

    def _open_pdb_file_in_pymol(self, path: Path):
        profile = self._collect_cluster_profile()
        executable = str(profile.get("pymol_executable") or "pymol").strip()
        try:
            subprocess.Popen([executable, str(path)])
            self._log(f"Opened PDB in PyMOL: {path}")
        except FileNotFoundError:
            messagebox.showerror(
                "PyMOL not found",
                "Could not find PyMOL executable.\n\n"
                "Add `pymol_executable` to cluster_profile.local.json, for example:\n"
                "\"pymol_executable\": \"C:/Program Files/PyMOL/PyMOLWin.exe\"",
            )
        except Exception as exc:
            self._show_error("Failed to open PyMOL", exc)

    def _apply_result_filter(self):
        if not hasattr(self, "results_tree"):
            return
        query = self.result_search_var.get().strip().lower()
        for item in self.results_tree.get_children():
            self.results_tree.delete(item)

        indexed_entries = list(enumerate(self.result_entries))
        indexed_entries.sort(
            key=lambda pair: self._result_sort_key(pair[1], self.result_sort_column),
            reverse=self.result_sort_reverse,
        )
        for original_index, entry in indexed_entries:
            if not self._entry_matches_result_filters(entry):
                continue
            haystack = " ".join(str(entry.get(key, "")) for key in ("target", "pilot", "stage", "shard", "type", "name", "path")).lower()
            if query and query not in haystack:
                continue
            self.results_tree.insert(
                "",
                END,
                iid=str(original_index),
                values=(
                    entry.get("target", ""),
                    entry.get("pilot", ""),
                    entry.get("stage", ""),
                    entry.get("shard", ""),
                    entry.get("type", ""),
                    entry.get("size", 0),
                    entry.get("mtime", ""),
                    entry.get("name", ""),
                    entry.get("path", ""),
                ),
            )
        self._update_result_status_display()

    def _sort_results(self, column: str):
        if self.result_sort_column == column:
            self.result_sort_reverse = not self.result_sort_reverse
        else:
            self.result_sort_column = column
            self.result_sort_reverse = column in {"mtime", "size"}
        self._apply_result_filter()

    def _result_sort_key(self, entry: dict, column: str):
        value = entry.get(column, "")
        if column == "size":
            return int(value or 0)
        return str(value).lower()

    def _schedule_auto_refresh(self):
        if not self.auto_refresh_var.get():
            return

        def tick():
            if not self.auto_refresh_var.get():
                return
            self.refresh_jobs()
            self.root.after(10000, tick)

        self.root.after(100, tick)

    def _log(self, message: str):
        self.log_lines.append(message)
        if len(self.log_lines) > 5000:
            self.log_lines = self.log_lines[-5000:]
        if self.log_text is not None and self.log_text.winfo_exists():
            self.log_text.insert(END, message + "\n")
            self.log_text.see(END)
        self._append_debug_log(f"[LOG] {message}\n")

    def _record_command_result(self, description: str, result: cluster_ops.CommandResult):
        timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        command = " ".join(str(part) for part in result.command)
        block = (
            f"\n===== {timestamp} | {description} | exit {result.returncode} =====\n"
            f"$ {command}\n"
            f"--- stdout ---\n{result.stdout.rstrip()}\n"
            f"--- stderr ---\n{result.stderr.rstrip()}\n"
        )
        self.raw_output_blocks.append(block)
        if len(self.raw_output_blocks) > 1000:
            self.raw_output_blocks = self.raw_output_blocks[-1000:]
        if self.raw_output_text is not None and self.raw_output_text.winfo_exists():
            self.raw_output_text.insert(END, block)
            self.raw_output_text.see(END)
        self._append_debug_log(block)

    def _append_debug_log(self, text: str):
        try:
            self.debug_log_path.parent.mkdir(parents=True, exist_ok=True)
            with self.debug_log_path.open("a", encoding="utf-8") as handle:
                handle.write(text)
        except Exception:
            pass

    def clear_raw_output(self):
        self.raw_output_blocks.clear()
        if self.raw_output_text is not None and self.raw_output_text.winfo_exists():
            self.raw_output_text.delete("1.0", END)

    def open_debug_log_folder(self):
        try:
            import os
            os.startfile(self.debug_log_path.parent)
        except Exception as exc:
            self._show_error("Failed to open debug log folder", exc)

    def _show_error(self, title: str, exc: Exception):
        detail = "".join(traceback.format_exception_only(type(exc), exc)).strip()
        self._set_status(f"Error: {title}")
        self._log(f"{title}: {detail}")
        messagebox.showerror(title, detail)

    def _set_status(self, message: str):
        if hasattr(self, "status_var"):
            self.status_var.set(message)
        self._append_debug_log(f"[STATUS] {message}\n")


def main():
    root = Tk()
    WorkflowGui(root)
    root.mainloop()


if __name__ == "__main__":
    main()
