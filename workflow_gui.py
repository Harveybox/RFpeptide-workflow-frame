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
            ("input_pdb", "Input PDB path", str),
            ("target_chain", "Target chain", str),
            ("binder_chain", "Binder chain", str),
            ("contigs", "RFdiffusion contigs", str),
            ("binder_length", "Binder length note", str),
        ],
    ),
    (
        "RFDiffusion",
        [
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
            ("proteinmpnn.relax_cycles", "Relax cycles", int),
            ("proteinmpnn.seqs_per_struct", "Seqs per structure", int),
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
            ("afcyc.rmsd_script", "RMSD script", str),
            ("afcyc.merge_script", "Merge script", str),
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
            ("pyrosetta.merge_script", "Merge script", str),
            ("pyrosetta.n_shards", "Shards", int),
            ("pyrosetta.ncpu", "CPU cores", int),
            ("pyrosetta.ptile", "PTILE", int),
            ("pyrosetta.pack_input", "Pack input", bool),
            ("pyrosetta.pack_separated", "Pack separated", bool),
            ("pyrosetta.packstat", "Packstat", bool),
        ],
    ),
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
        self.root.title("Binder workflow generator")
        self.root.geometry("980x760")
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

        self._build()
        self._load_values(self.config_data)
        self._log(f"Loaded config: {self.config_path}")

    def _build(self):
        top = Frame(self.root)
        top.pack(fill="x", padx=12, pady=10)

        self.config_path_var = StringVar(value=str(self.config_path))
        self.output_dir_var = StringVar(value=str(self.output_dir))
        self.top_local_pdb_var = StringVar()
        self.current_project_var = StringVar(value="Current project: not loaded")

        project_banner = Frame(top, bg="#17324d", bd=1, relief="solid")
        project_banner.grid(row=0, column=0, columnspan=5, sticky="ew", pady=(0, 8))
        Label(project_banner, text="CURRENT PROJECT CONFIG", bg="#17324d", fg="#ffffff").pack(side="left", padx=(10, 8), pady=7)
        Label(project_banner, textvariable=self.current_project_var, bg="#17324d", fg="#d7f2ff", anchor="w").pack(side="left", fill="x", expand=True, padx=(0, 10), pady=7)

        Label(top, text="Config").grid(row=1, column=0, sticky="w")
        Entry(top, textvariable=self.config_path_var).grid(row=1, column=1, sticky="ew", padx=6)
        Button(top, text="Open", command=self.open_config).grid(row=1, column=2, padx=3)
        Button(top, text="Save", command=self.save_config).grid(row=1, column=3, padx=3)
        Button(top, text="Save as", command=self.save_config_as).grid(row=1, column=4, padx=3)

        Label(top, text="Local PDB").grid(row=2, column=0, sticky="w", pady=(8, 0))
        Entry(top, textvariable=self.top_local_pdb_var).grid(row=2, column=1, sticky="ew", padx=6, pady=(8, 0))
        Button(top, text="Choose PDB", command=self.choose_local_input_pdb).grid(row=2, column=2, padx=3, pady=(8, 0))
        Button(top, text="Preview/Clean", command=self.preview_top_local_pdb).grid(row=2, column=3, padx=3, pady=(8, 0))

        Label(top, text="Output").grid(row=3, column=0, sticky="w", pady=(8, 0))
        Entry(top, textvariable=self.output_dir_var).grid(row=3, column=1, sticky="ew", padx=6, pady=(8, 0))
        Button(top, text="Choose", command=self.choose_output_dir).grid(row=3, column=2, padx=3, pady=(8, 0))
        Button(top, text="Generate workflow", command=self.generate).grid(row=3, column=3, columnspan=2, sticky="ew", padx=3, pady=(8, 0))
        top.columnconfigure(1, weight=1)

        notebook = ttk.Notebook(self.root)
        self.notebook = notebook
        notebook.pack(fill="both", expand=True, padx=12, pady=(0, 10))

        self.pdb_preprocess_frame = Frame(notebook)
        notebook.add(self.pdb_preprocess_frame, text="PDB Preprocess")
        self._build_pdb_preprocess(self.pdb_preprocess_frame)

        for group_name, fields in FIELD_GROUPS:
            frame = Frame(notebook)
            notebook.add(frame, text=group_name)
            self._build_group(frame, group_name, fields)
        for traced_key in ("target", "pilot", "scratch_date", "project_dir_name"):
            if traced_key in self.variables:
                self.variables[traced_key][0].trace_add("write", lambda *_: self._update_current_project_banner())
        cluster_settings_frame = Frame(notebook)
        notebook.add(cluster_settings_frame, text="Cluster Settings")
        self.cluster_settings_frame = cluster_settings_frame
        self._build_cluster_settings(cluster_settings_frame)
        cluster_frame = Frame(notebook)
        notebook.add(cluster_frame, text="Cluster Dashboard")
        self._build_cluster_dashboard(cluster_frame)
        results_frame = Frame(notebook)
        notebook.add(results_frame, text="Results Browser")
        self._build_results_browser(results_frame)
        data_frame = Frame(notebook)
        notebook.add(data_frame, text="Data Processing")
        self._build_data_processing(data_frame)
        scratch_frame = Frame(notebook)
        notebook.add(scratch_frame, text="Scratch Safety")
        self._build_scratch_safety(scratch_frame)

        log_frame = LabelFrame(self.root, text="Log")
        log_frame.pack(fill="both", expand=False, padx=12, pady=(0, 12))
        self.log_text = Text(log_frame, height=6, wrap="word")
        self.log_text.pack(fill="both", expand=True, padx=6, pady=6)

        raw_frame = LabelFrame(self.root, text="Raw SSH/LSF output")
        raw_frame.pack(fill="both", expand=False, padx=12, pady=(0, 12))
        raw_toolbar = Frame(raw_frame)
        raw_toolbar.pack(fill="x", padx=6, pady=(6, 0))
        Button(raw_toolbar, text="Clear raw output", command=self.clear_raw_output).pack(side="left", padx=(0, 4))
        Button(raw_toolbar, text="Open debug log folder", command=self.open_debug_log_folder).pack(side="left", padx=4)
        self.raw_output_text = Text(raw_frame, height=8, wrap="none")
        self.raw_output_text.pack(fill="both", expand=True, padx=6, pady=6)

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

        self.pdb_raw_path_var = StringVar()
        self.pdb_clean_output_var = StringVar()
        self.pdb_keep_protein_only_var = BooleanVar(value=True)
        self.pdb_renumber_atoms_var = BooleanVar(value=True)
        self.pdb_renumber_residues_var = BooleanVar(value=False)
        self.pdb_set_as_input_var = BooleanVar(value=True)
        self.pdb_preprocess_status_var = StringVar(value="Select a PDB and preview it.")

        input_frame = LabelFrame(body, text="Raw PDB")
        input_frame.pack(fill="x", pady=(0, 8))
        Label(input_frame, text="PDB file").grid(row=0, column=0, sticky="w", padx=6, pady=6)
        Entry(input_frame, textvariable=self.pdb_raw_path_var).grid(row=0, column=1, sticky="ew", padx=6, pady=6)
        Button(input_frame, text="Browse", command=self.choose_preprocess_pdb).grid(row=0, column=2, padx=3, pady=6)
        Button(input_frame, text="Preview sequence", command=self.preview_preprocess_pdb).grid(row=0, column=3, padx=3, pady=6)
        Button(input_frame, text="Use as workflow input", command=self.use_preprocess_pdb_as_input).grid(row=0, column=4, padx=3, pady=6)
        input_frame.columnconfigure(1, weight=1)

        options_frame = LabelFrame(body, text="Cleanup options")
        options_frame.pack(fill="x", pady=(0, 8))
        Label(options_frame, text="Cleaned output").grid(row=0, column=0, sticky="w", padx=6, pady=6)
        Entry(options_frame, textvariable=self.pdb_clean_output_var).grid(row=0, column=1, sticky="ew", padx=6, pady=6)
        Button(options_frame, text="Choose output", command=self.choose_clean_pdb_output).grid(row=0, column=2, padx=3, pady=6)
        Checkbutton(options_frame, text="Remove non-protein components", variable=self.pdb_keep_protein_only_var).grid(row=1, column=0, sticky="w", padx=6, pady=(0, 6))
        Checkbutton(options_frame, text="Renumber atom serials", variable=self.pdb_renumber_atoms_var).grid(row=1, column=1, sticky="w", padx=6, pady=(0, 6))
        Checkbutton(options_frame, text="Renumber residues continuously by chain", variable=self.pdb_renumber_residues_var).grid(row=1, column=2, sticky="w", padx=6, pady=(0, 6))
        Checkbutton(options_frame, text="Set cleaned PDB as workflow input", variable=self.pdb_set_as_input_var).grid(row=2, column=0, columnspan=2, sticky="w", padx=6, pady=(0, 6))
        Button(options_frame, text="Clean PDB", command=self.clean_preprocess_pdb).grid(row=2, column=2, sticky="ew", padx=6, pady=(0, 6))
        self._build_status_bar(options_frame, self.pdb_preprocess_status_var, "PDB").grid(row=2, column=3, sticky="ew", padx=6, pady=(0, 6))
        options_frame.columnconfigure(1, weight=1)
        options_frame.columnconfigure(3, weight=1)

        preview_frame = Frame(body)
        preview_frame.pack(fill="both", expand=True, pady=(0, 8))

        summary_frame = LabelFrame(preview_frame, text="Sequence summary")
        summary_frame.pack(side="left", fill="both", expand=True, padx=(0, 4))
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
        table_frame.pack(side="left", fill="both", expand=True, padx=(4, 0))
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
        outer = Frame(parent)
        outer.pack(fill="both", expand=True)
        scroll_canvas = Canvas(outer, highlightthickness=0)
        scroll_y = ttk.Scrollbar(outer, orient="vertical", command=scroll_canvas.yview)
        scroll_canvas.configure(yscrollcommand=scroll_y.set)
        scroll_y.pack(side="right", fill="y")
        scroll_canvas.pack(side="left", fill="both", expand=True)
        body = Frame(scroll_canvas)
        body_window = scroll_canvas.create_window((0, 0), window=body, anchor="nw")

        def update_scroll_region(_event=None):
            scroll_canvas.configure(scrollregion=scroll_canvas.bbox("all"))

        def fit_body_width(event):
            scroll_canvas.itemconfigure(body_window, width=event.width)

        def enable_mousewheel(_event=None):
            scroll_canvas.bind_all("<MouseWheel>", on_mousewheel)

        def disable_mousewheel(_event=None):
            scroll_canvas.unbind_all("<MouseWheel>")

        def on_mousewheel(event):
            scroll_canvas.yview_scroll(int(-1 * (event.delta / 120)), "units")

        body.bind("<Configure>", update_scroll_region)
        scroll_canvas.bind("<Configure>", fit_body_width)
        body.bind("<Enter>", enable_mousewheel)
        body.bind("<Leave>", disable_mousewheel)
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
        Button(actions, text="Submit stage", command=self.submit_selected_stage).grid(row=0, column=3, padx=4, pady=6, sticky="ew")
        self._build_status_bar(actions, self.status_var, "Cluster").grid(row=1, column=0, columnspan=4, sticky="ew", padx=6, pady=(0, 6))
        for column in range(4):
            actions.columnconfigure(column, weight=1)

        jobs_frame = LabelFrame(body, text="Jobs")
        jobs_frame.pack(fill="both", expand=True)
        jobs_toolbar = Frame(jobs_frame)
        jobs_toolbar.pack(fill="x", padx=6, pady=6)
        Button(jobs_toolbar, text="Refresh jobs", command=self.refresh_jobs).pack(side="left", padx=(0, 4))
        Button(jobs_toolbar, text="View bjobs -l", command=self.view_selected_job_details).pack(side="left", padx=4)
        Button(jobs_toolbar, text="Kill selected", command=self.kill_selected_jobs).pack(side="left", padx=4)
        Button(jobs_toolbar, text="Kill all visible", command=self.kill_all_visible_jobs).pack(side="left", padx=4)
        Checkbutton(jobs_toolbar, text="Auto refresh every 10s", variable=self.auto_refresh_var, command=self._schedule_auto_refresh).pack(side="left", padx=12)

        columns = ("jobid", "stat", "queue", "job_name", "exec_host", "user")
        self.jobs_tree = ttk.Treeview(jobs_frame, columns=columns, show="headings", selectmode="extended", height=10)
        for column in columns:
            self.jobs_tree.heading(column, text=column)
            width = 90 if column not in {"job_name", "exec_host"} else 240
            self.jobs_tree.column(column, width=width, anchor="w")
        self.jobs_tree.pack(fill="both", expand=True, padx=6, pady=(0, 6))

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
        outer = Frame(parent)
        outer.pack(fill="both", expand=True)
        scroll_canvas = Canvas(outer, highlightthickness=0)
        scroll_y = ttk.Scrollbar(outer, orient="vertical", command=scroll_canvas.yview)
        scroll_canvas.configure(yscrollcommand=scroll_y.set)
        scroll_y.pack(side="right", fill="y")
        scroll_canvas.pack(side="left", fill="both", expand=True)
        body = Frame(scroll_canvas)
        body_window = scroll_canvas.create_window((0, 0), window=body, anchor="nw")

        def update_scroll_region(_event=None):
            scroll_canvas.configure(scrollregion=scroll_canvas.bbox("all"))

        def fit_body_width(event):
            scroll_canvas.itemconfigure(body_window, width=event.width)

        def enable_mousewheel(_event=None):
            scroll_canvas.bind_all("<MouseWheel>", on_mousewheel)

        def disable_mousewheel(_event=None):
            scroll_canvas.unbind_all("<MouseWheel>")

        def on_mousewheel(event):
            scroll_canvas.yview_scroll(int(-1 * (event.delta / 120)), "units")

        body.bind("<Configure>", update_scroll_region)
        scroll_canvas.bind("<Configure>", fit_body_width)
        body.bind("<Enter>", enable_mousewheel)
        body.bind("<Leave>", disable_mousewheel)
        body.configure(padx=12, pady=12)

        self.data_csv_entries = []
        self.data_loaded_csv_path = None
        self.data_loaded_headers = []
        self.data_metric_var = StringVar()
        self.data_percentile_var = StringVar(value="90")
        self.data_threshold_var = StringVar()
        self.data_bins_var = StringVar(value="40")
        self.data_hist_min_var = StringVar()
        self.data_hist_max_var = StringVar()
        self.data_robust_hist_var = BooleanVar(value=True)
        self.data_status_var = StringVar(value="No merged CSV loaded")

        scan_frame = LabelFrame(body, text="Merged CSV scan")
        scan_frame.pack(fill="x", pady=(0, 8))
        Label(scan_frame, text="Scan root").grid(row=0, column=0, sticky="w", padx=6, pady=6)
        Label(scan_frame, textvariable=self.result_scan_status_var).grid(row=0, column=1, sticky="w", padx=6, pady=6)
        Button(scan_frame, text="Use current target", command=self.use_current_target_result_scan_settings).grid(row=0, column=2, padx=3, pady=6)
        Button(scan_frame, text="Scan merged CSV", command=self.scan_data_merged_csvs).grid(row=0, column=3, padx=3, pady=6)
        Button(scan_frame, text="Preview selected CSV", command=self.preview_data_csv).grid(row=0, column=4, padx=3, pady=6)
        Button(scan_frame, text="Load for analysis", command=self.load_data_csv_for_analysis).grid(row=0, column=5, padx=3, pady=6)
        self._build_status_bar(scan_frame, self.data_status_var, "Data").grid(row=1, column=0, columnspan=6, sticky="ew", padx=6, pady=(0, 6))
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
        Label(analysis_frame, text="Percentile").grid(row=0, column=2, sticky="w", padx=6, pady=6)
        Entry(analysis_frame, textvariable=self.data_percentile_var, width=8).grid(row=0, column=3, sticky="w", padx=6, pady=6)
        Label(analysis_frame, text="Custom threshold").grid(row=0, column=4, sticky="w", padx=6, pady=6)
        Entry(analysis_frame, textvariable=self.data_threshold_var, width=12).grid(row=0, column=5, sticky="w", padx=6, pady=6)
        Label(analysis_frame, text="Bins").grid(row=0, column=6, sticky="w", padx=6, pady=6)
        Entry(analysis_frame, textvariable=self.data_bins_var, width=8).grid(row=0, column=7, sticky="w", padx=6, pady=6)
        Label(analysis_frame, text="Hist min").grid(row=1, column=0, sticky="w", padx=6, pady=(0, 6))
        Entry(analysis_frame, textvariable=self.data_hist_min_var, width=12).grid(row=1, column=1, sticky="w", padx=6, pady=(0, 6))
        Label(analysis_frame, text="Hist max").grid(row=1, column=2, sticky="w", padx=6, pady=(0, 6))
        Entry(analysis_frame, textvariable=self.data_hist_max_var, width=12).grid(row=1, column=3, sticky="w", padx=6, pady=(0, 6))
        Checkbutton(analysis_frame, text="Auto robust range 1-99% when min/max blank", variable=self.data_robust_hist_var).grid(row=1, column=4, columnspan=4, sticky="w", padx=6, pady=(0, 6))
        Button(analysis_frame, text="Analyze metric", command=self.analyze_data_metric).grid(row=2, column=0, padx=6, pady=(0, 6))
        self._build_status_bar(analysis_frame, self.data_status_var, "Data").grid(row=2, column=1, columnspan=7, sticky="ew", padx=6, pady=(0, 6))
        analysis_frame.columnconfigure(1, weight=1)

        plot_frame = LabelFrame(body, text="Metric plots")
        plot_frame.pack(fill="both", expand=True, pady=(0, 8))
        hist_frame = LabelFrame(plot_frame, text="Histogram")
        hist_frame.pack(side="left", fill="both", expand=True, padx=(6, 3), pady=6)
        self.data_hist_canvas = Canvas(hist_frame, height=260, bg="white")
        self.data_hist_canvas.pack(fill="both", expand=True, padx=6, pady=6)
        scatter_frame = LabelFrame(plot_frame, text="Scatter by CSV order")
        scatter_frame.pack(side="left", fill="both", expand=True, padx=(3, 6), pady=6)
        self.data_scatter_canvas = Canvas(scatter_frame, height=260, bg="white")
        self.data_scatter_canvas.pack(fill="both", expand=True, padx=6, pady=6)

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
        if hasattr(self, "cluster_settings_frame"):
            self.notebook.select(self.cluster_settings_frame)

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

    def _load_values(self, config: dict):
        for key, (var, value_type) in self.variables.items():
            value = get_nested(config, key)
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
                try:
                    value = int(str(raw_value).strip())
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
        self._set_workflow_local_pdb(Path(path))
        self.pdb_raw_path_var.set(str(path))
        self.pdb_clean_output_var.set(str(self._default_clean_pdb_output_path(Path(path))))
        self._log(f"Selected local PDB: {path}")

    def preview_top_local_pdb(self):
        path_text = self.top_local_pdb_var.get().strip()
        if not path_text:
            self.choose_local_input_pdb()
            path_text = self.top_local_pdb_var.get().strip()
            if not path_text:
                return
        self.pdb_raw_path_var.set(path_text)
        self.notebook.select(self.pdb_preprocess_frame)
        self.preview_preprocess_pdb()

    def choose_preprocess_pdb(self):
        path = filedialog.askopenfilename(
            title="Choose raw PDB for preprocessing",
            initialdir=str(SCRIPT_DIR),
            filetypes=[("PDB files", "*.pdb *.ent"), ("All files", "*.*")],
        )
        if not path:
            return
        self.pdb_raw_path_var.set(path)
        self.pdb_clean_output_var.set(str(self._default_clean_pdb_output_path(Path(path))))
        self.preview_preprocess_pdb()

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
            )
            self.pdb_clean_output_var.set(str(output_path))
            summary = pdb_preprocess.parse_pdb(output_path)
            self._populate_pdb_preview(summary)
            if self.pdb_set_as_input_var.get():
                self._set_workflow_local_pdb(output_path)
            self.pdb_preprocess_status_var.set(
                f"Cleaned PDB: kept {result['kept_atoms']} atoms, removed {result['removed_atoms']} atoms -> {output_path.name}"
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
        relax_cycles = data.get("proteinmpnn", {}).get("relax_cycles", 4)
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
                f"Scratch: /scratch/{data.get('scratch_date')}/{data.get('cluster_user')}/{data.get('project_dir_name')}",
                f"RFDiffusion: {data['rfdiffusion'].get('n_shards')} shards × {data['rfdiffusion'].get('designs_per_shard')} designs, GPU queue={data['rfdiffusion'].get('queue')}, GPU cores={data['rfdiffusion'].get('gpu_ncpu')}",
                f"ProteinMPNN: {data['proteinmpnn'].get('n_shards')} shards, relax={data['proteinmpnn'].get('relax_cycles')}, queue={data['proteinmpnn'].get('queue')}",
                f"AfCycDesign: {data['afcyc'].get('n_shards')} shards, GPU queue={data['afcyc'].get('gpu_queue')}, GPU cores={data['afcyc'].get('gpu_ncpu')}, CPU queue={data['afcyc'].get('cpu_queue')}",
                f"PyRosetta: {data['pyrosetta'].get('n_shards')} shards, queue={data['pyrosetta'].get('queue')}",
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
        profile = self._collect_cluster_profile()
        script = profile.get("stage_scripts", {}).get(stage_name, "")
        self._log(f"Submit target workflow: {profile.get('remote_workflow_dir')} :: {stage_name} -> {script}")

        def on_success(_result):
            self.refresh_jobs()

        self._run_background(f"Submit {stage_name}", lambda: cluster_ops.submit_stage(profile, stage_name), on_success)

    def refresh_jobs(self):
        def worker():
            return cluster_ops.list_jobs(self._collect_cluster_profile())

        def on_success(result):
            command_result, jobs = result
            self._log(f"Refresh jobs exit {command_result.returncode}; {len(jobs)} jobs visible.")
            if not jobs:
                detail = (command_result.stdout.strip() or command_result.stderr.strip())
                if detail:
                    self._log("No jobs were parsed from bjobs output; check Raw SSH/LSF output.")
                else:
                    self._log("bjobs returned no visible jobs for the configured user.")
            self._populate_jobs(jobs)

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
            self.data_csv_entries = entries
            self._populate_data_csv_tree(entries)
            self.data_status_var.set(f"Found {len(entries)} merged CSV files")
            self._log(f"Merged CSV scan exit {command_result.returncode}; {len(entries)} files found.")

        self._run_background("Scan merged CSV files", worker, on_success)

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
        if local_path.exists():
            self.data_status_var.set(f"Using local merged CSV for analysis: {local_path}")
            self._load_local_data_csv(local_path)
            return

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

        self.data_status_var.set(f"Downloading merged CSV for analysis to {local_path.parent} ...")
        self._run_background(
            "Download merged CSV for analysis",
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
        self.data_status_var.set(f"Loaded {path.name}; {len(numeric_headers)} numeric columns from {sampled_rows} sampled rows")

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
            candidates, total_passing = self._read_candidate_rows(self.data_loaded_csv_path, metric, candidate_threshold, low, high)
            self._draw_metric_histogram(values, bins, percentile_value, custom_threshold, percentile, low, high, outside_count, range_note)
            self._draw_metric_scatter(values, percentile_value, custom_threshold, percentile, low, high, outside_count, range_note)
            self._populate_candidate_rows(candidates)
            threshold_label = "custom threshold" if custom_threshold is not None else f"P{percentile:g}"
            self.data_status_var.set(
                f"{metric}: n={len(values)}, displayed={len(display_values)}, P{percentile:g} within range={percentile_value:.4g}, "
                f"{threshold_label}>={candidate_threshold:.4g}, range=[{low:.4g}, {high:.4g}], "
                f"outside={outside_count}, candidates={total_passing} (showing {len(candidates)})"
            )
        except Exception as exc:
            self._show_error("Failed to analyze metric", exc)

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
                if value < threshold:
                    continue
                total += 1
                candidate_id = self._candidate_identifier(row)
                summary = self._candidate_summary(row, metric)
                item = (value, data_row_index, candidate_id, summary)
                if len(heap) < limit:
                    heapq.heappush(heap, item)
                elif value > heap[0][0]:
                    heapq.heapreplace(heap, item)
        candidates = sorted(heap, key=lambda item: item[0], reverse=True)
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
        for item in self.jobs_tree.get_children():
            self.jobs_tree.delete(item)
        for job in jobs:
            self.jobs_tree.insert(
                "",
                END,
                iid=job["jobid"],
                values=(job["jobid"], job["stat"], job["queue"], job["job_name"], job["exec_host"], job["user"]),
            )

    def _selected_job_ids(self) -> list[str]:
        return list(self.jobs_tree.selection())

    def kill_selected_jobs(self):
        job_ids = self._selected_job_ids()
        if not job_ids:
            messagebox.showinfo("No jobs selected", "Select one or more jobs first.")
            return
        if not messagebox.askyesno("Confirm bkill", f"Kill selected jobs?\n\n{', '.join(job_ids)}"):
            return
        self._run_background("Kill selected jobs", lambda: cluster_ops.kill_jobs(self._collect_cluster_profile(), job_ids), lambda _: self.refresh_jobs())

    def kill_all_visible_jobs(self):
        job_ids = list(self.jobs_tree.get_children())
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
