#!/usr/bin/env python3
import json
import sys
import traceback
from copy import deepcopy
from pathlib import Path
from tkinter import (
    BooleanVar,
    Button,
    Checkbutton,
    END,
    Entry,
    Frame,
    Label,
    LabelFrame,
    StringVar,
    Tk,
    Text,
    filedialog,
    messagebox,
    ttk,
)

SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

import generate_workflow


FIELD_GROUPS = [
    (
        "Project",
        [
            ("target", "Target name", str),
            ("pilot", "Pilot name", str),
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
        self.variables = {}

        self._build()
        self._load_values(self.config_data)
        self._log(f"Loaded config: {self.config_path}")

    def _build(self):
        top = Frame(self.root)
        top.pack(fill="x", padx=12, pady=10)

        self.config_path_var = StringVar(value=str(self.config_path))
        self.output_dir_var = StringVar(value=str(self.output_dir))

        Label(top, text="Config").grid(row=0, column=0, sticky="w")
        Entry(top, textvariable=self.config_path_var).grid(row=0, column=1, sticky="ew", padx=6)
        Button(top, text="Open", command=self.open_config).grid(row=0, column=2, padx=3)
        Button(top, text="Save", command=self.save_config).grid(row=0, column=3, padx=3)
        Button(top, text="Save as", command=self.save_config_as).grid(row=0, column=4, padx=3)

        Label(top, text="Output").grid(row=1, column=0, sticky="w", pady=(8, 0))
        Entry(top, textvariable=self.output_dir_var).grid(row=1, column=1, sticky="ew", padx=6, pady=(8, 0))
        Button(top, text="Choose", command=self.choose_output_dir).grid(row=1, column=2, padx=3, pady=(8, 0))
        Button(top, text="Generate workflow", command=self.generate).grid(row=1, column=3, columnspan=2, sticky="ew", padx=3, pady=(8, 0))
        top.columnconfigure(1, weight=1)

        notebook = ttk.Notebook(self.root)
        notebook.pack(fill="both", expand=True, padx=12, pady=(0, 10))

        for group_name, fields in FIELD_GROUPS:
            frame = Frame(notebook)
            notebook.add(frame, text=group_name)
            self._build_group(frame, group_name, fields)

        log_frame = LabelFrame(self.root, text="Log")
        log_frame.pack(fill="both", expand=False, padx=12, pady=(0, 12))
        self.log_text = Text(log_frame, height=8, wrap="word")
        self.log_text.pack(fill="both", expand=True, padx=6, pady=6)

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

    def generate(self):
        try:
            data = self._collect_config()
            output_dir = Path(self.output_dir_var.get()).expanduser()
            written = generate_workflow.generate(data, output_dir, force=True)
            self.config_data = data
            self.output_dir = output_dir
            self._log("Generated workflow files:")
            for path in written:
                self._log(f"  {path}")
            messagebox.showinfo("Workflow generated", f"Generated {len(written)} submit scripts.")
        except Exception as exc:
            self._show_error("Failed to generate workflow", exc)

    def _log(self, message: str):
        self.log_text.insert(END, message + "\n")
        self.log_text.see(END)

    def _show_error(self, title: str, exc: Exception):
        detail = "".join(traceback.format_exception_only(type(exc), exc)).strip()
        self._log(f"{title}: {detail}")
        messagebox.showerror(title, detail)


def main():
    root = Tk()
    WorkflowGui(root)
    root.mainloop()


if __name__ == "__main__":
    main()
