#!/usr/bin/env python3
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path


AA3_TO_1 = {
    "ALA": "A",
    "ARG": "R",
    "ASN": "N",
    "ASP": "D",
    "CYS": "C",
    "GLN": "Q",
    "GLU": "E",
    "GLY": "G",
    "HIS": "H",
    "ILE": "I",
    "LEU": "L",
    "LYS": "K",
    "MET": "M",
    "PHE": "F",
    "PRO": "P",
    "SER": "S",
    "THR": "T",
    "TRP": "W",
    "TYR": "Y",
    "VAL": "V",
    "SEC": "U",
    "PYL": "O",
    "MSE": "M",
}


@dataclass(frozen=True)
class ResidueRecord:
    chain: str
    resseq: str
    icode: str
    resname: str
    aa: str
    atom_count: int
    is_protein: bool


def _line_field(line: str, start: int, end: int) -> str:
    return line.rstrip("\n").ljust(80)[start:end]


def _residue_key(line: str) -> tuple[str, str, str, str]:
    return (
        _line_field(line, 21, 22).strip() or "_",
        _line_field(line, 22, 26).strip(),
        _line_field(line, 26, 27).strip(),
        _line_field(line, 17, 20).strip().upper(),
    )


def _is_atom_record(line: str) -> bool:
    return line.startswith("ATOM  ") or line.startswith("HETATM")


def parse_pdb(path: Path) -> dict:
    path = Path(path)
    residues: list[ResidueRecord] = []
    residue_order: list[tuple[str, str, str, str]] = []
    atom_counts: dict[tuple[str, str, str, str], int] = {}
    total_atoms = 0
    protein_atoms = 0
    nonprotein_atoms = 0
    nonprotein_residues: dict[str, int] = {}

    with path.open("r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if not _is_atom_record(line):
                continue
            total_atoms += 1
            key = _residue_key(line)
            if key not in atom_counts:
                residue_order.append(key)
                atom_counts[key] = 0
            atom_counts[key] += 1
            resname = key[3]
            if resname in AA3_TO_1:
                protein_atoms += 1
            else:
                nonprotein_atoms += 1
                nonprotein_residues[resname] = nonprotein_residues.get(resname, 0) + 1

    for key in residue_order:
        chain, resseq, icode, resname = key
        is_protein = resname in AA3_TO_1
        residues.append(
            ResidueRecord(
                chain=chain,
                resseq=resseq,
                icode=icode,
                resname=resname,
                aa=AA3_TO_1.get(resname, "X"),
                atom_count=atom_counts[key],
                is_protein=is_protein,
            )
        )

    protein_residues = [residue for residue in residues if residue.is_protein]
    return {
        "path": str(path),
        "total_atoms": total_atoms,
        "protein_atoms": protein_atoms,
        "nonprotein_atoms": nonprotein_atoms,
        "residues": residues,
        "protein_residues": protein_residues,
        "nonprotein_residues": nonprotein_residues,
        "gaps": residue_number_gaps(protein_residues),
    }


def residue_number_gaps(residues: list[ResidueRecord]) -> list[str]:
    gaps = []
    previous_by_chain: dict[str, int] = {}
    for residue in residues:
        if not residue.resseq.lstrip("-").isdigit() or residue.icode:
            previous_by_chain.pop(residue.chain, None)
            continue
        current = int(residue.resseq)
        previous = previous_by_chain.get(residue.chain)
        if previous is not None and current > previous + 1:
            gaps.append(f"chain {residue.chain}: {previous} -> {current}")
        previous_by_chain[residue.chain] = current
    return gaps


def sequence_blocks(residues: list[ResidueRecord], line_width: int = 60) -> str:
    chains: dict[str, list[ResidueRecord]] = {}
    for residue in residues:
        chains.setdefault(residue.chain, []).append(residue)

    lines = []
    for chain, chain_residues in chains.items():
        lines.append(f"Chain {chain}: {len(chain_residues)} protein residues")
        sequence = "".join(residue.aa for residue in chain_residues)
        for offset in range(0, len(sequence), line_width):
            chunk = sequence[offset:offset + line_width]
            start_residue = chain_residues[offset]
            end_residue = chain_residues[min(offset + len(chunk) - 1, len(chain_residues) - 1)]
            start_label = f"{start_residue.resseq}{start_residue.icode}".strip()
            end_label = f"{end_residue.resseq}{end_residue.icode}".strip()
            lines.append(f"  {start_label:>6} {chunk} {end_label:<6}")
        lines.append("")
    return "\n".join(lines).rstrip()


def summary_text(summary: dict) -> str:
    residues = summary["protein_residues"]
    nonprotein = summary["nonprotein_residues"]
    lines = [
        f"File: {summary['path']}",
        f"Atoms: total={summary['total_atoms']}, protein={summary['protein_atoms']}, non-protein={summary['nonprotein_atoms']}",
        f"Protein residues: {len(residues)}",
    ]
    if nonprotein:
        components = ", ".join(f"{name}:{count}" for name, count in sorted(nonprotein.items()))
        lines.append(f"Non-protein atom counts: {components}")
    else:
        lines.append("Non-protein atom counts: none")
    gaps = summary["gaps"]
    if gaps:
        lines.append(f"Residue number gaps: {len(gaps)} ({'; '.join(gaps[:8])}{' ...' if len(gaps) > 8 else ''})")
    else:
        lines.append("Residue number gaps: none detected")
    lines.append("")
    lines.append(sequence_blocks(residues))
    return "\n".join(lines).rstrip()


def clean_pdb(
    input_path: Path,
    output_path: Path,
    keep_protein_only: bool = True,
    renumber_atoms: bool = True,
    renumber_residues: bool = False,
) -> dict:
    input_path = Path(input_path)
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    atom_serial = 1
    kept_atoms = 0
    removed_atoms = 0
    residue_map: dict[tuple[str, str, str, str], int] = {}
    residue_counter_by_chain: dict[str, int] = {}
    last_chain = None
    wrote_any_atom = False

    with input_path.open("r", encoding="utf-8", errors="replace") as source, output_path.open("w", encoding="utf-8", newline="\n") as destination:
        destination.write(f"REMARK Cleaned from {input_path.name} by RFpeptide workflow frame\n")
        destination.write(f"REMARK protein_only={keep_protein_only} renumber_atoms={renumber_atoms} renumber_residues={renumber_residues}\n")
        for line in source:
            if not _is_atom_record(line):
                continue
            key = _residue_key(line)
            chain, _resseq, _icode, resname = key
            if keep_protein_only and resname not in AA3_TO_1:
                removed_atoms += 1
                continue

            if wrote_any_atom and chain != last_chain:
                destination.write("TER\n")

            padded = line.rstrip("\n").ljust(80)
            record = "ATOM  " if resname in AA3_TO_1 else padded[:6]
            if renumber_atoms:
                new_line = f"{record}{atom_serial:5d}{padded[11:]}"
            else:
                new_line = f"{record}{padded[6:]}"

            if renumber_residues:
                if key not in residue_map:
                    residue_counter_by_chain[chain] = residue_counter_by_chain.get(chain, 0) + 1
                    residue_map[key] = residue_counter_by_chain[chain]
                new_resseq = residue_map[key]
                new_line = f"{new_line[:21]}{chain if chain != '_' else ' '}{new_resseq:4d} {new_line[27:]}"

            destination.write(new_line.rstrip() + "\n")
            kept_atoms += 1
            atom_serial += 1
            last_chain = chain
            wrote_any_atom = True

        if wrote_any_atom:
            destination.write("TER\n")
        destination.write("END\n")

    return {
        "input_path": str(input_path),
        "output_path": str(output_path),
        "kept_atoms": kept_atoms,
        "removed_atoms": removed_atoms,
        "renumbered_residues": len(residue_map),
    }
