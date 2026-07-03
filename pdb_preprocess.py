#!/usr/bin/env python3
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import re
import shutil
from urllib.request import Request, urlopen


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


def normalize_pdb_id(pdb_id: str) -> str:
    value = str(pdb_id or "").strip().upper()
    if not re.fullmatch(r"[A-Z0-9]{4}", value):
        raise ValueError("PDB ID must contain exactly four letters/numbers")
    return value


def fetch_rcsb_pdb(pdb_id: str, destination_dir: Path, timeout: int = 60) -> dict:
    normalized_id = normalize_pdb_id(pdb_id)
    destination_dir = Path(destination_dir)
    destination_dir.mkdir(parents=True, exist_ok=True)
    output_path = destination_dir / f"{normalized_id}.pdb"
    temporary_path = output_path.with_suffix(".pdb.download")
    url = f"https://files.rcsb.org/download/{normalized_id}.pdb"
    request = Request(url, headers={"User-Agent": "RFpeptide-workflow-frame/experimental-pdb-fetch"})
    try:
        with urlopen(request, timeout=timeout) as response, temporary_path.open("wb") as handle:
            shutil.copyfileobj(response, handle)
        if temporary_path.stat().st_size < 100:
            raise ValueError(f"RCSB returned an unexpectedly small file for {normalized_id}")
        temporary_path.replace(output_path)
    except Exception:
        if temporary_path.exists():
            temporary_path.unlink()
        raise
    return {
        "pdb_id": normalized_id,
        "url": url,
        "path": str(output_path),
        "size": output_path.stat().st_size,
    }


def snapshot_original_pdb(path: Path, destination_dir: Path) -> Path:
    source = Path(path).resolve()
    destination_dir = Path(destination_dir)
    destination_dir.mkdir(parents=True, exist_ok=True)
    destination = destination_dir / f"{source.stem}_original{source.suffix or '.pdb'}"
    if source != destination.resolve():
        shutil.copy2(source, destination)
    return destination.resolve()


def _header_text(lines: list[str], record: str, start: int = 10) -> str:
    values = []
    for line in lines:
        if line.startswith(record):
            text = line[start:].strip()
            if text:
                values.append(text)
    return " ".join(values).strip()


def _parse_compound_chains(lines: list[str]) -> dict[str, str]:
    text = _header_text(lines, "COMPND")
    chain_descriptions: dict[str, str] = {}
    current_molecule = ""
    current_chains: list[str] = []
    for field in text.split(";"):
        key, separator, value = field.partition(":")
        if not separator:
            continue
        key = key.strip().upper()
        value = value.strip()
        if key == "MOLECULE":
            current_molecule = value
        elif key == "CHAIN":
            current_chains = [chain.strip() for chain in value.split(",") if chain.strip()]
            for chain in current_chains:
                chain_descriptions[chain] = current_molecule or "Unspecified molecule"
        elif key == "MOL_ID":
            current_molecule = ""
            current_chains = []
    return chain_descriptions


def _parse_dbrefs(lines: list[str]) -> dict[str, list[str]]:
    references: dict[str, list[str]] = {}
    for line in lines:
        if not line.startswith("DBREF "):
            continue
        padded = line.rstrip("\n").ljust(80)
        chain = padded[12:13].strip() or "_"
        database = padded[26:32].strip()
        accession = padded[33:41].strip()
        identifier = padded[42:54].strip()
        label = " ".join(part for part in (database, accession, identifier) if part)
        if label and label not in references.setdefault(chain, []):
            references[chain].append(label)
    return references


def _parse_hetnam(lines: list[str]) -> dict[str, str]:
    names: dict[str, list[str]] = {}
    for line in lines:
        if not line.startswith("HETNAM"):
            continue
        padded = line.rstrip("\n").ljust(80)
        component = padded[11:14].strip()
        description = padded[15:70].strip()
        if component and description:
            names.setdefault(component, []).append(description)
    return {component: " ".join(parts) for component, parts in names.items()}


def _parse_assemblies(lines: list[str]) -> list[dict]:
    assemblies = []
    current_ids: list[str] = []
    descriptions: list[str] = []

    def flush():
        nonlocal current_ids, descriptions
        if current_ids:
            description = " ".join(descriptions).strip() or "Biological assembly"
            for assembly_id in current_ids:
                assemblies.append({"id": assembly_id, "description": description})
        current_ids = []
        descriptions = []

    for line in lines:
        if not line.startswith("REMARK 350"):
            continue
        text = line[10:].strip()
        if text.startswith("BIOMOLECULE:"):
            flush()
            current_ids = [value.strip() for value in text.split(":", 1)[1].split(",") if value.strip()]
        elif current_ids and any(
            text.startswith(prefix)
            for prefix in (
                "AUTHOR DETERMINED BIOLOGICAL UNIT:",
                "SOFTWARE DETERMINED QUATERNARY STRUCTURE:",
                "SOFTWARE USED:",
            )
        ):
            descriptions.append(text)
    flush()
    return assemblies


def parse_pdb(path: Path) -> dict:
    path = Path(path)
    residues: list[ResidueRecord] = []
    residue_order: list[tuple[str, str, str, str]] = []
    atom_counts: dict[tuple[str, str, str, str], int] = {}
    total_atoms = 0
    protein_atoms = 0
    nonprotein_atoms = 0
    nonprotein_residues: dict[str, int] = {}
    component_residues: dict[str, set[tuple[str, str, str]]] = {}
    component_chains: dict[str, set[str]] = {}
    altloc_atoms = 0
    altloc_ids: set[str] = set()
    model_ids: set[str] = set()
    header_lines = []
    current_model = "1"

    with path.open("r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.startswith(("TITLE ", "COMPND", "DBREF ", "HETNAM", "REMARK 350")):
                header_lines.append(line)
            if line.startswith("MODEL "):
                current_model = line[10:14].strip() or str(len(model_ids) + 1)
                model_ids.add(current_model)
                continue
            if not _is_atom_record(line):
                continue
            model_ids.add(current_model)
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
                component_residues.setdefault(resname, set()).add((key[0], key[1], key[2]))
                component_chains.setdefault(resname, set()).add(key[0])
            altloc = _line_field(line, 16, 17).strip()
            if altloc:
                altloc_atoms += 1
                altloc_ids.add(altloc)

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
    chain_descriptions = _parse_compound_chains(header_lines)
    chain_dbrefs = _parse_dbrefs(header_lines)
    hetnam = _parse_hetnam(header_lines)
    chains = {}
    for residue in residues:
        chain = residue.chain
        info = chains.setdefault(
            chain,
            {
                "chain": chain,
                "description": chain_descriptions.get(chain, "Unspecified molecule"),
                "dbrefs": chain_dbrefs.get(chain, []),
                "protein_residues": 0,
                "nonprotein_residues": 0,
                "atoms": 0,
            },
        )
        info["atoms"] += residue.atom_count
        if residue.is_protein:
            info["protein_residues"] += 1
        else:
            info["nonprotein_residues"] += 1
    components = []
    for component, atoms in sorted(nonprotein_residues.items()):
        components.append(
            {
                "id": component,
                "description": hetnam.get(component, "No HETNAM description"),
                "atoms": atoms,
                "residue_instances": len(component_residues.get(component, set())),
                "chains": sorted(component_chains.get(component, set())),
            }
        )
    return {
        "path": str(path),
        "title": _header_text(header_lines, "TITLE "),
        "total_atoms": total_atoms,
        "protein_atoms": protein_atoms,
        "nonprotein_atoms": nonprotein_atoms,
        "residues": residues,
        "protein_residues": protein_residues,
        "nonprotein_residues": nonprotein_residues,
        "components": components,
        "chains": chains,
        "assemblies": _parse_assemblies(header_lines),
        "altloc_atoms": altloc_atoms,
        "altloc_ids": sorted(altloc_ids),
        "model_count": len(model_ids) or 1,
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


def component_summary_text(summary: dict) -> str:
    lines = [
        f"Title: {summary.get('title') or 'Not provided in PDB header'}",
        f"Models: {summary.get('model_count', 1)}",
        "",
        "CHAIN INFORMATION",
    ]
    chains = summary.get("chains", {})
    if chains:
        for chain in sorted(chains):
            info = chains[chain]
            references = ", ".join(info.get("dbrefs", [])) or "none"
            lines.extend(
                [
                    f"[Chain {chain}] {info.get('description') or 'Unspecified molecule'}",
                    f"  Protein residues: {info.get('protein_residues', 0)} | Other residues: {info.get('nonprotein_residues', 0)} | Atoms: {info.get('atoms', 0)}",
                    f"  Database references: {references}",
                ]
            )
    else:
        lines.append("No chains found.")

    lines.extend(["", "NON-STANDARD COMPONENTS"])
    components = summary.get("components", [])
    if components:
        for component in components:
            chains_text = ", ".join(component.get("chains", [])) or "unassigned"
            lines.extend(
                [
                    f"[{component.get('id')}] {component.get('description')}",
                    f"  Residue instances: {component.get('residue_instances', 0)} | Atoms: {component.get('atoms', 0)} | Chains: {chains_text}",
                ]
            )
    else:
        lines.append("None detected.")

    lines.extend(["", "BIOLOGICAL ASSEMBLIES"])
    assemblies = summary.get("assemblies", [])
    if assemblies:
        for assembly in assemblies:
            lines.append(f"[Assembly {assembly.get('id')}] {assembly.get('description')}")
    else:
        lines.append("No REMARK 350 assembly records found.")

    lines.extend(["", "ALTERNATE LOCATIONS"])
    altloc_atoms = summary.get("altloc_atoms", 0)
    if altloc_atoms:
        identifiers = ", ".join(summary.get("altloc_ids", [])) or "unspecified"
        lines.append(f"{altloc_atoms} atoms have alternate locations (IDs: {identifiers}).")
    else:
        lines.append("No alternate atom locations detected.")
    return "\n".join(lines).rstrip()


def _write_wrapped_pdb_record(destination, record: str, text: str, width: int = 69) -> None:
    words = str(text or "").split()
    if not words:
        return
    lines = []
    current = ""
    for word in words:
        candidate = f"{current} {word}".strip()
        if current and len(candidate) > width:
            lines.append(current)
            current = word
        else:
            current = candidate
    if current:
        lines.append(current)
    for index, value in enumerate(lines, start=1):
        continuation = "  " if index == 1 else f"{index:2d}"
        destination.write(f"{record}{continuation}  {value}\n")


def _write_cleaned_metadata(
    destination,
    summary: dict,
    source_lines: list[str],
    chain_map: dict[str, str],
    kept_components: set[str],
) -> None:
    header = next((line.rstrip("\n") for line in source_lines if line.startswith("HEADER")), "")
    if header:
        destination.write(header + "\n")
    _write_wrapped_pdb_record(destination, "TITLE ", summary.get("title", ""))

    molecule_index = 1
    for old_chain, new_chain in chain_map.items():
        info = summary.get("chains", {}).get(old_chain, {})
        description = info.get("description") or "Unspecified molecule"
        fields = f"MOL_ID: {molecule_index}; MOLECULE: {description}; CHAIN: {new_chain};"
        _write_wrapped_pdb_record(destination, "COMPND", fields)
        molecule_index += 1

    for line in source_lines:
        if not line.startswith("DBREF "):
            continue
        padded = line.rstrip("\n").ljust(80)
        old_chain = padded[12:13].strip() or "_"
        if old_chain not in chain_map:
            continue
        new_chain = chain_map[old_chain]
        destination.write(f"{padded[:12]}{new_chain if new_chain != '_' else ' '}{padded[13:]}".rstrip() + "\n")

    for line in source_lines:
        if not line.startswith(("HETNAM", "FORMUL")):
            continue
        padded = line.rstrip("\n").ljust(80)
        component = padded[11:14].strip() if line.startswith("HETNAM") else padded[12:15].strip()
        if component in kept_components:
            destination.write(line.rstrip("\n") + "\n")


def clean_pdb(
    input_path: Path,
    output_path: Path,
    keep_protein_only: bool = True,
    renumber_atoms: bool = True,
    renumber_residues: bool = False,
    keep_chains: set[str] | None = None,
    renumber_chains: bool = False,
) -> dict:
    input_path = Path(input_path)
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    source_lines = input_path.read_text(encoding="utf-8", errors="replace").splitlines(keepends=True)
    summary = parse_pdb(input_path)
    chain_order = list(summary.get("chains", {}))
    selected_chain_order = [
        chain
        for chain in chain_order
        if keep_chains is None or chain in keep_chains
    ]
    if renumber_chains and len(selected_chain_order) > 26:
        raise ValueError("PDB chain renumbering supports at most 26 retained chains")
    chain_map = {
        chain: chr(ord("A") + index) if renumber_chains else chain
        for index, chain in enumerate(selected_chain_order)
    }
    if not chain_map:
        raise ValueError("No chains remain after applying the chain selection")

    atom_serial = 1
    kept_atoms = 0
    removed_atoms = 0
    removed_chain_atoms = 0
    residue_map: dict[tuple[str, str, str, str], int] = {}
    residue_counter_by_chain: dict[str, int] = {}
    last_chain = None
    wrote_any_atom = False
    kept_components: set[str] = set()

    with output_path.open("w", encoding="utf-8", newline="\n") as destination:
        destination.write(f"REMARK Cleaned from {input_path.name} by RFpeptide workflow frame\n")
        destination.write(
            f"REMARK protein_only={keep_protein_only} renumber_atoms={renumber_atoms} "
            f"renumber_residues={renumber_residues} renumber_chains={renumber_chains}\n"
        )
        _write_cleaned_metadata(destination, summary, source_lines, chain_map, kept_components=set())
        atom_lines = []
        for line in source_lines:
            if not _is_atom_record(line):
                continue
            key = _residue_key(line)
            original_chain, _resseq, _icode, resname = key
            if original_chain not in chain_map:
                removed_atoms += 1
                removed_chain_atoms += 1
                continue
            if keep_protein_only and resname not in AA3_TO_1:
                removed_atoms += 1
                continue
            if resname not in AA3_TO_1:
                kept_components.add(resname)
            output_chain = chain_map[original_chain]

            if wrote_any_atom and output_chain != last_chain:
                atom_lines.append("TER\n")

            padded = line.rstrip("\n").ljust(80)
            record = "ATOM  " if resname in AA3_TO_1 else padded[:6]
            if renumber_atoms:
                new_line = f"{record}{atom_serial:5d}{padded[11:]}"
            else:
                new_line = f"{record}{padded[6:]}"
            new_line = f"{new_line[:21]}{output_chain if output_chain != '_' else ' '}{new_line[22:]}"

            if renumber_residues:
                if key not in residue_map:
                    residue_counter_by_chain[output_chain] = residue_counter_by_chain.get(output_chain, 0) + 1
                    residue_map[key] = residue_counter_by_chain[output_chain]
                new_resseq = residue_map[key]
                new_line = f"{new_line[:22]}{new_resseq:4d} {new_line[27:]}"

            atom_lines.append(new_line.rstrip() + "\n")
            kept_atoms += 1
            atom_serial += 1
            last_chain = output_chain
            wrote_any_atom = True

        if kept_components:
            destination.seek(0)
            destination.truncate()
            destination.write(f"REMARK Cleaned from {input_path.name} by RFpeptide workflow frame\n")
            destination.write(
                f"REMARK protein_only={keep_protein_only} renumber_atoms={renumber_atoms} "
                f"renumber_residues={renumber_residues} renumber_chains={renumber_chains}\n"
            )
            _write_cleaned_metadata(destination, summary, source_lines, chain_map, kept_components)
        destination.writelines(atom_lines)
        if wrote_any_atom:
            destination.write("TER\n")
        destination.write("END\n")

    if not wrote_any_atom:
        output_path.unlink(missing_ok=True)
        raise ValueError("No atoms remained after applying chain/component filters")

    return {
        "input_path": str(input_path),
        "output_path": str(output_path),
        "kept_atoms": kept_atoms,
        "removed_atoms": removed_atoms,
        "removed_chain_atoms": removed_chain_atoms,
        "renumbered_residues": len(residue_map),
        "kept_chains": selected_chain_order,
        "chain_map": chain_map,
    }
