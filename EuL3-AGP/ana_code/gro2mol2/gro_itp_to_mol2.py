#!/usr/bin/env python3
"""Convert GROMACS GRO + ITP files into Tripos MOL2 format."""

from __future__ import annotations

import argparse
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Tuple

# ------------------------------
# User-facing defaults and flags
# ------------------------------
DEFAULT_OUTPUT_SUFFIX = ".mol2"
DEFAULT_FILE_ENCODING = "utf-8"
DEFAULT_COORDINATE_SCALE = 10.0  # GRO uses nm, MOL2 uses angstrom

DEFAULT_MOL2_MOLECULE_TYPE = "SMALL"
DEFAULT_MOL2_CHARGE_TYPE = "USER_CHARGES"
DEFAULT_MOL2_NUM_FEATURES = 0
DEFAULT_MOL2_NUM_SETS = 0

DEFAULT_MOL2_BOND_TYPE = "1"
DEFAULT_MOL2_AROMATIC_BOND_TYPE = "ar"
DEFAULT_ENABLE_AROMATIC_BONDS = True

DEFAULT_ENABLE_UNIQUE_ATOM_NAMES = True
DEFAULT_UNIQUE_ATOM_NAME_SEPARATOR = ""

DEFAULT_FALLBACK_ATOM_NAME = "X"
DEFAULT_FALLBACK_RESIDUE_NAME = "MOL"
DEFAULT_FALLBACK_RESIDUE_ID = 1
DEFAULT_FALLBACK_MOLECULE_NAME = "Molecule"

DEFAULT_USE_ITP_TYPES_DIRECTLY = False
DEFAULT_FALLBACK_MOL2_ATOM_TYPE = "Du"

# ------------------------------
# GRO parsing constants
# ------------------------------
GRO_TITLE_LINE_INDEX = 0
GRO_ATOM_COUNT_LINE_INDEX = 1
GRO_ATOM_BLOCK_START_INDEX = 2

GRO_RESIDUE_ID_SLICE = slice(0, 5)
GRO_RESIDUE_NAME_SLICE = slice(5, 10)
GRO_ATOM_NAME_SLICE = slice(10, 15)
GRO_ATOM_ID_SLICE = slice(15, 20)
GRO_X_SLICE = slice(20, 28)
GRO_Y_SLICE = slice(28, 36)
GRO_Z_SLICE = slice(36, 44)

GRO_FALLBACK_MIN_TOKENS = 6
GRO_RESIDUE_TOKEN_PATTERN = re.compile(r"^(\d+)(\S+)$")

# ------------------------------
# ITP parsing constants
# ------------------------------
ITP_COMMENT_CHAR = ";"
ITP_PREPROCESSOR_PREFIX = "#"
ITP_SECTION_PATTERN = re.compile(r"^\[\s*([^\]]+)\s*\]$")

ITP_SECTION_MOLECULETYPE = "moleculetype"
ITP_SECTION_ATOMS = "atoms"
ITP_SECTION_BONDS = "bonds"

ITP_MIN_ATOM_FIELDS = 7
ITP_MIN_BOND_FIELDS = 2
ITP_DEFAULT_BOND_FUNCTYPE = "1"

# ------------------------------
# MOL2 formatting constants
# ------------------------------
MOL2_SECTION_MOLECULE = "@<TRIPOS>MOLECULE"
MOL2_SECTION_ATOM = "@<TRIPOS>ATOM"
MOL2_SECTION_BOND = "@<TRIPOS>BOND"
MOL2_SECTION_SUBSTRUCTURE = "@<TRIPOS>SUBSTRUCTURE"

MOL2_SUBSTRUCTURE_TYPE = "GROUP"
MOL2_SUBSTRUCTURE_DICT_TYPE = 0
MOL2_SUBSTRUCTURE_CHAIN = "****"
MOL2_SUBSTRUCTURE_INTER_BONDS = 0
MOL2_SUBSTRUCTURE_STATUS = "ROOT"

# ------------------------------
# Atom type mapping controls
# ------------------------------
ITP_TO_MOL2_TYPE_MAP: Dict[str, str] = {
    "c": "C.3",
    "c1": "C.1",
    "c2": "C.2",
    "c3": "C.3",
    "ca": "C.ar",
    "cp": "C.ar",
    "cq": "C.ar",
    "n": "N.3",
    "n1": "N.1",
    "n2": "N.2",
    "n3": "N.3",
    "n4": "N.3",
    "na": "N.ar",
    "nb": "N.ar",
    "n7": "N.am",
    "o": "O.2",
    "os": "O.2",
    "oh": "O.3",
    "s": "S.3",
    "h": "H",
    "h1": "H",
    "h4": "H",
    "ha": "H",
    "hc": "H",
    "hn": "H",
    "hx": "H",
}

ITP_TYPE_PREFIX_MAP: Dict[str, str] = {
    "eu": "Eu",
}

BOND_FUNCTYPE_TO_MOL2_TYPE_MAP: Dict[str, str] = {
    "1": "1",
}

AROMATIC_MOL2_TYPES = {"C.ar", "N.ar"}

# ------------------------------
# CLI flags
# ------------------------------
CLI_FLAG_GRO = "--gro"
CLI_FLAG_ITP = "--itp"
CLI_FLAG_OUT = "--out"
CLI_FLAG_NAME = "--name"
CLI_FLAG_KEEP_ITP_TYPES = "--keep-itp-types"
CLI_FLAG_NO_AROMATIC = "--no-aromatic-bonds"
CLI_FLAG_NO_UNIQUE_NAMES = "--no-unique-atom-names"


@dataclass(frozen=True)
class GroAtom:
    index: int
    residue_id: int
    residue_name: str
    atom_name: str
    x: float
    y: float
    z: float


@dataclass(frozen=True)
class ItpAtom:
    index: int
    atom_type: str
    residue_id: int
    residue_name: str
    atom_name: str
    charge: float


@dataclass(frozen=True)
class ItpBond:
    atom_i: int
    atom_j: int
    functype: str


@dataclass(frozen=True)
class Mol2Atom:
    index: int
    atom_name: str
    x: float
    y: float
    z: float
    atom_type: str
    substructure_id: int
    substructure_name: str
    charge: float


@dataclass(frozen=True)
class Mol2Bond:
    index: int
    atom_i: int
    atom_j: int
    bond_type: str


@dataclass(frozen=True)
class Mol2Substructure:
    index: int
    name: str
    root_atom: int


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Convert GRO + ITP to MOL2 using coordinates, charges, and bonds."
    )
    parser.add_argument(CLI_FLAG_GRO, required=True, type=Path, help="Input GRO file path")
    parser.add_argument(CLI_FLAG_ITP, required=True, type=Path, help="Input ITP file path")
    parser.add_argument(CLI_FLAG_OUT, type=Path, help="Output MOL2 file path")
    parser.add_argument(
        CLI_FLAG_NAME,
        default="",
        help="Optional molecule name override for MOL2 header",
    )
    parser.add_argument(
        CLI_FLAG_KEEP_ITP_TYPES,
        action="store_true",
        default=DEFAULT_USE_ITP_TYPES_DIRECTLY,
        help="Use ITP atom types directly in MOL2 instead of Sybyl-style mapping",
    )
    parser.add_argument(
        CLI_FLAG_NO_AROMATIC,
        action="store_true",
        default=not DEFAULT_ENABLE_AROMATIC_BONDS,
        help="Disable aromatic bond inference",
    )
    parser.add_argument(
        CLI_FLAG_NO_UNIQUE_NAMES,
        action="store_true",
        default=not DEFAULT_ENABLE_UNIQUE_ATOM_NAMES,
        help="Keep original atom names without appending atom index",
    )
    return parser.parse_args()


def _strip_itp_comment(line: str) -> str:
    return line.split(ITP_COMMENT_CHAR, 1)[0].strip()


def _sanitize_token(value: str, fallback: str) -> str:
    cleaned = "".join(value.split())
    return cleaned if cleaned else fallback


def _parse_gro_atom_line(line: str, line_number: int, coordinate_scale: float) -> GroAtom:
    try:
        residue_id = int(line[GRO_RESIDUE_ID_SLICE])
        residue_name = line[GRO_RESIDUE_NAME_SLICE].strip()
        atom_name = line[GRO_ATOM_NAME_SLICE].strip()
        atom_id = int(line[GRO_ATOM_ID_SLICE])
        x_coord = float(line[GRO_X_SLICE]) * coordinate_scale
        y_coord = float(line[GRO_Y_SLICE]) * coordinate_scale
        z_coord = float(line[GRO_Z_SLICE]) * coordinate_scale
        return GroAtom(
            index=atom_id,
            residue_id=residue_id,
            residue_name=residue_name or DEFAULT_FALLBACK_RESIDUE_NAME,
            atom_name=atom_name or DEFAULT_FALLBACK_ATOM_NAME,
            x=x_coord,
            y=y_coord,
            z=z_coord,
        )
    except ValueError:
        pass

    tokens = line.split()
    if len(tokens) < GRO_FALLBACK_MIN_TOKENS:
        raise ValueError(f"Invalid GRO atom line at {line_number}: {line.rstrip()}")

    residue_token = tokens[0]
    residue_match = GRO_RESIDUE_TOKEN_PATTERN.match(residue_token)
    if residue_match:
        residue_id = int(residue_match.group(1))
        residue_name = residue_match.group(2)
    else:
        residue_id = DEFAULT_FALLBACK_RESIDUE_ID
        residue_name = DEFAULT_FALLBACK_RESIDUE_NAME

    atom_name = tokens[1]
    atom_id = int(tokens[2])
    x_coord = float(tokens[3]) * coordinate_scale
    y_coord = float(tokens[4]) * coordinate_scale
    z_coord = float(tokens[5]) * coordinate_scale

    return GroAtom(
        index=atom_id,
        residue_id=residue_id,
        residue_name=residue_name or DEFAULT_FALLBACK_RESIDUE_NAME,
        atom_name=atom_name or DEFAULT_FALLBACK_ATOM_NAME,
        x=x_coord,
        y=y_coord,
        z=z_coord,
    )


def parse_gro_file(path: Path, coordinate_scale: float) -> Tuple[str, List[GroAtom]]:
    lines = path.read_text(encoding=DEFAULT_FILE_ENCODING).splitlines()
    minimum_required_lines = GRO_ATOM_BLOCK_START_INDEX + 1
    if len(lines) < minimum_required_lines:
        raise ValueError(f"GRO file is too short: {path}")

    title = lines[GRO_TITLE_LINE_INDEX].strip()
    atom_count_text = lines[GRO_ATOM_COUNT_LINE_INDEX].strip()
    if not atom_count_text:
        raise ValueError(f"Missing atom count in GRO file: {path}")

    atom_count = int(atom_count_text.split()[0])
    atom_start = GRO_ATOM_BLOCK_START_INDEX
    atom_end = atom_start + atom_count
    atom_lines = lines[atom_start:atom_end]

    if len(atom_lines) != atom_count:
        raise ValueError(
            f"GRO atom count mismatch in {path}: expected {atom_count}, found {len(atom_lines)}"
        )

    atoms: List[GroAtom] = []
    for offset, atom_line in enumerate(atom_lines):
        line_number = atom_start + offset + 1
        atoms.append(_parse_gro_atom_line(atom_line, line_number, coordinate_scale))

    return title, atoms


def parse_itp_file(path: Path) -> Tuple[str, List[ItpAtom], List[ItpBond]]:
    lines = path.read_text(encoding=DEFAULT_FILE_ENCODING).splitlines()

    current_section = ""
    molecule_name = ""
    atoms: List[ItpAtom] = []
    bonds: List[ItpBond] = []

    for raw_line in lines:
        stripped_line = raw_line.strip()
        if not stripped_line:
            continue
        if stripped_line.startswith(ITP_COMMENT_CHAR):
            continue
        if stripped_line.startswith(ITP_PREPROCESSOR_PREFIX):
            continue

        section_match = ITP_SECTION_PATTERN.match(stripped_line)
        if section_match:
            current_section = section_match.group(1).strip().lower()
            continue

        clean_line = _strip_itp_comment(raw_line)
        if not clean_line:
            continue

        fields = clean_line.split()

        if current_section == ITP_SECTION_MOLECULETYPE:
            if not molecule_name and fields:
                molecule_name = fields[0]
            continue

        if current_section == ITP_SECTION_ATOMS:
            if len(fields) < ITP_MIN_ATOM_FIELDS:
                continue
            atom_index = int(fields[0])
            atom_type = fields[1]
            residue_id = int(fields[2])
            residue_name = fields[3]
            atom_name = fields[4]
            charge = float(fields[6])
            atoms.append(
                ItpAtom(
                    index=atom_index,
                    atom_type=atom_type,
                    residue_id=residue_id,
                    residue_name=residue_name,
                    atom_name=atom_name,
                    charge=charge,
                )
            )
            continue

        if current_section == ITP_SECTION_BONDS:
            if len(fields) < ITP_MIN_BOND_FIELDS:
                continue
            atom_i = int(fields[0])
            atom_j = int(fields[1])
            functype = fields[2] if len(fields) > 2 else ITP_DEFAULT_BOND_FUNCTYPE
            bonds.append(ItpBond(atom_i=atom_i, atom_j=atom_j, functype=functype))

    return molecule_name, atoms, bonds


def _infer_element_symbol(atom_name: str, atom_type: str) -> str:
    atom_name_letters = "".join(ch for ch in atom_name if ch.isalpha())
    atom_type_letters = "".join(ch for ch in atom_type if ch.isalpha())
    candidate = atom_name_letters or atom_type_letters
    if not candidate:
        return "X"

    candidate_upper = candidate.upper()
    if len(candidate_upper) >= 2 and candidate_upper[:2] in {
        "CL",
        "BR",
        "NA",
        "MG",
        "ZN",
        "FE",
        "CA",
        "EU",
    }:
        return candidate_upper[:2].title()

    return candidate_upper[0]


def map_itp_atom_type_to_mol2(
    itp_type: str,
    atom_name: str,
    keep_itp_types: bool,
) -> str:
    normalized_type = itp_type.strip()
    normalized_type_lower = normalized_type.lower()

    if keep_itp_types:
        return _sanitize_token(normalized_type, DEFAULT_FALLBACK_MOL2_ATOM_TYPE)

    if normalized_type_lower in ITP_TO_MOL2_TYPE_MAP:
        return ITP_TO_MOL2_TYPE_MAP[normalized_type_lower]

    for prefix, mapped_type in ITP_TYPE_PREFIX_MAP.items():
        if normalized_type_lower.startswith(prefix):
            return mapped_type

    inferred_element = _infer_element_symbol(atom_name=atom_name, atom_type=itp_type)
    if inferred_element == "X":
        return DEFAULT_FALLBACK_MOL2_ATOM_TYPE
    return inferred_element


def make_atom_name(base_name: str, atom_index: int, use_unique_names: bool) -> str:
    sanitized_base = _sanitize_token(base_name, DEFAULT_FALLBACK_ATOM_NAME)
    if not use_unique_names:
        return sanitized_base
    return f"{sanitized_base}{DEFAULT_UNIQUE_ATOM_NAME_SEPARATOR}{atom_index}"


def determine_bond_type(
    bond: ItpBond,
    atom_type_by_index: Dict[int, str],
    enable_aromatic_bonds: bool,
) -> str:
    atom_i_type = atom_type_by_index.get(bond.atom_i, "")
    atom_j_type = atom_type_by_index.get(bond.atom_j, "")

    if (
        enable_aromatic_bonds
        and atom_i_type in AROMATIC_MOL2_TYPES
        and atom_j_type in AROMATIC_MOL2_TYPES
    ):
        return DEFAULT_MOL2_AROMATIC_BOND_TYPE

    return BOND_FUNCTYPE_TO_MOL2_TYPE_MAP.get(str(bond.functype), DEFAULT_MOL2_BOND_TYPE)


def _build_substructures(
    gro_atoms: List[GroAtom],
) -> Tuple[Dict[int, int], List[Mol2Substructure], Dict[int, str]]:
    residue_key_to_sub_id: Dict[Tuple[int, str], int] = {}
    atom_to_sub_id: Dict[int, int] = {}
    substructures: List[Mol2Substructure] = []
    sub_name_by_id: Dict[int, str] = {}

    for atom in gro_atoms:
        residue_key = (atom.residue_id, atom.residue_name or DEFAULT_FALLBACK_RESIDUE_NAME)
        if residue_key not in residue_key_to_sub_id:
            sub_id = len(residue_key_to_sub_id) + 1
            residue_key_to_sub_id[residue_key] = sub_id
            sub_name = _sanitize_token(residue_key[1], DEFAULT_FALLBACK_RESIDUE_NAME)
            substructures.append(Mol2Substructure(index=sub_id, name=sub_name, root_atom=atom.index))
            sub_name_by_id[sub_id] = sub_name
        atom_to_sub_id[atom.index] = residue_key_to_sub_id[residue_key]

    return atom_to_sub_id, substructures, sub_name_by_id


def build_mol2_records(
    gro_atoms: List[GroAtom],
    itp_atoms: List[ItpAtom],
    itp_bonds: List[ItpBond],
    keep_itp_types: bool,
    enable_aromatic_bonds: bool,
    use_unique_atom_names: bool,
) -> Tuple[List[Mol2Atom], List[Mol2Bond], List[Mol2Substructure]]:
    itp_atom_by_index = {atom.index: atom for atom in itp_atoms}

    gro_indices = {atom.index for atom in gro_atoms}
    itp_indices = set(itp_atom_by_index.keys())
    missing_in_itp = sorted(gro_indices - itp_indices)
    missing_in_gro = sorted(itp_indices - gro_indices)
    if missing_in_itp or missing_in_gro:
        raise ValueError(
            "Atom index mismatch between GRO and ITP. "
            f"Missing in ITP: {missing_in_itp[:10]}, Missing in GRO: {missing_in_gro[:10]}"
        )

    atom_to_sub_id, substructures, sub_name_by_id = _build_substructures(gro_atoms)

    mol2_atoms: List[Mol2Atom] = []
    mol2_atom_type_by_index: Dict[int, str] = {}

    for gro_atom in gro_atoms:
        itp_atom = itp_atom_by_index[gro_atom.index]
        mol2_atom_type = map_itp_atom_type_to_mol2(
            itp_type=itp_atom.atom_type,
            atom_name=gro_atom.atom_name,
            keep_itp_types=keep_itp_types,
        )
        sub_id = atom_to_sub_id[gro_atom.index]
        atom_name_source = itp_atom.atom_name or gro_atom.atom_name
        mol2_atom = Mol2Atom(
            index=gro_atom.index,
            atom_name=make_atom_name(
                base_name=atom_name_source,
                atom_index=gro_atom.index,
                use_unique_names=use_unique_atom_names,
            ),
            x=gro_atom.x,
            y=gro_atom.y,
            z=gro_atom.z,
            atom_type=mol2_atom_type,
            substructure_id=sub_id,
            substructure_name=sub_name_by_id[sub_id],
            charge=itp_atom.charge,
        )
        mol2_atoms.append(mol2_atom)
        mol2_atom_type_by_index[mol2_atom.index] = mol2_atom_type

    mol2_bonds: List[Mol2Bond] = []
    atom_index_set = {atom.index for atom in mol2_atoms}
    for bond_index, itp_bond in enumerate(itp_bonds, start=1):
        if itp_bond.atom_i not in atom_index_set or itp_bond.atom_j not in atom_index_set:
            raise ValueError(
                f"Bond references unknown atom index: {itp_bond.atom_i}-{itp_bond.atom_j}"
            )
        bond_type = determine_bond_type(
            bond=itp_bond,
            atom_type_by_index=mol2_atom_type_by_index,
            enable_aromatic_bonds=enable_aromatic_bonds,
        )
        mol2_bonds.append(
            Mol2Bond(
                index=bond_index,
                atom_i=itp_bond.atom_i,
                atom_j=itp_bond.atom_j,
                bond_type=bond_type,
            )
        )

    return mol2_atoms, mol2_bonds, substructures


def write_mol2_file(
    output_path: Path,
    molecule_name: str,
    atoms: List[Mol2Atom],
    bonds: List[Mol2Bond],
    substructures: List[Mol2Substructure],
) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)

    lines: List[str] = []
    lines.append(MOL2_SECTION_MOLECULE)
    lines.append(_sanitize_token(molecule_name, DEFAULT_FALLBACK_MOLECULE_NAME))
    lines.append(
        f"{len(atoms):>5d} {len(bonds):>5d} {len(substructures):>5d} "
        f"{DEFAULT_MOL2_NUM_FEATURES:>5d} {DEFAULT_MOL2_NUM_SETS:>5d}"
    )
    lines.append(DEFAULT_MOL2_MOLECULE_TYPE)
    lines.append(DEFAULT_MOL2_CHARGE_TYPE)
    lines.append("")

    lines.append(MOL2_SECTION_ATOM)
    for atom in atoms:
        lines.append(
            f"{atom.index:>7d} "
            f"{atom.atom_name:<8s} "
            f"{atom.x:>10.4f} {atom.y:>10.4f} {atom.z:>10.4f} "
            f"{atom.atom_type:<8s} "
            f"{atom.substructure_id:>4d} "
            f"{atom.substructure_name:<8s} "
            f"{atom.charge:>10.6f}"
        )

    lines.append(MOL2_SECTION_BOND)
    for bond in bonds:
        lines.append(
            f"{bond.index:>6d} {bond.atom_i:>5d} {bond.atom_j:>5d} {bond.bond_type:<4s}"
        )

    lines.append(MOL2_SECTION_SUBSTRUCTURE)
    for sub in substructures:
        lines.append(
            f"{sub.index:>5d} "
            f"{sub.name:<8s} "
            f"{sub.root_atom:>5d} "
            f"{MOL2_SUBSTRUCTURE_TYPE:<8s} "
            f"{MOL2_SUBSTRUCTURE_DICT_TYPE:>3d} "
            f"{MOL2_SUBSTRUCTURE_CHAIN:<4s} "
            f"{MOL2_SUBSTRUCTURE_INTER_BONDS:>3d} "
            f"{MOL2_SUBSTRUCTURE_STATUS}"
        )

    output_path.write_text("\n".join(lines) + "\n", encoding=DEFAULT_FILE_ENCODING)


def resolve_output_path(gro_path: Path, output_path: Path | None) -> Path:
    if output_path is not None:
        return output_path
    return gro_path.with_suffix(DEFAULT_OUTPUT_SUFFIX)


def main() -> None:
    args = parse_arguments()

    gro_path: Path = args.gro
    itp_path: Path = args.itp
    output_path: Path = resolve_output_path(gro_path=gro_path, output_path=args.out)

    gro_title, gro_atoms = parse_gro_file(
        path=gro_path,
        coordinate_scale=DEFAULT_COORDINATE_SCALE,
    )
    itp_name, itp_atoms, itp_bonds = parse_itp_file(path=itp_path)

    if not itp_atoms:
        raise ValueError(f"No [ atoms ] records found in ITP file: {itp_path}")
    if not itp_bonds:
        raise ValueError(f"No [ bonds ] records found in ITP file: {itp_path}")

    molecule_name = args.name or itp_name or gro_title or DEFAULT_FALLBACK_MOLECULE_NAME
    enable_aromatic_bonds = DEFAULT_ENABLE_AROMATIC_BONDS and (not args.no_aromatic_bonds)
    use_unique_names = DEFAULT_ENABLE_UNIQUE_ATOM_NAMES and (not args.no_unique_atom_names)

    mol2_atoms, mol2_bonds, mol2_substructures = build_mol2_records(
        gro_atoms=gro_atoms,
        itp_atoms=itp_atoms,
        itp_bonds=itp_bonds,
        keep_itp_types=args.keep_itp_types,
        enable_aromatic_bonds=enable_aromatic_bonds,
        use_unique_atom_names=use_unique_names,
    )

    write_mol2_file(
        output_path=output_path,
        molecule_name=molecule_name,
        atoms=mol2_atoms,
        bonds=mol2_bonds,
        substructures=mol2_substructures,
    )

    print(f"Wrote MOL2: {output_path}")
    print(
        f"Atoms: {len(mol2_atoms)} | Bonds: {len(mol2_bonds)} | "
        f"Substructures: {len(mol2_substructures)}"
    )


if __name__ == "__main__":
    main()
