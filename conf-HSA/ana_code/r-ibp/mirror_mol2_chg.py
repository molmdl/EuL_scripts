import sys
import os
import argparse


def mirror_mol2(input_path, output_path, axis="x"):
    axis_idx = {"x": 0, "y": 1, "z": 2}[axis.lower()]
    in_atom = False
    with open(input_path) as fin, open(output_path, "w") as fout:
        for line in fin:
            stripped = line.rstrip("\n")
            if stripped.startswith("@<TRIPOS>"):
                fout.write(stripped + "\n")
                in_atom = stripped == "@<TRIPOS>ATOM"
                continue
            if in_atom and stripped == "":
                in_atom = False
            if in_atom and stripped:
                parts = stripped.split()
                parts[2 + axis_idx] = f"{-float(parts[2 + axis_idx]):.4f}"
                fout.write("  " + "  ".join(parts) + "\n")
            else:
                fout.write(stripped + "\n")


def mirror_chg(input_path, output_path, axis="x"):
    axis_idx = {"x": 0, "y": 1, "z": 2}[axis.lower()]
    with open(input_path) as fin, open(output_path, "w") as fout:
        for line in fin:
            stripped = line.rstrip("\n")
            if not stripped:
                fout.write("\n")
                continue
            parts = stripped.split()
            parts[1 + axis_idx] = f"{-float(parts[1 + axis_idx]):12.6f}"
            fout.write("  ".join(parts) + "\n")


def main():
    parser = argparse.ArgumentParser(
        description="Mirror molecule coordinates in mol2 and chg files to generate enantiomer."
    )
    parser.add_argument("input_prefix", help="Input file prefix (e.g. ibp)")
    parser.add_argument("output_prefix", help="Output file prefix (e.g. ibp_r)")
    parser.add_argument("--axis", default="x", choices=["x", "y", "z"],
                        help="Axis to mirror (negate). Default: x")
    parser.add_argument("--dir", default=".", help="Directory containing input files. Default: .")
    args = parser.parse_args()

    d = args.dir
    mol2_in = os.path.join(d, f"{args.input_prefix}.mol2")
    chg_in = os.path.join(d, f"{args.input_prefix}.chg")
    mol2_out = os.path.join(d, f"{args.output_prefix}.mol2")
    chg_out = os.path.join(d, f"{args.output_prefix}.chg")

    if os.path.isfile(mol2_in):
        mirror_mol2(mol2_in, mol2_out, args.axis)
        print(f"Written: {mol2_out}")
    else:
        print(f"Skip mol2: {mol2_in} not found", file=sys.stderr)

    if os.path.isfile(chg_in):
        mirror_chg(chg_in, chg_out, args.axis)
        print(f"Written: {chg_out}")
    else:
        print(f"Skip chg: {chg_in} not found", file=sys.stderr)


if __name__ == "__main__":
    main()
