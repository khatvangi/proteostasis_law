#!/usr/bin/env python

import os
import argparse
import csv
from math import sqrt
from Bio.PDB.MMCIFParser import MMCIFParser
from Bio.PDB.NeighborSearch import NeighborSearch

# Reasonable default list; you can tweak later
DEFAULT_METALS = {
    "ZN", "FE", "CU", "MN", "MG", "CA",
    "CO", "NI", "CD", "HG", "K", "NA"
}

def get_element(atom):
    """
    Robustly infer element symbol from Bio.PDB Atom.
    """
    elem = atom.element
    if elem is None or elem.strip() == "":
        name = atom.get_id().strip()
        # crude but works for most atoms (C, N, O, S, FE, ZN, etc.)
        if len(name) >= 2 and name[0].isalpha() and name[1].isalpha():
            elem = name[0:2].upper()
        elif len(name) >= 1:
            elem = name[0].upper()
        else:
            elem = ""
    return elem.upper()


def is_metal(atom, metals_set):
    elem = get_element(atom)
    return elem in metals_set


def distance(a, b):
    v = a.coord - b.coord
    return sqrt((v * v).sum())


def process_structure(cif_path, metals_set, cutoff, writer):
    parser = MMCIFParser(QUIET=True)
    pdb_id = os.path.basename(cif_path).split(".")[0].upper()

    try:
        structure = parser.get_structure(pdb_id, cif_path)
    except Exception as e:
        print(f"[WARN] Failed to parse {cif_path}: {e}")
        return

    atoms = list(structure.get_atoms())

    # Separate metal atoms and all atoms for NeighborSearch
    metal_atoms = [a for a in atoms if is_metal(a, metals_set)]
    if not metal_atoms:
        return

    ns = NeighborSearch(atoms)

    seen_pairs = set()

    for m in metal_atoms:
        m_elem = get_element(m)
        m_serial = m.get_serial_number()
        m_coord = m.coord

        # Find neighbors within cutoff
        neighbors = ns.search(m_coord, cutoff, level="A")  # atom level

        for n in neighbors:
            # Skip the metal itself
            if n is m:
                continue

            res = n.get_parent()
            chain = res.get_parent()

            # Only standard protein residues (hetflag == " ")
            hetflag, resseq, icode = res.id
            if hetflag.strip() != "":
                continue

            resname = res.get_resname().strip()
            if len(resname) == 0:
                continue

            # Only amino acid residues (3-letter; we can refine later if needed)
            # You can also restrict to H/C/D/E later in analysis.
            if len(resname) != 3:
                continue

            chain_id = chain.id

            # De-duplicate residue–metal pairs
            key = (pdb_id, chain_id, resseq, icode, m_elem)
            if key in seen_pairs:
                continue
            seen_pairs.add(key)

            d = distance(m, n)

            writer.writerow({
                "pdb_id": pdb_id,
                "chain_id": chain_id,
                "resseq": resseq,
                "icode": icode if icode != " " else "",
                "resname": resname,
                "atom_name": n.get_id().strip(),
                "metal_element": m_elem,
                "metal_serial": m_serial,
                "distance": f"{d:.3f}",
                "cif_path": cif_path
            })


def main():
    ap = argparse.ArgumentParser(
        description="Scan mmCIF files for metal-coordinating residues."
    )
    ap.add_argument("--pdb_dir", required=True,
                    help="Directory with .cif files (E. coli structures).")
    ap.add_argument("--out_csv", required=True,
                    help="Output CSV file with metal-contacting residues.")
    ap.add_argument("--cutoff", type=float, default=3.0,
                    help="Distance cutoff (Å) for metal–atom contact.")
    ap.add_argument("--metals", nargs="*",
                    help="Override list of metal element symbols (e.g. ZN FE CU).")

    args = ap.parse_args()
    metals_set = set(args.metals) if args.metals else DEFAULT_METALS

    cif_files = [
        os.path.join(args.pdb_dir, f)
        for f in os.listdir(args.pdb_dir)
        if f.lower().endswith(".cif")
    ]
    cif_files.sort()

    with open(args.out_csv, "w", newline="") as f:
        fieldnames = [
            "pdb_id", "chain_id", "resseq", "icode", "resname", "atom_name",
            "metal_element", "metal_serial", "distance", "cif_path"
        ]
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()

        for i, cif_path in enumerate(cif_files, 1):
            print(f"[{i}/{len(cif_files)}] {os.path.basename(cif_path)}")
            process_structure(cif_path, metals_set, args.cutoff, writer)


if __name__ == "__main__":
    main()

