"""
Multi-candidate atom-inventory sweep for the 1JNX structure prep.

You're not sure which saved file corresponds to "structural completeness
confirmed by atom-inventory check" in the Methods paragraph (i.e. the state
right after PDBFixer, before any later FoldX RepairPDB / capping steps).
Rather than guess, this checks every candidate you have at once and prints
them side by side, so the differences settle it:

  - Residue count and numbering gaps, per candidate.
  - Whether ACE/NME cap residues are present (a capped file will show these;
    a pre-capping file will not - this is the main signal for telling
    file 1 and file 3 apart).
  - Target residues 1646-1648, 1694, 1817-1819: present/absent in the RAW
    file vs each candidate.
  - PDBFixer's own findMissingAtoms() re-run against each candidate.

Edit RAW_PDB_PATH and the CANDIDATES list below, then run as one cell/block.
"""

RAW_PDB_PATH = "/content/drive/MyDrive/PATH/TO/1JNX.pdb"  # <-- EDIT: raw, unprocessed RCSB download

CANDIDATES = {
    "1JNX_fixed_v2_Repair_1.pdb": "/content/drive/MyDrive/brca1_md_v2/inputs/1JNX_fixed_v2_Repair_1.pdb",
    "1JNX_capped.pdb":            "/content/drive/MyDrive/brca1_md_v2/inputs/1JNX_capped.pdb",
    "WT_1JNX_fixed_v2_Repair_1.pdb": "/content/drive/MyDrive/brca1_md_v2/inputs/WT_1JNX_fixed_v2_Repair_1.pdb",
    "WT_1JNX_capped.pdb":         "/content/drive/MyDrive/brca1_md_v2/inputs/WT_1JNX_capped.pdb",
}

import os

try:
    from pdbfixer import PDBFixer
except ImportError:
    raise SystemExit("pdbfixer not installed. Run: !pip install pdbfixer -q")

try:
    from Bio.PDB import PDBParser
except ImportError:
    raise SystemExit("biopython not installed. Run: !pip install biopython -q")

if not os.path.exists(RAW_PDB_PATH):
    raise SystemExit(f"RAW_PDB_PATH does not exist: {RAW_PDB_PATH}\nEdit the path at the top and re-run.")

missing_candidates = {name: path for name, path in CANDIDATES.items() if not os.path.exists(path)}
if missing_candidates:
    print("WARNING - these candidate paths were not found and will be skipped:")
    for name, path in missing_candidates.items():
        print(f"  {name}: {path}")
CANDIDATES = {name: path for name, path in CANDIDATES.items() if os.path.exists(path)}
if not CANDIDATES:
    raise SystemExit("None of the candidate paths exist. Check the CANDIDATES dict and re-run.")

TARGET_RESIDUES = [1646, 1647, 1648, 1649, 1694, 1817, 1818, 1819]
CAP_RESNAMES = {"ACE", "NME", "NHE", "NH2"}

parser = PDBParser(QUIET=True)


def inspect(pdb_path, label):
    structure = parser.get_structure(label, pdb_path)
    model = structure[0]
    chain_data = {}
    cap_residues_found = []
    for chain in model.get_chains():
        standard_residues = [r for r in chain if r.id[0] == " "]
        all_residues = list(chain)
        for r in all_residues:
            resname = r.get_resname().strip()
            if resname in CAP_RESNAMES:
                cap_residues_found.append((chain.id, resname, r.id[1]))
        if not standard_residues:
            continue
        nums = sorted(r.id[1] for r in standard_residues)
        chain_data[chain.id] = nums
    return chain_data, cap_residues_found


print("=" * 78)
print("RAW FILE")
print("=" * 78)
raw_chain_data, raw_caps = inspect(RAW_PDB_PATH, "raw")
raw_all = set()
for nums in raw_chain_data.values():
    raw_all.update(nums)
for chain_id, nums in raw_chain_data.items():
    print(f"Chain {chain_id}: {len(nums)} residues, range {nums[0]}-{nums[-1]}")

print("\n")

for name, path in CANDIDATES.items():
    print("=" * 78)
    print(f"CANDIDATE: {name}")
    print("=" * 78)

    chain_data, caps_found = inspect(path, name)
    for chain_id, nums in chain_data.items():
        print(f"\nChain {chain_id}: {len(nums)} residues")
        print(f"  Residue number range: {nums[0]} to {nums[-1]}")
        gaps = [(a, b) for a, b in zip(nums, nums[1:]) if b != a + 1]
        print(f"  Internal numbering gaps: {gaps if gaps else 'none, numerically contiguous'}")

    if caps_found:
        print(f"\n  ACE/NME/capping residues found: {caps_found}")
        print("  -> This file HAS capping applied.")
    else:
        print("\n  No ACE/NME/capping residues found.")
        print("  -> This file does NOT have capping applied (pre-capping stage).")

    all_nums = set()
    for nums in chain_data.values():
        all_nums.update(nums)

    print(f"\n  {'Residue':>8}  {'Raw':>8}  {'This file':>10}")
    for target in TARGET_RESIDUES:
        raw_status = "present" if target in raw_all else "ABSENT"
        this_status = "present" if target in all_nums else "ABSENT"
        print(f"  {target:>8}  {raw_status:>8}  {this_status:>10}")

    fixer = PDBFixer(filename=path)
    fixer.findMissingResidues()
    fixer.findMissingAtoms()
    if not fixer.missingAtoms:
        print("\n  PDBFixer re-check: zero missing heavy atoms.")
    else:
        print(f"\n  PDBFixer re-check: MISSING ATOMS STILL PRESENT: {fixer.missingAtoms}")

    print("\n")

print("=" * 78)
print("HOW TO READ THIS")
print("=" * 78)
print("The Methods sentence describes the state right after PDBFixer, before")
print("any later capping step. So: prefer whichever candidate (a) has NO")
print("capping residues, (b) shows zero missing heavy atoms, and (c) is the")
print("earliest-stage file in your pipeline order. If two candidates both")
print("qualify, they should report the SAME residue count - use that number.")
print("If a candidate shows capping residues, its residue count is almost")
print("certainly NOT what the manuscript sentence means, since caps are")
print("added after the step this sentence describes.")
