"""
Atom-inventory check for the PDBFixer-completed 1JNX structure.

Confirms, directly against the files (not by inference):
  1. Total residue count per chain in the fixed structure.
  2. Zero missing heavy atoms, by re-running PDBFixer's own detector against
     the fixed file (should return empty if the fix genuinely completed).
  3. Whether residues 1646-1648, 1694, 1817-1819 are PRESENT or ABSENT in
     the fixed file, cross-checked against their status in the raw 1JNX file.
     This answers the "filled vs left as a break" question directly.
  4. First and last modelled residue number, and any remaining internal
     numbering gaps.

Edit the two paths below, then run the whole script as one cell/block.
"""

RAW_PDB_PATH = "/content/drive/MyDrive/PATH/TO/1JNX.pdb"        # <-- EDIT THIS: original downloaded 1JNX
FIXED_PDB_PATH = "/content/drive/MyDrive/PATH/TO/1JNX_fixed.pdb"  # <-- EDIT THIS: PDBFixer output actually used for production

import os
import sys

# ---- dependencies ----
try:
    from pdbfixer import PDBFixer
except ImportError:
    raise SystemExit("pdbfixer not installed. Run: !pip install pdbfixer -q   (restart runtime if it was just installed, then re-run this cell)")

try:
    from Bio.PDB import PDBParser
except ImportError:
    raise SystemExit("biopython not installed. Run: !pip install biopython -q")

for label, path in [("RAW_PDB_PATH", RAW_PDB_PATH), ("FIXED_PDB_PATH", FIXED_PDB_PATH)]:
    if not os.path.exists(path):
        raise SystemExit(f"{label} does not exist: {path}\nEdit the path at the top of the script and re-run.")

TARGET_RESIDUES = [1646, 1647, 1648, 1649, 1694, 1817, 1818, 1819]

parser = PDBParser(QUIET=True)


def residue_numbers_by_chain(pdb_path, label):
    structure = parser.get_structure(label, pdb_path)
    model = structure[0]
    out = {}
    for chain in model.get_chains():
        residues = [r for r in chain if r.id[0] == " "]  # standard residues only; excludes waters/hetero/ions
        if not residues:
            continue
        out[chain.id] = sorted(r.id[1] for r in residues)
    return out


print("=" * 70)
print("STEP 1: Residue inventory, fixed structure")
print("=" * 70)

fixed_by_chain = residue_numbers_by_chain(FIXED_PDB_PATH, "fixed")
for chain_id, nums in fixed_by_chain.items():
    print(f"\nChain {chain_id}: {len(nums)} residues")
    print(f"  Residue number range: {nums[0]} to {nums[-1]}")
    gaps = [(a, b) for a, b in zip(nums, nums[1:]) if b != a + 1]
    if gaps:
        print(f"  Internal numbering gaps remaining: {gaps}")
    else:
        print("  No internal numbering gaps: chain is numerically contiguous.")

print("\n" + "=" * 70)
print("STEP 2: Raw vs fixed, target residues (1646-1648, 1694, 1817-1819)")
print("=" * 70)

raw_by_chain = residue_numbers_by_chain(RAW_PDB_PATH, "raw")
# Flatten to a single set per file for the presence check (report which chain if relevant)
raw_all = set()
for nums in raw_by_chain.values():
    raw_all.update(nums)
fixed_all = set()
for nums in fixed_by_chain.values():
    fixed_all.update(nums)

print(f"\n{'Residue':>8}  {'Raw 1JNX':>10}  {'Fixed file':>10}")
for target in TARGET_RESIDUES:
    raw_status = "present" if target in raw_all else "ABSENT"
    fixed_status = "present" if target in fixed_all else "ABSENT"
    print(f"{target:>8}  {raw_status:>10}  {fixed_status:>10}")

print("\nRead this table as: any target residue ABSENT in both columns was left as a break.")
print("Any residue ABSENT in raw but present in fixed was filled by PDBFixer.")

print("\n" + "=" * 70)
print("STEP 3: Re-running PDBFixer's own missing-atom detector on the FIXED file")
print("=" * 70)

fixer = PDBFixer(filename=FIXED_PDB_PATH)
fixer.findMissingResidues()
fixer.findMissingAtoms()

print(f"\nMissing residues detected (depends on SEQRES being present; informational only): {fixer.missingResidues}")
print(f"Missing heavy atoms detected: {fixer.missingAtoms}")
print(f"Missing terminal atoms detected: {fixer.missingTerminals}")

if not fixer.missingAtoms:
    print("\nCONFIRMED: zero missing heavy atoms in the fixed structure, per PDBFixer's own check.")
else:
    print("\nWARNING: PDBFixer still detects missing heavy atoms in this file. The fix did not fully complete, investigate before using this residue count.")

print("\n" + "=" * 70)
print("SUMMARY")
print("=" * 70)
for chain_id, nums in fixed_by_chain.items():
    print(f"Chain {chain_id} residue count for the manuscript's atom-inventory sentence: {len(nums)}")
