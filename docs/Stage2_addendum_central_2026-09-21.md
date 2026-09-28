# Central addendum (2026-09-21): consolidating Addenda 2, 3 and 4

Companion to `Stage2_rerun_criteria_and_addendum_2026-09-21.md` (the plan file, not edited) and to Addenda 2, 3 and 4, which this document summarises but does not replace. Where this summary and an individual addendum differ, the individual addendum is authoritative.

**Timeline, for context, not part of the scientific record above:** the six rerun production replicates were run over approximately four to five days, finishing 2026-09-20. The decisions below were made during analysis of those results on 2026-09-21, in the order given in Addenda 2, 3 and 4.

**Note on numbering:** the pre-analysis build and design decisions (peptide-bond, caps, equilibration length, cutoffs, software version) function as the first item in this addendum series but were recorded inline in `Stage2_rerun_criteria_and_addendum_2026-09-21.md` §2 rather than extracted as a separately numbered file. Numbering below begins at Addendum 2 for that reason.

## 1. Periodic-boundary handling (Addendum 2)

Raw analysis of all six stitched trajectories showed frame-to-frame backbone RMSD jumps of up to 4.3 nm and radius-of-gyration excursions to 3.2 nm, exceeding the plan file's 0.5 nm threshold — not plausible as real motion for a folded 206-residue domain. `gmx trjconv -pbc mol` was tried first and did not resolve this in five of six replicates. `gmx trjconv -pbc nojump` (System group) was applied identically to all six replicates instead, after which the plan file's analysis pipeline was run unchanged. All six then showed a maximum jump of 0.103 nm or less. Raw and `-pbc mol` outputs are not used in any reported result.

## 2. Table 5 cross-check: metric per hypothesis (Addendum 3)

Table 5 reports its fourteen original hypothesis residues by two different metrics: RMSF for five residues (43, 44, 45, 165, 174) and DSSP secondary-structure fraction for the other nine (sheet fraction at 28, 46, 47, 88, 89; helix fraction at 137, 162, 163, 164). The criterion-3 comparison is therefore run per residue using the metric and direction Table 5 itself specifies for that residue, not a single metric applied uniformly across all fourteen. Residues 171 to 173 remain a consistency observation carried over from Stage 1, not a criterion-3 test.

## 3. Residue 133 contact definition (Addendum 4)

Fixed before any contact result was computed or seen. A residue is in contact with residue 133 in a given frame if any heavy atom of that residue lies within 0.45 nm of any heavy atom of residue 133, computed per replicate across all 2001 frames of each 20 ns trajectory. This is exploratory throughout (criterion 4), not a criterion-3 test, and the definition was not altered after any result was seen.

## 4. Final criterion-3 result

Applying the Table 5 metric per residue (section 2) to the `nojump`-corrected trajectories (section 1): **2 of 14 original hypothesis residues confirmed, 12 not confirmed, 0 contradicted.** Confirmed: residue 28 (DSSP sheet gain, narrow margin) and residue 165 (RMSF, mutant lower). Both confirmed residues lie adjacent to a chain discontinuity in the model rather than in an uninterrupted region of the fold.

## 5. Where this belongs in the manuscript

This document is not quoted in the manuscript body. It is referenced from:
- **Methods**, where the `nojump` deviation from the plan file's original pipeline is reported as a stated deviation (section 1 above), and the per-residue metric assignment is reported as how the criterion-3 test was applied (section 2).
- **Results / Discussion**, where the final count (section 4) and the residue 133 contact definition (section 3) are reported as findings.
- **Data availability / Supplementary materials**, where this file and Addenda 2 to 4 are deposited alongside the plan file, so a reviewer can trace each deviation to the point in analysis where it was made.
