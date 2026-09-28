# Stage 2 rerun: confirmation criteria and addendum

Written 21 September 2026. This file is dated by the day it was written and is not backdated.

## 0. Provenance and status (read first)

- An earlier file named `Stage2_rerun_confirmation_criteria_2026-09-19.md` is referred to in the project notes, but no copy could be found. This is not that file.
- Criteria 1 to 6 in section 1 are transcribed unchanged from section 8 of `Paper2_MASTER_handoff_brief_2026-09-20.md` ("Six criteria (approved)"). That brief is committed alongside this file.
- This file was written after production finished for all six replicates and before any analysis of them. Sam confirmed on 2026-09-21 that no RMSD, RMSF, Rg, DSSP or interface metric has been computed on any Stage 2 replicate.
- Stage 1 (old-run) results were seen before these criteria were recorded, so this is a confirmation design, not a blind one.
- Because the 19 September file cannot be produced, the criteria should be described in the manuscript as specified before the Stage 2 trajectories were analysed, not as pre-registered.

## 1. The six criteria (unchanged from Master Brief section 8)

1. New runs primary, old runs sensitivity comparison only.
2. Hypotheses = original Table 5 residues 28, 43 to 45, 46, 47, 88 to 89, 137, 162 to 164, 165, 174, with 171 to 173 as a consistency observation.
3. Confirmed only if direction consistent and ranges non-overlapping, all 3 mutant vs all 3 WT replicates.
4. Failures reported as not confirmed, new signals labelled exploratory.
5. Residue 133, RMSD and Rg re-reported whatever they show.
6. Report everything, including contradictions.

## 2. Addendum: build and design decisions (dated 2026-09-21)

- **44/45 gap.** Decision confirmed: keep the peptide bond at the 44/45 gap in both WT and mutant. Reasoning: the gap (3.5 A) is too narrow for two caps (about 5 A needed) without clashing; keeping the bond is a smaller artefact, identical in both systems, and keeps residues 43 to 47 comparable to the old data. Caveat: residues 43 to 47 and residue 28 remain a "read with caution" area.
- **Caps.** ACE and NME at the 166/167 gap. Chain A = residues 1 to 166 + NME (207); chain B = ACE (208) + residues 167 to 206. All analysis uses the 206 real residues only; the caps are excluded.
- **NPT.** Equilibration extended from 100 ps (old runs) to 200 ps, equally for both systems because density had not flattened by 100 ps; production starts from the 200 ps frame. 
- **Cutoffs.** Plain 1.0 nm cutoffs with potential-shift and no dispersion correction, as in the old runs, to keep the comparison with them clean. Force-switch with 1.2 nm cutoffs (common CHARMM36 guidance) was noted and not adopted.
- **Software.** GROMACS 2025.4-conda_forge with CUDA on Tesla T4 GPUs.
- **Production.** Three WT and three mutant replicates, 10,000,000 steps (20 ns, 2 fs), compressed frames every 10 ps (2001 frames), new velocities per replicate, run in chunks with checkpoint restarts. WT rep 2 and mutant rep 3 each have three parts; the other four have two. Boundary frames were removed when parts were stitched.
- **PME tuning and drift (observation, not a criterion).** The PME grid chosen at startup differed between parts (80 / 1.000, 60 / 1.283, 72 / 1.067), and the conserved-energy drift tracked that grid (about 3.0e-05, 6.6e-05 and 7.6e-05 kJ/mol/ps per atom), in both WT and mutant and within replicates. No exclusion rule is attached to it. It is reported per part in the Methods or supplement.

## 3. Data fixed before analysis (stitched trajectories, 2001 frames, 0 to 20,000 ps)

| Replicate | File (under `brca1_md_v2/`) | Size (bytes) | sha256 |
|---|---|---|---|
| WT rep 1 | `wt_v2/stitched/wt_prod_rep1_0-20000ps.xtc` | 585,955,672 | `6406dd61e20531b1601987c942a5855c9d11519af0b4a29bd1c510248a439d1d` |
| WT rep 2 | `wt_v2/stitched/wt_prod_rep2_0-20000ps.xtc` | 585,953,948 | `ac881ce1167c2928b2bdd6dd4601770c27ad9b3eaa28d57f07f0bcb9ce63908c` |
| WT rep 3 | `wt_v2/stitched/wt_prod_rep3_0-20000ps.xtc` | 585,949,692 | `fac78dcd5cd64142655cc8f0fa77157ffddbd81f5bb6d97c9c0ebaf2165037f9` |
| Mutant rep 1 | `mut_v2/stitched/mut_prod_rep1_0-20000ps.xtc` | 585,995,072 | `2f8eca63810e6cfd4922a2aa2ed39c9cb93914ac0b8502d51d178407a7ea7d6c` |
| Mutant rep 2 | `mut_v2/stitched/mut_prod_rep2_0-20000ps.xtc` | 585,999,556 | `bfdda53b08b9bc261befd9962999b9edaaffe8c1af2fb7877ec9d13af725b45e` |
| Mutant rep 3 | `mut_v2/stitched/mut_prod_rep3_0-20000ps.xtc` | 585,990,996 | `abe8c18743511e0d8c3fd423bd4b9688d3370438b6cf4a58bc45ea75fd220879` |

## 4. Analysis statement (details fixed 2026-09-21, before any analysis)

- **Stage 1 commands, read from the `.xvg` headers in `brca1_md/analysis` (GROMACS 2026.3):** `gmx rms -s prod.tpr -f prod_full.xtc -tu ns` (backbone, least-squares fitted to backbone); `gmx rmsf -s prod.tpr -f prod_full.xtc -res`; `gmx gyrate -s prod.tpr -f prod_full.xtc`. DSSP output (`ss_*.dat`) had 206 characters per frame for WT and 207 for mutant; the DSSP command is not recorded, and the extra mutant character is inferred (from the old 166/167 difference) to be a chain-break marker.
- **Not recorded in Stage 1:** the atom group used for RMSF and Rg, the DSSP command, and any periodic-boundary processing.
- **Stage 2 pipeline (identical for all six replicates):** the full stitched trajectory (2001 frames, 0 to 20 ns) against the production `.tpr` starting structure. RMSD: Backbone restricted to residues 1 to 206, fitted to itself. RMSF (`-res`) and Rg: Protein group restricted to residues 1 to 206 (caps excluded). DSSP: `gmx dssp` on residues 1 to 206, chain-break markers removed, 206 characters per frame. Group sizes are checked against the built systems (3290 protein atoms WT, 3287 mutant, 618 backbone atoms). If a frame-to-frame RMSD jump above 0.5 nm indicates a periodic-boundary artefact, whole-molecule correction is applied identically to all six replicates and reported as a deviation.
- **Direction of the Stage 1 effect** for each criterion-2 residue is taken from the old RMSF files (mutant mean minus WT mean over three replicates) and cross-checked against Table 5 when its values are available (not yet seen when this was written).
- Results are judged strictly against section 1. Any deviation is reported as a deviation.

## 5. Additional analysis items registered (2026-09-21, before any analysis)

These are additions to the six criteria in section 1, not changes to them. They were added before any Stage 2 trajectory was analysed.

- **Mutation site.** Residue 133 (local numbering; Met1783 in UniProt numbering) is reported explicitly, for RMSF, DSSP and its local contacts, whatever it shows (this extends criterion 5).
- **Distal versus mutation site.** The Stage 1 hypothesis residues (criterion 2) lie away from residue 133 in the sequence. The 3D distance from residue 133 to each criterion-2 residue (minimum heavy-atom distance in the WT starting structure) is computed once and reported in a table. Results are described as "distal" or "near" only by reference to those reported distances. No distance threshold is chosen after seeing the results.
- **No target outcome.** The analysis is not aimed at confirming or changing the Stage 1 pattern. Confirmation, non-confirmation and any new signals (including at residue 133 itself, labelled exploratory) are reported under criteria 3, 4 and 6.
