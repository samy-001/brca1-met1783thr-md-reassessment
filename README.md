# BRCA1 p.Met1783Thr: Structural and Molecular Dynamics Reassessment

Analysis scripts, structural inputs, and derived results supporting the manuscript reassessing the BRCA1 variant of uncertain significance (VUS) p.Met1783Thr, one of 13 missense VUS originally reported in Onyia et al. (2025), using a confidence-gated FoldX and molecular dynamics pipeline.

Raw molecular dynamics trajectory files are too large for GitHub and are archived separately on Zenodo (see **Data availability** below).

## Repository structure

```
scripts/     Python analysis scripts (RMSD/RMSF/Rg/DSSP comparison, contact analysis,
             jump detection, figure generation)
data/
  mdp/                 GROMACS parameter files (energy minimisation, ion addition,
                        NVT, NPT, production)
  topology/wt_v2/      Wild-type topology files (rerun/primary system)
  topology/mut_v2/     Met1783Thr mutant topology files (rerun/primary system)
  individual_list_met133thr.txt   FoldX mutation specification (local numbering)
  *_capped.pdb         Final capped structures used to build the MD systems
results/
  old_run/     Analysis outputs (RMSD, RMSF, Rg, DSSP) and comparison plots from
               the initial exploratory run. Sensitivity comparison only; not the
               primary reported result.
  new_run/     Analysis outputs from the matched confirmatory rerun (primary
               result), including the *_nj_* files (periodic-boundary corrected,
               see Methods and docs/), manuscript Figures 2, 3 and 5, and
               equilibration_qc/ (NVT/NPT diagnostics supporting the extended
               equilibration decision documented in docs/).
docs/        Confirmation criteria, build/design decisions, and addenda
             documenting deviations from the original analysis plan, with
             traceability for reviewers.
```

## Data availability

- **Full molecular dynamics trajectories** (WT and Met1783Thr mutant, initial run and matched rerun): archived on Zenodo, DOI: [10.5281/zenodo.22985272](https://doi.org/10.5281/zenodo.22985272). Files are currently access-restricted pending peer review outcome; the record itself is public and citable.
- **This repository**: scripts, structural inputs, and derived analysis outputs (not the raw trajectories).

## Structural inputs

Starting structures were built from PDB entry 1JNX (the apo/unbound BRCT tandem repeat, selected over the peptide-bound structure 1T15 to provide an unbound baseline; see manuscript Methods) and AlphaFold DB models, using FoldX for mutant generation and PDBFixer/GROMACS `pdb2gmx` for system preparation. The CHARMM36m force field (March 2019 release) was used but is not redistributed here; it is freely available from the [MacKerell lab CHARMM force field page](http://mackerell.umaryland.edu/charmm_ff.shtml).

## Methodology notes and deviations

See `docs/` for the confirmation criteria fixed before the rerun was analysed, build/design decisions (chain-gap handling, capping, equilibration length, cutoffs), and documented deviations from the original analysis plan (periodic-boundary correction, per-residue metric assignment, residue 133 contact definition), each traceable to the point in analysis where the decision was made.

## Citation

If you use this code or data, please cite the associated manuscript (citation to be added upon publication) and the Zenodo dataset DOI above. See `CITATION.cff` for structured citation metadata.

## License

This repository, including code and data files, is released under CC0 1.0 Universal (Public Domain Dedication) — see `LICENSE`. This matches the license used for the accompanying Zenodo dataset.

## Contact

Questions about this repository or the underlying analysis can be directed to the corresponding author: adekunlesamuelsa@gmail.com
