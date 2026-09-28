#!/usr/bin/env python3
"""
find_jump_times.py

Pure-Python pass over already-computed rmsd.xvg files (raw and, where it
exists, whole-molecule-corrected) to find WHEN each replicate's big
frame-to-frame jump happens, and to dig into the mut rep 2 anomaly
(identical raw/corrected max jump, but very different last-frame value).

No GROMACS call. Just reads .xvg text.

EDIT THE CONFIG BLOCK BELOW to match your actual file paths on Drive.
Corrected-file naming is a guess (`_rmsd_corrected.xvg`) since that
convention wasn't fixed anywhere in the addendum -- change it to
whatever you actually named the post-`trjconv -pbc mol` reruns. If a
replicate has no corrected file yet, leave that entry as None and the
script will just report the raw side for it.
"""

import sys

# ---------------------------------------------------------------- CONFIG --
# One entry per replicate. "raw" = the O + N + "_rmsd.xvg" file written by
# stage2_analysis_an.py. "corrected" = the rmsd.xvg from a whole-molecule
# corrected rerun (trjconv -pbc mol -> gmx rms), or None if not made yet.
BASE = "/content/drive/MyDrive/brca1_md_v2/analysis2/"

REPLICATES = {
    "wt_prod_rep1":  {"raw": BASE + "wt_prod_rep1_rmsd.xvg",  "corrected": None},
    "wt_prod_rep2":  {"raw": BASE + "wt_prod_rep2_rmsd.xvg",  "corrected": BASE + "wt_prod_rep2_rmsd_corrected.xvg"},
    "wt_prod_rep3":  {"raw": BASE + "wt_prod_rep3_rmsd.xvg",  "corrected": BASE + "wt_prod_rep3_rmsd_corrected.xvg"},
    "mut_prod_rep1": {"raw": BASE + "mut_prod_rep1_rmsd.xvg", "corrected": BASE + "mut_prod_rep1_rmsd_corrected.xvg"},
    "mut_prod_rep2": {"raw": BASE + "mut_prod_rep2_rmsd.xvg", "corrected": BASE + "mut_prod_rep2_rmsd_corrected.xvg"},
    "mut_prod_rep3": {"raw": BASE + "mut_prod_rep3_rmsd.xvg", "corrected": BASE + "mut_prod_rep3_rmsd_corrected.xvg"},
}

JUMP_THRESHOLD = 0.5   # nm, same threshold used in stage2_analysis_an.py
TOP_N = 3              # how many biggest jumps to report per file
TAIL_N = 15            # how many trailing frames to dump for the mut-rep2 dig-in
NO_CHANGE_TOL = 0.05   # relative reduction below this = "correction did nothing, real motion"
OUT_CSV = BASE + "jump_decision_summary.csv"   # set to None to skip writing a CSV
# ---------------------------------------------------------------------- --


def load_xvg(path):
    """Return list of (time_ns, rmsd_nm) from a GROMACS .xvg, skipping headers."""
    pts = []
    with open(path) as fh:
        for line in fh:
            if not line.strip() or line[0] in "#@":
                continue
            parts = line.split()
            pts.append((float(parts[0]), float(parts[1])))
    return pts


def jumps(pts, threshold=JUMP_THRESHOLD):
    """
    Return every frame-to-frame jump |delta rmsd| as
    (time_from_ns, time_to_ns, delta_nm), sorted largest first.
    """
    out = []
    for (t0, y0), (t1, y1) in zip(pts, pts[1:]):
        d = abs(y1 - y0)
        out.append((t0, t1, d))
    out.sort(key=lambda r: -r[2])
    return out


def classify(raw_max, corr_max, threshold=JUMP_THRESHOLD, no_change_tol=NO_CHANGE_TOL, anomaly=False):
    """
    Decide raw vs corrected vs manual-review for one replicate, using an
    explicit, citable rule:

      - no corrected file at all             -> USE RAW
      - corrected max jump <= threshold      -> USE CORRECTED (artefact resolved)
      - correction reduced the jump by less
        than NO_CHANGE_TOL (relative)        -> USE RAW (correction did ~nothing -> real motion)
      - anything else (partial reduction,
        still above threshold)               -> MANUAL REVIEW
      - anomaly flag (same max jump, but
        last-frame value moved a lot) always
        overrides to MANUAL REVIEW, since a
        single max-jump comparison misses it
    """
    if anomaly:
        return "MANUAL REVIEW", "anomaly: max jump unchanged but a later frame diverged"
    if corr_max is None:
        tail = " (still exceeds threshold, uncorrected)" if raw_max > threshold else ""
        return "USE RAW", "no corrected file made" + tail
    if corr_max <= threshold:
        return "USE CORRECTED", f"artefact resolved: jump now {corr_max:.3f} <= {threshold} nm"
    rel_reduction = (raw_max - corr_max) / raw_max if raw_max else 0
    if rel_reduction < no_change_tol:
        return "USE RAW", f"correction changed jump by only {rel_reduction*100:.1f}% -> real motion, not PBC"
    return "MANUAL REVIEW", f"partial reduction ({rel_reduction*100:.1f}%), still {corr_max:.3f} > {threshold} nm"


def report_file(label, path, threshold=JUMP_THRESHOLD, top_n=TOP_N):
    try:
        pts = load_xvg(path)
    except FileNotFoundError:
        print(f"  {label}: FILE NOT FOUND ({path})")
        return None
    j = jumps(pts)
    print(f"  {label}: {len(pts)} frames, last t={pts[-1][0]:.1f} ns rmsd={pts[-1][1]:.3f} nm, "
          f"max rmsd={max(y for _, y in pts):.3f} nm")
    for t0, t1, d in j[:top_n]:
        flag = "  <-- >threshold" if d > threshold else ""
        print(f"    jump {d:.3f} nm between t={t0:.2f} -> t={t1:.2f} ns{flag}")
    return pts, j


def main():
    print(f"Jump threshold: {JUMP_THRESHOLD} nm, top {TOP_N} jumps per file, "
          f"no-change tolerance: {NO_CHANGE_TOL*100:.0f}%\n")

    summary = []  # (rep, raw_max, corr_max_or_None, decision, reason)

    for rep, files in REPLICATES.items():
        print(f"=== {rep} ===")
        raw_res = report_file("raw", files["raw"])
        corr_res = None
        if files["corrected"]:
            corr_res = report_file("corrected", files["corrected"])
        else:
            print("  corrected: (none configured / not yet made)")
        print()

        anomaly = False

        # Extra dig-in for any replicate whose raw and corrected max jump
        # are (near-)identical but whose last-frame value differs a lot --
        # this is exactly the mut rep 2 symptom, checked generically so it
        # will also catch it if another replicate does the same thing.
        if raw_res and corr_res:
            raw_pts, raw_j = raw_res
            corr_pts, corr_j = corr_res
            same_max = abs(raw_j[0][2] - corr_j[0][2]) < 0.01
            last_diff = abs(raw_pts[-1][1] - corr_pts[-1][1])
            if same_max and last_diff > 0.5:
                anomaly = True
                print(f"  ** anomaly flag: max jump unchanged by correction "
                      f"({raw_j[0][2]:.3f} nm at t={raw_j[0][0]:.2f} ns) but "
                      f"last-frame RMSD moved {raw_pts[-1][1]:.3f} -> {corr_pts[-1][1]:.3f} nm.")
                print(f"  ** this means the correction fixed a DIFFERENT, later jump than the "
                      f"single largest one. Trailing {TAIL_N} frames, raw vs corrected:")
                print(f"  {'t (ns)':>8}  {'raw rmsd':>9}  {'corr rmsd':>9}  {'raw d':>7}  {'corr d':>7}")
                raw_tail = raw_pts[-TAIL_N:]
                corr_tail = corr_pts[-TAIL_N:]
                for i in range(1, len(raw_tail)):
                    t = raw_tail[i][0]
                    ry, cy = raw_tail[i][1], corr_tail[i][1]
                    rd = raw_tail[i][1] - raw_tail[i - 1][1]
                    cd = corr_tail[i][1] - corr_tail[i - 1][1]
                    mark = "  <--" if abs(rd) > 0.3 or abs(cd) > 0.3 else ""
                    print(f"  {t:8.2f}  {ry:9.3f}  {cy:9.3f}  {rd:7.3f}  {cd:7.3f}{mark}")
                print()

        if raw_res:
            raw_max = raw_res[1][0][2]
            corr_max = corr_res[1][0][2] if corr_res else None
            decision, reason = classify(raw_max, corr_max, anomaly=anomaly)
            summary.append((rep, raw_max, corr_max, decision, reason))

    # ---- final decision table, citable in Methods -------------------------
    print("=" * 78)
    print("DECISION SUMMARY (rule: corrected if jump <= threshold after correction;")
    print("raw if correction changes jump by <%.0f%%; else flag for manual review;" % (NO_CHANGE_TOL * 100))
    print("anomaly cases always -> manual review)")
    print("=" * 78)
    header = f"{'replicate':<16}{'raw max':>9}{'corr max':>10}  {'decision':<16}reason"
    print(header)
    print("-" * len(header))
    for rep, raw_max, corr_max, decision, reason in summary:
        cm = f"{corr_max:.3f}" if corr_max is not None else "  n/a"
        print(f"{rep:<16}{raw_max:>9.3f}{cm:>10}  {decision:<16}{reason}")

    if OUT_CSV:
        try:
            with open(OUT_CSV, "w") as fh:
                fh.write("replicate,raw_max_jump_nm,corrected_max_jump_nm,decision,reason\n")
                for rep, raw_max, corr_max, decision, reason in summary:
                    cm = f"{corr_max:.3f}" if corr_max is not None else ""
                    fh.write(f"{rep},{raw_max:.3f},{cm},{decision},\"{reason}\"\n")
            print(f"\nWrote {OUT_CSV}")
        except OSError as e:
            print(f"\nCould not write {OUT_CSV}: {e}")

    print("\nMANUAL REVIEW entries are judgment calls, not automated -- decide those by eye")
    print("(e.g. inspect the trajectory in VMD/PyMOL around the flagged jump time) before")
    print("locking which trajectory (raw or corrected) is 'official' for that replicate.")
    print("\nIf a replicate shows 'FILE NOT FOUND', fix the path in the CONFIG block and rerun --")
    print("no GROMACS needed, this only parses .xvg text.")


if __name__ == "__main__":
    main()
