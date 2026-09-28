import os,sys,statistics as st
from collections import Counter
R=(sys.argv[1] if len(sys.argv)>1 else "/content/drive/MyDrive/brca1_md_v2/analysis2/")
WT=["wt_prod_rep%d"%i for i in(1,2,3)];MU=["mut_prod_rep%d"%i for i in(1,2,3)]
SH,HX="EB","HGI"
# residue, metric chars, Table 5 direction (+1 mutant higher), label, Table 5 Stage 1 ranges (WT / mutant), for reference only
T=[(28,SH,1,"sheet","0.44-0.90 / 0.91-0.99"),(46,SH,-1,"sheet","0.29-0.49 / 0.02-0.06"),(47,SH,-1,"sheet","0.63-0.98 / 0.02-0.33"),
(88,SH,-1,"sheet","0.40-0.98 / 0.00-0.33 (88-89)"),(89,SH,-1,"sheet","0.40-0.98 / 0.00-0.33 (88-89)"),
(137,HX,-1,"helix","0.54-0.60 / 0.15-0.48"),(162,HX,1,"helix","<0.01 / 0.26-0.61 (162-164)"),(163,HX,1,"helix","<0.01 / 0.26-0.61 (162-164)"),(164,HX,1,"helix","<0.01 / 0.26-0.61 (162-164)")]
C=[(171,HX,1,"helix","WT 0.87-0.88 in 1 of 3, 0.04-0.21 in 2 of 3; mut 0.76-0.88 (171-173)"),(172,HX,1,"helix",""),(173,HX,1,"helix","")]
def stop(m):raise SystemExit("STOP: "+m)
S={};seen=Counter()
for n in WT+MU:
    p=R+n+"_nj_ss206.dat"
    if not os.path.exists(p):stop("missing "+p)
    L=[l.rstrip("\n") for l in open(p) if l.strip()]
    if len(L)!=2001 or {len(s) for s in L}!={206}:stop(n+" ss206 shape wrong")
    S[n]=L
    for s in L:seen.update(s)
print("DSSP characters seen:",dict(seen))
def fr(n,r,ch):
    L=S[n];return sum(s[r-1] in ch for s in L)/len(L)
def rng(v):return "%.2f-%.2f"%(min(v),max(v))
def verdict(d1,w,m):
    if min(m)>max(w):d2=1
    elif max(m)<min(w):d2=-1
    else:return "NOT CONFIRMED (overlap)"
    return "CONFIRMED" if d2==d1 else "CONTRADICTED (opposite direction)"
rows=[]
print("\nDSSP fraction of frames over 20 ns. Dir: + = mutant higher in Table 5. Ranges over 3 reps.")
print("%-5s %-6s %-4s %-11s %-11s %-34s %s"%("res","type","dir","WT range","Mut range","verdict","Stage 1 (Table 5) WT / mut"))
for lab,L in(("criterion 2 (Table 5 metric)",T),("consistency observation, not a criterion-3 test",C)):
    print("--",lab)
    for r,ch,d,ty,ref in L:
        w=[fr(n,r,ch) for n in WT];m=[fr(n,r,ch) for n in MU]
        v=verdict(d,w,m)
        print("%-5d %-6s %-4s %-11s %-11s %-34s %s"%(r,ty,"+" if d>0 else "-",rng(w),rng(m),v,ref))
        rows.append([r,ty,d,rng(w),rng(m),v])
import csv
with open(R+"comparison_dssp_nj.csv","w",newline="") as f:
    w=csv.writer(f);w.writerow(["res","type","table5_dir","WT_range","Mut_range","verdict"]);w.writerows(rows)
print("\nEND OK")
