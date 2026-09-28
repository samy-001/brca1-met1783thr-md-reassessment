import os,sys,statistics as st
from collections import Counter
R=(sys.argv[1] if len(sys.argv)>1 else "/content/drive/MyDrive/brca1_md_v2/analysis2/")
OLD=(sys.argv[2] if len(sys.argv)>2 else "/content/drive/MyDrive/brca1_md/analysis/")
WT=["wt_prod_rep%d"%i for i in(1,2,3)];MU=["mut_prod_rep%d"%i for i in(1,2,3)]
OW=[OLD+"rmsf_rep%d.xvg"%i for i in(1,2,3)];OM=[OLD+"rmsf_mut_rep%d.xvg"%i for i in(1,2,3)]
RES=[28,43,44,45,46,47,88,89,137,162,163,164,165,174];CONS=[171,172,173]
def stop(m):raise SystemExit("STOP: "+m)
def xv(p):
    if not os.path.exists(p):stop("missing "+p)
    return [[float(x) for x in l.split()] for l in open(p) if l.strip() and l[0] not in "#@"]
bad=[];old={}
for p in OW+OM:
    d=xv(p);k=[int(a[0]) for a in d]
    print("old",os.path.basename(p),"residues",len(k),"first",k[0],"last",k[-1])
    if k!=list(range(1,207)):bad.append(os.path.basename(p))
    old[p]={int(a[0]):a[1] for a in d}
if bad:stop("old RMSF numbering is not 1..206 in: "+", ".join(bad)+". Paste the lines above; do not compare until this is understood.")
new={}
for n in WT+MU:
    d=xv(R+n+"_nj_rmsf.xvg");k=[int(a[0]) for a in d]
    if k!=list(range(1,207)):stop(n+" new RMSF not 1..206")
    new[n]={int(a[0]):a[1] for a in d}
def rng(v):return "%.3f-%.3f"%(min(v),max(v))
def ov(a,b):return not(min(a)>max(b) or max(a)<min(b))
def verdict(d1,w,m):
    if min(m)>max(w):d2=1
    elif max(m)<min(w):d2=-1
    else:return "NOT CONFIRMED (overlap)"
    return "CONFIRMED" if d2==d1 else "CONTRADICTED (opposite direction)"
rows=[]
print("\nRMSF (nm). S1 dir: + = mutant more flexible in old runs (mean mut - mean WT). Ranges over 3 reps.")
print("%-6s %-4s %-8s %-13s %-13s %s"%("res","S1","S1 ovlp","WT range","Mut range","criterion 3 verdict"))
for lab,L in(("criterion 2",RES),("consistency obs (171-173), not a criterion-3 test",CONS)):
    print("--",lab)
    for r in L:
        ow=[old[p][r] for p in OW];om=[old[p][r] for p in OM]
        d1=1 if st.mean(om)>st.mean(ow) else -1
        w=[new[n][r] for n in WT];m=[new[n][r] for n in MU]
        v=verdict(d1,w,m);o="yes" if ov(ow,om) else "no"
        print("%-6d %-4s %-8s %-13s %-13s %s"%(r,"+" if d1>0 else "-",o,rng(w),rng(m),v))
        rows.append([r,lab[:11],d1,o,rng(w),rng(m),v])
print("\nResidue 133 RMSF (nm): old WT",[round(old[p][133],3) for p in OW],"old mut",[round(old[p][133],3) for p in OM])
print("                       new WT",[round(new[n][133],3) for n in WT],"new mut",[round(new[n][133],3) for n in MU])
print("\nRMSD / Rg (new, nm): rmsd mean sd last max | Rg mean sd min max")
for n in WT+MU:
    y=[a[1] for a in xv(R+n+"_nj_rmsd.xvg")];g=[a[1] for a in xv(R+n+"_nj_gyrate.xvg")]
    print("%-14s %.3f %.3f %.3f %.3f | %.3f %.3f %.3f %.3f"%(n,st.mean(y),st.pstdev(y),y[-1],max(y),st.mean(g),st.pstdev(g),min(g),max(g)))
print("\nDSSP: residue 133 class fractions, and mean helix (H,G,I) fraction over residues 1-206")
for n in WT+MU:
    S=[l.rstrip("\n") for l in open(R+n+"_nj_ss206.dat") if l.strip()]
    if len(S)!=2001 or {len(s) for s in S}!={206}:stop(n+" ss206 shape wrong")
    c=Counter(s[132] for s in S);hx=st.mean(sum(ch in"HGI" for ch in s)/206 for s in S)
    print("%-14s res133 %s | helix %.3f"%(n," ".join("%s:%.2f"%(k,v/len(S)) for k,v in c.most_common()),hx))
import csv
with open(R+"comparison_rmsf_nj.csv","w",newline="") as f:
    w=csv.writer(f);w.writerow(["res","group","S1dir","S1overlap","WT_range","Mut_range","verdict"]);w.writerows(rows)
print("\nEND OK")
