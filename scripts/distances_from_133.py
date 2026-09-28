import os,sys,shutil,csv,subprocess as sp
R="/content/drive/MyDrive/brca1_md_v2/";W="/content/an/dist/";os.makedirs(W,exist_ok=True)
CRIT=[28,43,44,45,46,47,88,89,137,162,163,164,165,174];CONS=[171,172,173];SITE=133
WATER={"SOL","TIP3","HOH","WAT","NA","CL","SOD","CLA","K"}
def stop(m):raise SystemExit("STOP: "+m)
def to_gro(t,o):
    if shutil.which("gmx") is None:stop("gmx not found in this session (install GROMACS, or run with --gro file.gro)")
    if not os.path.exists(t):stop("missing "+t)
    r=sp.run(["gmx","editconf","-f",t,"-o",o],capture_output=True,text=True)
    if r.returncode or not os.path.exists(o):stop("editconf failed:\n"+(r.stdout+r.stderr)[-1200:])
def parse(g):
    L=open(g).read().split("\n");n=int(L[1]);at=[]
    for l in L[2:2+n]:
        rn=l[5:10].strip()
        if rn in WATER:break
        at.append((int(l[0:5]),rn,l[10:15].strip(),float(l[20:28]),float(l[28:36]),float(l[36:44])))
    b=[float(x) for x in L[2+n].split()]
    if len(b)>3 and any(abs(x)>1e-6 for x in b[3:]):stop("box is not orthorhombic")
    return at,b[:3]
def heavy(nm):return not nm.lstrip("0123456789").startswith("H")
def mind(A,B,box):
    best=(9e9,"","")
    for a in A:
        for b in B:
            d2=0.0
            for k in(0,1,2):
                d=a[3+k]-b[3+k];d-=box[k]*round(d/box[k]);d2+=d*d
            if d2<best[0]:best=(d2,a[2],b[2])
    return (best[0]**0.5,best[1],best[2])
def table(g):
    at,box=parse(g);res={}
    for a in at:
        if heavy(a[2]):res.setdefault(a[0],[]).append(a)
    if SITE not in res:stop("residue 133 not found")
    rn=res[SITE][0][1];print("  residue 133 is",rn,"with",len(res[SITE]),"heavy atoms; box %.3f %.3f %.3f nm; protein atoms read %d"%(box[0],box[1],box[2],len(at)))
    if rn!="MET":stop("residue 133 is not MET in this WT structure")
    out={}
    for r in CRIT+CONS:
        if r not in res:stop("residue %d not found"%r)
        out[r]=(res[r][0][1],)+mind(res[SITE],res[r],box)
    return out
G=[a for a in sys.argv[1:] if a!="--gro"]
if "--gro" in sys.argv:srcs=[("given",G[0])]
else:
    srcs=[]
    for i in(1,2,3):
        o=W+"wt_prod_rep%d_start.gro"%i;to_gro(R+"wt_v2/wt_prod_rep%d.tpr"%i,o);srcs.append(("wt rep %d"%i,o))
T=[]
for nm,g in srcs:
    print(nm);T.append(table(g))
mx=max(abs(t[r][1]-T[0][r][1]) for t in T for r in T[0]) if len(T)>1 else 0.0
print("\nMax difference in any distance between the WT starting structures: %.4f nm"%mx)
print("\nMinimum heavy-atom distance from residue 133 (Met, UniProt Met1783) in the WT production starting structure (minimum image, nm)")
print("%-5s %-5s %-8s %-9s %s"%("res","name","seq sep","dist nm","atom pair (133 - target)"))
rows=[]
for lab,L in(("criterion 2",CRIT),("consistency observation",CONS)):
    print("--",lab)
    for r in L:
        nm,d,a1,a2=T[0][r];print("%-5d %-5s %-8d %-9.3f %s - %s"%(r,nm,abs(r-SITE),d,a1,a2));rows.append([r,nm,abs(r-SITE),"%.3f"%d,a1,a2,lab])
out=R+"analysis2/distances_from_133.csv"
try:
    with open(out,"w",newline="") as f:
        w=csv.writer(f);w.writerow(["res","resname","seq_sep","min_heavy_dist_nm","atom133","atom_target","group"]);w.writerows(rows)
    print("\nWrote",out)
except OSError as e:print("\nCould not write CSV:",e)
print("END OK")
