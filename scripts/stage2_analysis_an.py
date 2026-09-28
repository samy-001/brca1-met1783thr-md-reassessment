import sys,os,re,shutil,subprocess as sp
N,F=sys.argv[1],sys.argv[2]
R="/content/drive/MyDrive/brca1_md_v2/"
D=R+F+"/";W="/content/an/"+N+"/";O=R+"analysis2/"
os.makedirs(W,exist_ok=True);os.makedirs(O,exist_ok=True)
def stop(m):raise SystemExit("STOP: "+m)
def gm(a,i=""):
    r=sp.run(["gmx"]+a,input=i,capture_output=True,text=True)
    if r.returncode:stop(" ".join(a[:2])+" failed:\n"+(r.stdout+r.stderr)[-1500:])
T0=D+N+".tpr";X0=D+"stitched/"+N+"_0-20000ps.xtc"
for p in(T0,X0):
    if not os.path.exists(p):stop("missing "+p)
    shutil.copy(p,W)
T=W+N+".tpr";X=W+os.path.basename(X0);I=W+"i.ndx"
gm(["make_ndx","-f",T,"-o",I],'r 1-206\n"Protein" & "r_1-206"\n"Backbone" & "r_1-206"\nq\n')
G={};k=None
for l in open(I):
    m=re.match(r"\[ (.+) \]",l)
    if m:k=m.group(1);G[k]=[]
    elif k:G[k]+=l.split()
P,B="Protein_&_r_1-206","Backbone_&_r_1-206"
if P not in G or B not in G:stop("index groups missing: "+str(list(G)))
nm=list(G);e=3290 if N.startswith("wt") else 3287
if len(G[P])!=e or len(G[B])!=618:stop("group sizes %d/%d, expected %d/618"%(len(G[P]),len(G[B]),e))
print("groups OK: protein",len(G[P]),"backbone",len(G[B]))
def run(t,a,i):gm([t,"-s",T,"-f",X,"-n",I]+a,i)
p,b=nm.index(P),nm.index(B)
run("rms",["-o",W+"rmsd.xvg","-tu","ns"],f"{b}\n{b}\n")
run("rmsf",["-o",W+"rmsf.xvg","-res"],f"{p}\n{p}\n")
run("gyrate",["-o",W+"gyrate.xvg"],f"{p}\n")
gm(["dssp","-s",T,"-f",X,"-o",W+"ss.dat","-sel","resnr 1 to 206"])
S=[l.strip().replace("=","") for l in open(W+"ss.dat")]
if {len(s) for s in S}!={206}:stop("ss line lengths "+str({len(s) for s in S}))
open(W+"ss206.dat","w").write("\n".join(S)+"\n")
for s,d in(("rmsd.xvg","rmsd.xvg"),("rmsf.xvg","rmsf.xvg"),("gyrate.xvg","gyrate.xvg"),("ss206.dat","ss206.dat")):shutil.copy(W+s,O+N+"_"+d)
def xv(q):return [[float(x) for x in l.split()] for l in open(q) if l[:1] not in "#@" and l.strip()]
r,g,f=xv(W+"rmsd.xvg"),xv(W+"gyrate.xvg"),dict((int(a),b) for a,b in xv(W+"rmsf.xvg"))
j=max(abs(y[1]-x[1]) for x,y in zip(r,r[1:]))
print(N,"| frames",len(r),"| rmsd last %.3f max %.3f max jump %.3f"%(r[-1][1],max(y[1] for y in r),j))
print("Rg min %.3f max %.3f | rmsf residues %d, res133 %.4f | ss frames %d"%(min(y[1] for y in g),max(y[1] for y in g),len(f),f[133],len(S)))
if j>0.5:print("WARNING: possible PBC jump")
print("END OK")
