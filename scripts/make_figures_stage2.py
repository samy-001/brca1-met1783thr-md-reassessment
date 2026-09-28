import os,sys,glob,csv
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib.lines import Line2D
_c=["/content/drive/MyDrive/brca1_md_v2/"]+[q+"/" for q in glob.glob("/content/drive/.shortcut-targets-by-id/*/brca1_md_v2")]
A=sys.argv[1] if len(sys.argv)>1 else next((q+"analysis2/" for q in _c if os.path.isdir(q)),None)
if A is None:raise SystemExit("STOP: brca1_md_v2 not found (add the shortcut to My Drive and mount Drive)")
O=sys.argv[2] if len(sys.argv)>2 else A+"figures/"
os.makedirs(O,exist_ok=True)
WT=["wt_prod_rep%d"%i for i in(1,2,3)];MU=["mut_prod_rep%d"%i for i in(1,2,3)]
CW,CM="#1f77b4","#ff7f0e";LS=["-","--",":"]
HYP=[28,43,44,45,46,47,88,89,137,162,163,164,165,174];CONF=[28,165]
def xv(p):
    if not os.path.exists(p):raise SystemExit("STOP: missing "+p)
    return np.array([[float(x) for x in l.split()[:2]] for l in open(p) if l.strip() and l[0] not in "#@"])
def rm(y,w=50):return np.convolve(y,np.ones(w)/w,mode="valid")
def marks(ax,res=True):
    for r in HYP:ax.axvspan(r-.5,r+.5,color="#b0b0b0" if r not in CONF else "#e6550d",alpha=.25 if r not in CONF else .35,lw=0)
    ax.axvline(133,color="k",ls="--",lw=.8)
    for g in(44.5,166.5):ax.axvline(g,color="k",ls=":",lw=.8)
def runs(pat,fn):
    return [(n,fn(A+n+pat),CW if n.startswith("wt") else CM,LS[int(n[-1])-1]) for n in WT+MU]
def leg(ax,extra=()):
    h=[Line2D([],[],color=CW,label="WT (3 replicates)"),Line2D([],[],color=CM,label="Met1783Thr (3 replicates)")]+list(extra)
    ax.legend(handles=h,loc="lower left",bbox_to_anchor=(0,1.02),ncol=3,fontsize=8,frameon=False)
# Figure A: RMSD and Rg
fig,ax=plt.subplots(2,1,figsize=(7,6),sharex=True)
for n,d,c,l in runs("_nj_rmsd.xvg",xv):ax[0].plot(d[24:len(d)-25,0],rm(d[:,1]),color=c,ls=l,lw=1)
for n,d,c,l in runs("_nj_gyrate.xvg",xv):ax[1].plot(d[24:len(d)-25,0]/1000,rm(d[:,1]),color=c,ls=l,lw=1)
for a,t in zip(ax,("Backbone RMSD (nm)","Radius of gyration (nm)")):a.set_ylabel(t);a.axvspan(10,20,color="#dddddd",alpha=.5,lw=0)
ax[1].set_xlabel("Time (ns)");leg(ax[0])
fig.tight_layout();fig.savefig(O+"fig_rerun_global.png",dpi=300);plt.close(fig)
# Figure B: RMSF, helix, sheet
def dssp(p):
    S=[l.rstrip("\n") for l in open(p) if l.strip()]
    h=np.array([[ch in"HGI" for ch in s] for s in S]).mean(0);e=np.array([[ch in"EB" for ch in s] for s in S]).mean(0)
    return h,e
D={n:dssp(A+n+"_nj_ss206.dat") for n in WT+MU}
fig,ax=plt.subplots(3,1,figsize=(10,9),sharex=True)
for n,d,c,l in runs("_nj_rmsf.xvg",xv):ax[0].plot(d[:,0],d[:,1],color=c,ls=l,lw=1)
for n in WT+MU:
    c=CW if n.startswith("wt") else CM;l=LS[int(n[-1])-1];x=np.arange(1,207)
    ax[1].plot(x,D[n][0],color=c,ls=l,lw=1);ax[2].plot(x,D[n][1],color=c,ls=l,lw=1)
for a,t in zip(ax,("RMSF (nm)","Helix fraction (H, G, I)","Sheet fraction (E, B)")):a.set_ylabel(t);marks(a)
ax[2].set_xlabel("Residue (local numbering)");ax[2].set_xlim(1,206)
leg(ax[0],[Patch(color="#b0b0b0",alpha=.4,label="Table 5 hypothesis residue"),Patch(color="#e6550d",alpha=.4,label="Confirmed (28, 165)"),Line2D([],[],color="k",ls="--",lw=.8,label="Residue 133"),Line2D([],[],color="k",ls=":",lw=.8,label="Chain gaps (44/45, 166/167)")])
fig.tight_layout();fig.savefig(O+"fig_rerun_per_residue.png",dpi=300);plt.close(fig)
# Figure C: contacts with residue 133
p=A+"contacts_133_all.csv"
if os.path.exists(p):
    rows=list(csv.DictReader(open(p)))
    fig,ax=plt.subplots(figsize=(10,3.5))
    for n in WT+MU:
        x=[int(r["res"]) for r in rows];y=[float(r[n]) for r in rows]
        ax.plot(x,y,color=CW if n.startswith("wt") else CM,ls=LS[int(n[-1])-1],lw=1)
    marks(ax);ax.set_xlim(1,206);ax.set_ylim(0,1.02);ax.set_xlabel("Residue (local numbering)");ax.set_ylabel("Contact fraction with residue 133")
    leg(ax)
    fig.tight_layout();fig.savefig(O+"fig_rerun_contacts133.png",dpi=300);plt.close(fig)
else:print("contacts csv not found, skipped figure C")
print("Wrote:",sorted(os.listdir(O)));print("END OK")
