import os,statistics as st
R="/content/drive/MyDrive/brca1_md_v2/analysis2/"
WT=["wt_prod_rep%d"%i for i in(1,2,3)];MU=["mut_prod_rep%d"%i for i in(1,2,3)]
def stop(m):raise SystemExit("STOP: "+m)
def xv(p):
    if not os.path.exists(p):stop("missing "+p)
    return [[float(x) for x in l.split()] for l in open(p) if l.strip() and l[0] not in "#@"]
for n in WT+MU:
    r=xv(R+n+"_nj_rmsd.xvg");g=xv(R+n+"_nj_gyrate.xvg")
    rl=[a for a in r if a[0]>=10.0]          # rmsd.xvg time is ns (an_nojump.py used -tu ns)
    gl=[a for a in g if a[0]>=10000.0]       # gyrate.xvg time is ps (no -tu flag was passed)
    if not rl:stop(n+" rmsd: no frames >=10ns, check time column/units, got max %.3f"%max(a[0] for a in r))
    if not gl:stop(n+" gyrate: no frames >=10000ps, check time column/units, got max %.3f"%max(a[0] for a in g))
    ry=[a[1] for a in rl];gy=[a[1] for a in gl]
    print("%-14s rmsd mean %.4f sd %.4f last %.4f max %.4f | rg mean %.4f sd %.4f min %.4f max %.4f  (n_rmsd=%d n_rg=%d)"%(n,st.mean(ry),st.pstdev(ry),ry[-1],max(ry),st.mean(gy),st.pstdev(gy),min(gy),max(gy),len(rl),len(gl)))
print("END OK")
