import glob,os
R="/content/drive/MyDrive/brca1_md_v2/"
def xv(q):return [[float(x) for x in l.split()] for l in open(q) if l[:1] not in "#@" and l.strip()]
reps=["wt_prod_rep1","wt_prod_rep2","wt_prod_rep3","mut_prod_rep1","mut_prod_rep2","mut_prod_rep3"]
for n in reps:
    for tag,folder in (("raw","analysis2"),("pbc","analysis2_pbc")):
        p=R+folder+"/"+n+"_rmsd.xvg"
        if not os.path.exists(p):
            print(n,tag,"MISSING",p);continue
        d=xv(p)
        jumps=[(abs(d[i+1][1]-d[i][1]),d[i][0],d[i+1][0]) for i in range(len(d)-1)]
        jumps.sort(reverse=True)
        top=jumps[:3]
        print(n,tag,"| top jumps (nm, at ns->ns):",["%.3f at %.2f->%.2f"%(j,a,b) for j,a,b in top])
print("END OK")
