import os,sys,glob,shutil,csv,subprocess as sp
_c=["/content/drive/MyDrive/brca1_md_v2/"]+[q+"/" for q in glob.glob("/content/drive/.shortcut-targets-by-id/*/brca1_md_v2")]
R=next((q for q in _c if os.path.isdir(q)),None)
def stop(m):raise SystemExit("STOP: "+m)
THR=0.45;SHOW=0.10;SITE=133;NFR=2001
REP=[("wt_prod_rep%d"%i,"wt_v2") for i in(1,2,3)]+[("mut_prod_rep%d"%i,"mut_v2") for i in(1,2,3)]
RES=[r for r in range(1,207) if r!=SITE]
CRIT=[28,43,44,45,46,47,88,89,137,162,163,164,165,174];CONS=[171,172,173]
CHK={137:0.302,89:1.139,47:1.366,174:2.335}
def xv(p):return [[float(x) for x in l.split()] for l in open(p) if l.strip() and l[0] not in "#@"]
def rng(v):return "%.2f-%.2f"%(min(v),max(v))
def one(N,F):
    W="/content/an/c133_"+N+"/";os.makedirs(W,exist_ok=True)
    D=R+F+"/"
    for p in(D+N+".tpr",D+"stitched/"+N+"_0-20000ps.xtc"):
        if not os.path.exists(p):stop("missing "+p)
        shutil.copy(p,W)
    T=W+N+".tpr";X=W+N+"_0-20000ps.xtc";O=W+"pd.xvg"
    ref='group "Protein" and resnr %d and not name "H*"'%SITE
    sel='group "Protein" and resnr 1 to 206 and not resnr %d and not name "H*"'%SITE
    r=sp.run(["gmx","pairdist","-s",T,"-f",X,"-ref",ref,"-sel",sel,"-refgrouping","all","-selgrouping","res","-type","min","-o",O],capture_output=True,text=True)
    if r.returncode or not os.path.exists(O):stop("pairdist failed:\n"+(r.stdout+r.stderr)[-1500:])
    d=xv(O)
    if len(d)!=NFR:stop("%s: %d frames, expected %d"%(N,len(d),NFR))
    if {len(x) for x in d}!={len(RES)+1}:stop("%s: columns %s, expected %d"%(N,sorted({len(x) for x in d}),len(RES)+1))
    if N.startswith("wt"):
        bad=[(k,round(d[0][RES.index(k)+1],3),v) for k,v in CHK.items() if abs(d[0][RES.index(k)+1]-v)>0.02]
        if bad:stop("%s: frame-0 check failed (res, got, expected): %s"%(N,bad))
        print(N,"frame-0 check OK",{k:round(d[0][RES.index(k)+1],3) for k in CHK})
    f={r:sum(x[i+1]<=THR for x in d)/len(d) for i,r in enumerate(RES)}
    shutil.rmtree(W,ignore_errors=True)
    return f
def summarize(F,out=None):
    WT=[F[n] for n,_ in REP[:3]];MU=[F[n] for n,_ in REP[3:]]
    print("\nContact fraction with residue 133 (heavy atoms within %.2f nm), fraction of %d frames. Ranges over 3 replicates. EXPLORATORY."%(THR,NFR))
    print("%-5s %-6s %-11s %-11s %s"%("res","list","WT range","Mut range","ranges"))
    rows=[]
    for r in RES:
        w=[x[r] for x in WT];m=[x[r] for x in MU]
        lab="crit2" if r in CRIT else ("171-3" if r in CONS else "")
        ov="overlap" if not(min(m)>max(w) or max(m)<min(w)) else ("mut higher" if min(m)>max(w) else "mut lower")
        rows.append([r,lab,rng(w),rng(m),ov]+["%.3f"%x for x in w+m])
        if max(w+m)>=SHOW or lab:print("%-5d %-6s %-11s %-11s %s"%(r,lab,rng(w),rng(m),ov))
    if out:
        with open(out,"w",newline="") as f:
            c=csv.writer(f);c.writerow(["res","list","WT_range","Mut_range","ranges"]+[n for n,_ in REP]);c.writerows(rows)
        print("\nWrote",out)
def main():
    if R is None:stop("brca1_md_v2 not found (add the shared folder shortcut to My Drive, mount Drive)")
    if shutil.which("gmx") is None:stop("gmx not found")
    print("Using folder:",R);F={}
    for N,fo in REP:
        F[N]=one(N,fo);print(N,"done")
    summarize(F,R+"analysis2/contacts_133_all.csv")
    print("END OK")
if __name__=="__main__":main()
