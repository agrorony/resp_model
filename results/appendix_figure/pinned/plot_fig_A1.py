import csv, numpy as np, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
D="results/appendix_figure/pinned/"
COL={"vertisol":"#6B4F35","sand":"#D97B4F"}; LAB={"vertisol":"Vertisol (Bnei Re'em)","sand":"Sand (Rehovot)"}
sig={r['metric']:r for r in csv.DictReader(open("results/stage12/stage12_significance.csv"))}['peak_frac']
fig,ax=plt.subplots(figsize=(6.2,3.6))
x=np.linspace(0,1,500)
for s in ["vertisol","sand"]:
    c=np.load(D+f"curves_{s}.npy"); m=c.mean(0); se=c.std(0,ddof=1)/np.sqrt(c.shape[0])
    ax.plot(x,m,color=COL[s],lw=2,label=f"{LAB[s]} (n = {c.shape[0]})"); ax.fill_between(x,m-se,m+se,color=COL[s],alpha=0.35,lw=0)
ax.spines[['top','right']].set_visible(False)
ax.set_xlabel("Fraction of observation window"); ax.set_ylabel("Simulated respiration R(t) (a.u.)")
ax.legend(frameon=False,fontsize=8.5,loc="upper right")
ax.text(0.98,0.62,f"Peak position: Mann–Whitney\nrank-biserial r = {abs(float(sig['rank_biserial_r'])):.1f}, p < 1e-14",transform=ax.transAxes,ha="right",fontsize=8.5,
        bbox=dict(boxstyle="round",fc="white",ec="0.6"))
fig.tight_layout(); fig.savefig(D+"fig_A1_stage12_mean_se.png",dpi=200); fig.savefig(D+"fig_A1_stage12_mean_se.svg"); print("saved")
