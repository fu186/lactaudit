import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt, json, glob, os
import numpy as np, pandas as pd

plt.rcParams['font.sans-serif']=['WenQuanYi Zen Hei','Noto Sans CJK JP','DejaVu Sans']
plt.rcParams['axes.unicode_minus']=False
RES='/root/audit/results'
files=sorted(glob.glob(os.path.join(RES,'TAA_TAA*.json')))
recs=[]
for f in files:
    d=json.load(open(f)); n=d['name'].replace('TAA_','').replace('_GSE155468','')
    recs.append({'p':n,'top':d.get('reproduced_top'),
                 'rho':d.get('T4',{}).get('rho_detect_vs_score'),
                 'ssr':d.get('T1',{}).get('SSR')})
df=pd.DataFrame(recs)

fig,axes=plt.subplots(1,2,figsize=(14,5.2),gridspec_kw={'width_ratios':[1,1.15]})
CMAP={'Macrophage':'#c0392b','Neutrophil':'#2980b9','Endothelial':'#27ae60'}

# left: top cell type flips
a=axes[0]
a.scatter(range(len(df)),[0]*len(df),s=0,c='white')
for i,r in df.reset_index(drop=True).iterrows():
    a.scatter(i,0,s=900,color=CMAP[r['top']],edgecolor='white',linewidth=2,zorder=3)
    a.text(i,0.28,r['top'][:6],ha='center',fontsize=9,fontweight='bold',color=CMAP[r['top']])
    a.text(i,-0.30,f"ρ={r['rho']:.2f}",ha='center',fontsize=8.5,color='#555')
a.set_xticks(range(len(df))); a.set_xticklabels(df['p'],fontsize=10)
a.set_xlabel('Thoracic aortic aneurysm patient samples (8 independent patients)',fontsize=11)
a.set_title('Same disease, same platform:\n"the cell type with the highest lactylation" flips repeatedly between patients',fontsize=12,fontweight='bold')
a.set_ylim(-0.6,0.6); a.set_yticks([])
for s in ['top','right','left']: a.spines[s].set_visible(False)
from matplotlib.lines import Line2D
a.legend(handles=[Line2D([],[],marker='o',ls='',ms=12,color=c,label=k) for k,c in CMAP.items()],
         loc='upper center',ncol=3,fontsize=9,frameon=False,bbox_to_anchor=(0.5,1.0))

# right: T4 rho vs SSR
a=axes[1]
x=np.arange(len(df)); w=0.38
a.bar(x-w/2,df['rho'],w,color='#c0392b',label='T4: ρ(detection rate, score)',edgecolor='white')
a.bar(x+w/2,df['ssr'],w,color='#2980b9',label='T1: SSR',edgecolor='white')
a.axhline(0.5,color='#c0392b',ls='--',lw=1.3,alpha=0.7)
a.axhline(1.0,color='#2980b9',ls='--',lw=1.3,alpha=0.7)
a.text(len(df)-0.4,0.52,'T4 pass line 0.5',fontsize=8,color='#c0392b',ha='right')
a.text(len(df)-0.4,1.03,'T1 pass line 1.0',fontsize=8,color='#2980b9',ha='right')
a.set_xticks(x); a.set_xticklabels(df['p'],fontsize=10)
a.set_ylabel('value',fontsize=11); a.set_ylim(0,1.35)
a.set_title('T4 exceeds the threshold throughout (7/7 > 0.7); T1 fails in half',fontsize=12,fontweight='bold')
a.legend(fontsize=9,frameon=False); a.grid(axis='y',alpha=0.3)
for s in ['top','right']: a.spines[s].set_visible(False)

plt.suptitle('Within-patient replication in GSE155468: the artefact reproduces stably across 8 independent patients',fontsize=13,fontweight='bold',y=1.03)
plt.tight_layout()
out='/workspace/审计图_TAA患者重复实验.png'
plt.savefig(out,dpi=170,bbox_inches='tight'); print('saved',out)

# also write a CSV
df.to_csv('/workspace/TAA八患者审计明细.csv',index=False,encoding='utf-8-sig')
print('saved TAA eight-patient audit table (.csv)')
