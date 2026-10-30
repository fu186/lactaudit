import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import re
import numpy as np, json, glob, os

plt.rcParams['font.sans-serif'] = ['WenQuanYi Zen Hei', 'Noto Sans CJK JP', 'DejaVu Sans']
plt.rcParams['axes.unicode_minus'] = False

RES='/root/audit/results'
rows=[]
for f in sorted(glob.glob(os.path.join(RES,'*.json'))):
    d=json.load(open(f)); t4=d.get('T4'); t1=d.get('T1')
    rows.append({'name':d['name'],'t4':t4.get('rho_detect_vs_score') if t4 else None,
                 'ssr':t1.get('SSR') if t1 else None,
                 'p':d.get('T2',{}).get('empirical_p') if d.get('T2') else None})

fig,axes=plt.subplots(1,3,figsize=(16,5.2))
RED='#c0392b'; GRN='#27ae60'; GRY='#95a5a6'

# Panel 1: T4 rho
a=axes[0]
sub=[r for r in rows if r['t4'] is not None]
sub.sort(key=lambda x:-x['t4'])
names=[re.sub(r'_GSE.*','',r['name']).replace('AS_mouse','AS(m)').replace('_ours','(ours)').replace('TAA_TAA','TAA') for r in sub]; vals=[r['t4'] for r in sub]
cols=[RED if v>0.5 else GRN for v in vals]
a.barh(range(len(vals)),vals,color=cols,edgecolor='white')
a.axvline(0.5,color='k',ls='--',lw=1.5,label='pass threshold ρ=0.5')
a.set_yticks(range(len(names))); a.set_yticklabels(names,fontsize=9)
a.set_xlabel('ρ (gene-set detection rate vs. lactylation score)',fontsize=10)
a.set_title('T4: does cell-type ranking just follow sequencing depth?',fontsize=11,fontweight='bold')
a.legend(fontsize=9); a.set_xlim(0,1.05); a.grid(axis='x',alpha=0.3)

# Panel 2: SSR
a=axes[1]
sub=[r for r in rows if r['ssr'] is not None]
sub.sort(key=lambda x:x['ssr'])
names=[re.sub(r'_GSE.*','',r['name']).replace('AS_mouse','AS(m)').replace('_ours','(ours)').replace('TAA_TAA','TAA') for r in sub]; vals=[r['ssr'] for r in sub]
cols=[GRN if v>1 else RED for v in vals]
a.barh(range(len(vals)),vals,color=cols,edgecolor='white')
a.axvline(1.0,color='k',ls='--',lw=1.5,label='pass threshold SSR=1')
a.set_yticks(range(len(names))); a.set_yticklabels(names,fontsize=9)
a.set_xlabel('SSR = H(lactylation genes) / H(housekeeping genes)',fontsize=10)
a.set_title('T1: are housekeeping genes more discriminative?',fontsize=11,fontweight='bold')
a.legend(fontsize=9); a.grid(axis='x',alpha=0.3)

# Panel 3: T2 p
a=axes[2]
sub=[r for r in rows if r['p'] is not None]
sub.sort(key=lambda x:x['p'])
names=[re.sub(r'_GSE.*','',r['name']).replace('AS_mouse','AS(m)').replace('_ours','(ours)').replace('TAA_TAA','TAA') for r in sub]; vals=[r['p'] for r in sub]
cols=[GRN if v<0.05 else RED for v in vals]
a.barh(range(len(vals)),vals,color=cols,edgecolor='white')
a.axvline(0.05,color='k',ls='--',lw=1.5,label='pass threshold p=0.05')
a.set_yticks(range(len(names))); a.set_yticklabels(names,fontsize=9)
a.set_xlabel('Empirical p (expression-matched null)',fontsize=10)
a.set_title('T2: can random gene sets reach the same effect?',fontsize=11,fontweight='bold')
a.legend(fontsize=9); a.grid(axis='x',alpha=0.3)

plt.suptitle('Reproduction audit of a bioinformatics lactylation paper: 16 datasets x 3 control tests',fontsize=13,fontweight='bold',y=1.02)
plt.tight_layout()
out='/workspace/审计结果图_T1T2T4.png'
plt.savefig(out,dpi=170,bbox_inches='tight')
print('saved',out)
