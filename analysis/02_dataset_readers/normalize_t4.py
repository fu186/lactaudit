#!/usr/bin/env python3
"""
normalize_t4.py -- unify the inconsistently named T4 fields in results/*.json into canonical
key names, and write out the single authoritative summary table used in the paper, paper_table.csv.

Canonical keys:
  rho_detect_vs_score  -- rank correlation (LRG detection rate, lactylation score)  i.e. the T4 primary metric
  rho_ngene_vs_score   -- rank correlation (mean gene count, lactylation score)
  depth_driven         -- True means T4 fails
Legacy keys accepted: rho_detect_vs_score / rho_ngene_vs_score / 'p' (for p_ngene only)
"""
import json, glob, os, re
import pandas as pd, numpy as np

RES = '/root/audit/results'
# Kept as the original directory name: make_bib_figs.py reads from this path.
OUT = '/workspace/论文数据表'

os.makedirs(OUT, exist_ok=True)

# Manually curated metadata (data source, species, disease)
META = {
    'AAA_GSE166676':      dict(accession='GSE166676', disease='Abdominal aortic aneurysm', species='Human', tissue='Aortic wall'),
    'AS_GSE159677_ours':  dict(accession='GSE159677', disease='Carotid atherosclerosis',  species='Human', tissue='Carotid plaque'),
    'AS_mouse_GSE214611': dict(accession='GSE214611', disease='Atherosclerosis',          species='Mouse', tissue='Aorta'),
    'CAD_GSE131778':      dict(accession='GSE131778', disease='Coronary artery disease',  species='Human', tissue='Coronary artery'),
    'DFU_GSE165816':      dict(accession='GSE165816', disease='Diabetic foot ulcer',      species='Human', tissue='Foot ulcer tissue'),
    'PDAC_GSE212966':     dict(accession='GSE212966', disease='Pancreatic ductal adenocarcinoma', species='Human', tissue='Tumour'),
    'ccRCC_GSE156632':    dict(accession='GSE156632', disease='Clear cell renal cell carcinoma',  species='Human', tissue='Tumour'),
    'ccRCC2_GSE171306':   dict(accession='GSE171306', disease='Clear cell renal cell carcinoma',  species='Human', tissue='Tumour'),
}

def _num(x):
    try:
        v = float(x)
        return v if v == v else None
    except (TypeError, ValueError):
        return None

def canon_t4(t4):
    if not t4:
        return None
    det = _num(t4.get('rho_detect_vs_score'))
    ng  = _num(t4.get('rho_ngene_vs_score'))
    p_det = _num(t4.get('p_detect', t4.get('p')))
    dd = t4.get('depth_driven')
    if dd is None:
        dd = bool((det or 0) > 0.5 or (ng or 0) > 0.5)
    return dict(rho_detect_vs_score=det, rho_ngene_vs_score=ng,
                p_detect=p_det, depth_driven=bool(dd))

def canon_t1(t1):
    if not t1: return None
    return dict(H_real=t1.get('H_real'), H_hk=t1.get('H_hk_mean', t1.get('H_hk')),
                SSR=t1.get('SSR'),
                pass_=bool(t1.get('SSR', 0) > 1) if t1.get('SSR') is not None else None)

def canon_t2(t2):
    if not t2: return None
    return dict(observed_effect=t2.get('observed_effect', t2.get('observed')),
                null_mean=t2.get('null_mean'),
                null_p95=t2.get('null_p95'),
                empirical_p=t2.get('empirical_p'),
                pass_=t2.get('pass'))

rows = []
patched = []
for f in sorted(glob.glob(os.path.join(RES, '*.json'))):
    d = json.load(open(f))
    name = d['name']
    t4 = canon_t4(d.get('T4'))
    # Write the canonicalised fields back
    if d.get('T4') and t4:
        d['T4'] = t4
        json.dump(d, open(f, 'w'), indent=2, ensure_ascii=False)
        patched.append(name)
    meta = META.get(name, {})
    if name == 'TAA_GSE155468':
        meta = dict(accession='GSE155468', disease='Thoracic aortic aneurysm',
                    species='Human', tissue='Ascending aorta', patient='TAA1')
    elif name.startswith('TAA_TAA'):
        meta = dict(accession='GSE155468', disease='Thoracic aortic aneurysm',
                    species='Human', tissue='Ascending aorta',
                    patient=name.split('_')[1])
    rows.append(dict(
        dataset=name,
        accession=d.get('accession', meta.get('accession', '')),
        disease=meta.get('disease'), species=meta.get('species'), tissue=meta.get('tissue'),
        patient=meta.get('patient'),
        n_cells=d.get('n_cells'),
        n_lrg=d.get('n_lrg_avail'),
        top=d.get('reproduced_top'),
        T1_SSR=(d.get('T1') or {}).get('SSR'),
        T1_H_real=(d.get('T1') or {}).get('H_real'),
        T1_H_hk=(d.get('T1') or {}).get('H_hk_mean', (d.get('T1') or {}).get('H_hk')),
        T2_p=(d.get('T2') or {}).get('empirical_p'),
        T3_rho_ngene=(d.get('T3') or {}).get('rho_ngene', (d.get('T3') or {}).get('rho_hk')) and (d.get('T3') or {}).get('rho_ngene'),
        T4_rho=(d.get('T4') or {}).get('rho_detect_vs_score'),
        T4_rho_ngene=(d.get('T4') or {}).get('rho_ngene_vs_score'),
        T4_pass=(not (d.get('T4') or {}).get('depth_driven')) if d.get('T4') else None,
        overall_pass=(d.get('audit') or {}).get('overall_pass'),
    ))

df = pd.DataFrame(rows)
df.to_csv(f'{OUT}/paper_table.csv', index=False)
print(df.to_string(index=False))
print(f"\n[normalized] {len(patched)} json files rewritten")
print(f"[saved] {OUT}/paper_table.csv")

# Summary statistics
print("\n=== Summary statistics ===")
print(f"Number of datasets: {len(df)}")
print(f"Surviving (overall_pass=True): {int((df.overall_pass == True).sum())}")
print(f"T1 failed (SSR<=1): {int((df.T1_SSR <= 1).sum())} / {df.T1_SSR.notna().sum()}")
print(f"T2 failed (p>=0.05): {int((df.T2_p >= 0.05).sum())} / {df.T2_p.notna().sum()}")
t4 = df.T4_rho.dropna()
print(f"T4 computable: {len(t4)}  mean={t4.mean():+.4f}  median={t4.median():+.4f}  min={t4.min():+.4f}  max={t4.max():+.4f}")
print(f"T4 >0.5: {int((t4>0.5).sum())}/{len(t4)}   T4 >0.7: {int((t4>0.7).sum())}/{len(t4)}")
