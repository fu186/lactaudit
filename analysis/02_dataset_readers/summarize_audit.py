#!/usr/bin/env python3
"""
summarize_audit.py -- aggregate all audit JSONs under /root/audit/results and build the summary table
"""
import json, os, glob, sys
import pandas as pd

RES = '/root/audit/results'

META = {
    'AAA_GSE166676':  ('Abdominal aortic aneurysm', 'GSE166676', 'human aortic tissue scRNA-seq', 'highest in neutrophils'),
    'ccRCC_GSE156632':('Clear cell renal cell carcinoma', 'GSE156632', 'human ccRCC tissue scRNA-seq', 'highest in VSMC/plasma cells'),
    'AS_GSE159677_ours':('Atherosclerosis (this study)', 'GSE159677', 'human carotid plaque scRNA-seq', 'highest in macrophages'),
    'AS_mouse_GSE214611': ('Atherosclerosis (mouse)', 'GSE214611', 'mouse aorta scRNA-seq', 'highest in VSMC'),
    'CAD_GSE131778':  ('Coronary artery disease', 'GSE131778', 'human coronary artery scRNA-seq', 'highest in DC'),
    'DFU_GSE165816':  ('Diabetic foot ulcer', 'GSE165816', 'human DFU tissue scRNA-seq', 'highest in neutrophils'),
    'PDAC_GSE212966': ('Pancreatic ductal adenocarcinoma', 'GSE212966', 'human PDAC tissue scRNA-seq', 'highest in macrophages'),
    'TAA_GSE155468':  ('Thoracic aortic aneurysm (sample TAA1)', 'GSE155468', 'human TAA tissue scRNA-seq', 'highest in neutrophils'),
    'TAA_TAA2_GSE155468': ('Thoracic aortic aneurysm (sample TAA2)', 'GSE155468', 'human TAA tissue scRNA-seq', 'highest in neutrophils'),
    'TAA_TAA3_GSE155468': ('Thoracic aortic aneurysm (sample TAA3)', 'GSE155468', 'human TAA tissue scRNA-seq', 'highest in macrophages'),
    'TAA_TAA4_GSE155468': ('Thoracic aortic aneurysm (sample TAA4)', 'GSE155468', 'human TAA tissue scRNA-seq', 'highest in macrophages'),
    'TAA_TAA5_GSE155468': ('Thoracic aortic aneurysm (sample TAA5)', 'GSE155468', 'human TAA tissue scRNA-seq', 'highest in neutrophils'),
    'TAA_TAA6_GSE155468': ('Thoracic aortic aneurysm (sample TAA6)', 'GSE155468', 'human TAA tissue scRNA-seq', 'highest in macrophages'),
    'TAA_TAA7_GSE155468': ('Thoracic aortic aneurysm (sample TAA7)', 'GSE155468', 'human TAA tissue scRNA-seq', 'highest in macrophages'),
    'TAA_TAA8_GSE155468': ('Thoracic aortic aneurysm (sample TAA8)', 'GSE155468', 'human TAA tissue scRNA-seq', 'highest in endothelial cells'),
    'ccRCC2_GSE171306':('Clear cell renal cell carcinoma (2)', 'GSE171306', 'human ccRCC tissue scRNA-seq', '—'),
    'GBM_GSE162631':  ('Glioblastoma', 'GSE162631', 'human GBM scRNA-seq', '—'),
}


def main():
    rows = []
    for f in sorted(glob.glob(os.path.join(RES, '*.json'))):
        d = json.load(open(f))
        name = d['name']
        a = d.get('audit', {})
        t1, t2, t3, t4 = d.get('T1'), d.get('T2'), d.get('T3'), d.get('T4')
        rows.append({
            'dataset': name,
            'disease': META.get(name, ('?',))[0],
            'GEO': META.get(name, ('?', '?'))[1],
            'n_cells': d.get('n_cells'),
            'reproduced_top': str(d.get('reproduced_top'))[:14],
            'SSR': round(t1['SSR'], 3) if t1 and 'SSR' in t1 else None,
            'T1': ('PASS' if t1 and t1.get('SSR', 0) > 1 else 'FAIL') if t1 else '—',
            'T2p': round(t2.get('empirical_p'), 4) if t2 and 'empirical_p' in t2 else None,
            'T2': ('PASS' if t2 and t2.get('pass') else 'FAIL') if t2 else '—',
            'rho(det. rate, score)': round(t4['rho_detect_vs_score'], 3) if t4 and 'rho_detect_vs_score' in t4 else None,
            'T4': ('PASS' if t4 and not t4.get('depth_driven') else 'FAIL') if t4 else '—',
            'verdict': 'PASS' if a.get('overall_pass') else 'FAIL',
        })
    df = pd.DataFrame(rows)
    out = os.path.join(RES, 'audit_summary.csv')
    df.to_csv(out, index=False)
    print(df.to_string(index=False))
    n = len(df)
    npass = (df['verdict'] == 'PASS').sum()
    print(f"\naudited {n} papers/datasets in total; original conclusions pass in {npass} and fail in {n - npass}, "
          f"survival rate {npass / n * 100:.1f}%")
    print(f"[saved] {out}")


if __name__ == '__main__':
    main()
