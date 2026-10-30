# -*- coding: utf-8 -*-
"""Path resolution for the submission package -- shared by every script in the package.

The same code must run under two layouts: the working directory and the submission package.

  Working-directory layout
    <WORK>/投稿版_BMP/Manuscript_BMP_MethodsArticle.md        (directory names as-is)
    <WORK>/投稿版_3分/Manuscript_BMC_Bioinformatics_ResearchArticle.md
    <WORK>/投稿版_BIB/Manuscript_BIB_Article.md
    <WORK>/复算/t4r_bootstrap_corrected.csv
    <TOOL>/results/*.json              TOOL = the audit toolkit directory
    <TOOL>/lactcheck.py

  Submission-package layout
    <PKG>/01_manuscripts/A_BMP_methods_article/Manuscript_BMP_MethodsArticle.md
    <PKG>/01_manuscripts/B_BMC_research_article/Manuscript_BMC_Bioinformatics_ResearchArticle.md
    <PKG>/01_manuscripts/Z_archived_BIB/Manuscript_BIB_Article.md
    <PKG>/04_code_and_data/01_audit_library/lactcheck.py
    <PKG>/04_code_and_data/03_recompute/t4r_bootstrap_corrected.csv
    <PKG>/04_code_and_data/04_simulation/sim_t4_t4r.json
    <PKG>/04_code_and_data/07_result_records/*.json
    <PKG>/02_figures/*.png

Usage:
    from audit_paths import ROOT, RES, TOOL, manuscript, find, suppdir

    ROOT   root directory holding the manuscripts and all materials
    RES    per-dataset result-record directory (holds *.json and *_T4.csv)
    TOOL   directory of the audit library (holds lactcheck.py)
    manuscript('BMC'|'BMP'|'BIB') -> Path
    find('filename') -> Path
    suppdir() -> supplementary-material directory

The environment variables AUDIT_ROOT / AUDIT_DATA can override these.
"""
import os
from pathlib import Path

HERE = Path(__file__).resolve().parent

# manuscript filename -> tag
MS_FILES = {
    'BMP': 'Manuscript_BMP_MethodsArticle.md',
    'BMC': 'Manuscript_BMC_Bioinformatics_ResearchArticle.md',
    'BIB': 'Manuscript_BIB_Article.md',
}

# candidate names for the result-record directory
RESNAMES = ('07_result_records', 'results')

# candidate names for the supplementary-material directory
SUPNAMES = ('补充材料', '03_supplementary_tables')


def _candidates():
    """Build candidate search roots: several parent levels of the script directory,
    plus the first-level subdirectories of those levels."""
    out = []
    env = os.environ.get('AUDIT_ROOT')
    if env:
        out.append(Path(env))
    # go up 4 levels
    p = HERE
    for _ in range(5):
        out.append(p)
        p = p.parent
    # first- and second-level subdirectories of each candidate root
    for base in list(out):
        if base.is_dir():
            try:
                for sub in sorted(base.iterdir()):
                    if sub.is_dir():
                        out.append(sub)
            except PermissionError:
                pass
    # de-duplicate while preserving order
    seen, uniq = set(), []
    for p in out:
        k = str(p)
        if k not in seen:
            seen.add(k)
            uniq.append(p)
    return uniq


CANDS = _candidates()


def _first_dir_with(pred, limit=None):
    for d in CANDS[:limit]:
        try:
            if d.is_dir() and pred(d):
                return d
        except (PermissionError, OSError):
            continue
    return None


# ---- ROOT: a directory containing any of the named manuscripts -----------------
def _root():
    for d in CANDS:
        try:
            if not d.is_dir():
                continue
            for fn in MS_FILES.values():
                if list(d.glob(f'**/{fn}')):
                    return d
        except (PermissionError, OSError):
            continue
    raise SystemExit('Cannot locate the manuscript root directory; please set the AUDIT_ROOT environment variable.')


ROOT = _root()


def find(name):
    """Locate a file recursively under ROOT and return the first hit."""
    hits = sorted(ROOT.glob(f'**/{name}'))
    if not hits:
        raise SystemExit(f'{name} not found (search root: {ROOT}); please set AUDIT_ROOT.')
    return hits[0]


def manuscript(tag):
    return find(MS_FILES[tag])


def manuscripts():
    """Return {tag: Path}; tags whose file is missing are omitted."""
    out = {}
    for tag, fn in MS_FILES.items():
        hits = sorted(ROOT.glob(f'**/{fn}'))
        if hits:
            out[tag] = hits[0]
    return out


# ---- RES: per-dataset result-record directory --------------------------------
def _res():
    for d in CANDS:
        try:
            if not d.is_dir():
                continue
            for rn in RESNAMES:
                c = d / rn
                if c.is_dir() and list(c.glob('*.json')):
                    return c
        except (PermissionError, OSError):
            continue
    # fall back: any directory containing AAA_GSE166676.json
    for d in CANDS:
        try:
            if (d / 'AAA_GSE166676.json').exists():
                return d
        except (PermissionError, OSError):
            continue
    raise SystemExit('Cannot locate the result-record directory; please set AUDIT_ROOT.')


RES = _res()


# ---- TOOL: directory of the audit library -----------------------------------
def _tool():
    for d in CANDS:
        try:
            if d.is_dir() and (d / 'lactcheck.py').exists():
                return d
        except (PermissionError, OSError):
            continue
    # submission-package layout: 01_audit_library contains lactcheck.py
    for d in CANDS:
        try:
            if d.is_dir():
                hits = list(d.glob('**/lactcheck.py'))
                if hits:
                    return hits[0].parent
        except (PermissionError, OSError):
            continue
    return ROOT


TOOL = _tool()


# ---- supplementary-material directory ---------------------------------------
def suppdir():
    for d in CANDS:
        try:
            if not d.is_dir():
                continue
            for sn in SUPNAMES:
                c = d / sn
                if c.is_dir() and list(c.glob('Table_S*.csv')):
                    return c
        except (PermissionError, OSError):
            continue
    hits = sorted(ROOT.glob('**/Table_S1_dataset_metadata.csv'))
    if hits:
        return hits[0].parent
    return ROOT


# ---- raw-data directory (needed only by the estimator scripts) --------------
def datadir():
    env = os.environ.get('AUDIT_DATA')
    cands = [Path(env)] if env else []
    cands += [ROOT / 'data', ROOT.parent / 'data', Path('/root/audit/data')]
    for d in cands:
        try:
            if d.is_dir() and any((d / x).is_dir() for x in ('gse155468', 'gse166676')):
                return d
        except (PermissionError, OSError):
            continue
    return Path('/root/audit/data')


if __name__ == '__main__':
    print(f'ROOT = {ROOT}')
    print(f'RES  = {RES}')
    print(f'TOOL = {TOOL}')
    print(f'SUPP = {suppdir()}')
    print(f'DATA = {datadir()}')
    for t, p in manuscripts().items():
        print(f'  {t}: {p.name}')
