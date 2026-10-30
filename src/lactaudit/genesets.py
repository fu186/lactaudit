"""
Predefined gene sets.

All gene symbols are uppercase HGNC-style.  Mouse gene symbols in GEO matrices
are usually also uppercase (Actb -> ACTB), so the same sets are used for both
species after uppercasing the input gene names.

Provenance of LACTYLATION_SET
-----------------------------
The 116-gene lactylation-related set (LRG) was assembled from the union of
genes reported in lactylation-focused single-cell studies:

    Davis FM et al. J Exp Med 2021;218:e20201839
    Li Y et al.    Circulation 2020;142:1374-1388
    Rho H et al.   Cell Metab 2023;35:1406-1423
    Yang Z et al.  Transl Lung Cancer Res 2025;14:5447-5464

It has three functional blocks:
  (1) chromatin writers / erasers  (EP300, CREBBP, HDAC1-11, SIRT1-7, KAT*)
  (2) lactate production & transport (LDHA/B, SLC16A1-14, PKM, HK1-3, GPI, ENO*)
  (3) repeatedly reported lactylation-responsive genes
      (NLRP3, GSDMD, HMGB1, MECP2, METTL3, YTHDF2, SLC7A11, ...)

The set is intentionally broad.  That is precisely why it fails the depth
audit: many of its members are ubiquitously expressed, so a summed score over
them tracks sequencing depth more than cell identity.
"""

# ---------------------------------------------------------------------------
# Lactylation-related gene set (LRG) — 116 genes
# ---------------------------------------------------------------------------
LACTYLATION_SET = [
    # chromatin writers / erasers
    "EP300", "CREBBP", "KAT2A", "KAT2B", "KAT5", "KAT6A", "KAT6B", "KAT7", "KAT8",
    "HDAC1", "HDAC2", "HDAC3", "HDAC4", "HDAC5", "HDAC6", "HDAC7", "HDAC8",
    "HDAC9", "HDAC10", "HDAC11",
    "SIRT1", "SIRT2", "SIRT3", "SIRT4", "SIRT5", "SIRT6", "SIRT7",
    "AARS1", "AARS2", "TRIM33", "DPF2", "NEDD4",
    # lactate production & transport
    "LDHA", "LDHB", "LDHC",
    "SLC16A1", "SLC16A2", "SLC16A3", "SLC16A4", "SLC16A5", "SLC16A6", "SLC16A7",
    "SLC16A8", "SLC16A9", "SLC16A10", "SLC16A11", "SLC16A12", "SLC16A13", "SLC16A14",
    "PKM", "PKM2", "PKLR", "PFKP", "PFKL", "PFKM",
    "HK1", "HK2", "HK3", "GPI", "ENO1", "ENO2", "ENO3",
    "PDK1", "PDK2", "PDK3", "PDK4",
    "HIF1A", "HIF1AN", "MPC1", "MPC2",
    "PDHA1", "PDHB", "PDHX", "DLAT", "DLD",
    # reported lactylation-responsive genes
    "IL10", "ARG1", "IDH2", "MRC1", "CD163", "TGFB1",
    "SDHA", "FH", "ACLY", "OGDH", "DLST", "SUCLG1", "MDH2",
    "CDKN1A", "CDKN2A", "SERPINE1", "MMP9", "VEGFA", "LRG1",
    "YTHDF2", "METTL3", "SLC7A11", "NLRP3", "GSDMD", "CASP1",
    "HMGB1", "MECP2", "RUNX1", "SOX10", "NR4A3", "TRAP1", "CHI3L1",
    "SNAI1", "SCD", "VCAM1", "IL6", "IL1B", "TNF", "NFKB1", "RELA", "STAT3",
]

# ---------------------------------------------------------------------------
# Housekeeping panel (HK) — 40 genes
# ---------------------------------------------------------------------------
# Ribosomal proteins, cytoskeletal genes and universal metabolic enzymes.
# These are the genes a summed score *should* be blind to, because every cell
# needs them.  Any apparent cell-type signal in this panel is an artefact.
HOUSEKEEPING_SET = [
    "ACTB", "GAPDH", "TUBB", "PPIA", "RPL13A", "RPS18", "B2M", "HPRT1", "TBP",
    "YWHAZ", "PGK1", "EEF1A1", "RPLP0", "UBC", "RPS3", "RPL8", "RPSA",
    "RPL10", "RPL11", "RPL12", "RPS5", "RPS6", "RPS7", "RPS9", "RPL4", "RPL6",
    "RPS2", "RPS4X", "RPL3", "RPL5", "ACTG1", "TUBB4B", "TUBA1B", "PPIB",
    "CCT3", "CCT5", "PSMA1", "PSMB1", "RAN", "NAA10",
]

# ---------------------------------------------------------------------------
# Expanded housekeeping pool (139 genes)
# ---------------------------------------------------------------------------
# Used for the *upsizing* arm of the size-matched null: if the target set has
# 110 genes and the HK panel only has 40, one can either shrink the target or
# grow the control.  This pool makes the second option possible.
HOUSEKEEPING_POOL = sorted(set(HOUSEKEEPING_SET) | {
    "ACTB2", "TUBA1A", "RPS10", "RPS11", "RPS12", "RPS13", "RPS14", "RPS15",
    "RPS15A", "RPS16", "RPS17", "RPS19", "RPS20", "RPS21", "RPS23", "RPS24",
    "RPS25", "RPS26", "RPS27", "RPS27A", "RPS28", "RPS29",
    "RPL7", "RPL7A", "RPL9", "RPL10A", "RPLP1", "RPLP2",
    "RPL14", "RPL15", "RPL17", "RPL18", "RPL18A", "RPL19", "RPL21", "RPL22",
    "RPL23", "RPL23A", "RPL24", "RPL26", "RPL27", "RPL27A", "RPL28", "RPL29",
    "RPL30", "RPL31", "RPL32", "RPL34", "RPL35", "RPL35A", "RPL36", "RPL36A",
    "RPL37", "RPL37A", "RPL38", "RPL39", "RPL41",
    "EEF1B2", "EEF2", "EIF4B",
    "ALDOA", "ENO1", "PGAM1", "TPI1", "SLC25A3",
    "ATP5F1A", "ATP5F1B", "UQCRB", "NDUFA1", "NDUFB1",
    "COX4I1", "COX5A", "COX6C", "COX7A2", "COX7C", "SDHB", "CS", "OGDH",
    "ACTN4", "FLNA", "VIM", "CLTA", "CLTB", "ARF1", "RAB7A", "CAPZB",
    "CFL1", "GDI2", "YWHAB", "YWHAG", "YWHAE", "CALM1", "CALM2",
    "PSMA2", "PSMA3", "PSMA4", "PSMB2", "PSMB3", "PSMB4", "PSMB5", "PSMB6",
})

# ---------------------------------------------------------------------------
# Marker panel for ten major cell types
# ---------------------------------------------------------------------------
# Deliberately canonical and small (5-10 genes per type).  The goal is not
# fine-grained annotation but a stable, comparable label across datasets run
# by different groups with different pipelines.
CELL_TYPE_MARKERS = {
    "Endothelial": ["PECAM1", "VWF", "CDH5", "CLDN5", "EGFL7", "RAMP2", "EMCN", "PLVAP"],
    "VSMC": ["ACTA2", "MYH11", "TAGLN", "CNN1", "DES", "LMOD1", "MYL9"],
    "Fibroblast": ["DCN", "LUM", "COL1A1", "COL1A2", "PDGFRA", "COL3A1"],
    "Macrophage": ["CD68", "CD14", "LYZ", "CSF1R", "CD163", "MRC1", "C1QA", "C1QB", "ITGAM", "FCGR3A"],
    "Tcell": ["CD3D", "CD3E", "CD2", "CD3G", "TRAC", "IL7R", "CD8A", "CD4"],
    "Bcell": ["CD79A", "CD79B", "MS4A1", "CD19", "BANK1", "IGHM"],
    "Plasma": ["MZB1", "JCHAIN", "IGHG1", "DERL3", "SDC1"],
    "Mast": ["TPSAB1", "TPSB2", "CPA3", "MS4A2", "KIT"],
    "Neutrophil": ["FCGR3B", "CSF3R", "CXCR2", "S100A8", "S100A9"],
    "DC": ["CD1C", "CLEC9A", "LILRA4", "IRF8", "FCER1A"],
}

# ---------------------------------------------------------------------------
# Non-metabolic reference signatures
# ---------------------------------------------------------------------------
# Used to test whether the depth failure is specific to metabolism.  If a
# hypoxia or stemness signature fails just as badly, the problem is the
# estimator, not the biology being scored.
REFERENCE_SIGNATURES = {
    "Hypoxia": [
        "VEGFA", "SLC2A1", "PGK1", "PDK1", "LDHA", "BNIP3", "BNIP3L", "ADM",
        "NDRG1", "HILPDA", "ANGPTL4", "SERPINE1", "CA9", "P4HA1", "SLC2A3",
        "ENO1", "HK2", "PFKFB3", "EGLN3", "DDIT4", "MAFF", "PLIN2", "AK4",
        "INSIG2", "MIF", "TFRC", "CXCR4", "GAPDH",
    ],
    "Senescence": [
        "CDKN2A", "CDKN1A", "TP53", "SERPINE1", "IGFBP3", "IGFBP7", "IL6",
        "CXCL8", "CXCL1", "MMP3", "MMP1", "TGFB1", "SPARC", "TIMP1", "PLAU",
        "PLAT", "LMNB1", "HMGA1", "HMGA2", "TWIST1", "SERPINB2", "THBS1",
        "ICAM1", "CCL2", "CSF1", "AXL", "GDF15", "MKI67",
    ],
    "Stemness": [
        "SOX2", "NANOG", "POU5F1", "KLF4", "MYC", "LIN28A", "LIN28B",
        "ALDH1A1", "PROM1", "THY1", "ENG", "NES", "NOTCH1", "WNT3A", "LEF1",
        "TCF7", "BMI1", "EZH2", "KIT", "LGR5", "EPCAM", "ITGA6", "CD44",
        "CXCR4", "MME", "ABCG2", "TWIST1", "SNAI1",
    ],
    # Negative control: a ubiquitous programme with no expected cell-type
    # specificity.  If it shows a strong apparent signal, the estimator is at
    # fault by construction.
    "Apoptosis_NegCtrl": [
        "BAX", "BAK1", "CASP3", "CASP8", "CASP9", "BCL2", "BCL2L1", "MCL1",
        "APAF1", "CYCS", "DIABLO", "BID", "BAD", "PMAIP1", "FAS",
        "TNFRSF10B", "XIAP", "BIRC5",
    ],
}
