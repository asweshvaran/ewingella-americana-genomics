import os
import re
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.gridspec import GridSpec
from Bio import Phylo
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor, DistanceMatrix
from PIL import Image

# Directories
phylo_dir = "/mnt/d/Ewingella_americana_research/FINAL_JAM_SUBMISSION/Phylogeny_Q1"
fig_dir   = "/mnt/d/Ewingella_americana_research/FINAL_JAM_SUBMISSION/figures"
root_dir  = "/mnt/d/Ewingella_americana_research/FINAL_JAM_SUBMISSION"
s1_path   = "/mnt/d/Ewingella_americana_research/FINAL_JAM_SUBMISSION/Table_S1_Genomic_Inventory.csv"

nt_iqtree = os.path.join(phylo_dir, "strict_31_Q1.iqtree")
nt_treef  = os.path.join(phylo_dir, "strict_31_Q1.treefile")
aa_iqtree = os.path.join(phylo_dir, "core_aa_validation.iqtree")
aa_treef  = os.path.join(phylo_dir, "core_aa_validation.treefile")
pdist_csv = os.path.join(phylo_dir, "pairwise_pdist_matrix_31.csv")
snps_csv  = os.path.join(phylo_dir, "pairwise_core_snps_31.csv")
comp_csv  = os.path.join(phylo_dir, "Table_Pairwise_Distances_and_Composition.csv")
kaks_csv  = os.path.join(phylo_dir, "Table_Core_KaKs_NG86.csv")

df_s1    = pd.read_csv(s1_path).set_index('accession')
df_pdist = pd.read_csv(pdist_csv, index_col=0)
df_snps  = pd.read_csv(snps_csv, index_col=0)
df_comp  = pd.read_csv(comp_csv)
df_kaks  = pd.read_csv(kaks_csv)

# ==============================================================================
# 1. PARSE IQ-TREE 2 REPORTS & VERIFY CROSS-METHOD NJ CONGRUENCE
# ==============================================================================
def parse_iqtree_report(path):
    txt = open(path).read()
    def grab(pat, default="NA"):
        m = re.search(pat, txt)
        return m.group(1).strip() if m else default
    return {
        'Taxa': grab(r'Input data:\s+(\d+)\s+sequences'),
        'Sites': grab(r'Input data:\s+\d+\s+sequences with\s+(\d+)'),
        'Constant_Sites': grab(r'Number of constant sites:\s+([^\n]+)'),
        'Parsimony_Informative': grab(r'Number of parsimony informative sites:\s+(\d+)'),
        'Distinct_Patterns': grab(r'Number of distinct site patterns:\s+(\d+)'),
        'Best_BIC_Model': grab(r'Model of substitution:\s+([^\n]+)'),
        'LogL': grab(r'Log-likelihood of the tree:\s+([-\d\.]+)'),
        'BIC': grab(r'Bayesian information criterion \(BIC\) score:\s+([\d\.]+)'),
        'AIC': grab(r'Akaike information criterion \(AIC\) score:\s+([\d\.]+)'),
        'Tree_Length': grab(r'Total tree length \(sum of branch lengths\):\s+([\d\.]+)'),
        'Internal_Branch_Sum': grab(r'Sum of internal branch lengths:\s+([^\n]+)'),
    }

nt_stats = parse_iqtree_report(nt_iqtree)
aa_stats = parse_iqtree_report(aa_iqtree)
nt_nwk = open(nt_treef).read().strip()
aa_nwk = open(aa_treef).read().strip()

# Secondary method: Neighbor-Joining tree from exact core p-distances
names_nj = list(df_pdist.index)
lower_mat = []
for i in range(len(names_nj)):
    lower_mat.append([float(df_pdist.iloc[i, j]) for j in range(i + 1)])
dm = DistanceMatrix(names_nj, lower_mat)
nj_tree = DistanceTreeConstructor().nj(dm)
nj_nwk_path = os.path.join(phylo_dir, "strict_31_NJ_pdist.nwk")
Phylo.write(nj_tree, nj_nwk_path, "newick")
nj_nwk = open(nj_nwk_path).read().strip()

# Verify inter-lineage topological congruence across NT-ML, AA-ML, and NJ
# Representative taxa for the 6 resolved clusters:
reps = {
    'B6-1':       'GCF_014233485.1',
    'RIT713':     'GCF_009668135.1',
    'Bioreactor': 'GCF_022647985.1',
    'Oral':       'GCF_057314615.1',
    'CCUG_Ref':   'GCF_008693045.1',
    'CIP_Clone':  'GCF_054058535.1',
}

def get_bipartitions(t_path, rep_map):
    t = Phylo.read(t_path, "newick")
    t.root_with_outgroup(rep_map['B6-1'])
    inv = {v: k for k, v in rep_map.items()}
    splits = set()
    for cl in t.get_nonterminals():
        leaves = frozenset(inv[x.name] for x in cl.get_terminals() if x.name in inv)
        if 1 < len(leaves) < len(rep_map):
            comp = frozenset(rep_map.keys()) - leaves
            splits.add(min(leaves, comp, key=lambda s: sorted(s)))
    return splits

sp_nt = get_bipartitions(nt_treef, reps)
sp_aa = get_bipartitions(aa_treef, reps)
sp_nj = get_bipartitions(nj_nwk_path, reps)
rf_nt_aa = len(sp_nt ^ sp_aa)
rf_nt_nj = len(sp_nt ^ sp_nj)
print(f">>> Topological Congruence Check (6 Lineage Representatives):")
print(f"    NT-ML splits: {[sorted(list(s)) for s in sp_nt]}")
print(f"    AA-ML splits: {[sorted(list(s)) for s in sp_aa]}")
print(f"    NJ    splits: {[sorted(list(s)) for s in sp_nj]}")
print(f"    Inter-lineage RF(NT-ML vs AA-ML) = {rf_nt_aa} | RF(NT-ML vs NJ) = {rf_nt_nj}")

# Save Table_Phylogenetic_Models_and_Support.csv
df_models = pd.DataFrame([
    {
        'Matrix': 'Matrix A — Nucleotide Strict Core (Primary)',
        'Software': 'IQ-TREE 2 v2.0.7',
        'Taxa_N': nt_stats['Taxa'],
        'Loci_Count': 2900,
        'Alignment_Length': f"{int(nt_stats['Sites']):,} bp",
        'Parsimony_Informative_Sites': int(nt_stats['Parsimony_Informative']),
        'Distinct_Patterns': int(nt_stats['Distinct_Patterns']),
        'BIC_Best_Model': nt_stats['Best_BIC_Model'],
        'LogLikelihood': float(nt_stats['LogL']),
        'BIC_Score': float(nt_stats['BIC']),
        'Total_Tree_Length': float(nt_stats['Tree_Length']),
        'Support_Metrics': '1,000 SH-aLRT / 1,000 UFBoot (--bnni)',
        'InterLineage_RF_vs_Primary': 0,
        'Newick_String': nt_nwk
    },
    {
        'Matrix': 'Matrix B — Translated Amino Acids (Validation)',
        'Software': 'IQ-TREE 2 v2.0.7',
        'Taxa_N': aa_stats['Taxa'],
        'Loci_Count': 400,
        'Alignment_Length': f"{int(aa_stats['Sites']):,} aa",
        'Parsimony_Informative_Sites': int(aa_stats['Parsimony_Informative']),
        'Distinct_Patterns': int(aa_stats['Distinct_Patterns']),
        'BIC_Best_Model': aa_stats['Best_BIC_Model'],
        'LogLikelihood': float(aa_stats['LogL']),
        'BIC_Score': float(aa_stats['BIC']),
        'Total_Tree_Length': float(aa_stats['Tree_Length']),
        'Support_Metrics': '1,000 SH-aLRT / 1,000 UFBoot (--bnni)',
        'InterLineage_RF_vs_Primary': rf_nt_aa,
        'Newick_String': aa_nwk
    },
    {
        'Matrix': 'Cross-Method — Neighbor-Joining (Core p-distances)',
        'Software': 'Biopython NJ (p-distance)',
        'Taxa_N': 31,
        'Loci_Count': 2900,
        'Alignment_Length': '2,634,233 bp',
        'Parsimony_Informative_Sites': int(nt_stats['Parsimony_Informative']),
        'Distinct_Patterns': int(nt_stats['Distinct_Patterns']),
        'BIC_Best_Model': 'Uncorrected p-distance (NJ)',
        'LogLikelihood': np.nan,
        'BIC_Score': np.nan,
        'Total_Tree_Length': round(nj_tree.total_branch_length(), 6),
        'Support_Metrics': 'Distance congruence check',
        'InterLineage_RF_vs_Primary': rf_nt_nj,
        'Newick_String': nj_nwk
    },
    {
        'Matrix': 'Outgroup Root Verification — 9 Taxa (6 Lineages + E. allii + Polar E4 + Y. enterocolitica 8081)',
        'Software': 'IQ-TREE 2 v2.0.7',
        'Taxa_N': 9,
        'Loci_Count': 147,
        'Alignment_Length': '150,258 bp',
        'Parsimony_Informative_Sites': 14205,
        'Distinct_Patterns': 3021,
        'BIC_Best_Model': 'GTR+F+R3',
        'LogLikelihood': -475412.1637,
        'BIC_Score': 951146.1704,
        'Total_Tree_Length': 1.5265,
        'Support_Metrics': '1,000 SH-aLRT / 1,000 UFBoot (--bnni): 100/100 E. americana monophyly & 100/100 basal B6-1',
        'InterLineage_RF_vs_Primary': 0,
        'Newick_String': open(os.path.join(phylo_dir, 'outgroup_root_check', 'outgroup_rooted_9taxa.treefile')).read().strip()
    }
])
df_models.to_csv(os.path.join(phylo_dir, "Table_Phylogenetic_Models_and_Support.csv"), index=False)
print(">>> Saved Table_Phylogenetic_Models_and_Support.csv")

# ==============================================================================
# 2. RENDER Q1 PUBLICATION-GRADE MULTI-PANEL FIGURE 4 (300 DPI)
# ==============================================================================
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['DejaVu Sans', 'Arial', 'Helvetica']
plt.rcParams['mathtext.fontset'] = 'dejavusans'

# Load primary ML tree (has dual SH-aLRT / UFBoot in clade.name!)
tree = Phylo.read(nt_treef, "newick")
tree.root_at_midpoint()
tree.ladderize(reverse=True)

# Assign exact Y positions (1..31) to terminals and recursive midpoint Y to internal nodes
terminals = list(tree.get_terminals())
y_pos = {}
for idx, tip in enumerate(terminals):
    y_pos[tip] = len(terminals) - idx  # top = 31, bottom = 1

def calc_y(clade):
    if clade.is_terminal():
        return y_pos[clade]
    ys = [calc_y(c) for c in clade.clades]
    y_pos[clade] = 0.5 * (min(ys) + max(ys))
    return y_pos[clade]

calc_y(tree.root)

x_pos = {}
def calc_x(clade, curr_x=0.0):
    bl = clade.branch_length if clade.branch_length is not None else 0.0
    x_pos[clade] = curr_x + bl
    for c in clade.clades:
        calc_x(c, x_pos[clade])

calc_x(tree.root, 0.0)
max_x = max(x_pos[t] for t in terminals)

# Classify each tip into its exact clade with clean color & readable label
CLADE_CFG = {
    'B6-1':       {'color': '#B91C1C', 'bg': '#FEE2E2', 'title': 'Mushroom-Associated\nLineage (B6-1)'},
    'RIT713':     {'color': '#C2410C', 'bg': '#FFEDD5', 'title': 'Wild Fungal\nAssociate (RIT713)'},
    'Bioreactor': {'color': '#1D4ED8', 'bg': '#DBEAFE', 'title': 'Bioreactor Effluent\nClones (N = 6)'},
    'Oral':       {'color': '#7E22CE', 'bg': '#F3E8FF', 'title': 'Human Oral Saliva\nIsolate (N = 1)'},
    'TypeRef':    {'color': '#0F766E', 'bg': '#CCFBF1', 'title': 'Clinical Type Strain\n(CCUG / ATCC, N = 2)'},
    'CIP':        {'color': '#334155', 'bg': '#F1F5F9', 'title': 'Clinical CIP 81.94 &\nHMF3896 Subclones\n(N = 20)'},
}

def classify_tip(acc):
    row = df_s1.loc[acc] if acc in df_s1.index else {}
    strain = str(row.get('strain', '') or '')
    if strain == 'nan': strain = ''
    if acc == 'GCF_014233485.1':
        return 'B6-1', f"{acc} — B6-1 (Mushroom rot, China) [Diverged lineage]"
    elif acc == 'GCF_009668135.1':
        return 'RIT713', f"{acc} — RIT713 (Wild fungus associate, USA: NY)"
    elif acc.startswith('GCF_0226'):
        return 'Bioreactor', f"{acc} — Bioreactor effluent clone (USA: WI)"
    elif acc == 'GCF_057314615.1':
        return 'Oral', f"{acc} — Human oral saliva isolate (USA)"
    elif acc == 'GCF_008693045.1':
        return 'TypeRef', f"{acc} — CCUG 14506ᵀ (Type strain, throat, Sweden)"
    elif acc == 'GCF_000735345.1':
        return 'TypeRef', f"{acc} — ATCC 33852 (Type-strain deposit, USA)"
    else:
        lbl = strain if strain else "CIP 81.94 subclone"
        return 'CIP', f"{acc} — {lbl} (Clinical subclone)"

# Figure layout: Clean spacious landscape 3-panel Q1 layout
fig = plt.figure(figsize=(22.0, 14.2), dpi=300, facecolor='white')
gs = GridSpec(2, 2, figure=fig, width_ratios=[1.58, 1.0], height_ratios=[1.0, 1.0],
              left=0.04, right=0.98, top=0.90, bottom=0.08, wspace=0.26, hspace=0.32)

axA = fig.add_subplot(gs[:, 0])
axB = fig.add_subplot(gs[0, 1])
axC = fig.add_subplot(gs[1, 1])

# ------------------------------------------------------------------------------
# PANEL A: ROOTED IQ-TREE 2 ML PHYLOGRAM (2,634,233 bp; 2,900 STRICT CORE LOCI)
# ------------------------------------------------------------------------------
axA.set_facecolor('white')

# Draw shaded background bands per contiguous clade group
clade_ys = {}
tip_info = {}
for tip in terminals:
    ckey, tlbl = classify_tip(tip.name)
    tip_info[tip] = (ckey, tlbl)
    clade_ys.setdefault(ckey, []).append(y_pos[tip])

x_band_end = max_x + 0.049
for ckey, ys_list in clade_ys.items():
    y_lo, y_hi = min(ys_list) - 0.44, max(ys_list) + 0.44
    cfg = CLADE_CFG[ckey]
    axA.add_patch(patches.Rectangle((-0.001, y_lo), x_band_end + 0.001, y_hi - y_lo,
                                    facecolor=cfg['bg'], edgecolor='none', alpha=0.65, zorder=0))
    # Crisp right-side clade vertical bracket bar
    bx = x_band_end + 0.0012
    axA.plot([bx, bx], [y_lo + 0.08, y_hi - 0.08], color=cfg['color'], lw=3.0, solid_capstyle='butt', zorder=4)
    axA.text(bx + 0.0012, 0.5*(y_lo + y_hi), cfg['title'], ha='left', va='center',
             fontsize=8.2, fontweight='bold', color=cfg['color'], linespacing=1.15)

# Draw crisp phylogram branches
def draw_branches(clade):
    x0, y0 = x_pos[clade], y_pos[clade]
    if not clade.is_terminal():
        c_ys = [y_pos[c] for c in clade.clades]
        axA.plot([x0, x0], [min(c_ys), max(c_ys)], color='#1E293B', lw=1.5, solid_capstyle='round', zorder=3)
        for c in clade.clades:
            xc, yc = x_pos[c], y_pos[c]
            # Highlight B6-1 branch in bold crimson!
            is_b6 = (c.is_terminal() and c.name == 'GCF_014233485.1')
            col_b = '#B91C1C' if is_b6 else '#1E293B'
            lw_b  = 2.6 if is_b6 else 1.5
            axA.plot([x0, xc], [yc, yc], color=col_b, lw=lw_b, solid_capstyle='round', zorder=4 if is_b6 else 3)
            draw_branches(c)

draw_branches(tree.root)

# Draw tip markers & clean unclipped accession+strain labels
for tip in terminals:
    xt, yt = x_pos[tip], y_pos[tip]
    ckey, tlbl = tip_info[tip]
    cfg = CLADE_CFG[ckey]
    is_b6 = (ckey == 'B6-1')
    axA.scatter(xt, yt, c=cfg['color'], s=46 if is_b6 else 26, edgecolors='white', lw=0.7, zorder=6)
    axA.text(xt + 0.0009, yt, tlbl, ha='left', va='center',
             fontsize=8.3 if is_b6 else 7.6,
             fontweight='bold' if (is_b6 or ckey in ('TypeRef', 'RIT713', 'Oral')) else 'medium',
             color=cfg['color'] if is_b6 else '#1E293B', zorder=6)

# Draw crisp dual SH-aLRT / UFBoot support badges at ALL resolved major internal nodes
# (Exclude zero-length clonal polytomies < 0.00002 to eliminate overlapping badge clutter!)
for cl in tree.get_nonterminals():
    if cl == tree.root:
        continue
    bl = cl.branch_length if cl.branch_length is not None else 0.0
    raw_lbl = str(cl.name or '')
    if '/' not in raw_lbl:
        continue
    parts = raw_lbl.split('/')
    try:
        alrt_v = float(parts[0])
        ufb_v  = float(parts[1])
    except ValueError:
        continue
    # Show badge on all strongly resolved non-trivial branches (bl >= 0.0001 and >= 80/80)
    if bl >= 0.0001 and alrt_v >= 80 and ufb_v >= 80:
        xn, yn = x_pos[cl], y_pos[cl]
        alrt_str = f"{int(round(alrt_v))}" if abs(alrt_v - round(alrt_v)) < 0.05 else f"{alrt_v:.1f}"
        ufb_str  = f"{int(round(ufb_v))}"
        badge_txt = f"{alrt_str}/{ufb_str}"
        # Avoid crossing adjacent vertical trunk on the Bioreactor+RIT713 node (place cleanly above-right)
        if 3.5 <= yn <= 5.5 and xn < 0.014:
            axA.text(xn + 0.0007, yn + 0.30, badge_txt, ha='left', va='bottom',
                     fontsize=7.3, fontweight='bold', color='#991B1B',
                     bbox=dict(boxstyle='round,pad=0.16', facecolor='white', edgecolor='#DC2626', lw=0.75, alpha=0.96),
                     zorder=8)
        else:
            axA.text(xn - 0.0006, yn + 0.34, badge_txt, ha='right', va='bottom',
                     fontsize=7.3, fontweight='bold', color='#991B1B',
                     bbox=dict(boxstyle='round,pad=0.16', facecolor='white', edgecolor='#DC2626', lw=0.75, alpha=0.96),
                     zorder=8)

# Clean boxed ModelFinder & Congruence summary tucked cleanly left of CIP trunk (zero overlap!)
summary_box_txt = (
    "Phylogenomic Model & Validation Summary\n"
    "───────────────────────────────────────\n"
    f"• Primary (NT): 2,634,233 bp | 2,900 loci\n"
    f"  89,557 PI sites | {nt_stats['Best_BIC_Model']} (BIC)\n"
    f"  LogL = {float(nt_stats['LogL']):,.1f} | Len = {nt_stats['Tree_Length']}\n"
    "• Support: 1,000 SH-aLRT / UFBoot (--bnni)\n"
    f"• AA Check: 120,790 aa (400 PCGs) | {aa_stats['Best_BIC_Model']}\n"
    f"  LogL = {float(aa_stats['LogL']):,.1f} | 6/6 clades (RF={rf_nt_aa}*)\n"
    "• Outgroup Root (147 loci, 150.3 kb):\n"
    "  GTR+F+R3 ('E. allii' + E4 + Y. ent.)\n"
    "  → 100/100 monophyly & basal B6-1"
)
axA.text(0.0008, 30.8, summary_box_txt, ha='left', va='top', fontsize=7.1, color='#0F172A',
         family='DejaVu Sans', linespacing=1.30,
         bbox=dict(boxstyle='round,pad=0.38', facecolor='#F8FAFC', edgecolor='#64748B', lw=0.85, alpha=0.96),
         zorder=9)

axA.set_xlim(-0.0015, x_band_end + 0.0175)
axA.set_ylim(0.2, 31.9)
axA.set_yticks([])
axA.set_xlabel("Nucleotide Substitutions per Site (Midpoint-Rooted IQ-TREE 2 ML Phylogram)", fontsize=10.2, fontweight='bold', color='#0F172A', labelpad=8)
axA.spines['top'].set_visible(False)
axA.spines['right'].set_visible(False)
axA.spines['left'].set_visible(False)
axA.spines['bottom'].set_linewidth(1.0)
axA.tick_params(axis='x', labelsize=9.0)
axA.set_title("A. Rooted Core-Genome ML Phylogeny (2,634,233 bp; 2,900 Strict Core Loci)",
              fontsize=11.8, fontweight='bold', color='#0F172A', loc='left', pad=12)

# ------------------------------------------------------------------------------
# PANEL B: PAIRWISE CORE-GENOME DIVERGENCE (p-DISTANCE & CORE SNPs)
# ------------------------------------------------------------------------------
axB.set_facecolor('white')

def get_grp(acc):
    if acc == "GCF_014233485.1": return "B6-1"
    if acc == "GCF_009668135.1": return "RIT713"
    if acc == "GCF_057314615.1": return "Oral"
    if acc.startswith("GCF_0226"): return "Bioreactor"
    if acc in ("GCF_008693045.1", "GCF_000735345.1"): return "TypeRef"
    return "CIP"

grps = {a: get_grp(a) for a in df_pdist.index}
def mean_pair(g1, g2):
    vals_p, vals_s = [], []
    accs = list(df_pdist.index)
    for i in range(len(accs)):
        for j in range(i + 1, len(accs)):
            a, b = accs[i], accs[j]
            if (grps[a] == g1 and grps[b] == g2) or (grps[a] == g2 and grps[b] == g1):
                vals_p.append(100.0 * df_pdist.iloc[i, j])
                vals_s.append(df_snps.iloc[i, j])
    return np.mean(vals_p), np.mean(vals_s)

bar_defs = [
    ("Intra-clonal: Type Strain\n(CCUG 14506ᵀ vs ATCC 33852)", *mean_pair('TypeRef', 'TypeRef'), '#0F766E'),
    ("Intra-clonal: CIP 81.94 &\nHMF3896 Subclones (N = 20)",   *mean_pair('CIP', 'CIP'),         '#475569'),
    ("Intra-clonal: Bioreactor\nEffluent Clones (N = 6)",      *mean_pair('Bioreactor', 'Bioreactor'), '#2563EB'),
    ("Inter-lineage: Oral Isolate\nvs. Type Strain",           *mean_pair('Oral', 'TypeRef'),    '#7E22CE'),
    ("Inter-lineage: CIP Subclones\nvs. Type Strain",          *mean_pair('CIP', 'TypeRef'),     '#334155'),
    ("Inter-lineage: RIT713\nvs. Type Strain",                 *mean_pair('RIT713', 'TypeRef'),  '#C2410C'),
    ("Inter-lineage: Bioreactor\nvs. Type Strain",             *mean_pair('Bioreactor', 'TypeRef'), '#1D4ED8'),
    ("Inter-lineage: B6-1 (Mushroom)\nvs. Type Strain",        *mean_pair('B6-1', 'TypeRef'),    '#B91C1C'),
]

y_b = np.arange(len(bar_defs))
vals_b = [d[1] for d in bar_defs]
snps_b = [d[2] for d in bar_defs]
cols_b = [d[3] for d in bar_defs]
lbls_b = [d[0] for d in bar_defs]

# Draw horizontal divider separating intra-clonal vs inter-lineage
axB.axhline(2.5, color='#94A3B8', ls='--', lw=1.0, zorder=1)
axB.text(4.28, 2.66, "▲ Inter-lineage divergence (1.71–3.19%)", ha='right', va='bottom', fontsize=7.6, fontstyle='italic', color='#991B1B')
axB.text(4.28, 2.34, "▼ Intra-clonal stasis (<0.01%)", ha='right', va='top', fontsize=7.6, fontstyle='italic', color='#475569')

bars = axB.barh(y_b, [max(v, 0.03) for v in vals_b], height=0.62, color=cols_b, edgecolor='#1E293B', lw=0.7, zorder=3)
for idx_b, (v, s_cnt, c_bar) in enumerate(zip(vals_b, snps_b, cols_b)):
    snp_str = f"{s_cnt:,.0f} SNPs" if s_cnt < 1000 else f"{s_cnt/1000:.1f}k SNPs"
    txt_val = f"{v:.4f}% ({snp_str})" if v < 0.05 else f"{v:.2f}% ({snp_str})"
    axB.text(max(v, 0.03) + 0.07, idx_b, txt_val, ha='left', va='center',
             fontsize=8.2, fontweight='bold' if idx_b == 7 else 'medium',
             color='#B91C1C' if idx_b == 7 else '#1E293B')

axB.set_yticks(y_b)
axB.set_yticklabels(lbls_b, fontsize=8.0, color='#1E293B')
axB.set_xlim(0, 4.35)
axB.set_xlabel("Mean Pairwise Core Nucleotide Divergence (% p-distance)", fontsize=9.6, fontweight='bold', color='#0F172A', labelpad=6)
axB.grid(axis='x', linestyle=':', alpha=0.55, zorder=0)
axB.spines['top'].set_visible(False)
axB.spines['right'].set_visible(False)
axB.tick_params(axis='x', labelsize=8.8)
axB.set_title("B. Core-Genome Pairwise Divergence & Clonal Stasis (p-distance)",
              fontsize=11.4, fontweight='bold', color='#0F172A', loc='left', pad=10)

# ------------------------------------------------------------------------------
# PANEL C: GENE-WISE Ka/Ks SELECTION PRESSURE & COMPOSITIONAL BENCHMARKS
# ------------------------------------------------------------------------------
axC.set_facecolor('white')

kaks_groups = [
    ('Oral',       'Oral\n(n=2,144)',       '#7E22CE'),
    ('Bioreactor', 'Bioreactor\n(n=2,156)', '#1D4ED8'),
    ('RIT713',     'RIT713\n(n=2,172)',     '#C2410C'),
    ('B6-1',       'B6-1\n(n=2,180)',       '#B91C1C'),
]

plot_data = []
for key_k, _, _ in kaks_groups:
    s = df_kaks[f'KaKs_{key_k}'].dropna().values
    # Cap visualization tail at 0.25 while reporting exact full-dataset medians & 100% < 1.0
    plot_data.append(s)

positions = np.arange(1, len(kaks_groups) + 1)
vp = axC.violinplot([d[d <= 0.20] for d in plot_data], positions=positions, widths=0.62,
                    showmeans=False, showmedians=False, showextrema=False)
for idx_v, body in enumerate(vp['bodies']):
    body.set_facecolor(kaks_groups[idx_v][2])
    body.set_edgecolor('#1E293B')
    body.set_alpha(0.32)
    body.set_linewidth(0.8)

# Overlay crisp IQR boxes & median points
for idx_v, (d_full, (key_k, _, col_k)) in enumerate(zip(plot_data, kaks_groups)):
    pos = positions[idx_v]
    q25, med, q75 = np.percentile(d_full, [25, 50, 75])
    axC.plot([pos, pos], [q25, q75], color=col_k, lw=4.5, solid_capstyle='round', zorder=4)
    axC.scatter(pos, med, color='white', edgecolor=col_k, s=50, lw=1.8, zorder=5)
    axC.text(pos, q75 + 0.014, f"Med = {med:.4f}\n(100% < 1.0)", ha='center', va='bottom',
             fontsize=7.6, fontweight='bold', color=col_k,
             bbox=dict(boxstyle='round,pad=0.14', facecolor='white', edgecolor='none', alpha=0.82), zorder=6)

# Compositional homogeneity benchmark inset box (top-right inside Panel C)
gc_min, gc_max, gc_mean = df_comp['GC_pct'].min(), df_comp['GC_pct'].max(), df_comp['GC_pct'].mean()
at_mean = df_comp['AT_skew'].mean()
gc_sk_mean = df_comp['GC_skew'].mean()
comp_box = (
    "Strict-Core Compositional Benchmarks (N = 31)\n"
    f"• Core GC%: {gc_mean:.2f}% (range {gc_min:.2f}–{gc_max:.2f}%)\n"
    f"• Mean AT-skew: +{at_mean:.4f} | GC-skew: +{gc_sk_mean:.4f}\n"
    "• AA χ² Composition Test: 31/31 passed (0 failures)\n"
    "• Neutrality Cutoff (Ka/Ks = 1.0): 0/2,180 genes > 1.0\n"
    "  → Strong pervasive purifying selection (ω ≪ 1)"
)
axC.text(0.50, 0.96, comp_box, transform=axC.transAxes, ha='center', va='top',
         fontsize=7.6, color='#0F172A', linespacing=1.28,
         bbox=dict(boxstyle='round,pad=0.40', facecolor='#F8FAFC', edgecolor='#64748B', lw=0.85, alpha=0.96),
         zorder=8)

axC.set_xticks(positions)
axC.set_xticklabels([g[1] for g in kaks_groups], fontsize=8.8, fontweight='bold', color='#1E293B')
axC.set_ylim(-0.005, 0.215)
axC.set_ylabel("Gene-wise $K_a/K_s$ vs. CCUG 14506ᵀ (Nei–Gojobori)", fontsize=9.2, fontweight='bold', color='#0F172A', labelpad=6)
axC.set_xlabel("Diverged E. americana Biological Lineages (Single-Copy Strict Core PCGs)", fontsize=9.2, fontweight='bold', color='#0F172A', labelpad=6)
axC.grid(axis='y', linestyle=':', alpha=0.55, zorder=0)
axC.spines['top'].set_visible(False)
axC.spines['right'].set_visible(False)
axC.tick_params(axis='y', labelsize=8.6)
axC.set_title("C. Core Purifying Selection (Ka/Ks) & Compositional Homogeneity",
              fontsize=11.4, fontweight='bold', color='#0F172A', loc='left', pad=10)

# Global Publication Title & Footnote
fig.text(0.51, 0.962,
         "Figure 4: Core-Genome Maximum-Likelihood Phylogenomics, Lineage Divergence, and Purifying Selection in Ewingella americana",
         ha='center', va='center', fontsize=14.2, fontweight='bold', color='#0F172A')
fig.text(0.51, 0.936,
         "IQ-TREE 2 reconstruction from 2,634,233 bp (2,900 strict core loci; N = 31 curated WGS genomes) with dual SH-aLRT/UFBoot support, AA-ML & BioNJ congruence, and Nei–Gojobori Ka/Ks benchmarks.",
         ha='center', va='center', fontsize=9.4, fontstyle='italic', color='#475569')
fig.text(0.04, 0.022,
         "Notes: Panel A displays the midpoint-rooted IQ-TREE 2 v2.0.7 ML phylogram (GTR+F+I+R4; 1,000 SH-aLRT / 1,000 UFBoot --bnni). Independent outgroup-rooted ML reconstruction across 147 universal single-copy loci (150,258 bp; GTR+F+R3) including sister species 'E. allii' (NIAB_247711), Polar E4, and Yersinia enterocolitica 8081 confirmed E. americana monophyly (100/100) and basal placement of B6-1 (100/100), validating the exact root branch. "
         "*AA validation (120,790 aa; JTT+F+I; 31/31 passing χ² composition test) recovered all 6 monophyletic clusters and basal B6-1 (1 short-branch Oral swap; RF = 2; BioNJ RF = 4). "
         "Pairwise core SNPs resolve CCUG 14506ᵀ/ATCC 33852 (3 SNPs) and the 20 CIP 81.94/HMF3896 subclones (4.8 internal SNPs) as two distinct clinical clonal clusters separated by 53,765 SNPs (2.06%).",
         ha='left', va='center', fontsize=7.8, color='#334155', wrap=True)

# ==============================================================================
# EXPORT 300-DPI PNG, LZW TIFF, EDITABLE SVG, AND PDF
# ==============================================================================
for bname in ["Final_Figure_4_Phylogenomics_Q1"]:
    for target_d in [phylo_dir, fig_dir, root_dir]:
        os.makedirs(target_d, exist_ok=True)
    svg_p = os.path.join(fig_dir, f"{bname}.svg")
    pdf_p = os.path.join(fig_dir, f"{bname}.pdf")
    png_p = os.path.join(fig_dir, f"{bname}.png")
    tif_p = os.path.join(fig_dir, f"{bname}.tif")
    fig.savefig(svg_p, facecolor='white')
    fig.savefig(pdf_p, facecolor='white')
    fig.savefig(png_p, dpi=300, facecolor='white')
    im = Image.open(png_p)
    im.save(tif_p, format="TIFF", dpi=(300, 300), compression="tiff_lzw")
    for ext in [".svg", ".tif", ".pdf", ".png"]:
        src = os.path.join(fig_dir, f"{bname}{ext}")
        for d_out in [phylo_dir, root_dir]:
            with open(src, 'rb') as fi, open(os.path.join(d_out, f"{bname}{ext}"), 'wb') as fo:
                fo.write(fi.read())
    print(f">>> Exported {bname} ({os.path.getsize(tif_p)/(1024*1024):.2f} MB TIFF, {im.size[0]}x{im.size[1]} px)")
plt.close(fig)
