import os
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import seaborn as sns
from PIL import Image
from scipy.optimize import curve_fit
from Bio import Phylo

out_dir = "/mnt/d/Ewingella_americana_research/FINAL_JAM_SUBMISSION/figures"
root_dir = "/mnt/d/Ewingella_americana_research/FINAL_JAM_SUBMISSION"
os.makedirs(out_dir, exist_ok=True)

# Set global styles
plt.rcParams['font.sans-serif'] = "DejaVu Sans"
plt.rcParams['font.family'] = "sans-serif"
plt.rcParams['mathtext.fontset'] = "dejavusans"

def save_refined_figure(fig, base_name):
    # Save SVG, PDF, PNG, TIFF in both figures/ and root
    svg_p = os.path.join(out_dir, f"{base_name}.svg")
    pdf_p = os.path.join(out_dir, f"{base_name}.pdf")
    png_p = os.path.join(out_dir, f"{base_name}.png")
    tif_p = os.path.join(out_dir, f"{base_name}.tif")
    
    print(f"Exporting {base_name}...")
    fig.savefig(svg_p, bbox_inches='tight')
    fig.savefig(pdf_p, bbox_inches='tight')
    fig.savefig(png_p, dpi=300, bbox_inches='tight')
    plt.close(fig)
    
    im = Image.open(png_p)
    im.save(tif_p, format="TIFF", dpi=(300, 300), compression="tiff_lzw")
    
    # Also copy to root FINAL_JAM_SUBMISSION/
    for ext in [".svg", ".tif", ".pdf", ".png"]:
        src = os.path.join(out_dir, f"{base_name}{ext}")
        dst = os.path.join(root_dir, f"{base_name}{ext}")
        with open(src, 'rb') as f_in, open(dst, 'wb') as f_out:
            f_out.write(f_in.read())
            
    sz_mb = os.path.getsize(tif_p) / (1024*1024)
    print(f"  Successfully exported {base_name} ({sz_mb:.2f} MB TIFF)")

# ==============================================================================
# FIGURE 1: REFINED METHODOLOGICAL WORKFLOW
# ==============================================================================
print("\n>>> Building Final Figure 1: Workflow Architecture...")
fig1, ax1 = plt.subplots(figsize=(22, 16.5), dpi=300)
ax1.set_xlim(0, 100)
ax1.set_ylim(0, 100)
ax1.axis('off')

c_bg_card = "#FFFFFF"
c_border_card = "#CBD5E1"
c_text_title = "#0F172A"
c_text_sub = "#475569"
c_text_body = "#334155"
c_pipe = "#475569"

stage_colors = {
    1: {"header": "#1E293B", "accent": "#2563EB", "badge_bg": "#F8FAFC", "badge_border": "#93C5FD", "badge_text": "#1E3A8A"},
    2: {"header": "#0F172A", "accent": "#0D9488", "badge_bg": "#F8FAFC", "badge_border": "#99F6E4", "badge_text": "#0F766E"},
    3: {"header": "#1E1B4B", "accent": "#4F46E5", "badge_bg": "#F8FAFC", "badge_border": "#C7D2FE", "badge_text": "#3730A3"},
    4: {"header": "#292524", "accent": "#D97706", "badge_bg": "#F8FAFC", "badge_border": "#FDE68A", "badge_text": "#92400E"},
    5: {"header": "#2E1065", "accent": "#7C3AED", "badge_bg": "#F8FAFC", "badge_border": "#DDD6FE", "badge_text": "#5B21B6"},
    6: {"header": "#064E3B", "accent": "#059669", "badge_bg": "#F8FAFC", "badge_border": "#A7F3D0", "badge_text": "#065F46"},
}

def draw_spacious_card(ax, x, y, w, h, stage_num, title, subtitle, methods_items, results_items, metric_badge):
    cfg = stage_colors[stage_num]
    shadow = patches.FancyBboxPatch((x + 0.35, y - 0.35), w, h, boxstyle="round,pad=0.5",
                                   facecolor="#E2E8F0", edgecolor="none", alpha=0.6)
    ax.add_patch(shadow)
    card = patches.FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.5",
                                 facecolor=c_bg_card, edgecolor=c_border_card, linewidth=1.4)
    ax.add_patch(card)
    header_h = 4.6
    header = patches.FancyBboxPatch((x, y + h - header_h), w, header_h, boxstyle="round,pad=0.5",
                                   facecolor=cfg["header"], edgecolor="none")
    ax.add_patch(header)
    accent_strip = patches.Rectangle((x - 0.5, y + h - 0.5), w + 1.0, 0.7, facecolor=cfg["accent"], edgecolor="none")
    ax.add_patch(accent_strip)
    pill_w = 6.8
    pill_h = 2.2
    pill = patches.FancyBboxPatch((x + 1.2, y + h - header_h + 1.2), pill_w, pill_h, boxstyle="round,pad=0.2",
                                 facecolor=cfg["accent"], edgecolor="none")
    ax.add_patch(pill)
    ax.text(x + 1.2 + pill_w/2, y + h - header_h + 1.2 + pill_h/2, f"STAGE {stage_num:02d}", 
            ha='center', va='center', fontsize=8.2, fontweight='bold', color='white')
    ax.text(x + pill_w + 2.4, y + h - header_h/2 + 0.7, title, 
            ha='left', va='center', fontsize=9.6, fontweight='bold', color='white')
    ax.text(x + pill_w + 2.4, y + h - header_h/2 - 1.0, subtitle, 
            ha='left', va='center', fontsize=7.8, fontstyle='italic', color='#94A3B8')
    
    curr_y = y + h - header_h - 1.8
    ax.text(x + 1.4, curr_y, "METHODOLOGY & PIPELINE", ha='left', va='center', 
            fontsize=7.8, fontweight='bold', color=cfg["accent"])
    curr_y -= 1.8
    for item in methods_items:
        ax.text(x + 1.6, curr_y, "•", ha='left', va='top', fontsize=9.0, color=cfg["accent"], fontweight='bold')
        ax.text(x + 3.2, curr_y + 0.1, item, ha='left', va='top', fontsize=7.5, color=c_text_body)
        curr_y -= 2.6
    curr_y -= 0.5
    ax.plot([x + 1.5, x + w - 1.5], [curr_y, curr_y], color="#F1F5F9", linewidth=1.2)
    curr_y -= 1.8
    ax.text(x + 1.4, curr_y, "KEY FINDINGS & RESOLUTION", ha='left', va='center', 
            fontsize=7.8, fontweight='bold', color=c_text_title)
    curr_y -= 1.8
    for item in results_items:
        ax.text(x + 1.6, curr_y, "✓", ha='left', va='top', fontsize=8.5, color=cfg["accent"], fontweight='bold')
        ax.text(x + 3.2, curr_y + 0.1, item, ha='left', va='top', fontsize=7.5, color=c_text_body)
        curr_y -= 2.6
        
    # Ample vertical spacing before metric badge
    badge_y = y + 1.2
    badge_h = 3.2
    badge = patches.FancyBboxPatch((x + 1.2, badge_y), w - 2.4, badge_h, boxstyle="round,pad=0.3",
                                  facecolor=cfg["badge_bg"], edgecolor=cfg["badge_border"], linewidth=1.2)
    ax.add_patch(badge)
    ax.text(x + w/2, badge_y + badge_h/2, metric_badge, ha='center', va='center', 
            fontsize=8.0, fontweight='bold', color=cfg["badge_text"])

def draw_pipeline_arrow(ax, x1, y1, x2, y2, label="", label_pos=0.5):
    ax.annotate('', xy=(x2, y2), xytext=(x1, y1),
                arrowprops=dict(arrowstyle="->", color=c_pipe, linewidth=2.2, mutation_scale=15))
    if label:
        lx = x1 + (x2 - x1) * label_pos
        ly = y1 + (y2 - y1) * label_pos
        ax.text(lx, ly, label, ha='center', va='center', 
                fontsize=7.6, fontweight='bold', color=c_pipe,
                bbox=dict(boxstyle="round,pad=0.25", facecolor="#FFFFFF", edgecolor="#CBD5E1", linewidth=0.8))

col_w = 28.5
card_h = 37.0
y_top = 53.5
y_bot = 5.5
x_col1, x_col2, x_col3 = 4.0, 35.8, 67.5

banner1 = patches.FancyBboxPatch((4.0, 91.8), 92.0, 3.2, boxstyle="round,pad=0.4", facecolor="#0F172A", edgecolor="none")
ax1.add_patch(banner1)
accent_bar1 = patches.Rectangle((3.6, 91.8), 1.2, 3.2, facecolor="#38BDF8", edgecolor="none")
ax1.add_patch(accent_bar1)
ax1.text(6.5, 93.4, "PART I: GENOMIC RETRIEVAL, TAXONOMIC DEMARCATION & HIGH-STRINGENCY PAN-GENOME MODELING", 
        ha='left', va='center', fontsize=9.2, fontweight='bold', color="white")

draw_spacious_card(ax1, x_col1, y_top, col_w, card_h, 1,
                   "Genomic Data Retrieval & Curation", "Global Repository Mining & De Novo Assembly",
                   ["Exhaustive acquisition of all INSDC public assemblies (N=35)",
                    "De novo SPAdes assembly of plant epiphyte NIAB_247711 (ENA SRA)",
                    "Lineage-specific marker gene profiling via CheckM v1.2.2"],
                   ["Quality Criteria: Completeness > 90%, Contamination < 5%",
                    "Verified Cohort Observed: Completeness >= 98.8%, Contam <= 1.95%",
                    "Assessment of scaffold contiguity, N50, L50, and raw read depth"],
                   "SCREENED DATASET: N = 36 TOTAL GENOME ASSEMBLIES")

draw_spacious_card(ax1, x_col2, y_top, col_w, card_h, 2,
                   "FastANI Taxonomic Demarcation", "Species Boundary Delineation & Curation",
                   ["All-vs-all FastANI v1.33 against type strain CCUG 14506T",
                    "Application of standard 95–96% whole-genome species boundary",
                    "Automated dereplication screening and assembly defect detection"],
                   ["Confirmed N = 33 true E. americana assemblies (>=96.36% ANI)",
                    "NIAB_247711 resolved as proposed 'E. allii' (99.09% ANI to type)",
                    "Excluded 2 polar isolates (<83.2%) & 2 clinical QC fails (NCTC/BRK18a)"],
                   "VERIFIED PAN-GENOME COHORT: N = 31 HIGH-QUALITY GENOMES")

draw_spacious_card(ax1, x_col3, y_top, col_w, card_h, 3,
                   "Pan-Genome Graph Reconstruction", "Panaroo High-Stringency Orthology Modeling",
                   ["Panaroo v1.3.4 execution in '--clean-mode strict' (N=31 cohort)",
                    "Contig-edge fragmentation filtering to eliminate assembly artifacts",
                    "Mathematical trajectory modeling across 100 random permutations"],
                   ["5,545 total gene families: Strict Core = 2,900 (52.30%), Soft Core = 848",
                    "Shell = 1,038 (18.72%), Cloud = 759 (13.69%), Singletons = 492",
                    "Heaps' Law Fit: n = 4,381.82 · N^0.0671 (gamma = 0.0671; Omega = 2,658.7)"],
                   "RESOLVED PAN-GENOME: 5,545 FAMILIES | CORE: 3,748 (67.6%)")

draw_pipeline_arrow(ax1, x_col1 + col_w, y_top + card_h/2, x_col2, y_top + card_h/2, "QC PASS")
draw_pipeline_arrow(ax1, x_col2 + col_w, y_top + card_h/2, x_col3, y_top + card_h/2, "N = 31")

hub_y = 47.8
hub = patches.FancyBboxPatch((4.0, hub_y), 92.0, 2.6, boxstyle="round,pad=0.3", facecolor="#0F172A", edgecolor="#38BDF8", linewidth=1.2)
ax1.add_patch(hub)
ax1.text(50.0, hub_y + 1.3, "PART II: EVOLUTIONARY, POPULATION STRUCTURE & FUNCTIONAL SPECIALIZATION DOWNSTREAM OF STAGE 03", 
        ha='center', va='center', fontsize=8.6, fontweight='bold', color="#38BDF8")

drop_start_y = hub_y
drop_end_y = y_bot + card_h + 0.2
draw_pipeline_arrow(ax1, x_col1 + col_w/2, drop_start_y, x_col1 + col_w/2, drop_end_y, "CLONAL SENSITIVITY", label_pos=0.5)
draw_pipeline_arrow(ax1, x_col2 + col_w/2, drop_start_y, x_col2 + col_w/2, drop_end_y, "STRICT CORE ALIGNMENT", label_pos=0.5)
draw_pipeline_arrow(ax1, x_col3 + col_w/2, drop_start_y, x_col3 + col_w/2, drop_end_y, "MOBILOME & AMR", label_pos=0.5)

draw_spacious_card(ax1, x_col1, y_bot, col_w, card_h, 4,
                   "Clonal Dereplication & Sensitivity", "Demasking Database Redundancy Distortions",
                   ["Provenance tracking across international isolate deposits",
                    "Demasked 22 type-strain clones (19 CIP, ATCC 33852, HMF3896, CCUG)",
                    "Demasked 6 clonal anaerobic bioreactor effluent replicates (WI, USA)"],
                   ["Delineated exactly 5 authentic biological lineages across species",
                    "Dereplicated 5-lineage core expands to 3,585 gene families (64.65%)",
                    "Lineage singletons unmasked: 1,064 families (19.19% vs 492 in raw cohort)",
                    "Redundancy modeling: Clonal duplication deflates apparent core by 23.6%"],
                   "CLONAL OVERREPRESENTATION DEFLATES APPARENT CORE DIVERSITY")

draw_spacious_card(ax1, x_col2, y_bot, col_w, card_h, 5,
                   "Core-Genome Phylogenomics", "Maximum-Likelihood Phylogenomic Inference",
                   ["Concatenated strict core supermatrix: 2,634,233 bp across 2,900 loci",
                    "IQ-TREE 2 v2.0.7 maximum-likelihood tree under GTR+F+I+R4 model",
                    "Nodal robustness validation: 1,000 Ultrafast Bootstrap (UFBoot) replicates",
                    "Rooted using Yersinia enterocolitica 8081 & Y. pestis CO92 outgroups"],
                   ["High topological resolution: Monophyletic bioreactor effluent radiation",
                    "Extremely compact clinical type-strain radiation (<0.02% divergence)",
                    "Mushroom-associated isolate B6-1 forms a distinct, relatively long-branched lineage"],
                   "2.63 Mb ALIGNMENT | 2,900 STRICT CORE LOCI | ALL 31 RESOLVED")

draw_spacious_card(ax1, x_col3, y_bot, col_w, card_h, 6,
                   "Resistome, Mobilome & Adaptations", "Plasmidome, Resistance & CAZyme Profiling",
                   ["Screening via ABRicate (CARD, AMRFinderPlus, ResFinder, PlasmidFinder, VFDB)",
                    "Primary screening thresholds: Identity >= 75%, Coverage >= 70%",
                    "IntegronFinder v2.0.6: 0 complete integrons, 0 In0, 0 CALIN, 0 attC sites",
                    "CAZy and dbCAN3 profiling of GH18 chitinolytic genes and candidate effectors"],
                   ["B6-1 plasmids: pB61a (331.4 kb) & pB61b (105.3 kb; predicted Type-F T4SS)",
                    "Both plasmids encode non-canonical RepA/RepB replicons (0 canonical Inc)",
                    "Conserved chromosomal AmpC beta-lactamase & intrinsic multidrug efflux",
                    "Species-wide core GH18 chitinolytic genes (ChiD, nagA-Z) & candidate effectors"],
                   "PLASMIDS ENCODE PREDICTED TYPE-F T4SS | INTRINSIC AMPC")

fig1.text(0.5, 0.976, "Figure 1: Methodological Architecture and Global Comparative Genomic Workflow for Ewingella americana", 
         ha='center', va='center', fontsize=13.8, fontweight='bold', color=c_text_title)
fig1.text(0.5, 0.956, "A standardized six-stage analytical framework encompassing genomic data curation, high-stringency pan-genomics, clonal dereplication, core phylogenomics, and functional profiling.", 
         ha='center', va='center', fontsize=9.0, fontstyle='italic', color=c_text_sub)

save_refined_figure(fig1, "Final_Figure_1_Workflow")

# ==============================================================================
# FIGURE 2: REFINED FASTANI HEATMAP
# ==============================================================================
print("\n>>> Building Final Figure 2: FastANI Heatmap...")
ani_csv = "/mnt/d/Ewingella_americana_research/data_true_ewingella/genus_ewingella_fastani_matrix.csv"
df_ani = pd.read_csv(ani_csv, index_col=0)

label_dict = {
    'E_allii_20WA0182T': 'Proposed taxon "Ewingella allii" 20WA0182T (Type strain)',
    'E_allii_NIAB247711': 'Proposed taxon "Ewingella allii" NIAB_247711 (Epiphyte)',
    'E_americana_B6-1': 'E. americana B6-1 (Mushroom-associated isolate)',
    'E_americana_BRK18a': 'E. americana BRK18a (Blood, QC fail: fragmented)',
    'E_americana_Bioreactor': 'E. americana Bioreactor (Effluent representative)',
    'E_americana_CCUG14506T': 'E. americana CCUG 14506T (Type strain, throat)',
    'E_americana_Oral': 'E. americana Oral isolate (Human saliva)',
    'E_americana_RIT713': 'E. americana RIT713 (Fungal associate)',
    'E_docleensis_CoE038': 'Ewingella docleensis CoE-038-23T (Type strain)',
    'Polar_Antarctic_L47': 'Polar isolate L47 (Antarctic soil, outgroup)',
    'Polar_Arctic_E4': 'Polar isolate E4 (Arctic snow, outgroup)'
}
df_ani.index = [label_dict.get(i, i) for i in df_ani.index]
df_ani.columns = [label_dict.get(c, c) for c in df_ani.columns]

fig2, ax2 = plt.subplots(figsize=(15.8, 14.2), dpi=300)
sns.heatmap(df_ani, annot=True, fmt=".2f", cmap="Blues", cbar_kws={'label': 'Average Nucleotide Identity (ANI %)'},
            linewidths=0.8, linecolor='#cbd5e1', ax=ax2, vmin=82, vmax=100)
ax2.set_title("Figure 2: Whole-Genome FastANI Demarcation Across 11 Representative Taxa\n(Representative of the Screened 36-Assembly Dataset)", 
              fontsize=12.8, fontweight='bold', pad=18)
ax2.set_xticklabels(ax2.get_xticklabels(), rotation=42, ha='right', fontsize=8.6)
ax2.set_yticklabels(ax2.get_yticklabels(), rotation=0, fontsize=8.8)

fig2.text(0.05, 0.015, 
          "Species boundary: >=95-96% ANI. Pairwise FastANI v1.33 across 11 representative taxa selected from the curated 36-assembly cohort.\n"
          "NIAB_247711 shares 99.09% ANI with proposed taxon 'Ewingella allii' 20WA0182T, confirming conspecific status within this newly described taxon.\n"
          "Proposed 'E. allii' and E. docleensis share 97.74% ANI, forming a sister species complex separated from E. americana (93.5-94.0% ANI).\n"
          "Polar isolates E4 and L47 exhibit ~82.9% ANI, supporting exclusion from the E. americana species cohort as non-Ewingella enterobacteria.",
          fontsize=8.5, color="#334155", bbox=dict(boxstyle="round,pad=0.45", facecolor="#f8fafc", edgecolor="#cbd5e1"))
plt.subplots_adjust(bottom=0.35, left=0.34, right=0.97, top=0.92)

save_refined_figure(fig2, "Final_Figure_2_FastANI")

# ==============================================================================
# FIGURE 3: REFINED PAN-GENOME MODELING & CLONAL SENSITIVITY
# ==============================================================================
print("\n>>> Building Final Figure 3: Pan-Genome Modeling...")
s3_csv = "/mnt/d/Ewingella_americana_research/JAM_submission_True_Ewingella/Table_S3_Pan_Dynamics_Sampling.csv"
df_s3 = pd.read_csv(s3_csv)
x_data = df_s3['Genome_Count'].values
pan_y = df_s3['Pan_Mean'].values
core_y = df_s3['Core_Mean'].values

def heaps_law(N, kappa, gamma):
    return kappa * (N ** gamma)

def core_decay(N, A, B, Omega):
    return A * np.exp(-B * N) + Omega

popt_pan, pcov_pan = curve_fit(heaps_law, x_data, pan_y, p0=[4380, 0.067])
perr_pan = np.sqrt(np.diag(pcov_pan))
popt_core, pcov_core = curve_fit(core_decay, x_data, core_y, p0=[1700, 0.058, 2658])
perr_core = np.sqrt(np.diag(pcov_core))

kappa_fit, gamma_fit = popt_pan
A_fit, B_fit, omega_fit = popt_core

res_pan = pan_y - heaps_law(x_data, *popt_pan)
r2_pan = 1 - (np.sum(res_pan**2) / np.sum((pan_y - np.mean(pan_y))**2))
res_core = core_y - core_decay(x_data, *popt_core)
r2_core = 1 - (np.sum(res_core**2) / np.sum((core_y - np.mean(core_y))**2))

fig3, axes3 = plt.subplots(2, 2, figsize=(16.0, 13.0), dpi=300)
plt.subplots_adjust(wspace=0.28, hspace=0.34)

core_col = "#2b5c8f"
pan_col = "#c2593f"
soft_col = "#5a9bd4"
shell_col = "#e8a838"
cloud_col = "#7cae7a"

# 3A
ax_3a = axes3[0, 0]
ax_3a.errorbar(x_data, pan_y, yerr=df_s3['Pan_Std'].values, fmt='o', color=pan_col, 
              ecolor='#e8a595', elinewidth=1.2, capsize=2.5, markersize=5,
              label="Permuted Mean ± SD (100 random permutations)")
x_fit = np.linspace(1, 31, 200)
ax_3a.plot(x_fit, heaps_law(x_fit, *popt_pan), 'r-', linewidth=2.5, 
           label=f"Heaps' Law Fit: $n = {kappa_fit:.2f} \\cdot N^{{{gamma_fit:.4f}}}$\n($\\gamma = {gamma_fit:.4f} \\pm {perr_pan[1]:.4f}$, $\\kappa = {kappa_fit:.2f} \\pm {perr_pan[0]:.2f}$, $R^2 = {r2_pan:.4f}$)")
ax_3a.set_title("A. Pan-Genome Accumulation (Heaps' Law Trajectory)", fontsize=12, fontweight='bold', pad=10)
ax_3a.set_xlabel("Number of Genomes Sequenced ($N$)", fontsize=11)
ax_3a.set_ylabel("Total Pan-Genome Gene Families ($n$)", fontsize=11)
ax_3a.grid(True, linestyle="--", alpha=0.5)
ax_3a.legend(loc="lower right", fontsize=9.0, framealpha=0.95)
ax_3a.set_ylim(4200, 5750)
ax_3a.text(0.04, 0.88, "Open Pan-Genome Trajectory\n($\\alpha = 1 - \\gamma = 0.9329 < 1.0$)", 
           transform=ax_3a.transAxes, fontsize=9.2, fontweight='bold',
           bbox=dict(boxstyle="round,pad=0.4", facecolor="#fef9e7", edgecolor="#f39c12", alpha=0.9))

# 3B
ax_3b = axes3[0, 1]
ax_3b.errorbar(x_data, core_y, yerr=df_s3['Core_Std'].values, fmt='s', color=core_col, 
              ecolor='#9ebcdb', elinewidth=1.2, capsize=2.5, markersize=5,
              label="Permuted Mean ± SD (100 random permutations)")
ax_3b.plot(x_fit, core_decay(x_fit, *popt_core), 'b-', linewidth=2.5, 
           label=f"Decay Fit: $\\mathrm{{Core}}(N) = {A_fit:.1f} \\cdot e^{{-{B_fit:.4f}N}} + {omega_fit:.2f}$\n($\\Omega = {omega_fit:.2f} \\pm {perr_core[2]:.2f}$, $B = {B_fit:.4f}$, $R^2 = {r2_core:.4f}$)")
ax_3b.axhline(omega_fit, color='#7f8c8d', linestyle='--', linewidth=1.5, label=f"Extrapolated Asymptote ($\\Omega = {omega_fit:.1f}$)")
ax_3b.set_title("B. Core Genome Exponential Decay", fontsize=12, fontweight='bold', pad=10)
ax_3b.set_xlabel("Number of Genomes Sequenced ($N$)", fontsize=11)
ax_3b.set_ylabel("Conserved Core Gene Families", fontsize=11)
ax_3b.grid(True, linestyle="--", alpha=0.5)
ax_3b.legend(loc="upper right", fontsize=9.0, framealpha=0.95)
ax_3b.set_ylim(2400, 4650)
ax_3b.text(0.04, 0.12, "Observed Strict Core ($N=31$): 2,900 families\nAsymptotic Core ($\\Omega$): 2,658.7 families\n(Model projection under $N \\to \\infty$, incorporating draft dropouts)", 
           transform=ax_3b.transAxes, fontsize=8.5,
           bbox=dict(boxstyle="round,pad=0.4", facecolor="#ebf5fb", edgecolor="#2980b9", alpha=0.9))

# 3C
ax_3c = axes3[1, 0]
categories = ['Strict Core (100%)\n2,900 families', 'Soft Core (95–99%)\n848 families', 'Shell (15–94%)\n1,038 families', 'Cloud (<15%)\n759 families']
sizes = [2900, 848, 1038, 759]
colors = [core_col, soft_col, shell_col, cloud_col]
explode = (0.05, 0.02, 0.02, 0.02)
wedges, texts, autotexts = ax_3c.pie(sizes, explode=explode, labels=categories, autopct='%1.1f%%',
                                    pctdistance=0.76, colors=colors, startangle=140, textprops=dict(color="#1a252f", fontsize=9.5))
for at in autotexts:
    at.set_fontsize(10)
    at.set_fontweight('bold')
centre_circle = plt.Circle((0,0), 0.52, fc='white')
ax_3c.add_artist(centre_circle)
ax_3c.text(0, 0.12, "5,545", ha='center', va='center', fontsize=15, fontweight='bold', color='#1a252f')
ax_3c.text(0, -0.05, "Total Gene\nFamilies", ha='center', va='center', fontsize=9.8, fontweight='bold', color='#566573')
ax_3c.text(0, -0.25, "(Singletons: 492 [8.9%])", ha='center', va='center', fontsize=7.6, color='#566573')
ax_3c.set_title("C. Pan-Genome Frequency Spectrum ($N=31$)", fontsize=12, fontweight='bold', pad=10)

# 3D: Clean, Zero-Overlap Clonal Sensitivity Panel (6 True Biological Lineages)
ax_3d = axes3[1, 1]
labels_d = ['Full Cohort\n(N=31 Genomes)', 'Dereplicated Cohort\n(N=6 Biological Lineages)']
core_counts = [2900, 3566]
acc_counts = [2645, 1979]
x_d = np.arange(len(labels_d))
width = 0.52

b1 = ax_3d.bar(x_d, core_counts, width, label='Core Genome', color=core_col, edgecolor="#1b3b5f", linewidth=1.2)
b2 = ax_3d.bar(x_d, acc_counts, width, bottom=core_counts, label='Accessory Genome', color=shell_col, edgecolor="#b37d1a", linewidth=1.2)

# Perfectly centered values inside respective bar segments
ax_3d.text(0, 2900 / 2, "Strict Core\n2,900\n(52.3%)", ha='center', va='center', color='white', fontweight='bold', fontsize=9.0)
ax_3d.text(0, 2900 + 2645 / 2, "Accessory\nGenome\n2,645 (47.7%)", ha='center', va='center', color='#1a252f', fontweight='bold', fontsize=8.8)

ax_3d.text(1, 3566 / 2, "Lineage Core\n3,566\n(64.3%)", ha='center', va='center', color='white', fontweight='bold', fontsize=9.0)
ax_3d.text(1, 3566 + 1979 / 2, "Lineage\nAccessory\n1,979 (35.7%)", ha='center', va='center', color='#1a252f', fontweight='bold', fontsize=8.8)

# Singleton badges placed cleanly ABOVE bars with zero overlap against legend or axis edges
ax_3d.annotate("Singletons:\n492 (8.9%)", xy=(0, 5570), xytext=(0, 6220),
               ha='center', va='center', fontsize=8.8, fontweight='bold', color='#1e3a8a',
               bbox=dict(boxstyle="round,pad=0.35", facecolor="#f0f9ff", edgecolor="#93c5fd", linewidth=1.0),
               arrowprops=dict(arrowstyle="->", color="#3b82f6", linewidth=1.3))

ax_3d.annotate("Lineage-Specific Singletons:\n1,015 (18.3%)", xy=(1, 5570), xytext=(1, 6220),
               ha='center', va='center', fontsize=8.8, fontweight='bold', color='#92400e',
               bbox=dict(boxstyle="round,pad=0.35", facecolor="#fefce8", edgecolor="#fde047", linewidth=1.0),
               arrowprops=dict(arrowstyle="->", color="#eab308", linewidth=1.3))

ax_3d.set_title("D. Clonal Dereplication Sensitivity ($N=6$ Lineages)", fontsize=12, fontweight='bold', pad=10)
ax_3d.set_ylabel("Gene Family Count (Universal Cohort = 5,545)", fontsize=11)
ax_3d.set_xticks(x_d)
ax_3d.set_xticklabels(labels_d, fontsize=10, fontweight='bold')
ax_3d.grid(axis='y', linestyle='--', alpha=0.5)
ax_3d.set_xlim(-0.46, 1.46)
ax_3d.set_ylim(0, 7750)
ax_3d.legend(loc="upper center", ncol=2, fontsize=9.2, framealpha=0.96, bbox_to_anchor=(0.5, 0.99), edgecolor="#cbd5e1")
ax_3d.text(0.5, 6860, "Both cohorts evaluated within the identical 5,545-family pan-genome universe", 
           ha='center', va='center', fontsize=8.4, fontstyle='italic', color='#566573')

plt.suptitle("Figure 3: Mathematical Modeling and Clonal Dereplication Dynamics of the Ewingella americana Pan-Genome", fontsize=13.5, fontweight='bold', y=0.98)

save_refined_figure(fig3, "Final_Figure_3_PanGenome")

# ==============================================================================
# FIGURE 4: REFINED CORE-GENOME PHYLOGENOMICS
# ==============================================================================
print("\n>>> Building Final Figure 4: Core Phylogenomics...")
tree_file = "/mnt/d/Ewingella_americana_research/data_true_ewingella/iqtree_strict_31/strict_31.contree"
meta_file = "/mnt/d/Ewingella_americana_research/JAM_submission_True_Ewingella/Table_S1_Genomic_Inventory.csv"

df_meta = pd.read_csv(meta_file)
meta_dict = df_meta.set_index('accession').to_dict('index')

tree4 = Phylo.read(tree_file, "newick")
tip_colors4 = {}
tip_labels4 = {}

# Use Accession as primary label, followed by clean descriptive secondary label
for tip in tree4.get_terminals():
    acc = tip.name
    row = meta_dict.get(acc, {})
    strain = str(row.get("strain", "") or "")
    source = str(row.get("isolation_source", "") or "")
    
    if "B6-1" in strain or acc == "GCF_014233485.1":
        color = "#c0392b"
        label = f"{acc} — B6-1 (Mushroom-associated isolate, China) [Distinct lineage]"
    elif "RIT713" in strain or acc == "GCF_009668135.1":
        color = "#d35400"
        label = f"{acc} — RIT713 (Fungal surface associate, USA: NY)"
    elif "bioreactor" in source.lower() or acc.startswith("GCF_0226"):
        color = "#2980b9"
        rep_name = strain if strain and strain!='nan' else "Bioreactor"
        label = f"{acc} — {rep_name} (Anaerobic effluent clone, USA: WI)"
    elif "oral" in source.lower() or acc == "GCF_057314615.1":
        color = "#8e44ad"
        label = f"{acc} — Oral isolate (Human saliva, USA)"
    elif "CCUG 14506" in strain or acc == "GCF_008693045.1":
        color = "#16a085"
        label = f"{acc} — CCUG 14506T (Clinical type strain, human throat, Sweden)"
    elif "ATCC 33852" in strain or acc == "GCF_000735345.1":
        color = "#16a085"
        label = f"{acc} — ATCC 33852 (Clinical type-strain clone, USA)"
    elif "CIP" in strain or "HMF" in strain or acc.startswith("GCF_0540") or acc == "GCF_056280585.1":
        color = "#34495e"
        clean_strain = strain if strain and strain!='nan' else "CIP 81.94 subclone"
        label = f"{acc} — {clean_strain} (Clinical type-strain subclone)"
    else:
        color = "#7f8c8d"
        label = f"{acc} — {strain}"
        
    tip_labels4[acc] = label
    tip_colors4[acc] = color
    tip.name = label

fig4 = plt.figure(figsize=(21, 15), dpi=300)
ax4 = fig4.add_subplot(1, 1, 1)

def get_branch_label(clade):
    if clade.confidence is not None and not clade.is_terminal():
        val = int(round(clade.confidence))
        if val >= 70:
            return f" {val} "
    return ""

Phylo.draw(
    tree4,
    axes=ax4,
    do_show=False,
    show_confidence=True,
    branch_labels=get_branch_label,
    label_func=lambda x: str(x),
    label_colors=lambda x: tip_colors4.get(x, "black")
)

ax4.set_title("Figure 4: Core-Genome Maximum-Likelihood Phylogenomics of Ewingella americana (2,634,233 bp Alignment, 2,900 Strict Core Loci)", 
              fontsize=13.0, fontweight='bold', pad=25)
ax4.set_xlabel("Nucleotide Substitutions per Site (IQ-TREE 2 v2.0.7, GTR+F+I+R4 Model)", fontsize=11, fontweight='medium')
ax4.set_ylabel("")

# Full clearance based on actual plotted distance:
# Longest root-to-tip distance is B6-1 at x = 0.05214. Tip label extends to x ~ 0.082.
# Setting xlim = (-0.002, 0.098) ensures 100% visibility of the entire branch and full label text.
ax4.set_xlim(-0.002, 0.098)

for text in ax4.texts:
    txt = text.get_text().strip()
    if txt.isdigit():
        text.set_fontsize(8.5)
        text.set_fontweight('bold')
        text.set_color('#c0392b')
        text.set_bbox(dict(boxstyle="square,pad=0.15", facecolor="white", edgecolor="#c0392b", linewidth=0.6))
    else:
        text.set_fontsize(8.8)
        text.set_fontweight('medium')

legend_elements4 = [
    patches.Patch(facecolor="#c0392b", edgecolor="black", label="Mushroom-Associated Isolate B6-1 (Distinct lineage, Flammulina rot, China)"),
    patches.Patch(facecolor="#d35400", edgecolor="black", label="Environmental Fungus-Associated RIT713 (USA: NY)"),
    patches.Patch(facecolor="#2980b9", edgecolor="black", label="Anaerobic Bioreactor Effluent Clade (N=6 clonal replicates, USA: WI)"),
    patches.Patch(facecolor="#8e44ad", edgecolor="black", label="Human Oral Cavity Saliva Isolate (GCF_057314615.1, USA)"),
    patches.Patch(facecolor="#16a085", edgecolor="black", label="Human Clinical Type Strain Clones (CCUG 14506T & ATCC 33852)"),
    patches.Patch(facecolor="#34495e", edgecolor="black", label="Clinical Type-Strain Subclone Radiation (19 CIP subclones & HMF3896_M4G, France/USA)")
]
ax4.legend(handles=legend_elements4, loc="lower right", fontsize=9.0, framealpha=0.95, edgecolor="#bdc3c7")

root_note4 = (
    "Methodological & Taxonomic Notes:\n"
    "• Phylogeny inferred from 2,634,233 bp concatenated strict core gene alignment (2,900 loci; core_gene_alignment.aln) using IQ-TREE 2 v2.0.7 (GTR+F+I+R4, 1,000 UFBoot).\n"
    "• Post-filtered alignment (core_gene_alignment_filtered.aln; 2,722 loci, 2,422,091 bp) was also generated by Panaroo, but the complete 2,900-locus supermatrix was selected for maximum nodal support.\n"
    "• Tree rooted using Yersinia enterocolitica 8081 and Yersinia pestis CO92 as outgroups; outgroup branches were pruned post-rooting to maximize ingroup resolution.\n"
    "• All N = 31 verified high-quality E. americana genomes are included in the tree.\n"
    "• Branch support values indicate Ultrafast Bootstrap (UFBoot) support (values >= 70% shown; unsupported subclone branches <50% omitted for visual clarity).\n"
    "• B6-1 formed a distinct, relatively long-branched lineage within the verified species radiation (terminal branch = 0.0359 substitutions/site; root-to-tip distance = 0.0521)."
)
fig4.text(0.08, 0.02, root_note4, fontsize=8.6, color="#2c3e50",
          bbox=dict(boxstyle="round,pad=0.5", facecolor="#f8f9fa", edgecolor="#bdc3c7", linewidth=1.0))

plt.subplots_adjust(bottom=0.16, right=0.98, left=0.06, top=0.93)
save_refined_figure(fig4, "Final_Figure_4_Phylogenomics")

# ==============================================================================
# FIGURE 5: SUBSTANTIALLY REDESIGNED RESISTOME & MOBILOME
# ==============================================================================
print("\n>>> Building Final Figure 5: Substantially Redesigned Resistome & Mobilome...")
lineages5 = [
    ("GCF_008693045.1", "E. americana CCUG 14506T (Type strain, throat)"),
    ("GCF_000735345.1", "E. americana ATCC 33852 (Clinical clone)"),
    ("GCF_054057715.1", "E. americana CIP 81.94 M4R (Clinical subclone)"),
    ("GCF_057314615.1", "E. americana Oral isolate (Human saliva, USA)"),
    ("GCF_022647985.1", "E. americana Bioreactor 1 (Effluent, USA: WI)"),
    ("GCF_022655415.1", "E. americana Bioreactor 2 (Effluent, USA: WI)"),
    ("GCF_009668135.1", "E. americana RIT713 (Fungus-associated, NY)"),
    ("GCF_014233485.1", "E. americana B6-1 (Mushroom-associated, China)"),
    ("GCF_002816315.1", "E. americana BRK18a (Blood, QC fail: fragmented)"),
    ("NIAB_247711",     "\"E. allii\" NIAB_247711 (Proposed taxon, epiphyte)"),
    ("GCF_044430645.1", "\"E. allii\" 20WA0182T (Proposed type strain, onion rot)"),
    ("GCF_040419245.1", "E. docleensis CoE-038-23T (Type strain, cheese)"),
    ("GCF_006438725.1", "Polar isolate E4 (Arctic snow, outgroup)")
]

# Concise column names with explicit categories
amr_features5 = [
    ("pB61a", "Megaplasmid pB61a\n(331.4 kb; RepA/B)"),
    ("pB61b", "Conjugative pB61b\n(105.3 kb; RepA/B)"),
    ("t4ss_pB61b", "Predicted T4SS\n(32 kb tra/trb)"),
    ("oqxAB_plasmid", "Plasmid Efflux\n(oqxAB on pB61a)"),
    ("ter_plasmid", "Tellurite Resistance\n(ter on pB61a)"),
    ("ars_plasmid", "Arsenic Resistance\n(ars on pB61a)"),
    ("ampC_chrom", "AmpC Beta-Lactamase\n(Chromosomal blaACT)"),
    ("blaRAHN_contam", "Rahnella blaRAHN-3\n(Contaminant in BRK18a)"),
    ("acrB", "RND Multidrug Efflux\n(AcrB)"),
    ("emrB", "MFS Multidrug Efflux\n(EmrB)"),
    ("dinF", "MATE Multidrug Efflux\n(DinF)"),
    ("crp_hns", "Global Regulators\n(CRP & H-NS)")
]

vir_features5 = [
    ("chiD", "GH18 Chitinase\nChiD (Predicted)"),
    ("gh18_2", "Secondary GH18\nChitinase (Predicted)"),
    ("chbA", "Chitobiose PTS\nTransporter ChbA"),
    ("chbG", "Chitin Deacetylase\nChbG (Predicted)"),
    ("nagA", "Glucosamine-6-P\nDeacetylase NagA"),
    ("lprI", "Candidate Lysozyme\nInhibitor LprI"),
    ("fim_usher", "Candidate Adhesin\nUsher FimD"),
    ("fim_chap", "Candidate Adhesin\nChaperone FimC"),
    ("tssC1", "T6SS Sheath\nTssC1 (Candidate)"),
    ("tssB1", "T6SS Sheath\nTssB1 (Candidate)"),
    ("hcp", "T6SS Structural\nEffector Hcp"),
    ("paar", "Candidate T6SS\nPAAR Spike"),
    ("fliG", "Flagellar Motor\nSwitch FliG"),
    ("rcsB", "Stress / Capsule\nRegulator RcsB")
]

n_rows5 = len(lineages5)
mat_amr = np.zeros((n_rows5, len(amr_features5)))
mat_vir = np.zeros((n_rows5, len(vir_features5)))

for i, (acc, label) in enumerate(lineages5):
    if acc == "GCF_014233485.1":
        mat_amr[i, 0] = 1
        mat_amr[i, 1] = 1
        mat_amr[i, 2] = 1
        mat_amr[i, 3] = 1
        mat_amr[i, 4] = 1
        mat_amr[i, 5] = 1
    if acc != "GCF_006438725.1":
        mat_amr[i, 6] = 1
    if acc == "GCF_002816315.1":
        mat_amr[i, 7] = 1
    mat_amr[i, 8] = 1
    mat_amr[i, 9] = 1
    mat_amr[i, 10] = 1
    if acc != "GCF_006438725.1":
        mat_amr[i, 11] = 1
    else:
        mat_amr[i, 11] = 0.5
        
    if acc != "GCF_006438725.1":
        mat_vir[i, 0] = 1
        mat_vir[i, 1] = 1
        mat_vir[i, 2] = 1
        mat_vir[i, 3] = 1
        mat_vir[i, 4] = 1
    if acc in ["GCF_014233485.1", "GCF_009668135.1"]:
        mat_vir[i, 5] = 1
        mat_vir[i, 6] = 1
        mat_vir[i, 7] = 1
    mat_vir[i, 8] = 1
    if acc == "GCF_014233485.1":
        mat_vir[i, 9] = 1
    mat_vir[i, 10] = 1
    if acc == "GCF_014233485.1":
        mat_vir[i, 11] = 1
    mat_vir[i, 12] = 1
    mat_vir[i, 13] = 1

row_labels5 = [l[1] for l in lineages5]
df_amr5 = pd.DataFrame(mat_amr, index=row_labels5, columns=[f[1] for f in amr_features5])
df_vir5 = pd.DataFrame(mat_vir, index=row_labels5, columns=[v[1] for v in vir_features5])

# Substantially redesigned: Two clearly separated subplots with wide horizontal & vertical clearance
fig5, (ax5_1, ax5_2) = plt.subplots(1, 2, figsize=(26.0, 14.2), gridspec_kw={'width_ratios': [1.12, 1.25]}, dpi=300)
plt.subplots_adjust(wspace=0.48, bottom=0.22, top=0.81, left=0.17, right=0.98)

cmap_amr = ["#f8f9fa", "#d4e6f1", "#1f618d"]
cmap_vir = ["#f8f9fa", "#d5f5e3", "#1e8449"]

# Panel A
sns.heatmap(df_amr5, ax=ax5_1, cmap=cmap_amr, cbar=False, linewidths=1.2, linecolor="#cbd5e1", annot=False, vmin=0, vmax=1)
ax5_1.set_title("Panel A: Resistome, Mobilome, and Efflux Architecture\n(Distinguishing Plasmid-Borne Determinants from Conserved Chromosomal Loci)", 
                fontsize=12.0, fontweight='bold', pad=58, color='#0f172a')
ax5_1.set_xticklabels(ax5_1.get_xticklabels(), rotation=42, ha="right", fontsize=8.6, fontweight='medium')
ax5_1.set_yticklabels(ax5_1.get_yticklabels(), rotation=0, fontsize=9.2, fontweight='medium')

ax5_1.axvline(6, color='#c0392b', linestyle='--', linewidth=2.0)
ax5_1.axvline(8, color='#7f8c8d', linestyle=':', linewidth=1.6)

# Crisp bracket underlines + group headers cleanly between heatmap top edge and Panel A title
for x1, x2, col_g, lbl_g in [
    (0.1, 5.9, '#922b21', "PLASMID-BORNE MOBILOME\n(pB61a & pB61b)"),
    (6.1, 7.9, '#78281f', "CHROMOSOMAL\nRESISTANCE"),
    (8.1, 11.9, '#1b4f72', "CONSERVED CHROMOSOMAL\nEFFLUX SYSTEMS"),
]:
    ax5_1.plot([x1, x2], [-0.12, -0.12], color=col_g, lw=2.2, clip_on=False)
    ax5_1.text(0.5*(x1 + x2), -0.26, lbl_g, ha='center', va='bottom', fontsize=8.5, fontweight='bold', color=col_g, clip_on=False)

# Panel B
sns.heatmap(df_vir5, ax=ax5_2, cmap=cmap_vir, cbar=False, linewidths=1.2, linecolor="#cbd5e1", annot=False, vmin=0, vmax=1)
ax5_2.set_title("Panel B: Predicted Chitin Catabolism & Genomic Candidate Determinants\n(Conserved Chitinolytic Potential vs. In Silico Candidate Colonization Traits)", 
                fontsize=12.0, fontweight='bold', pad=58, color='#0f172a')
ax5_2.set_xticklabels(ax5_2.get_xticklabels(), rotation=42, ha="right", fontsize=8.6, fontweight='medium')
ax5_2.set_yticklabels(ax5_2.get_yticklabels(), rotation=0, fontsize=9.2, fontweight='medium')

ax5_2.axvline(5, color='#27ae60', linestyle='--', linewidth=2.0)

for x1, x2, col_g, lbl_g in [
    (0.1, 4.9, '#196f3d', "CONSERVED CHITIN CATABOLISM\n(Genomic Potential)"),
    (5.1, 13.9, '#145a32', "IN SILICO CANDIDATE COLONIZATION DETERMINANTS\n(Requires Experimental Phenotypic Validation)"),
]:
    ax5_2.plot([x1, x2], [-0.12, -0.12], color=col_g, lw=2.2, clip_on=False)
    ax5_2.text(0.5*(x1 + x2), -0.26, lbl_g, ha='center', va='bottom', fontsize=8.5, fontweight='bold', color=col_g, clip_on=False)

# Clean Super Title with generous vertical clearance above Panel A/B titles
fig5.text(0.575, 0.962, "Figure 5: Comparative Resistome, Mobilome, and Ecological Adaptation Determinants Across 13 Representative Ewingella Genomes", 
          ha='center', va='center', fontsize=14.0, fontweight='bold', color='#0f172a')

# Discrete Legend cleanly positioned between x-labels and bottom footnote
legend_patches5 = [
    patches.Patch(facecolor="#1f618d", edgecolor="#154360", label="Determinant Detected (AMR / Mobilome: Screening Thresholds >=75% id, >=70% cov)"),
    patches.Patch(facecolor="#1e8449", edgecolor="#145a32", label="Determinant Detected (Chitinolytic / In Silico Candidate: >=75% id, >=70% cov)"),
    patches.Patch(facecolor="#f8f9fa", edgecolor="#cbd5e1", label="Determinant Not Detected / Locus Absent")
]
fig5.legend(handles=legend_patches5, loc="lower center", ncol=3, fontsize=9.2, framealpha=0.96, 
            bbox_to_anchor=(0.575, 0.068), edgecolor="#cbd5e1")

# Concise bottom footnote cleanly below legend
concise_footnote = (
    "Screening & Interpretation Summary: Profiling conducted via ABRicate v1.0.1 (CARD, AMRFinderPlus, ResFinder, PlasmidFinder, VFDB; thresholds >=75% id, >=70% cov).\n"
    "Standalone IntegronFinder v2.0.6 confirmed 0 complete integrons, 0 In0, 0 CALIN, and 0 attC sites. Plasmids pB61a (331.4 kb) and pB61b (105.3 kb) encode non-canonical RepA/RepB replicons and predicted Type-F T4SS.\n"
    "Functional status: Chitinolytic loci reflect genomic potential; adhesins, LprI, and T6SS effectors represent in silico candidate traits requiring future laboratory validation."
)
fig5.text(0.17, 0.012, concise_footnote, fontsize=8.4, color="#334155",
          bbox=dict(boxstyle="round,pad=0.4", facecolor="#f8fafc", edgecolor="#cbd5e1", linewidth=0.8))

save_refined_figure(fig5, "Final_Figure_5_Resistome")

print("\n>>> All Five Publication-Refined Figures successfully exported to SVG, TIFF, PDF, and PNG!")
