import os
import glob
import math
import numpy as np
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq

# Paths
aln_path   = "/mnt/d/Ewingella_americana_research/data_true_ewingella/panaroo_strict_31/core_gene_alignment.aln"
genes_dir  = "/mnt/d/Ewingella_americana_research/data_true_ewingella/panaroo_strict_31/aligned_gene_sequences"
s1_path    = "/mnt/d/Ewingella_americana_research/FINAL_JAM_SUBMISSION/Table_S1_Genomic_Inventory.csv"
out_dir    = "/mnt/d/Ewingella_americana_research/FINAL_JAM_SUBMISSION/Phylogeny_Q1"
os.makedirs(out_dir, exist_ok=True)

df_s1 = pd.read_csv(s1_path).set_index('accession')

# ==============================================================================
# 1. NUCLEOTIDE COMPOSITION, SKEW & EXACT PAIRWISE p-DISTANCES (2,634,233 bp)
# ==============================================================================
print(">>> 1. Loading 2,634,233 bp strict core alignment...")
records = list(SeqIO.parse(aln_path, "fasta"))
taxa = [r.id for r in records]
n_taxa = len(taxa)
aln_len = len(records[0].seq)
print(f"    Loaded {n_taxa} taxa x {aln_len:,} bp")

# Convert to numpy uint8 array for blazing-fast exact vectorized calculation
seq_mat = np.zeros((n_taxa, aln_len), dtype=np.uint8)
comp_rows = []
for i, rec in enumerate(records):
    s = str(rec.seq).upper().encode('ascii')
    arr = np.frombuffer(s, dtype=np.uint8)
    seq_mat[i] = arr
    cA = int(np.sum(arr == ord('A')))
    cT = int(np.sum(arr == ord('T')))
    cG = int(np.sum(arr == ord('G')))
    cC = int(np.sum(arr == ord('C')))
    valid = cA + cT + cG + cC
    gaps = aln_len - valid
    gc_pct = 100.0 * (cG + cC) / valid if valid else 0.0
    at_skew = (cA - cT) / (cA + cT) if (cA + cT) else 0.0
    gc_skew = (cG - cC) / (cG + cC) if (cG + cC) else 0.0
    meta = df_s1.loc[rec.id] if rec.id in df_s1.index else {}
    comp_rows.append({
        'Accession': rec.id,
        'Strain': meta.get('strain', '') if isinstance(meta, pd.Series) else '',
        'Lineage': meta.get('Biological_Lineage_Assignment', '') if isinstance(meta, pd.Series) else '',
        'Aligned_bp': aln_len,
        'Valid_ACGT_bp': valid,
        'Gap_Ambig_pct': round(100.0 * gaps / aln_len, 4),
        'GC_pct': round(gc_pct, 4),
        'AT_skew': round(at_skew, 6),
        'GC_skew': round(gc_skew, 6)
    })

df_comp = pd.DataFrame(comp_rows)

# Compute exact pairwise p-distances & raw core SNP counts (only unambiguous A/C/G/T sites)
print(">>> 2. Computing exact pairwise nucleotide p-distances & core SNP counts...")
is_acgt = (seq_mat == ord('A')) | (seq_mat == ord('C')) | (seq_mat == ord('G')) | (seq_mat == ord('T'))

pdist_mat = np.zeros((n_taxa, n_taxa), dtype=float)
snp_mat   = np.zeros((n_taxa, n_taxa), dtype=int)
for i in range(n_taxa):
    for j in range(i + 1, n_taxa):
        both_valid = is_acgt[i] & is_acgt[j]
        n_valid = int(np.sum(both_valid))
        n_diff = int(np.sum((seq_mat[i] != seq_mat[j]) & both_valid))
        pd_val = n_diff / n_valid if n_valid > 0 else 0.0
        pdist_mat[i, j] = pdist_mat[j, i] = pd_val
        snp_mat[i, j]   = snp_mat[j, i]   = n_diff

df_pdist = pd.DataFrame(pdist_mat, index=taxa, columns=taxa)
df_snps  = pd.DataFrame(snp_mat, index=taxa, columns=taxa)
df_pdist.to_csv(os.path.join(out_dir, "pairwise_pdist_matrix_31.csv"))
df_snps.to_csv(os.path.join(out_dir, "pairwise_core_snps_31.csv"))

# Add p-distance & SNPs relative to Type Strain CCUG 14506T (GCF_008693045.1)
ref_acc = "GCF_008693045.1"
ref_idx = taxa.index(ref_acc)
df_comp['pDist_to_CCUG14506T'] = [round(pdist_mat[i, ref_idx], 6) for i in range(n_taxa)]
df_comp['Pct_Divergence_to_CCUG14506T'] = [round(100.0 * pdist_mat[i, ref_idx], 4) for i in range(n_taxa)]
df_comp['Core_SNPs_to_CCUG14506T'] = [int(snp_mat[i, ref_idx]) for i in range(n_taxa)]
df_comp.to_csv(os.path.join(out_dir, "Table_Pairwise_Distances_and_Composition.csv"), index=False)
print("    Saved Table_Pairwise_Distances_and_Composition.csv")

# Print quick summary of pairwise divergences by lineage
def get_short_lin(acc):
    if acc == "GCF_014233485.1": return "B6-1"
    if acc == "GCF_009668135.1": return "RIT713"
    if acc == "GCF_057314615.1": return "Oral"
    if acc.startswith("GCF_0226"): return "Bioreactor"
    if acc in ("GCF_008693045.1", "GCF_000735345.1"): return "TypeStrain_Ref"
    return "CIP_HMF_Subclones"

lin_labels = [get_short_lin(a) for a in taxa]
pairs_summary = {}
for i in range(n_taxa):
    for j in range(i + 1, n_taxa):
        k = tuple(sorted([lin_labels[i], lin_labels[j]]))
        pairs_summary.setdefault(k, []).append((pdist_mat[i, j], snp_mat[i, j]))

print("    Lineage Pairwise Summary:")
for k, vals in sorted(pairs_summary.items()):
    mean_p = 100.0 * np.mean([v[0] for v in vals])
    mean_s = np.mean([v[1] for v in vals])
    print(f"      {k[0]:18s} vs {k[1]:18s}: mean p-dist = {mean_p:.4f}% ({mean_s:,.1f} SNPs, n={len(vals)})")

# ==============================================================================
# 2. TRANSLATED AMINO ACID VALIDATION MATRIX + NEI-GOJOBORI Ka/Ks + RSCU
# ==============================================================================
print("\n>>> 3. Building Translated AA Validation Matrix + Nei-Gojobori Ka/Ks & RSCU...")

# Standard bacterial genetic code (Table 11)
CODON_TABLE = {}
bases = ['T', 'C', 'A', 'G']
aa_str = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
idx_c = 0
for b1 in bases:
    for b2 in bases:
        for b3 in bases:
            CODON_TABLE[b1 + b2 + b3] = aa_str[idx_c]
            idx_c += 1

STOP_CODONS = {'TAA', 'TAG', 'TGA'}

# Precompute Nei-Gojobori (1986) synonymous & nonsynonymous sites per sense codon
S_SITES = {}
N_SITES = {}
for cod, aa in CODON_TABLE.items():
    if aa == '*':
        continue
    s_cnt = 0.0
    for pos in range(3):
        for nb in bases:
            if nb == cod[pos]:
                continue
            mut = cod[:pos] + nb + cod[pos+1:]
            if CODON_TABLE[mut] == aa:
                s_cnt += 1.0 / 3.0
    S_SITES[cod] = s_cnt
    N_SITES[cod] = 3.0 - s_cnt

# Precompute exact NG86 pathway-averaged S_diff and N_diff for codon pairs
NG_PAIR_CACHE = {}
def ng_codon_diff(c1, c2):
    if c1 == c2:
        return (0.0, 0.0)
    key = (c1, c2)
    if key in NG_PAIR_CACHE:
        return NG_PAIR_CACHE[key]
    diffs = [p for p in range(3) if c1[p] != c2[p]]
    if len(diffs) == 1:
        res = (1.0, 0.0) if CODON_TABLE[c1] == CODON_TABLE[c2] else (0.0, 1.0)
    elif len(diffs) == 2:
        p1, p2 = diffs
        paths = [
            (c1, c1[:p1] + c2[p1] + c1[p1+1:], c2),
            (c1, c1[:p2] + c2[p2] + c1[p2+1:], c2),
        ]
        sd_tot, nd_tot, valid_p = 0.0, 0.0, 0
        for a, m, b in paths:
            if m in STOP_CODONS:
                continue
            valid_p += 1
            sd_tot += (1.0 if CODON_TABLE[a] == CODON_TABLE[m] else 0.0) + (1.0 if CODON_TABLE[m] == CODON_TABLE[b] else 0.0)
            nd_tot += (0.0 if CODON_TABLE[a] == CODON_TABLE[m] else 1.0) + (0.0 if CODON_TABLE[m] == CODON_TABLE[b] else 1.0)
        res = (sd_tot / valid_p, nd_tot / valid_p) if valid_p else (0.0, 2.0)
    else:
        # 3 differences: average over 6 paths avoiding stops
        import itertools
        sd_tot, nd_tot, valid_p = 0.0, 0.0, 0
        for perm in itertools.permutations([0, 1, 2]):
            m1 = list(c1); m1[perm[0]] = c2[perm[0]]; m1 = "".join(m1)
            m2 = list(m1); m2[perm[1]] = c2[perm[1]]; m2 = "".join(m2)
            if m1 in STOP_CODONS or m2 in STOP_CODONS:
                continue
            valid_p += 1
            for x, y in [(c1, m1), (m1, m2), (m2, c2)]:
                if CODON_TABLE[x] == CODON_TABLE[y]:
                    sd_tot += 1.0
                else:
                    nd_tot += 1.0
        res = (sd_tot / valid_p, nd_tot / valid_p) if valid_p else (0.0, 3.0)
    NG_PAIR_CACHE[key] = res
    return res

def compute_ng86_kaks(seq1, seq2):
    S_tot = N_tot = Sd_tot = Nd_tot = 0.0
    L = min(len(seq1), len(seq2))
    L -= L % 3
    for i in range(0, L, 3):
        c1 = seq1[i:i+3]
        c2 = seq2[i:i+3]
        if c1 not in S_SITES or c2 not in S_SITES:
            continue
        S_tot += 0.5 * (S_SITES[c1] + S_SITES[c2])
        N_tot += 0.5 * (N_SITES[c1] + N_SITES[c2])
        sd, nd = ng_codon_diff(c1, c2)
        Sd_tot += sd
        Nd_tot += nd
    if S_tot == 0 or N_tot == 0:
        return (np.nan, np.nan, np.nan)
    pS = Sd_tot / S_tot
    pN = Nd_tot / N_tot
    # Jukes-Cantor correction (Nei & Gojobori 1986)
    if pS >= 0.75 or pN >= 0.75:
        return (np.nan, np.nan, np.nan)
    dS = -0.75 * math.log(1.0 - (4.0 / 3.0) * pS) if pS > 0 else 0.0
    dN = -0.75 * math.log(1.0 - (4.0 / 3.0) * pN) if pN > 0 else 0.0
    kaks = (dN / dS) if dS > 1e-6 else np.nan
    return (dN, dS, kaks)

# Scan strict-core single-copy aligned genes in aligned_gene_sequences/
gene_files = sorted(glob.glob(os.path.join(genes_dir, "*.aln.fas")))
print(f"    Found {len(gene_files)} aligned gene files in aligned_gene_sequences/")

# Representative lineages for gene-wise Ka/Ks vs Type Strain CCUG 14506T
rep_targets = {
    'B6-1':       'GCF_014233485.1',
    'RIT713':     'GCF_009668135.1',
    'Bioreactor': 'GCF_022647985.1',
    'Oral':       'GCF_057314615.1',
}

aa_concat = {t: [] for t in taxa}
kaks_rows = []
codon_counts_ref = {c: 0 for c in S_SITES.keys()}
n_aa_genes = 0
n_aa_sites = 0

for gf in gene_files:
    gname = os.path.basename(gf).replace(".aln.fas", "")
    recs = list(SeqIO.parse(gf, "fasta"))
    if len(recs) != n_taxa:
        continue
    # Map headers back to genome accession (Panaroo prefixes header with accession;)
    gmap = {}
    for r in recs:
        acc_key = r.id.split(';')[0]
        if acc_key in aa_concat:
            gmap[acc_key] = str(r.seq).upper()
    if len(gmap) != n_taxa:
        continue
    g_len = len(next(iter(gmap.values())))
    if g_len < 300 or g_len % 3 != 0:
        continue
    # Verify clean reading frame (no internal stops in reference CCUG 14506T)
    ref_s = gmap[ref_acc]
    codons_ref = [ref_s[i:i+3] for i in range(0, g_len - 3, 3)]
    n_stops_internal = sum(1 for c in codons_ref if c in STOP_CODONS)
    if n_stops_internal > 0:
        continue

    # Accumulate RSCU codon counts on CCUG 14506T
    for c in codons_ref:
        if c in codon_counts_ref:
            codon_counts_ref[c] += 1

    # Translate gap-free codon columns for AA validation matrix (use clean core loci, cap at 500 loci or all if fast!)
    # Let's include ALL clean strict-core single-copy PCGs up to 600 evenly distributed, or ALL clean ones!
    # Actually translating 2,000 genes for 31 taxa is ~600k aa — IQ-TREE2 runs in ~60 sec on 300-500 clean loci (~150k aa)!
    # Let's collect Ka/Ks across ALL clean strict-core genes, and build a super-clean 400-locus (~130,000 aa) gap-free AA validation matrix!
    row_k = {'Gene': gname, 'Length_bp': g_len}
    for rlbl, racc in rep_targets.items():
        dN, dS, w = compute_ng86_kaks(ref_s, gmap[racc])
        row_k[f'Ka_{rlbl}'] = round(dN, 6) if not np.isnan(dN) else np.nan
        row_k[f'Ks_{rlbl}'] = round(dS, 6) if not np.isnan(dS) else np.nan
        row_k[f'KaKs_{rlbl}'] = round(w, 6) if not np.isnan(w) else np.nan
    kaks_rows.append(row_k)

    if n_aa_genes < 400:
        # Keep only gap-free codon columns across all 31 taxa for Matrix B
        good_cols = []
        for pos in range(0, g_len - 3, 3):
            ok = True
            for t in taxa:
                c = gmap[t][pos:pos+3]
                if c not in S_SITES:
                    ok = False
                    break
            if ok:
                good_cols.append(pos)
        if len(good_cols) >= 80:
            for t in taxa:
                aa_seq = "".join(CODON_TABLE[gmap[t][p:p+3]] for p in good_cols)
                aa_concat[t].append(aa_seq)
            n_aa_genes += 1
            n_aa_sites += len(good_cols)

df_kaks = pd.DataFrame(kaks_rows)
df_kaks.to_csv(os.path.join(out_dir, "Table_Core_KaKs_NG86.csv"), index=False)

# Compute RSCU table
aa_to_codons = {}
for c, a in CODON_TABLE.items():
    if a != '*':
        aa_to_codons.setdefault(a, []).append(c)

rscu_rows = []
for a, clist in sorted(aa_to_codons.items()):
    tot_a = sum(codon_counts_ref[c] for c in clist)
    k_syn = len(clist)
    for c in clist:
        cnt = codon_counts_ref[c]
        rscu = (cnt * k_syn / tot_a) if tot_a > 0 else 0.0
        rscu_rows.append({'AminoAcid': a, 'Codon': c, 'Count': cnt, 'RSCU': round(rscu, 4)})
df_rscu = pd.DataFrame(rscu_rows)
df_rscu.to_csv(os.path.join(out_dir, "Table_Core_RSCU.csv"), index=False)

# Write Translated Amino Acid Validation Matrix (Matrix B)
aa_fasta = os.path.join(out_dir, "core_aa_validation.faa")
with open(aa_fasta, 'w') as fh:
    for t in taxa:
        fh.write(f">{t}\n{''.join(aa_concat[t])}\n")

print(f"    Evaluated Ka/Ks across {len(df_kaks)} clean single-copy strict-core PCGs!")
for rlbl in rep_targets:
    s = df_kaks[f'KaKs_{rlbl}'].dropna()
    n_pur = int((s < 1.0).sum())
    print(f"      CCUG 14506T vs {rlbl:10s}: n={len(s):4d} genes, median Ka/Ks = {s.median():.4f}, mean = {s.mean():.4f}, purifying (<1) = {100.0*n_pur/len(s):.2f}%")
print(f"    Wrote Matrix B (Translated AA Validation): {n_aa_genes} single-copy PCGs, {n_aa_sites:,} gap-free aa sites -> {aa_fasta}")
