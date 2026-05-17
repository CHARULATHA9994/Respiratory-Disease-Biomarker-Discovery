# ============================================================
# QUANTITATIVE CROSS-OMICS PATHWAY INTEGRATION
# ============================================================
# Project  : Respiratory Disease Biomarker Discovery
# Author   : M. Charulatha
# Degree   : MSc Bioinformatics and Data Science
# Purpose  : Integrating Transcriptomics, Proteomics, and
#            Metabolomics data to identify shared biological
#            pathways across allergic respiratory diseases
#            (Dust Allergy, Sinusitis, Allergic Asthma)
# ============================================================

import pandas as pd

# ============================================================
# STEP 1 — LOAD FILES
# Loading RNA, Protein, and Metabolite data from Excel files
# Each file contains pathway enrichment results from sinusitis
# multi-omics datasets
# ============================================================

rna  = pd.read_excel("/content/RNA.SINUSITIS.CSV.xlsx")
prot = pd.read_excel("/content/PROTEIN.SINUSITIS.CSV.xlsx")
meta = pd.read_excel("/content/METABOLITE.SINUSITIS.CSV.xlsx")

# ============================================================
# STEP 2 — CLEAN COLUMN NAMES
# Removing extra spaces from column headers to avoid errors
# Renaming columns consistently across all three omics files
# ============================================================

for df in [rna, prot, meta]:
    df.columns = df.columns.str.strip()

rna  = rna.rename(columns={'Pathway': 'pathway', 'log10(q)': 'log10_q'})
prot = prot.rename(columns={'Pathway': 'pathway', 'log10(q)': 'log10_q'})
meta = meta.rename(columns={'Pathway': 'pathway', 'log10(q)': 'log10_q'})

# ============================================================
# STEP 3 — CONVERT LOG10(Q) TO ENRICHMENT SCORE
# Converting negative log10(q) values into positive scores
# Higher score = more statistically significant pathway
# Keeping only pathway name and score for each omics layer
# ============================================================

rna['RNA']           = -rna['log10_q']
prot['Protein']      = -prot['log10_q']
meta['Metabolite']   = -meta['log10_q']

rna  = rna[['pathway', 'RNA']]
prot = prot[['pathway', 'Protein']]
meta = meta[['pathway', 'Metabolite']]

# ============================================================
# STEP 4 — MERGE ALL THREE OMICS LAYERS
# Combining RNA, Protein, and Metabolite data into one table
# Using outer join so no pathway is lost from any layer
# Missing values filled with 0 (pathway absent in that layer)
# ============================================================

merged = rna.merge(prot, on='pathway', how='outer') \
            .merge(meta, on='pathway', how='outer') \
            .fillna(0)

# ============================================================
# STEP 5 — FUNCTIONAL CATEGORY CLASSIFICATION
# Each pathway is assigned to a biological theme category
# based on keywords in its name
# Categories: Cytokine signaling, Immune activation,
# Inflammation, Host-pathogen interaction,
# Metabolic adaptation, Cellular processes, Other
# ============================================================

def categorize(p):
    p = str(p).lower()

    # CYTOKINE SIGNALING
    if any(x in p for x in [
        'interleukin', 'cytokine', 'jak', 'stat', 'tyrobp',
        'immune signaling', 'immune effector', 'inflammatory signaling',
        'immunoregulatory'
    ]):
        return 'Cytokine signaling'

    # IMMUNE ACTIVATION / REGULATION
    elif any(x in p for x in [
        'adaptive immune', 'lymphocyte', 'leukocyte', 'chemotaxis',
        'neutrophil', 'phagocytosis', 'cell activation',
        'immune response', 'immune system process'
    ]):
        return 'Immune activation'

    # INFLAMMATION
    elif 'inflamm' in p:
        return 'Inflammation'

    # HOST-PATHOGEN / DISEASE
    elif any(x in p for x in [
        'bacter', 'infection', 'pathogen', 'lupus'
    ]):
        return 'Host-pathogen interaction'

    # METABOLIC ADAPTATION
    elif any(x in p for x in [
        'metabolism', 'biosynthesis', 'catabolic',
        'folate', 'amino', 'lipid', 'glutathione', 'nitrogen'
    ]):
        return 'Metabolic adaptation'

    # CELLULAR / STRUCTURAL
    elif any(x in p for x in [
        'organelle', 'locomotion', 'matrisome'
    ]):
        return 'Cellular processes'

    # OTHER
    else:
        return 'Other'

merged['Category'] = merged['pathway'].apply(categorize)

# ============================================================
# STEP 6 — CALCULATE COMBINED SCORE
# Summing RNA + Protein + Metabolite scores for each pathway
# Higher combined score = stronger cross-omics evidence
# ============================================================

merged['Combined_Score'] = merged[['RNA', 'Protein', 'Metabolite']].sum(axis=1)

# ============================================================
# STEP 7 — CONSISTENCY SCORE
# Checking how many omics layers support each pathway
# Consistency = 3 means all three omics agree (strongest)
# Consistency = 2 means two omics agree
# Consistency = 1 means only one omics layer detected it
# ============================================================

def consistency_info(row):
    present = []
    if row['RNA'] > 0:
        present.append('RNA')
    if row['Protein'] > 0:
        present.append('Protein')
    if row['Metabolite'] > 0:
        present.append('Metabolite')
    return len(present), " + ".join(present)

merged[['Consistency', 'Omics']] = merged.apply(
    lambda row: pd.Series(consistency_info(row)), axis=1
)

# ============================================================
# STEP 8 — GROUP BY FUNCTIONAL CATEGORY
# Aggregating all pathways by their biological theme
# Summing scores and finding maximum consistency per category
# ============================================================

grouped = merged.groupby('Category').agg({
    'RNA'           : 'sum',
    'Protein'       : 'sum',
    'Metabolite'    : 'sum',
    'Combined_Score': 'sum',
    'Consistency'   : 'max'
}).reset_index()

# ============================================================
# STEP 9 — CALCULATE OVERLAP PERCENTAGE
# Finding what percentage of pathways appear in 2+ omics
# True overlap = same pathway detected across multiple layers
# ============================================================

total           = len(grouped)
common          = len(grouped[grouped['Consistency'] >= 2])
overlap_percent = (common / total) * 100

merged['Pathway_Consistency'] = (
    (merged['RNA']        > 0).astype(int) +
    (merged['Protein']    > 0).astype(int) +
    (merged['Metabolite'] > 0).astype(int)
)

true_overlap = merged[merged['Pathway_Consistency'] >= 2]

# ============================================================
# STEP 10 — CHECK UNCATEGORISED PATHWAYS
# Reviewing pathways that did not fit any category above
# ============================================================

other_df = merged[merged['Category'] == 'Other']
print("\nOTHER PATHWAYS:\n")
print(other_df[['pathway']])

# ============================================================
# STEP 11 — PRINT FINAL RESULTS
# Displaying final ranked pathway categories by combined score
# ============================================================

print("\nFINAL RESULT:\n")
print(grouped.sort_values(by='Combined_Score', ascending=False))
print("\nOVERLAP %:", overlap_percent)

# ============================================================
# STEP 12 — SAVE OUTPUT FILES
# Saving results to CSV for further analysis and reporting
# ============================================================

grouped.to_csv("final___multiomics.csv",   index=False)
merged.to_csv("detailed___mapping.csv",    index=False)

print("\nFiles saved: final___multiomics.csv and detailed___mapping.csv")
