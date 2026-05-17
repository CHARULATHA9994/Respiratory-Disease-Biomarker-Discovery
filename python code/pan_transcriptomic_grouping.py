# ============================================================
# PAN-TRANSCRIPTOMIC PATHWAY GROUPING
# ============================================================
# Project  : Respiratory Disease Biomarker Discovery
# Author   : M. Charulatha
# Degree   : MSc Bioinformatics and Data Science
# Purpose  : Integrating enrichment results from all three
#            respiratory diseases (Dust Allergy, Sinusitis,
#            Allergic Asthma) into one pan-transcriptomic
#            pathway score per functional category
# ============================================================

import pandas as pd

# ============================================================
# STEP 1 — LOAD PAN-TRANSCRIPTOMIC DATA
# This file contains pathway enrichment results combined
# across all three diseases from the transcriptomic analysis
# ============================================================

pan = pd.read_excel("/content/pan transcr.xlsx")
pan.columns = pan.columns.str.strip()
pan = pan.rename(columns={'Pathway': 'pathway', 'Log10(q)': 'log10_q'})

# ============================================================
# STEP 2 — CONVERT TO PAN SCORE
# Negative log10(q) gives a positive enrichment score
# Higher score = more statistically significant pathway
# ============================================================

pan['Pan_Score'] = -pan['log10_q']

# ============================================================
# STEP 3 — CLASSIFY PATHWAYS INTO FUNCTIONAL CATEGORIES
# Each pathway is mapped to a biological theme based on
# keywords in the pathway name
# This allows comparison across diseases at a theme level
# ============================================================

def categorize(p):
    p = str(p).lower()

    # CYTOKINE SIGNALING
    # Pathways involving interleukin, cytokine, JAK-STAT signaling
    if any(x in p for x in [
        'interleukin', 'cytokine', 'jak', 'stat', 'tyrobp',
        'immune signaling'
    ]):
        return 'Cytokine signaling'

    # IMMUNE ACTIVATION / REGULATION
    # Pathways involving T cells, lymphocytes, immune responses
    elif any(x in p for x in [
        't cell', 'lymphocyte', 'immune response', 'cell activation'
    ]):
        return 'Immune activation'

    # INFLAMMATION
    # Pathways with inflammatory processes
    elif 'inflamm' in p:
        return 'Inflammation'

    # HOST-PATHOGEN INTERACTION
    # Pathways involving bacterial infection or pathogens
    elif any(x in p for x in [
        'bacter', 'infection', 'pathogen'
    ]):
        return 'Host-pathogen interaction'

    # METABOLIC ADAPTATION
    # Pathways involving metabolic processes and biosynthesis
    elif any(x in p for x in [
        'metabolism', 'biosynthesis', 'acid'
    ]):
        return 'Metabolic adaptation'

    # OTHER — does not fit any defined category
    else:
        return 'Other'

pan['Category'] = pan['pathway'].apply(categorize)

# ============================================================
# STEP 4 — GROUP BY FUNCTIONAL CATEGORY
# Sum Pan_Score for all pathways within each category
# This gives one total score per biological theme
# ============================================================

pan_grouped = pan.groupby('Category')['Pan_Score'].sum().reset_index()

# ============================================================
# STEP 5 — VIEW PATHWAY-CATEGORY MAPPING
# Printing how each pathway was assigned to a category
# Useful for verification and quality check
# ============================================================

print("\nPATHWAY TO CATEGORY MAPPING:\n")
print(pan[['pathway', 'Category']])

# ============================================================
# STEP 6 — PRINT FINAL GROUPED RESULT
# Displaying categories ranked by total Pan_Score
# Result showed: Host-pathogen (3.44), Immune activation (2.55),
# Inflammation (2.39), Cytokine signaling (2.36),
# Metabolic adaptation (1.17)
# ============================================================

print("\nPAN GROUPED RESULT:\n")
print(pan_grouped.sort_values(by='Pan_Score', ascending=False))

# ============================================================
# STEP 7 — SAVE OUTPUT
# Saving grouped result for use in the Bridge Step
# ============================================================

pan_grouped.to_csv("pan_grouped.csv", index=False)
print("\nSaved: pan_grouped.csv")
