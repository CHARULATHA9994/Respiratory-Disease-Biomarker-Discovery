# ============================================================
# BRIDGE STEP — PAN-TRANSCRIPTOMIC vs MULTI-OMICS INTEGRATION
# ============================================================
# Project  : Respiratory Disease Biomarker Discovery
# Author   : M. Charulatha
# Degree   : MSc Bioinformatics and Data Science
# Purpose  : Comparing pan-transcriptomic pathway scores with
#            multi-omics (RNA + Protein + Metabolite) scores
#            to identify pathways consistently supported
#            across all biological data layers
# ============================================================

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

# ============================================================
# STEP 1 — LOAD DATA
# pan_grouped.csv   = output from pan-transcriptomic step
# final_multiomics.csv = output from cross-omics integration
# Both files contain pathway categories and their scores
# ============================================================

pan   = pd.read_csv("/content/pan_grouped.csv")
multi = pd.read_csv("/content/final___multiomics.csv")

pan.columns   = pan.columns.str.strip()
multi.columns = multi.columns.str.strip()

# ============================================================
# STEP 2 — RENAME COLUMNS FOR CLARITY
# Renaming score columns to clearly identify their source
# Pan_Score  = score from transcriptomics only
# Multi_Score = score from RNA + Protein + Metabolite combined
# ============================================================

pan   = pan.rename(columns={'Pan_Score': 'Pan_Score'})
multi = multi.rename(columns={'Combined_Score': 'Multi_Score'})

# ============================================================
# STEP 3 — MERGE (BRIDGE STEP)
# Combining both datasets on the shared 'Category' column
# Outer join ensures no category is lost from either side
# Missing values filled with 0
# ============================================================

bridge = pd.merge(multi, pan, on='Category', how='outer').fillna(0)

# ============================================================
# STEP 4 — NORMALIZE SCORES
# Dividing each score by its maximum to bring both
# onto a 0-1 scale — this ensures fair comparison
# Without normalisation, one score could dominate unfairly
# ============================================================

bridge['Pan_norm']   = bridge['Pan_Score']   / bridge['Pan_Score'].max()
bridge['Multi_norm'] = bridge['Multi_Score'] / bridge['Multi_Score'].max()

# ============================================================
# STEP 5 — CALCULATE DIFFERENCE AND AGREEMENT
# Difference = how much Pan and Multi scores diverge
# Agreement  = 1 - Difference (higher = more consistent)
# Agreement close to 1.0 = both omics layers agree strongly
# ============================================================

bridge['Difference'] = bridge['Pan_norm'] - bridge['Multi_norm']
bridge['Agreement']  = 1 - abs(bridge['Difference'])

# ============================================================
# STEP 6 — CLASSIFY PATHWAY RELATIONSHIP
# Each category is labelled based on which omics support it:
# "Common core pathway"   = both Pan AND Multi strongly agree
# "Transcript-dominant"   = only transcriptomics shows it
# "Multi-layer supported" = multi-omics shows it more strongly
# "Low relevance"         = weak signal in both
# ============================================================

def interpretation(row):
    if row['Pan_norm'] > 0.3 and row['Multi_norm'] > 0.3:
        return "Common core pathway"
    elif row['Pan_norm'] > 0.3 and row['Multi_norm'] <= 0.3:
        return "Transcript-dominant"
    elif row['Pan_norm'] <= 0.3 and row['Multi_norm'] > 0.3:
        return "Multi-layer supported"
    else:
        return "Low relevance"

bridge['Interpretation'] = bridge.apply(interpretation, axis=1)

# ============================================================
# STEP 7 — SORT BY MULTI-OMICS SCORE
# Ranking categories by multi-omics evidence strength
# ============================================================

bridge = bridge.sort_values(by='Multi_norm', ascending=False)

# ============================================================
# STEP 8 — PRINT AND SAVE BRIDGE RESULTS
# ============================================================

print("\nBRIDGE ANALYSIS RESULT:\n")
print(bridge)

bridge.to_csv("bridge_analysis.csv", index=False)
print("\nSaved: bridge_analysis.csv")

# ============================================================
# VISUALISATION 1 — BAR CHART: Pan vs Multi-omics Comparison
# Comparing normalised scores side by side for each category
# Blue bars = Pan-transcriptomics score
# Orange bars = Multi-omics score
# ============================================================

categories = bridge['Category'].tolist()
pan_scores = bridge['Pan_norm'].tolist()
multi_scores = bridge['Multi_norm'].tolist()

plt.figure()
x = np.arange(len(categories))
plt.bar(x - 0.2, pan_scores,   width=0.4, label='Pan-transcriptomics')
plt.bar(x + 0.2, multi_scores, width=0.4, label='Multi-omics')
plt.xticks(x, categories, rotation=40, ha='right')
plt.ylabel("Normalized Score")
plt.title("Pan vs Multi-omics Comparison")
plt.legend()
plt.tight_layout()
plt.show()

# ============================================================
# VISUALISATION 2 — HEATMAP: Omics Contribution per Category
# Shows how much RNA, Protein, and Metabolite each contribute
# to each biological pathway category
# Colour scale: light yellow (low) to dark red (high)
# ============================================================

data = bridge[['RNA', 'Protein', 'Metabolite']].values

log_data  = np.log1p(data)
norm_data = log_data / log_data.max()

colors = ["#ffffb2", "#fecc5c", "#fd8d3c", "#f03b20", "#bd0026"]
cmap   = LinearSegmentedColormap.from_list("warm_strong", colors)

plt.figure(figsize=(6, 5))
im = plt.imshow(norm_data, cmap=cmap, aspect='auto')
plt.grid(which='both', color='white', linestyle='-', linewidth=0.4, alpha=0.2)
plt.tick_params(axis='both', which='both', length=0)
plt.colorbar(im, label='Normalized contribution')
plt.xticks([0, 1, 2], ['RNA', 'Protein', 'Metabolite'])
plt.yticks(range(len(bridge)), bridge['Category'])
plt.title("Omics Contribution Heatmap")
plt.tight_layout()
plt.show()

# ============================================================
# VISUALISATION 3 — BAR CHART: Pathway Agreement Summary
# Shows how many pathway categories fall into each
# interpretation class (Common core, Transcript-dominant, etc.)
# ============================================================

counts = bridge['Interpretation'].value_counts()

plt.figure(figsize=(5, 4))
bars = plt.bar(
    counts.index,
    counts.values,
    width=0.5,
    edgecolor='black',
    linewidth=0.6
)

plt.ylabel('Number of pathways', fontsize=11)
plt.title('Pathway Agreement Summary', fontsize=13, weight='bold')
plt.xticks(rotation=20, ha='right')
plt.gca().spines['top'].set_visible(False)
plt.gca().spines['right'].set_visible(False)

for i, v in enumerate(counts.values):
    plt.text(i, v + 0.05, str(v), ha='center', fontsize=9)

plt.tight_layout()
plt.show()
