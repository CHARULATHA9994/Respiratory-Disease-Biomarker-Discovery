# Respiratory Disease Biomarker Discovery
### MSc Bioinformatics and Data Science — Research Project

## What This Project Is About
Identified shared molecular biomarkers across three allergic 
respiratory diseases — dust allergy, chronic rhinosinusitis 
with nasal polyps (CRSwNP), and allergic asthma — using 
pan-transcriptomic and multi-omics analysis.

## Diseases Studied
- Dust Allergy (GSE9150)
- Chronic Rhinosinusitis with Nasal Polyps (GSE23552)
- Allergic Asthma (GSE65204)

## Key Findings
- Identified 4 hub biomarkers: CCL5, CCL20, IL7R, and IDO1
- CCL20 and IDO1 validated in 2 out of 3 external datasets
- Machine learning confirmed strong predictive power:
  - Random Forest AUC: 0.776
  - Logistic Regression AUC: 0.816

## Methods Used
- Pan-transcriptomic analysis using R (limma, GEOquery)
- Functional enrichment: Metascape, g:Profiler, Enrichr
- Protein-Protein Interaction network: STRING + Cytoscape
- Multi-omics integration using custom Python framework
- Machine Learning: Random Forest + Logistic Regression
- Evaluation: 5-fold cross-validation, AUC-ROC, Precision, Recall, F1-score

## Tools & Technologies
| Category | Tools Used |
|---|---|
| Programming | Python, R |
| Python Libraries | Pandas, Scikit-learn, Matplotlib |
| R Libraries | limma, GEOquery, Bioconductor |
| Network Analysis | STRING, Cytoscape, CytoHubba |
| Enrichment Analysis | Metascape, g:Profiler, Enrichr |
| Multi-omics | MetaboAnalyst 5.0, Reactome |

## Author
M. Charulatha
MSc Bioinformatics and Data Science
