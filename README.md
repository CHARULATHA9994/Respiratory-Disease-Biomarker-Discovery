# 🫁 Respiratory Disease Biomarker Discovery

**MSc Bioinformatics and Data Science — Research Project**

> Identifying shared molecular biomarkers across dust allergy, chronic rhinosinusitis, and allergic asthma using pan-transcriptomic and multi-omics analysis with machine learning.

---

## 📌 Project Objective

This project identifies **shared molecular biomarkers** across three allergic respiratory diseases using:
- Pan-transcriptomic differential expression analysis (R/limma)
- Functional enrichment and protein-protein interaction network analysis
- Multi-omics integration and machine learning-based biomarker validation

The goal is to uncover conserved disease mechanisms and candidate biomarkers that span multiple allergic conditions.

---

## 📂 Repository Structure

```
Respiratory-Disease-Biomarker-Discovery/
│
├── data/                          # Raw and processed expression data (not tracked — see below)
│
├── scripts/
│   ├── analysis.R                 # DEA pipeline: download, normalise, QC, limma
│   ├── bridge_step_analysis.py    # Bridge between R DEA outputs and Python pipeline
│   ├── cross_omics_integration.py # Multi-omics integration framework
│   ├── machine_learning_biomarker_evaluation.py  # ML training, validation, AUC
│   └── pan_transcriptomic_grouping.py            # Cross-dataset gene grouping
│
├── results/                       # Output CSVs: DEG lists, gene-level summaries
├── figures/                       # QC plots, PCA, PPI network visualisations
│
├── README.md
├── requirements.txt               # Python dependencies
└── Project_Thesis_FINAL.pdf       # Full MSc thesis
```

> ⚠️ Raw GEO expression data is not stored here. It is downloaded programmatically by `analysis.R` using GEOquery. See [How to Run](#-how-to-run) below.

---

## 🧬 Datasets Used

| GEO Accession | Disease | Platform |
|---|---|---|
| GSE9150 | Dust Mite Allergy vs Control | Affymetrix GPL570 |
| GSE23552 | Chronic Rhinosinusitis with Nasal Polyps (CRSwNP) vs Normal | Affymetrix HuEx-1_0-st |
| GSE65204 | Allergic Asthma vs Control | Agilent GPL14550 |

---

## 🔬 Analysis Pipeline

```
GEO Datasets (GSE9150, GSE23552, GSE65204)
        │
        ▼
[1] Download & Normalise (R: GEOquery, limma)
        │  Log2 transform, quantile normalisation
        ▼
[2] Quality Control
        │  Boxplots, PCA plots per dataset
        ▼
[3] Differential Expression Analysis (limma/eBayes)
        │  Strict list: adj.P.Val < 0.05
        │  Discovery list: adj.P.Val < 0.1 or p < 0.05 & |logFC| > 0.5
        ▼
[4] Probe → Gene Symbol Mapping
        │  GPL570 (GEOquery), HuEx DB (Bioconductor), GPL14550 (GEOquery)
        ▼
[5] Pan-Transcriptomic Integration (Python)
        │  Identify genes shared across all 3 datasets
        ▼
[6] Functional Enrichment
        │  Metascape, g:Profiler, Enrichr
        ▼
[7] PPI Network Analysis
        │  STRING database + Cytoscape + CytoHubba (hub gene identification)
        ▼
[8] Multi-Omics Integration (Python: cross_omics_integration.py)
        │  MetaboAnalyst 5.0, Reactome pathway mapping
        ▼
[9] Machine Learning Validation (Python: machine_learning_biomarker_evaluation.py)
        │  Random Forest + Logistic Regression
        │  5-fold cross-validation, AUC-ROC, Precision, Recall, F1
        ▼
[10] Hub Biomarker Identification
         CCL5, CCL20, IL7R, IDO1
```

---

## 🏆 Key Findings

### 4 Hub Biomarkers Identified

| Gene | Full Name | Role |
|---|---|---|
| **CCL5** | C-C Motif Chemokine Ligand 5 | T-cell and eosinophil recruitment |
| **CCL20** | C-C Motif Chemokine Ligand 20 | Dendritic cell attraction, mucosal immunity |
| **IL7R** | Interleukin-7 Receptor | T-cell development and survival |
| **IDO1** | Indoleamine 2,3-Dioxygenase 1 | Immune tolerance, tryptophan metabolism |

### Cross-Dataset Validation
- **CCL20** and **IDO1** validated in 2 out of 3 external datasets

### Machine Learning Performance

| Model | AUC-ROC |
|---|---|
| Logistic Regression | **0.816** |
| Random Forest | 0.776 |

Evaluation: 5-fold cross-validation with Precision, Recall, and F1-score metrics.

---

## 🛠 Tools & Technologies

| Category | Tools |
|---|---|
| Programming | Python 3, R |
| R Libraries | limma, GEOquery, Bioconductor |
| Python Libraries | Pandas, Scikit-learn, Matplotlib, NumPy |
| Network Analysis | STRING, Cytoscape, CytoHubba |
| Enrichment Analysis | Metascape, g:Profiler, Enrichr |
| Multi-omics | MetaboAnalyst 5.0, Reactome |

---

## ▶️ How to Run

### Prerequisites

**R (≥ 4.2)** with Bioconductor:
```r
if (!requireNamespace("BiocManager")) install.packages("BiocManager")
BiocManager::install(c("GEOquery", "limma", "huex10sttranscriptcluster.db"))
install.packages(c("ggplot2", "dplyr", "readr"))
```

**Python (≥ 3.9)**:
```bash
pip install -r requirements.txt
```

### Step-by-step

```bash
# Step 1: Clone the repository
git clone https://github.com/CHARULATHA994/Respiratory-Disease-Biomarker-Discovery.git
cd Respiratory-Disease-Biomarker-Discovery

# Step 2: Run DEA pipeline (downloads data automatically from GEO)
Rscript scripts/analysis.R
# Outputs: data/ and results/ CSVs, figures/ QC PDFs

# Step 3: Bridge R outputs to Python pipeline
python scripts/bridge_step_analysis.py

# Step 4: Pan-transcriptomic grouping across datasets
python scripts/pan_transcriptomic_grouping.py

# Step 5: Multi-omics integration
python scripts/cross_omics_integration.py

# Step 6: Machine learning evaluation
python scripts/machine_learning_biomarker_evaluation.py
```

> 💡 All GEO data is downloaded automatically — no manual data download required.

---

## 📄 Thesis

The full MSc thesis is available in this repository: [`Project_Thesis_FINAL.pdf`](./Project_Thesis_FINAL.pdf)

---

## 👩‍🔬 Author

**M. Charulatha**  
MSc Bioinformatics and Data Science  

---

## 📜 License

This project is for academic research purposes. Please cite appropriately if you use or build upon this work.

