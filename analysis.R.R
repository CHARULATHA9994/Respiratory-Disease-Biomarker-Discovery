# =============================================================================
# COMPLETE GEO MICROARRAY ANALYSIS PIPELINE
# Datasets: GSE9150 (Dust Mite Allergy vs Control)
#           GSE23552 (Asthma with Nasal Polyps vs Normal)
#           GSE65204 (Asthma vs Control)
# Steps: Download -> Normalize -> QC -> DEA -> Probe-to-Gene Mapping
# =============================================================================


# =============================================================================
# SECTION 1: GSE9150 - Dust Mite Allergy vs Control
# Download, Normalize, QC, and Differential Expression Analysis
# =============================================================================

# --- Step 0: Install and load required packages ---
if (!requireNamespace("GEOquery")) install.packages("GEOquery")
if (!requireNamespace("limma"))    install.packages("limma")
if (!requireNamespace("ggplot2"))  install.packages("ggplot2")

library(GEOquery)
library(limma)
library(ggplot2)

# --- Step 1: Download dataset from GEO ---
gset  <- getGEO("GSE9150", GSEMatrix = TRUE, getGPL = FALSE)[[1]]
expr  <- exprs(gset)      # Expression matrix (probes x samples)
pheno <- pData(gset)      # Phenotype/sample metadata

# --- Step 2: Define sample groups (Control vs DustMite) ---
pheno$Group <- ifelse(
  pheno$`characteristics_ch1` == "non allergic",
  "Control",
  "DustMite"
)
pheno$Group <- factor(pheno$Group, levels = c("Control", "DustMite"))

cat("\nGroup counts:\n")
print(table(pheno$Group, useNA = "ifany"))

# Save phenotype info to CSV
write.csv(pheno, "GSE9150_pheno.csv", row.names = TRUE)

# --- Step 3: Log2 transform expression data if not already transformed ---
qx   <- as.numeric(quantile(expr, c(0, .25, .5, .75, 1.0), na.rm = TRUE))
logC <- (qx[5] > 100) || (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2)
if (logC) {
  expr <- log2(expr + 1)
}

# Save normalized expression matrix
write.csv(expr, "GSE9150_normalized_expr.csv")

# --- Step 4: QC Plots (Boxplot + PCA) ---
pdf("GSE9150_QC.pdf", width = 10, height = 6)
par(mfrow = c(1, 2))

# Boxplot of normalized expression across samples
boxplot(expr, main = "Boxplot of Normalized Expression", las = 2, col = "lightblue")

# PCA plot to visualise sample clustering by group
pca <- prcomp(t(expr), scale. = TRUE)
plot(pca$x[,1], pca$x[,2],
     col  = pheno$Group,
     pch  = 19,
     xlab = "PC1", ylab = "PC2",
     main = "PCA Plot")
legend("topright", legend = levels(pheno$Group), col = 1:2, pch = 19)

dev.off()

# --- Step 5: Differential Expression Analysis (limma) ---
design <- model.matrix(~ Group, data = pheno)
fit    <- lmFit(expr, design)
fit    <- eBayes(fit)

# Extract full DEG table, sorted by p-value
deg <- topTable(
  fit,
  coef          = "GroupDustMite",
  number        = Inf,
  adjust.method = "BH",
  sort.by       = "P"
)

# Filter to significant DEGs: adj.P.Val < 0.05 AND |logFC| > 1
deg_filtered <- subset(deg, adj.P.Val < 0.05 & abs(logFC) > 1)

# Save results to CSV
write.csv(deg,          "GSE9150_DEG_FullList.csv")
write.csv(deg_filtered, "GSE9150_DEG_Significant.csv")

cat("\n✅ Analysis complete.\n")
cat("Full DEG list saved as: GSE9150_DEG_FullList.csv\n")
cat("Significant DEGs saved as: GSE9150_DEG_Significant.csv\n")
cat("QC plots saved as: GSE9150_QC.pdf\n")


# =============================================================================
# SECTION 2: GSE23552 - Asthma with Chronic Rhinosinusitis (CRS) vs Normal
# Download, Normalize, QC, and Differential Expression Analysis
# =============================================================================

# --- Load libraries (already loaded above, repeated for standalone use) ---
library(GEOquery)
library(limma)
library(ggplot2)

# --- Parameters ---
gse_id <- "GSE23552"

# --- Step 1: Download processed expression data ---
gset  <- getGEO(gse_id, GSEMatrix = TRUE, getGPL = FALSE)[[1]]

# --- Step 2: Extract phenotype data and define groups ---
pheno <- pData(gset)

# Extract disease state from the characteristics column
pheno$Group <- sub("disease state: ", "", pheno$`characteristics_ch1`)
pheno$Group <- make.names(pheno$Group)   # Replace spaces with dots for R compatibility

# Keep only Normal and CRS (Asthmatic Chronic Rhinosinusitis with Nasal Polyps) samples
pheno <- pheno[pheno$Group %in% c("Normal", "Asthmatic.Chronic.Rhinosinusitis.with.Nasal.Polyps"), ]
gset  <- gset[, rownames(pheno)]

# --- Step 3: Extract and optionally log2-transform expression data ---
expr <- exprs(gset)
qx   <- as.numeric(quantile(expr, c(0, 0.25, 0.5, 0.75, 1.0), na.rm = TRUE))
if (qx[5] > 100 || (qx[5] - qx[1] > 50 && qx[2] > 0)) {
  expr <- log2(expr + 1)
}

# --- Step 4: Quantile normalisation across arrays ---
expr <- normalizeBetweenArrays(expr)

# Save normalised expression and phenotype to CSV
write.csv(expr,  paste0(gse_id, "_normalized_expr.csv"))
write.csv(pheno, paste0(gse_id, "_pheno.csv"))

# --- Step 5: QC Plots ---
pdf(paste0(gse_id, "_QC.pdf"))

# Boxplot
boxplot(expr, las = 2,
        main = "Boxplot of normalized expression",
        col  = "lightblue")

# PCA using ggplot2
pca    <- prcomp(t(expr), scale. = TRUE)
pca_df <- data.frame(PC1 = pca$x[, 1], PC2 = pca$x[, 2], Group = pheno$Group)
print(
  ggplot(pca_df, aes(PC1, PC2, color = Group)) +
    geom_point(size = 3) +
    theme_minimal() +
    ggtitle("PCA Plot")
)

dev.off()

# --- Step 6: Differential Expression Analysis (limma with contrasts) ---
pheno$Group <- factor(
  pheno$Group,
  levels = c("Normal", "Asthmatic.Chronic.Rhinosinusitis.with.Nasal.Polyps")
)
design      <- model.matrix(~0 + pheno$Group)
colnames(design) <- levels(pheno$Group)

# Define contrast: CRS vs Normal
contrast <- makeContrasts(
  Asthmatic.Chronic.Rhinosinusitis.with.Nasal.Polyps - Normal,
  levels = design
)

fit  <- lmFit(expr, design)
fit2 <- contrasts.fit(fit, contrast)
fit2 <- eBayes(fit2)

# Extract DEG table
deg           <- topTable(fit2, coef = 1, number = Inf, adjust.method = "BH")
deg$ProbeID   <- rownames(deg)

# Save DEG table
write.csv(deg, paste0(gse_id, "_DEG_Control_vs_CRS.csv"), row.names = FALSE)
cat("✅ Done - QC PDF, normalized expr CSV, phenotype CSV, and DEG CSV saved for", gse_id, "\n")


# =============================================================================
# SECTION 3: GSE65204 - Asthma vs Control
# Download, Normalize, QC, and Differential Expression Analysis
# =============================================================================

# --- Load required packages ---
if (!requireNamespace("GEOquery")) install.packages("GEOquery")
if (!requireNamespace("limma"))    install.packages("limma")
if (!requireNamespace("ggplot2"))  install.packages("ggplot2")

library(GEOquery)
library(limma)
library(ggplot2)

# --- Step 1: Download GSE65204 dataset ---
gset  <- getGEO("GSE65204", GSEMatrix = TRUE, getGPL = FALSE)[[1]]
expr  <- exprs(gset)
pheno <- pData(gset)

# --- Step 2: Extract group information (Asthma vs Control) ---
pheno$Group <- ifelse(grepl("TRUE", pheno$characteristics_ch1), "Asthma", "Control")
pheno$Group <- factor(pheno$Group, levels = c("Control", "Asthma"))

# Save phenotype table
write.csv(pheno, "GSE65204_pheno.csv", row.names = TRUE)

# --- Step 3: Log2 transform if necessary ---
qx   <- as.numeric(quantile(expr, c(0, 0.25, 0.5, 0.75, 1.0), na.rm = TRUE))
logC <- (qx[5] > 100) || (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2)
if (logC) {
  expr <- log2(expr + 1)
}

# Save normalised expression
write.csv(expr, "GSE65204_normalized_expr.csv")

# --- Step 4: QC Plots ---
pdf("GSE65204_QC.pdf", width = 10, height = 6)
par(mfrow = c(1, 2))

boxplot(expr, main = "Boxplot of Normalized Expression", las = 2)

pca <- prcomp(t(expr), scale. = TRUE)
plot(pca$x[,1], pca$x[,2],
     col  = pheno$Group,
     pch  = 19,
     xlab = "PC1", ylab = "PC2",
     main = "PCA Plot")
legend("topright", legend = levels(pheno$Group), col = 1:2, pch = 19)

dev.off()

# --- Step 5: Differential Expression Analysis (limma) ---
design <- model.matrix(~ Group, data = pheno)
fit    <- lmFit(expr, design)
fit    <- eBayes(fit)
deg    <- topTable(fit, coef = "GroupAsthma", number = Inf, adjust.method = "BH")

# Save DEG results
write.csv(deg, "GSE65204_DEG_Control_vs_Asthma.csv")

cat("✅ Analysis complete. Files saved:\n")
cat("  - GSE65204_pheno.csv\n")
cat("  - GSE65204_normalized_expr.csv\n")
cat("  - GSE65204_QC.pdf\n")
cat("  - GSE65204_DEG_Control_vs_Asthma.csv\n")


# =============================================================================
# SECTION 4: GSE9150 - Probe ID to Gene Symbol Mapping
# Map Affymetrix GPL570 probe IDs -> Gene Symbols, then filter DEG lists
# =============================================================================

# --- Install required packages ---
if (!requireNamespace("GEOquery", quietly = TRUE)) install.packages("GEOquery")
if (!requireNamespace("dplyr",    quietly = TRUE)) install.packages("dplyr")
if (!requireNamespace("readr",    quietly = TRUE)) install.packages("readr")

library(GEOquery)
library(dplyr)
library(readr)

# --- Step 1: Load DEG file produced in Section 1 ---
deg_file <- "GSE9150_DEG_Significant.csv"   # Change path if needed
deg      <- read_csv(deg_file)

cat("First rows of DEG file:\n")
print(head(deg))

# --- Step 2: Download GPL570 annotation from GEO ---
gpl       <- getGEO("GPL570", AnnotGPL = TRUE)
gpl_table <- Table(gpl)

cat("Columns in GPL table:\n")
print(colnames(gpl_table))

# --- Step 3: Extract Probe ID and Gene Symbol columns ---
gene_symbol_col <- grep("Gene.?Symbol", colnames(gpl_table),
                        ignore.case = TRUE, value = TRUE)[1]

annot <- gpl_table %>%
  select(ID, all_of(gene_symbol_col)) %>%
  rename(ProbeID = ID, GeneSymbol = all_of(gene_symbol_col)) %>%
  mutate(GeneSymbol = ifelse(GeneSymbol == "" | GeneSymbol == "---", NA, GeneSymbol))

# Remove rows where gene symbol is missing
annot <- annot %>% filter(!is.na(GeneSymbol))

# --- Step 4: Merge annotation with DEG data ---
# Ensure the first column in DEG data is labelled "ProbeID"
colnames(deg)[1] <- "ProbeID"
deg_annot <- merge(deg, annot, by = "ProbeID")

# --- Step 5: Collapse multiple probes per gene (average/min statistics) ---
deg_gene <- deg_annot %>%
  group_by(GeneSymbol) %>%
  summarise(
    logFC    = mean(logFC,    na.rm = TRUE),
    AveExpr  = mean(AveExpr,  na.rm = TRUE),
    t        = mean(t,        na.rm = TRUE),
    P.Value  = min(P.Value,   na.rm = TRUE),
    adj.P.Val = min(adj.P.Val, na.rm = TRUE),
    B        = mean(B,        na.rm = TRUE)
  ) %>%
  ungroup()

# --- Step 6: Create Strict and Discovery gene lists ---
# Strict: adj.P.Val < 0.05 only
deg_strict <- deg_gene %>% filter(adj.P.Val < 0.05)

# Discovery: adj.P.Val < 0.1 OR (P.Value < 0.05 AND |logFC| > 0.5)
deg_discovery <- deg_gene %>%
  filter(adj.P.Val < 0.1 | (P.Value < 0.05 & abs(logFC) > 0.5))

# --- Step 7: Save outputs ---
write_csv(deg_gene,      "GSE9150_DEG_GeneLevel.csv")
write_csv(deg_strict,    "GSE9150_DEG_Strict.csv")
write_csv(deg_discovery, "GSE9150_DEG_Discovery.csv")

cat("✅ Mapping, collapsing, and filtering completed.\n")
cat("Files saved:\n")
cat("  - GSE9150_DEG_GeneLevel.csv (collapsed)\n")
cat("  - GSE9150_DEG_Strict.csv (adj.P.Val < 0.05)\n")
cat("  - GSE9150_DEG_Discovery.csv (adj.P.Val<0.1 or p<0.05 & |logFC|>0.5)\n")


# =============================================================================
# SECTION 5: GSE23552 - Probe ID to Gene Symbol Mapping
# Map Affymetrix HuEx-1_0-st (huex10sttranscriptcluster) probes -> Gene Symbols
# =============================================================================

# --- Install required packages ---
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

packages <- c("tidyverse", "annotate", "huex10sttranscriptcluster.db")
for (pkg in packages) {
  if (!requireNamespace(pkg, quietly = TRUE))
    BiocManager::install(pkg)
}

library(tidyverse)
library(annotate)
library(huex10sttranscriptcluster.db)

# --- Step 1: Load DEG file produced in Section 2 ---
deg <- read.csv("GSE23552_DEG_Control_vs_CRS.csv")

cat("Column names in DEG file:\n")
print(colnames(deg))

# Confirm ProbeID column exists
if (!"ProbeID" %in% colnames(deg)) {
  stop("Your file must have a column named 'ProbeID'.")
}

# --- Step 2: Map ProbeID -> Gene Symbol using Bioconductor annotation DB ---
deg$GeneSymbol <- mapIds(
  huex10sttranscriptcluster.db,
  keys     = as.character(deg$ProbeID),
  column   = "SYMBOL",
  keytype  = "PROBEID",
  multivals = "first"
)

# Remove entries without a gene symbol
deg <- deg %>% filter(!is.na(GeneSymbol))

# --- Step 3: Collapse duplicate probes per gene (average logFC) ---
deg_genelevel <- deg %>%
  group_by(GeneSymbol) %>%
  summarise(
    logFC     = mean(logFC,     na.rm = TRUE),
    AveExpr   = mean(AveExpr,   na.rm = TRUE),
    t         = mean(t,         na.rm = TRUE),
    P.Value   = min(P.Value,    na.rm = TRUE),
    adj.P.Val = min(adj.P.Val,  na.rm = TRUE),
    B         = mean(B,         na.rm = TRUE),
    .groups   = "drop"
  )

# --- Step 4: Create Strict and Discovery gene lists ---
deg_strict <- deg_genelevel %>%
  filter(adj.P.Val < 0.05)

deg_discovery <- deg_genelevel %>%
  filter(adj.P.Val < 0.1 | (P.Value < 0.05 & abs(logFC) > 0.5))

# --- Step 5: Save outputs ---
write.csv(deg_genelevel,  "GSE23552_genelevel.csv",  row.names = FALSE)
write.csv(deg_strict,     "GSE23552_strict.csv",     row.names = FALSE)
write.csv(deg_discovery,  "GSE23552_discovery.csv",  row.names = FALSE)

cat("✅ Files saved: GSE23552_genelevel.csv, GSE23552_strict.csv, GSE23552_discovery.csv\n")


# =============================================================================
# SECTION 6: GSE65204 - Probe ID to Gene Symbol Mapping
# Map Agilent GPL14550 probe IDs -> Gene Symbols, then filter DEG lists
# =============================================================================

# Script: Map Agilent GPL14550 probe IDs to Gene Symbols for GSE65204

# --- Install and load required packages ---
if (!requireNamespace("GEOquery", quietly = TRUE)) install.packages("GEOquery")
if (!requireNamespace("dplyr",    quietly = TRUE)) install.packages("dplyr")
if (!requireNamespace("readr",    quietly = TRUE)) install.packages("readr")

library(GEOquery)
library(dplyr)
library(readr)

# --- Step 1: Load DEG file produced in Section 3 ---
deg_file <- "GSE65204_DEG_Control_vs_Asthma.csv"   # Update path if needed
deg      <- read_csv(deg_file)

cat("First rows of DEG file:\n")
print(head(deg))

# --- Step 2: Download GPL14550 (Agilent) annotation from GEO ---
gpl       <- getGEO("GPL14550", AnnotGPL = TRUE)
gpl_table <- Table(gpl)

# Convert to tibble for safe dplyr operations
gpl_table <- as_tibble(gpl_table)

cat("Columns in GPL table:\n")
print(colnames(gpl_table))

# --- Step 3: Extract Probe ID and Gene Symbol columns ---
gene_symbol_col <- grep("Gene.?Symbol", colnames(gpl_table),
                        ignore.case = TRUE, value = TRUE)[1]

annot <- gpl_table %>%
  dplyr::select(ID, all_of(gene_symbol_col)) %>%
  dplyr::rename(ProbeID = ID, GeneSymbol = !!gene_symbol_col) %>%
  mutate(GeneSymbol = ifelse(GeneSymbol == "" | GeneSymbol == "---", NA, GeneSymbol))

# Remove probes with no gene symbol
annot <- annot %>% filter(!is.na(GeneSymbol))

# --- Step 4: Merge annotation with DEG data ---
colnames(deg)[1] <- "ProbeID"
deg_annot <- merge(deg, annot, by = "ProbeID")

# --- Step 5: Collapse multiple probes per gene ---
deg_gene <- deg_annot %>%
  group_by(GeneSymbol) %>%
  summarise(
    logFC     = mean(logFC,     na.rm = TRUE),
    AveExpr   = mean(AveExpr,   na.rm = TRUE),
    t         = mean(t,         na.rm = TRUE),
    P.Value   = min(P.Value,    na.rm = TRUE),
    adj.P.Val = min(adj.P.Val,  na.rm = TRUE),
    B         = mean(B,         na.rm = TRUE),
    .groups   = "drop"
  )

# --- Step 6: Create Strict and Discovery gene lists ---
deg_strict <- deg_gene %>% filter(adj.P.Val < 0.05)

deg_discovery <- deg_gene %>%
  filter(adj.P.Val < 0.1 | (P.Value < 0.05 & abs(logFC) > 0.5))

# --- Step 7: Save outputs ---
write_csv(deg_gene,      "GSE65204_DEG_GeneLevel.csv")
write_csv(deg_strict,    "GSE65204_DEG_Strict.csv")
write_csv(deg_discovery, "GSE65204_DEG_Discovery.csv")

cat("✅ Mapping, collapsing, and filtering completed.\n")
cat("Files saved:\n")
cat("  - GSE65204_DEG_GeneLevel.csv\n")
cat("  - GSE65204_DEG_Strict.csv\n")
cat("  - GSE65204_DEG_Discovery.csv\n")
