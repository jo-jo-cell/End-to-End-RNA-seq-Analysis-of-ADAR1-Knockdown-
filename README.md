# End-to-End-RNA-seq-Analysis-of-ADAR1-Knockdown
Transcriptomic re-analysis of a public RNA-seq dataset (n = 6) comparing ADAR1 knockdown to scrambled control. Covers QC, DESeq2 differential expression, volcano and heatmap visualisation, and GO/GSEA/KEGG functional enrichment.

## Aim
To characterise transcriptional changes following ADAR1 depletion and identify
the biological processes and pathways affected.

## Workflow
1. Count matrix import and low-count filtering (14,752 genes tested)
2. Differential expression with **DESeq2** (design: ~ treatment), LFC shrinkage
3. QC: PCA, p-value histogram, MA plot
4. Visualisation: volcano plot (**EnhancedVolcano**), heatmap of top DEGs (VST, row z-score)
5. Functional analysis with **clusterProfiler**:
   - GO BP over-representation of upregulated genes
   - GO BP GSEA (gseGO) on ranked genes
   - KEGG over-representation of significant DEGs

## Key results
- **QC:** PC1 (36%) and PC2 (19%) show clear structure between conditions.
- **Top DEGs:** include UGT2B17, UGT2B11, AZGP1, DDC, OPRK1, DOC2A, GRB10, NDUFB9.
- **GO (up in siADAR1):** regulation of hormone levels, axonogenesis, estrogen metabolic process.
- **GSEA GO BP:**
  - Suppressed: myoblast/cell–cell fusion, myotube differentiation, oxidative
    phosphorylation, aerobic respiration, cytoplasmic translation, striated muscle contraction
  - Activated: ferroptosis and its regulation, endosomal transport, cholesterol
    biosynthesis, protein dephosphorylation
- **KEGG:** ascorbate/aldarate metabolism, retinol metabolism, bile secretion,
  pentose/glucuronate interconversions (largely driven by UGT family genes),
  dilated cardiomyopathy.

## Figures
| | |
|---|---|
| PCA | ![](results/figures/pca.png) |
| MA plot | ![](results/figures/ma_plot.png) |
| P-value distribution | ![](results/figures/pvalue_hist.png) |
| Volcano | ![](results/figures/volcano.png) |
| Top DEG heatmap | ![](results/figures/heatmap.png) |
| GO BP (ORA) | ![](results/figures/go_bp_up.png) |
| GSEA GO BP | ![](results/figures/gsea_gobp.png) |
| KEGG | ![](results/figures/kegg.png) |

## Limitations
- Small sample size (n = 3 per group).
- Some top hits (e.g. PCDH11Y, a Y-linked gene; NPIPA7/NPIPB6, multi-copy
  genes) may reflect sample sex/donor differences or multi-mapping rather than
  ADAR1 biology; interpret with caution.
- Findings are descriptive therefore no experimental validation.

## Tools
R, DESeq2, EnhancedVolcano, clusterProfiler, pheatmap/ggplot2

## Author
Moses Mwangi Muiruri, Biotechnology, Technical University of Kenya
