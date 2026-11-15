# Transcriptomic re-analysis of ADAR1 knockdown in LNCaP prostate cancer cells

An RNA-seq re-analysis comparing ADAR1 knockdown (siADAR1) with scrambled siRNA
control in LNCaP cells, covering quality control, differential expression, and
GO, GSEA and KEGG functional enrichment.

## Background
ADAR1 is an RNA-editing enzyme that converts adenosine to inosine in
double-stranded RNA. This editing helps cells distinguish their own RNA from
viral RNA and limits innate immune sensing of endogenous transcripts. ADAR1 is
also implicated in several cancers, including prostate cancer. LNCaP is an
androgen-sensitive prostate cancer cell line, which makes it a useful model for
asking how ADAR1 depletion affects hormone-responsive transcriptional programmes.

## Aim
To characterise the transcriptional changes that follow ADAR1 knockdown in LNCaP
cells and identify the biological processes and pathways affected.

## Data source
- Public RNA-seq dataset: GSE327582, SRA runs SRR38030412, SRR38030416, SRR38030420, SRR38030424, SRR38030428 and SRR38030432 (n = 6)
- Original study: [citation]
- Design: three Scramble control samples and three siADAR1 samples

## Methods
Reads were quantified to gene-level counts ([aligner/quantifier and version]).
Low-count genes were filtered, leaving 14,752 genes for testing. Differential
expression was assessed with DESeq2 ([version]) using treatment as the design
factor, with log fold-change shrinkage applied for visualisation. Quality control
comprised principal component analysis, a p-value distribution check and an MA
plot. Results were visualised with a volcano plot (EnhancedVolcano) and a heatmap
of top differentially expressed genes (variance-stabilised counts, row z-scores).
Functional analysis used clusterProfiler ([version]): Gene Ontology Biological
Process over-representation of upregulated genes, GO Biological Process gene set
enrichment analysis (GSEA) on ranked genes, and KEGG over-representation of
significant genes. Significance was defined as adjusted p-value < [0.05] and
|log2 fold change| > [value].

## Results

### Quality control
Principal component analysis showed that PC1 explained 36% of variance and PC2
explained 19%, with samples separating by treatment. [Add one line on the
p-value distribution and MA plot, e.g. whether the p-value histogram shows an
enrichment of low p-values.]

### Differential expression
The most significant genes included UGT2B17, UGT2B11, AZGP1, DDC, OPRK1, DOC2A,
GRB10 and NDUFB9. [Add the number of significantly up and downregulated genes.]

### Functional enrichment
- **GO Biological Process (upregulated in siADAR1):** regulation of hormone
  levels, axonogenesis and estrogen metabolic process.
- **GSEA, suppressed:** myoblast and cell-cell fusion, myotube differentiation,
  oxidative phosphorylation, aerobic respiration, cytoplasmic translation and
  striated muscle contraction.
- **GSEA, activated:** ferroptosis and its regulation, endosomal transport,
  cholesterol biosynthesis and protein dephosphorylation.
- **KEGG:** ascorbate and aldarate metabolism, retinol metabolism, bile
  secretion, pentose and glucuronate interconversions, and dilated
  cardiomyopathy.

## Interpretation
Several of the top genes (UGT2B15, UGT2B17, UGT2B11) encode enzymes that
glucuronidate androgens, which is consistent with the enrichment of hormone and
steroid metabolism terms in an androgen-sensitive cell line. This suggests, as a
hypothesis, that ADAR1 loss may alter steroid handling in LNCaP cells. The
suppression of oxidative phosphorylation, aerobic respiration and translation
terms is a second theme that may point to a reduction in cellular metabolic and
biosynthetic activity. The activation of ferroptosis-related terms is also
notable given the role of lipid and cholesterol metabolism in that process.

## Limitations
- Small sample size (n = 3 per group). The findings are descriptive and have not
  been experimentally validated.
## Author
Moses Mwangi Muiruri
