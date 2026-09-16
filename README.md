# Rice Drought RNA-seq Analysis

## Project Overview

This project presents a computational RNA-seq analysis of **rice (*****Oryza sativa*****) under drought stress** using six paired-end RNA-seq samples: three Water controls and three Drought-treated samples.

The analysis covers quality control, read trimming, genome alignment, gene-level quantification, differential expression analysis, functional enrichment, protein interaction analysis, exploratory machine learning, co-expression analysis, and integrated candidate-gene prioritization.

The main objective was to identify genes and biological processes associated with the rice drought response and prioritize candidate genes for further investigation.

---

# Key Results at a Glance

| Analysis                       |                                                Result |
| ------------------------------ | ----------------------------------------------------: |
| RNA-seq samples                |                               6 (3 Water + 3 Drought) |
| Sequencing                     |                                   Paired-end, ~150 bp |
| HISAT2 alignment               |                                          93.43–93.83% |
| Annotated genes quantified     |                                                38,993 |
| Differentially expressed genes |                                                 2,991 |
| Upregulated genes              |                                                 1,466 |
| Downregulated genes            |                                                 1,525 |
| PCA                            |                                             PC1 ≈ 98% |
| GO background                  |                                          19,591 genes |
| Significant GO terms           |                                                    46 |
| Significant KEGG pathways      |                                                    28 |
| STRING/CytoHubba               |                              Top 10 MCC hubs analyzed |
| Random Forest                  |                       Exploratory gene prioritization |
| WGCNA                          | Exploratory drought-associated co-expression analysis |
| Final candidate genes          |                                                    15 |

---

# 1. Research Objective

The aim of this project was to investigate transcriptional changes in rice under drought stress and identify genes and biological pathways associated with the drought response.

The analysis was designed to move from raw RNA-seq data toward biologically relevant candidate genes through multiple complementary computational analyses.

The main comparison was:

**Drought vs Water**

---

# 2. Dataset

### Species

*Oryza sativa*

### Experimental groups

* Water: 3 biological replicates
* Drought: 3 biological replicates
* Total samples: 6

### SRA accessions

* SRR36685773
* SRR36685782
* SRR36685783
* SRR36685784
* SRR36685785
* SRR36685786

### Reference genome

**Oryza sativa IRGSP-1.0**

### Annotation

**Oryza_sativa.IRGSP-1.0.63.gtf**

### Protein annotation

The corresponding rice protein FASTA was used for downstream annotation.

---

# 3. Analysis Workflow

The analysis followed this general workflow:

**Raw RNA-seq reads → FastQC → fastp trimming → MultiQC → HISAT2 alignment → samtools processing → featureCounts → DESeq2 → Functional enrichment → STRING/CytoHubba → Random Forest → WGCNA → Candidate-gene prioritization**

---

# 4. Quality Control

Initial sequencing quality was assessed using **FastQC** and summarized using **MultiQC**.

The 12 FASTQ files represented six paired-end samples.

Reads were approximately 150 bp long, with GC content ranging from approximately 51–53%.

FastQC duplication levels ranged approximately from 22–35%. The R2 files generally showed higher duplication than R1 files.

The FastQC report contained some modules flagged as failed. These values represent **FastQC modules flagged by the program**, not the percentage of unusable reads.

Overall, the sequencing data were suitable for downstream RNA-seq processing.

### Post-trimming QC

FastQC and MultiQC were also used to evaluate the processed reads after trimming.

The post-trimming QC figures are available in:

`post trimming qc/`

---

# 5. Read Trimming

Adapter and low-quality sequence processing was performed using **fastp 1.3.6**.

The same paired-end processing strategy was applied to all six samples.

The processed reads were subsequently used for genome alignment.

---

# 6. Genome Alignment

Trimmed reads were aligned to the **Oryza sativa IRGSP-1.0 reference genome** using **HISAT2 2.2.3**.

Overall alignment rates were:

| Sample      | Condition | Overall alignment |
| ----------- | --------- | ----------------: |
| SRR36685773 | Water     |            93.74% |
| SRR36685782 | Water     |            93.73% |
| SRR36685783 | Water     |            93.77% |
| SRR36685784 | Drought   |            93.43% |
| SRR36685785 | Drought   |            93.58% |
| SRR36685786 | Drought   |            93.83% |

The alignment rates were consistently above 93% across all six samples.

---

# 7. Gene-Level Quantification

Gene-level read counting was performed using **featureCounts** with the rice GTF annotation.

The analysis used a paired-end, reverse-stranded configuration.

The resulting count matrix contained:

**38,993 annotated genes across 6 RNA-seq samples.**

The count matrix was used as input for DESeq2.

---

# 8. Differential Expression Analysis

Differential expression analysis was performed using **DESeq2**.

The experimental design was:

**~ condition**

with:

* Water as the reference condition
* Drought as the comparison condition

The contrast used was:

**Drought vs Water**

Genes were considered differentially expressed when they satisfied both criteria:

* adjusted p-value (`padj`) < 0.05
* absolute log2 fold change ≥ 1

Multiple-testing correction was applied through DESeq2-adjusted p-values.

### Differential expression results

A total of:

**2,991 DEGs**

were identified.

* **1,466 upregulated**
* **1,525 downregulated**

The complete DEG results are available in:

`results/DEGs/DEGs_all_with_geneIDs.csv`

Upregulated genes:

`results/DEGs/DEGs_upregulated_with_geneIDs.csv`

Downregulated genes:

`results/DEGs/DEGs_downregulated_with_geneIDs.csv`

---

# 9. Principal Component Analysis

Principal Component Analysis (PCA) was performed using DESeq2's `plotPCA()` function on variance-stabilized expression data (`vsd`), with samples grouped by experimental condition.

The analysis used the **top 500 most variable genes** selected by DESeq2's `plotPCA()` function.

The first principal component explained approximately **98% of the total variation**. Water and Drought samples showed clear separation in the PCA space, indicating a strong transcriptional difference between the two experimental conditions.

PCA was generated using:

```r
plotPCA(vsd, intgroup = "condition")
```

PCA figure:

`figures/DEGs/PCA.png.png`


---

# 10. DEG Visualization

A volcano plot was generated to visualize the relationship between statistical significance and fold change.

A heatmap was also generated to visualize expression patterns of differentially expressed genes.

Figures:

`figures/DEGs/Volcano_Plot.png.png`

`figures/DEGs/DEG_Heatmap.png.png`

---

# 11. Gene Ontology Enrichment

Gene Ontology (GO) enrichment analysis was performed in R using the clusterProfiler package separately for upregulated and downregulated genes.

The GO analysis used a background of:

**19,591 unique GO-annotated rice genes.**

Significant enrichment was evaluated using FDR-adjusted significance values.

A total of:

**46 significant GO terms**

were identified across Biological Process (BP), Cellular Component (CC), and Molecular Function (MF).

Examples of enriched biological processes included:

* response to oxidative stress
* defense response to fungus
* response to wounding
* response to water deprivation
* response to auxin
* response to abscisic acid
* cellular response to hypoxia
* response to light stimulus
* xenobiotic metabolic process

The GO results are available in:

`results/GO/`

GO figures are available in:

`figures/GO/`

---

# 12. KEGG Pathway Enrichment

KEGG pathway enrichment was performed in R using the clusterProfiler package using the same GO-annotated background of 19,591 genes.

A total of:

**28 significant KEGG pathways**

were identified.

Important enriched pathways included pathways related to:

* secondary metabolite biosynthesis
* starch and sucrose metabolism
* MAPK signaling in plants
* plant hormone signal transduction
* flavonoid biosynthesis
* diterpenoid biosynthesis
* alpha-linolenic acid metabolism
* nitrogen metabolism
* cutin, suberine and wax biosynthesis

The complete KEGG results are available in:

`results/KEGG/`

Figures:

`figures/KEGG/KEGG_Downregulated_Top10.png`

`figures/KEGG/KEGG_Upregulated_All_Significant.png`

---

# 13. Functional Candidate Prioritization

To move beyond a simple DEG list, genes were examined using their differential-expression results together with functional evidence from GO and KEGG analyses.

This allowed genes with stronger functional support within the drought-response analysis to be prioritized for further investigation.

The integrated candidate analysis combined:

* differential expression
* GO evidence
* KEGG pathway evidence

This should be considered **integrated transcriptomic candidate prioritization**, rather than experimental validation.

---

# 14. Final Candidate Genes

A final set of **15 candidate genes** was prioritized based on combined differential-expression and functional-enrichment evidence.

The complete candidate list is retained as a project result and can be found at:

`results/final_candidates/Final_Candidate_Genes.csv`

The corresponding integrated visualization is available at:

`figures/final_candidate_genes/Final_Candidate_Integration_NEW.png`

The final candidate table contains:

* gene ID
* preferred/annotation-derived name
* log2 fold change
* adjusted p-value
* GO evidence count
* KEGG evidence count

The candidate list is intended for **further biological investigation**, not as a set of experimentally validated drought biomarkers.

---

# 15. STRING and CytoHubba Network Analysis

STRING protein-protein interaction analysis was performed as an **additional network-level analysis** using the DEG-associated interaction network.

CytoHubba was then used with the MCC method to identify highly connected hub proteins.

The top 10 MCC-ranked hubs were examined and mapped to rice gene identifiers where possible.

Examples of identified hub genes included:

* BADH1 — Os04g0464000
* Q0DCC7_ORYSJ — Os06g0347100
* A0A0P0WAB9 — Os02g0764100
* A0A0P0WMZ2 — Os05g0447580
* 4CLL6 — Os01g0901600

The complete hub-gene result is available at:

`results/STRING_CytoHubba/hubgenes.csv`

Figures:

`figures/Hubgenes/STRING network.png`

`figures/Hubgenes/string_interactions.tsv_1_MCC_top10.png`

### Important interpretation

STRING/CytoHubba was treated as a **separate network analysis**.

Hub status was **not used as a direct criterion for the final 15 candidate-gene list**.

This keeps network connectivity and functional candidate prioritization as two complementary but separate analyses.

---

# 16. Random Forest Analysis

Random Forest was used as an **exploratory machine-learning approach for gene prioritization**.

Variance-stabilized expression data were used for the analysis.

Because the dataset contained only six samples and thousands of genes, feature reduction was performed before the exploratory modeling step.

An exploratory Leave-One-Out Cross-Validation (LOOCV) procedure classified all six samples correctly:

**6/6 samples classified correctly.**

However, this result should not be interpreted as an unbiased estimate of predictive performance.

Feature selection was performed before cross-validation, meaning information from the complete dataset influenced the selected features. Therefore, the Random Forest analysis is presented as **exploratory candidate-gene prioritization**, rather than as a validated predictive model.

The main Random Forest results are available in:

`results/RandomForest/RF_LOOCV_predictions.csv`

and:

`results/RandomForest/RF_top20_biomarkers_annotated.csv`

---

# 17. WGCNA Co-expression Analysis

Weighted Gene Co-expression Network Analysis (WGCNA) was performed as an additional exploratory analysis to investigate co-expression patterns associated with drought treatment.

The analysis was based on variance-stabilized expression data and focused on highly variable genes.

A drought-associated co-expression module was identified and its relationship with the drought condition was examined.

The WGCNA results are available in:

`results/WGCNA/`

including:

`results/WGCNA/WGCNA_Module2_Drought_Associated_Genes.csv`

`results/WGCNA/WGCNA_Module_Summary.csv`

### Important limitation

Only six samples were available for the analysis.

Therefore, the WGCNA results should be interpreted as an **exploratory co-expression analysis**, rather than robust network validation. The small sample size limits the reliability and generalizability of module-trait relationships.

---

# 18. Integrated Interpretation

The different analyses provide complementary evidence for understanding the rice drought response.

### Differential expression

Identified **2,991 genes** with significant expression changes between Drought and Water conditions.

### GO enrichment

Highlighted biological processes related to stress response, defense, water deprivation, hormones, oxidative stress, and other drought-associated functions.

### KEGG enrichment

Identified pathways involving plant signaling, secondary metabolism, carbohydrate metabolism, hormones, and other metabolic processes.

### STRING/CytoHubba

Provided a network-level view of highly connected proteins within the DEG-associated interaction network.

### Random Forest

Provided an exploratory data-driven approach for prioritizing genes based on their expression patterns.

### WGCNA

Provided an exploratory co-expression perspective on genes associated with drought treatment.

Together, these analyses provide multiple computational perspectives for prioritizing genes and pathways for future biological investigation.

---

# 19. Gene Annotation Note

Rice gene identifiers were retained as the primary identifiers throughout the analysis.

For example:

**Os12g0611100**

is a rice gene/locus identifier. Annotation-derived preferred names may originate from the functional annotation workflow and can represent transferred or orthology-based annotations.

Therefore, the **Os gene identifier should be treated as the primary gene identifier**, while preferred names are used as supporting annotation information.

---

# 20. Software and Resources

### Software

* FastQC 0.12.1
* fastp 1.3.6
* MultiQC 1.35
* HISAT2 2.2.3
* samtools 1.24
* featureCounts
* DESeq2
* WGCNA
* Random Forest
* STRING
* Cytoscape / CytoHubba

### Reference resources

* *Oryza sativa* IRGSP-1.0 reference genome
* IRGSP-1.0 gene annotation
* Gene Ontology
* KEGG
* eggNOG functional annotation
* STRING protein interaction database

---

# 21. Reproducibility

The repository contains the main analysis scripts and shell commands used throughout the RNA-seq workflow.

The raw sequencing files, SRA files, BAM files, and reference genome files are not included in the repository because of their large size.

The analysis can be followed from the provided scripts, processed result tables, figures, and documented reference information.

---

# 22. Limitations

Several limitations should be considered when interpreting the results:

* The study contains only six RNA-seq samples.
* Random Forest results are exploratory because feature selection was performed before LOOCV.
* WGCNA results are exploratory because of the small sample size.
* Computational candidate prioritization does not replace experimental validation.
* Experimental batch information was not available in the analyzed metadata, so batch effects could not be formally evaluated.
* Functional annotations may include transferred or orthology-based annotations and should therefore be interpreted together with the primary rice gene identifiers.

---

# 23. Summary

This project provides an end-to-end computational analysis of rice RNA-seq data under drought stress.

The analysis identified **2,991 differentially expressed genes**, including **1,466 upregulated and 1,525 downregulated genes**.

Functional enrichment identified **46 significant GO terms** and **28 significant KEGG pathways**, highlighting biological processes and pathways associated with stress response, signaling, metabolism, and drought-related responses.

Additional STRING/CytoHubba, Random Forest, and WGCNA analyses provided network, machine-learning, and co-expression perspectives.

Finally, an integrated DEG + GO + KEGG analysis produced a shortlist of **15 candidate genes** for further investigation.

The complete candidate list is available at:

`results/final_candidates/Final_Candidate_Genes.csv`

