# rice-drought-rnaseq
RNA-seq analysis of rice drought response using differential expression, functional enrichment, PPI, machine learning, and WGCNA
# Rice Drought Response RNA-seq Analysis

This project presents a computational RNA-seq analysis of the drought response in *Oryza sativa* (rice).

## Project Overview

The study compares rice samples under drought stress with well-watered control samples using RNA-seq data.

The analysis includes:

* Raw sequencing quality control
* Read trimming and filtering
* Post-trimming quality assessment
* MultiQC quality summary
* HISAT2 read alignment
* Samtools alignment processing
* FeatureCounts gene-level quantification
* PCA and expression visualization
* Differential expression analysis using DESeq2
* eggNOG functional annotation
* Gene Ontology (GO) enrichment
* KEGG pathway enrichment
* STRING protein-protein interaction analysis
* CytoHubba hub-gene analysis
* Random Forest-based candidate prioritization
* WGCNA co-expression analysis
* Integrated candidate-gene prioritization

## Dataset

Species: *Oryza sativa*

Experimental comparison:

**Drought vs Water**

Samples:

* SRR36685773
* SRR36685782
* SRR36685783
* SRR36685784
* SRR36685785
* SRR36685786

The raw sequencing data are not included in this repository. The NCBI SRA accession IDs are provided so that the dataset can be independently accessed.

## Main Results

Differential expression analysis identified:

* 2,991 differentially expressed genes
* 1,466 upregulated genes
* 1,525 downregulated genes

WGCNA identified a drought-associated module containing 614 genes.

The integrated analysis prioritized several candidate genes showing evidence from differential expression, functional/pathway analysis, machine learning, and/or co-expression analysis.

## Important Note

The computationally prioritized genes should be considered candidate drought-response genes rather than experimentally validated drought-tolerance genes.

## Tools

Major tools and packages used include:

* FastQC
* fastp
* MultiQC
* HISAT2
* Samtools
* FeatureCounts
* DESeq2
* eggNOG
* STRING
* Cytoscape / CytoHubba
* Random Forest
* WGCNA
* R

## Repository Structure

```text
rice-drought-rnaseq/
├── README.md
├── scripts/
├── bash/
├── qc/
├── results/
├── figures/
└── docs/
```
