# Genetic Architecture of Cutaneous Melanoma Across the Variant Spectrum

This repository contains main analysis code supporting the study:

**“The Genetic Architecture of Cutaneous Melanoma Across the Variant Spectrum”**

We systematically characterized the genetic architecture of cutaneous melanoma (CM) across the variant allele-frequency spectrum and functional genomic landscape, encompassing both coding and non-coding regions. We conducted a large-scale GWAS meta-analysis of 1,398,034 individuals. Spatially resolved mapping of heritability revealed that heritability enrichment in melanocytes is context-dependent, shaped by cellular functional states and their spatial organization within the tumor–immune microenvironment, where T cell infiltration and interactions contribute. We extended rare-variant association analyses beyond protein-coding genes to include noncoding RNAs and epidermal melanocyte–specific cCREs, and validated these findings across multiple independent cohorts (GEL, AoU, MGBB). We further demonstrate the contributions of rare variants to trait heritability, genetic risk prediction, and clinical relevance, including associations with earlier age at diagnosis. Finally, we constructed a biological network that contextualizes interactions among CM-associated genes prioritized from GWAS and rare-variant analyses, providing a systems-level view of melanoma biology that integrates genetic signals across the variant spectrum.

---

## 🧭Overview

This study integrates multiple analytical components to provide a unified view of the genetic architecture of cutaneous melanoma:

- **Part 1:** GWAS meta-analysis and post-GWAS characterization
- **Part 2:** Spatially resolved mapping of trait heritability within the tumor–immune microenvironment  
- **Part 3:** Identification and validation of rare-variant associations across protein-coding genes, noncoding RNAs, and cCREs  
- **Part 4:** Integrative analysis of common and rare variants to delineate their joint, as well as distinct and complementary, contributions to CM genetic risk and underlying biology

<p align="center">
  <img src="figures/Study.Overview.png" width="750">
</p>

<p align="center">
  <b>Figure 1.</b> Study Overview.
</p>

---

## ⚙️Execution Workflow

Within each analysis folder, scripts are organized and named according to their execution order.

---

## 🖥️Environment

Analyses were performed across:
- HPC clusters (SLURM-based)
- UK Biobank Research Analysis Platform (UKB-RAP)  
- Genomics England Research Environment (GEL)

*Paths in scripts may require adaptation.*

---

## 📊Data Availability

- UK Biobank data used in this study are available to approved researchers for health-related research through application and can be accessed via the cloud-based Research Analysis Platform (UKB-RAP).  

- Data from the National Genomic Research Library (NGRL) used in this study are available within the secure Genomics England Research Environment. Data used in this research include Aggregated Variant Calls (AggV2) derived from release 10 of the 100,000 Genomes Project, comprising 78,195 germline genomes, together with associated clinical and phenotype data. Relevant result files generated in this study are available within the Genomics England Research Environment.  

- Exome-wide rare-variant association summary statistics for cutaneous melanoma (PheCode 172.1) from the All of Us and Mass General Brigham Biobank cohorts are publicly available at https://hugeamp.org:8000/research.html?pageid=600_traits_app_home.  

- Summary statistics from the GWAS meta-analysis of cutaneous melanoma, functional annotations for fine-mapped candidate causal variants, and PRS-CS weight files for common-variant polygenic risk score construction have been deposited in the Zenodo repository (https://doi.org/10.5281/zenodo.19562796).  

---

## 📄 Manuscript

**The Genetic Architecture of Cutaneous Melanoma Across the Variant Spectrum**

*Manuscript in preparation.*

---

## 📬Contact

Shiyu Zhang (Shelley)  
Email: shiyuzhang0522@gmail.com 

---