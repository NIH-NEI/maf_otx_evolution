# maf_otx_evolution

Source code and analysis workflows accompanying the manuscript on the **sequence and expression evolution of the MAF and OTX gene families across Metazoa**, with a particular focus on their evolution and roles in retinal development and function.

---

## Overview

This repository contains the computational workflows, scripts, notebooks, and configuration files used to generate the analyses presented in the accompanying manuscript.

The project was developed primarily using **Snakemake**, together with **Python** and **R** scripts, and includes analyses spanning comparative genomics, phylogenetics, sequence evolution, synteny, bulk RNA-seq, single-cell RNA-seq, and protein domain coevolution.

The repository is organized as follows:

- **`workflow/`** contains the Snakemake workflow files (`.smk`).
- **`workflow/scripts/`** and its subdirectories contain the Python and R scripts called by the workflows.
- Several analyses and visualizations are implemented as Jupyter notebooks (`.ipynb`).

Unlike a software package, this repository represents the complete research workflow used for the manuscript. Although we have attempted to document the analyses as thoroughly as possible, the workflows are **not fully automated**. Many scripts contain hardcoded paths and system-specific configurations corresponding to the NIH Biowulf computing environment on which the analyses were originally performed.

Researchers wishing to reproduce or adapt these analyses should expect to modify file paths, software locations, and configuration parameters for their own computing environment.

---

# Repository Structure

## A. Proteome preparation

Preparation of proteomes downloaded from NCBI RefSeq, Ensembl, and custom genome annotations.

- `A0_prepare_species_sample_table.smk`
- `A1a_prepare_proteome_from_ncbi_refseq.smk`
- `A1b_prepare_proteome_from_ensembl.smk`
- `A1c_prepare_proteome_from_custom_annotation.smk`

---

## B. Identification and curation of MAF and OTX genes

Identification and curation of annotated and de novo predicted MAF and OTX proteins using OrthoFinder, BLAST, InterProScan, and protein language model embeddings.

- `B1_prepare_inputs_for_orthofinder.smk`
- `B2_extract_annotated_maf_otx_from_orthofinder.smk`
- `B3a_create_blastdb.smk`
- `B3b_run_blast.smk`
- `B4a_start_interproscan.smk`
- `B4b_consolidate_interproscan.smk`
- `B5a_combine_interpro_blast.smk`
- `B5b_prepare_mini_orthofinder.smk`
- `B5c_run_mini_orthofinder.smk`
- `B5d_annotate_denovo.smk`
- `B6_quantify_samples_to_transdecoder_cds.smk`
- `B7_make_first_grand_list.smk`
- `B8_present_gene_finds.smk`
- `B9a_generate_esm3_300M_embeddings.ipynb`
- `B9b_plot_esm3_300M_embeddings.ipynb`

---

## C. Phylogenetic analyses

Construction, refinement, and evaluation of MAF and OTX phylogenies.

- `C0_subset_grand_list.smk`
- `C1_hmm_correct_maf_otx.smk`
- `C2_draw_ortholog_tree_alignment.smk`
- `C3_export_nexus.smk`
- `C3_justify_tree_constraints.smk`

---

## D. Synteny analyses

Comparative synteny analyses in jawless and jawed vertebrates.

- `D1_start_orthofinder_for_jawless_synteny.smk`
- `D2_extract_cyclostomes_annotated_maf_otx_from_orthofinder.smk`
- `D3_jawless_synteny.smk`
- `D4_jawed_synteny.smk`

---

## E. Sequence evolution analyses

HyPhy-based analyses of coding sequence evolution.

- `E1_get_CDS_for_grand_list.smk`
- `E2_prepare_hyphy_inputs.smk`
- `E3_run_hyphy_absrel.smk`

---

## F. Bulk RNA-seq analyses

Orthology assignment, expression quantification, normalization, visualization, and GEO submission preparation.

- `F1_prepare_inputs_for_orthofinder_for_expression.smk`
- `F2_create_ortho_table_for_bulkrna.smk`
- `F3_quantify_bulkrna.smk`
- `F4_normalize_bulk_rna.smk`
- `F5_visualize_bulk_expr.smk`
- `F6_prepare_GEO_submission.smk`

---

## G. Single-cell RNA-seq analyses

Visualization and analysis of single-cell RNA-seq datasets.

- `G1_hahn_scrna_visualization.smk`

---

## H. Reanalysis of the NRL upstream regulatory element

- `H1_revisit_nrl_re.smk`

---

## I. Protein domain coevolution analyses

Identification of protein domains and analysis of domain coevolution within the MAF and OTX gene families.

- `I1_start_interpro_grand_list.smk`
- `I2_consolidate_interpro_grand_list.smk`
- `I3_mark_maf_otx_domains.smk`
- `I4_analyze_domain_coevolution.smk`

---

# Software Requirements

The workflows make use of a variety of third-party software packages. Major dependencies include:

- Snakemake
- Python (≥3.10 recommended)
- R (≥4.2 recommended)
- OrthoFinder
- BLAST+
- InterProScan
- HMMER
- HyPhy
- MAFFT
- IQ-TREE
- TransDecoder
- Salmon
- samtools
- bedtools

Additional Python and R package dependencies are documented within the individual scripts and notebooks.

---

# Data Availability

The analyses use publicly available genomic, proteomic, and transcriptomic datasets obtained from sources including:

- NCBI RefSeq
- Ensembl
- GEO
- Additional species-specific genome annotations

Because many datasets were downloaded and curated over an extended period, users attempting to reproduce the analyses may need to update download locations or accession numbers.

---

# Reproducibility

This repository is intended to accompany the manuscript and provide complete access to the computational analyses.

The workflows should be regarded as **research workflows** rather than a turnkey software package. Reproducing the analyses will generally require:

1. Installing all third-party software dependencies.
2. Downloading the required genome, proteome, and RNA-seq datasets.
3. Updating hardcoded file paths.
4. Modifying configuration files for the local computing environment.
5. Adapting resource requirements for the available computational infrastructure.

Many workflows were originally executed on the NIH Biowulf high-performance computing cluster and therefore include cluster-specific settings that may require modification before execution elsewhere.

---

# Citation

If you use this repository or build upon these analyses, please cite:

> Soumitra Pal, Noor D. White Carreiro, Zachary Batz, Duffie Judy, Anand Swaroop.  *Evolution of MAF and OTX families during the emergence and stabilization of rod photoreceptors across Metazoans.* *(Manuscript under review.)*

The citation will be updated upon publication.

---

# Contact

For questions regarding the analyses or repository, please open a GitHub Issue or contact soumitra.pal@nih.gov.

---

# License

Unless otherwise specified, the source code in this repository is released under the MIT License.

The datasets used in this study remain subject to the licenses and terms of their original providers.
