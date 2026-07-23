# WARNING: This Snakefile is not intended to be run as-is.
# It depends on OrthoFinder results with hard-coded output paths and/or filenames,
# which are specific to a particular run and may not be valid for your dataset.

rule all:
    input:
        "exports/curated_genes/draft_annotated_maf_otx.tsv",
        "scratch/annotated_proteins/annotated_maf_otx.faa",

rule extract_annotated_maf_otx:
    input:
        "scratch/orthofinder_annotated_mammals/results/Results_Dec01/Orthogroups/Orthogroups.tsv",
        "scratch/orthofinder_annotated_birds/results/Results_Dec01/Orthogroups/Orthogroups.tsv",
        "scratch/orthofinder_annotated_nonbird_sauropsids/results/Results_Dec01/Orthogroups/Orthogroups.tsv",
        "scratch/orthofinder_annotated_fish_amphibians/results/Results_Dec01/Orthogroups/Orthogroups.tsv",
        "scratch/orthofinder_annotated_teleosts/results/Results_Dec03/Orthogroups/Orthogroups.tsv",
        "scratch/orthofinder_annotated_invertebrates/results/Results_Nov27/Orthogroups/Orthogroups.tsv",
    output:
        "exports/curated_genes/draft_annotated_maf_otx.tsv",
    script:
        "scripts/annotate_genes/extract_maf_otx_from_orthofinder_on_annotated_proteins.py"

rule make_fasta_annotated_maf_otx:
    input:
        "exports/curated_genes/draft_annotated_maf_otx.tsv",
    output:
        "scratch/annotated_proteins/annotated_maf_otx.faa",
        "scratch/annotated_proteins/annotated_vertebrate_maf_otx.faa",
        "scratch/annotated_proteins/annotated_invertebrate_maf_otx.faa",
    script:
        "scripts/annotate_genes/generate_fasta_tsv_annotated_maf_otx.py"

