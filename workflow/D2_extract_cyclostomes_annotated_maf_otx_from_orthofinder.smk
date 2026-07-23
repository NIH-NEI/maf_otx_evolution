# WARNING: This Snakefile is not intended to be run as-is.
# It depends on OrthoFinder results with hard-coded output paths and/or filenames,
# which are specific to a particular run and may not be valid for your dataset.

rule all:
    input:
        "scratch/grand_list/extra_cyclostomes.tsv"

rule extract_annotated_maf_otx:
    input:
        "scratch/orthofinder_diverse/results/Results_Feb11/Orthogroups/Orthogroups.tsv",
    output:
        "exports/curated_genes/draft_annotated_maf_otx_longformat_cyclostomes.tsv",
    script:
        "scripts/annotate_genes/extract_cyclostomes_maf_otx_from_orthofinder_on_diverse_set.py"

rule make_to_add_to_grand_list:
    input:
        "exports/curated_genes/draft_annotated_maf_otx_longformat_cyclostomes.tsv",
    output:
        "scratch/grand_list/extra_cyclostomes.tsv"
    script:
        "scripts/annotate_genes/generate_fasta_tsv_annotated_maf_otx_cyclostomes.py"
