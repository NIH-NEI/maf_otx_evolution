import pandas as pd

ALL_SAMPLES_FILE = "configs/species107_rnaseq_samples.tsv"
ALL_SPECIES_FILE = "configs/species181_table.tsv"

rule all:
    input:
        "exports/bulk_heatmap/jawless_bulk_heatmap.pdf",
        "exports/bulk_heatmap/jawed_bulk_heatmap.pdf",
        "exports/bulk_heatmap/invertebrate_bulk_heatmap.pdf",


rule visualize_jawed_expression:
    input:
        norm = "exports/kallisto_normalized/expr_matrix__jawed.tsv",
        samples = ALL_SAMPLES_FILE,
        orgs = ALL_SPECIES_FILE,
    output:
        pdf = "exports/bulk_heatmap/jawed_bulk_heatmap.pdf",
    script:
        "scripts/quantify_rnaseq/show_mafs_otx_jawed.R"


rule visualize_jawless_expression:
    input:
        norm = "exports/kallisto_normalized/expr_matrix__jawless.tsv",
        samples = ALL_SAMPLES_FILE,
        orgs = ALL_SPECIES_FILE,
    output:
        pdf = "exports/bulk_heatmap/jawless_bulk_heatmap.pdf",
    script:
        "scripts/quantify_rnaseq/show_mafs_otx_jawless.R"


rule visualize_invertebrate_expression:
    input:
        norm = "exports/kallisto_normalized/expr_matrix__invertebrate.tsv",
        samples = ALL_SAMPLES_FILE,
        orgs = ALL_SPECIES_FILE,
    output:
        pdf = "exports/bulk_heatmap/invertebrate_bulk_heatmap.pdf",
    script:
        "scripts/quantify_rnaseq/show_mafs_otx_invertebrate.R"


