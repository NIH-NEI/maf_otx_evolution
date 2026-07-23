import pandas as pd

ALL_SAMPLES_FILE = "configs/species107_rnaseq_samples.tsv"
ALL_SPECIES_FILE = "configs/species181_table.tsv"

rule all:
    input:
        "exports/bulk_heatmap/jawless_bulk_heatmap_samplewise.pdf",
        "exports/bulk_heatmap/jawed_bulk_heatmap_samplewise.pdf",
        "exports/bulk_heatmap/invertebrate_bulk_heatmap_samplewise.pdf",
        "exports/bulk_heatmap/jawless_bulk_heatmap_orgwise.pdf",
        "exports/bulk_heatmap/jawed_bulk_heatmap_orgwise.pdf",
        "exports/bulk_heatmap/invertebrate_bulk_heatmap_orgwise.pdf",


rule visualize_jawed_samplewise_expression:
    input:
        norm = "exports/kallisto_normalized/expr_matrix__jawed.tsv",
        samples = ALL_SAMPLES_FILE,
        orgs = ALL_SPECIES_FILE,
    output:
        pdf = "exports/bulk_heatmap/jawed_bulk_heatmap_samplewise.pdf",
    script:
        "scripts/quantify_rnaseq/show_mafs_otx_jawed_samplewise.R"

rule visualize_jawed_orgwise_expression:
    input:
        norm = "exports/kallisto_normalized/expr_matrix__jawed.tsv",
        samples = ALL_SAMPLES_FILE,
        orgs = ALL_SPECIES_FILE,
    output:
        pdf = "exports/bulk_heatmap/jawed_bulk_heatmap_orgwise.pdf",
    script:
        "scripts/quantify_rnaseq/show_mafs_otx_jawed_orgwise.R"


rule visualize_jawless_samplewise_expression:
    input:
        norm = "exports/kallisto_normalized/expr_matrix__jawless.tsv",
        samples = ALL_SAMPLES_FILE,
        orgs = ALL_SPECIES_FILE,
    output:
        pdf = "exports/bulk_heatmap/jawless_bulk_heatmap_samplewise.pdf",
    script:
        "scripts/quantify_rnaseq/show_mafs_otx_jawless_samplewise.R"

rule visualize_jawless_orgwise_expression:
    input:
        norm = "exports/kallisto_normalized/expr_matrix__jawless.tsv",
        samples = ALL_SAMPLES_FILE,
        orgs = ALL_SPECIES_FILE,
    output:
        pdf = "exports/bulk_heatmap/jawless_bulk_heatmap_orgwise.pdf",
    script:
        "scripts/quantify_rnaseq/show_mafs_otx_jawless_orgwise.R"


rule visualize_invertebrate_samplewise_expression:
    input:
        norm = "exports/kallisto_normalized/expr_matrix__invertebrate.tsv",
        samples = ALL_SAMPLES_FILE,
        orgs = ALL_SPECIES_FILE,
    output:
        pdf = "exports/bulk_heatmap/invertebrate_bulk_heatmap_samplewise.pdf",
    script:
        "scripts/quantify_rnaseq/show_mafs_otx_invertebrate_samplewise.R"

rule visualize_invertebrate_orgwise_expression:
    input:
        norm = "exports/kallisto_normalized/expr_matrix__invertebrate.tsv",
        samples = ALL_SAMPLES_FILE,
        orgs = ALL_SPECIES_FILE,
    output:
        pdf = "exports/bulk_heatmap/invertebrate_bulk_heatmap_orgwise.pdf",
    script:
        "scripts/quantify_rnaseq/show_mafs_otx_invertebrate_orgwise.R"


