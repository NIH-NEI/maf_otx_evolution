import pandas as pd



rule all:
    input:
        "exports/gene_presence/gene_presence_heatmap_AllMaf.pdf",
        "exports/gene_presence/gene_presence_table_AllMafOtx.tsv",

rule draw_gene_presence_heatmap:
    input:
        "scratch/grand_list/grand_list.tsv"
    output:
        "exports/gene_presence/gene_presence_heatmap_{genefamily}.pdf"
    script:
        "scripts/present_gene_finds/draw_presence_heatmap.R"

rule create_gene_presence_table:
    input:
        tsv = "scratch/grand_list/grand_list.tsv"
    output:
        "exports/gene_presence/gene_presence_table_{genefamily}.tsv"
    script:
        "scripts/present_gene_finds/create_presence_table.py"
