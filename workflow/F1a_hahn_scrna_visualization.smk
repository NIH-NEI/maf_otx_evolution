import pandas as pd

SHEKHAR_TOP = "/data/pals2/celltype-evolution/data/ShekharData"

keep_species = [
#    "seaLamprey",
    "zebrafish",
    "brownAnole",
    "chicken",
    "opossum",
    "cow",
    "sheep",
    "ferret",
    "treeShrew",
    "groundSquirrel",
    "deerMouse",
    "grassMouse",
    "labMouse",
    "marmoset",
    "crabMacaque",
    "human",
]

samples = (
    pd.read_csv("configs/hahnscrna17.tsv", sep = "\t")
    .merge(pd.read_csv("configs/species_table.tsv", sep = "\t"), how = "left")
    .query("OrganismShortName in @keep_species")
)

print(samples)

my2hahn = samples.set_index("OrganismShortName")["HahnName"].to_dict()
print(my2hahn)


rule all:
    input:
        "scratch/scrna_dotplots/allorgs_dotplot.pdf",
        "scratch/scrna_cell_proportions/allorgs_cell_proportions.pdf",


rule draw_dot_plots:
    input:
        seu = lambda wc: SHEKHAR_TOP + f"/SeuratObjects/{my2hahn[wc.org]}_initial.rds",
    output:
        "scratch/scrna_dotplots/{org}_dotplot.pdf"
    script:
        "scripts/visualize_scrna/draw_hahn_dotplots_maf_otx.R"

rule extract_dot_plot_data:
    input:
        seu = lambda wc: SHEKHAR_TOP + f"/SeuratObjects/{my2hahn[wc.org]}_initial.rds",
    output:
        "scratch/scrna_dotplots/{org}_dotplot.tsv"
    script:
        "scripts/visualize_scrna/extract_dotplot_data_maf_otx.R"

rule combine_dotplots_data:
    input:
        expand("scratch/scrna_dotplots/{org}_dotplot.tsv", org = samples.OrganismShortName),
    output:
        "scratch/scrna_dotplots/allorgs_dotplot.tsv"
    run:
        import pandas as pd
        df = pd.concat([
            pd.read_csv(infile, sep = "\t")
            for infile in input
        ])
        df.to_csv(output[0], sep = "\t")

rule create_combine_dotplot:
    input:
        "scratch/scrna_dotplots/allorgs_dotplot.tsv"
    output:
        "scratch/scrna_dotplots/allorgs_dotplot.pdf"
    script:
        "scripts/visualize_scrna/draw_hahn_combined_orgs_dotplots_maf_otx.R"



rule extract_cell_proportions:
    input:
        seu = lambda wc: SHEKHAR_TOP + f"/SeuratObjects/{my2hahn[wc.org]}_initial.rds",
    output:
        "scratch/scrna_cell_proportions/{org}_cell_proportions.tsv"
    script:
        "scripts/visualize_scrna/extract_cell_proportions.R"

rule combine_cell_proportions_data:
    input:
        expand("scratch/scrna_cell_proportions/{org}_cell_proportions.tsv", org = samples.OrganismShortName),
    output:
        "scratch/scrna_cell_proportions/allorgs_cell_proportions.tsv"
    run:
        import pandas as pd
        df = pd.concat([
            pd.read_csv(infile, sep = "\t")
            for infile in input
        ])
        df.to_csv(output[0], sep = "\t")

rule draw_cell_proportions:
    input:
        "scratch/scrna_cell_proportions/allorgs_cell_proportions.tsv"
    output:
        "scratch/scrna_cell_proportions/allorgs_cell_proportions.pdf"
    script:
        "scripts/visualize_scrna/draw_hahn_combined_orgs_cell_proportions.R"
