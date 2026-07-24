import pandas as pd

SHEKHAR_TOP = "/data/pals2/celltype-evolution/data/ShekharData"

keep_species = [
#    "seaLamprey",
    "zebrafish",
    "brownAnole",
    "chicken",
    "opossum",
    "pig",
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
    .merge(pd.read_csv("configs/species181_table.tsv", sep = "\t"), how = "left")
    .query("OrganismShortName in @keep_species")
)

print(samples)

my2hahn = samples.set_index("OrganismShortName")["HahnName"].to_dict()
print(my2hahn)


rule all:
    input:
        "scratch/scrna_dotplots/allorgs_avg_expr_normalized.tsv",
        "scratch/scrna_dotplots/allorgs_avg_expr_normalized.pdf",


rule extract_average_expression:
    input:
        seu = lambda wc: SHEKHAR_TOP + f"/SeuratObjects/{my2hahn[wc.org]}_initial.rds",
    output:
        "scratch/scrna_dotplots/{org}_avg_expr_normalized.tsv"
    script:
        "scripts/visualize_scrna/get_normalized_average_expression.R"

rule combine_average_expression:
    input:
        expand("scratch/scrna_dotplots/{org}_avg_expr_normalized.tsv", org = samples.OrganismShortName),
    output:
        "scratch/scrna_dotplots/allorgs_avg_expr_normalized.tsv"
    run:
        import pandas as pd
        df = pd.concat([
            pd.read_csv(infile, sep = "\t")
            for infile in input
        ])
        df.to_csv(output[0], sep = "\t")

rule create_combine_avg_expr:
    input:
        "scratch/scrna_dotplots/allorgs_avg_expr_normalized.tsv",
        "scratch/scrna_cell_proportions/allorgs_cell_proportions.tsv",
    output:
        "scratch/scrna_dotplots/allorgs_avg_expr_normalized.pdf"
    script:
        "scripts/visualize_scrna/draw_hahn_combined_orgs_norm_expr_maf_otx.R"


