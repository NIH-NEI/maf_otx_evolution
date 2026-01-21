import pandas as pd
import numpy as np

grand_tsv = snakemake.input["tsv"]
gene_family = snakemake.wildcards["genefamily"]

family_list = {
    "AllMafOtx": ["LargeMAF", "SmallMAF", "AllOTX"],
    "AllMaf": ["LargeMAF", "SmallMAF"],
    "AllOtx": ["AllOTX"],
    "SmallMaf": ["SmallMAF"],
    "LargeMaf": ["LargeMAF"],
}[gene_family]

print(gene_family)
print(family_list)


orgs = (
    pd.read_csv("configs/species173_table.tsv", sep = "\t")
    [["OrganismShortName", "OrganismID", "OrganismColor"]]
    .merge(
        (
            pd.read_csv("configs/species173_annotations.tsv", sep = "\t")
            [["OrganismShortName", "AnnotationSource"]]
        ),
        how  ="left",
    )
)

print(orgs)

df = (
    pd.read_csv(grand_tsv, sep = "\t")
    .query("GeneFamily in @family_list")
    .query("GeneGroup != 'Ignore'")
    .assign(GeneGroup = lambda tdf: np.where(tdf.GeneGroup == "MAFL", "CMAF", tdf.GeneGroup))
    .assign(GeneGroup = lambda tdf: np.where(tdf.GeneGroup == "MAFS", "MAFF", tdf.GeneGroup))
    .assign(GeneGroup = lambda tdf: np.where(tdf.GeneGroup == "OTX", "OTX1", tdf.GeneGroup))
    .pivot_table(
        values = "ProteinID",
        index = ["OrganismShortName", "ProteinSource"],
        columns = "GeneGroup",
        aggfunc= ','.join,
    )
    .fillna("")
)

df = df[[x
    for x in ["CMAF", "NRL", "MAFA", "MAFB", "MAFF", "MAFG", "MAFK", "OTX1", "OTX2", "CRX"]
    if x in df.columns
]]

#df["nMAF"] = (df != "").sum(axis=1)

print(df)
orgs = orgs.merge(df.reset_index(), how = "left")
print(orgs)
print(orgs.query("OrganismColor == 'Birds'"))
orgs.to_csv(snakemake.output[0], sep = "\t", index=False)
