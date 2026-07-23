import pandas as pd
import numpy as np

atmost1gene = snakemake.input["atmost1gene"]
grand_list = snakemake.input["grand_list"]
maf_otx_focused = snakemake.output["maf_otx_focused"]

ogs = (
    pd.read_csv(atmost1gene, sep = "\t", dtype = str, low_memory = False)
)
print(ogs)

organisms = ogs.columns
print(organisms)

maf_otx_info = (
    pd.read_csv(grand_list, sep = "\t")
    .query("GeneFamily in ['AllOTX', 'LargeMAF', 'SmallMAF']")
    .query("GeneGroup not in ['Ignore', 'Outgroup']")
    .query("OrganismShortName in @organisms")
    .query("Status == 'Final'")
    .query("ProteinSource == 'Annotated'")
)

df = (
    ogs
    .melt(id_vars = "Orthogroup", var_name = "OrganismShortName", value_name = "ProteinID")
    .query("ProteinID.notna()")
    .assign(ProteinID = lambda tdf: tdf.ProteinID.apply(lambda x: x.split(", ")))
    .explode("ProteinID")
)
print(df)

maf_otx_ogs = df.query("ProteinID in @maf_otx_info.ProteinID").Orthogroup.unique().tolist()
print(maf_otx_ogs)
others = (
    ogs
    .query("Orthogroup not in @maf_otx_ogs")
    .melt(id_vars = "Orthogroup", var_name = "OrganismShortName", value_name = "ProteinID")
    .query("ProteinID.notna()")
    .assign(ProteinID = lambda tdf: tdf.ProteinID.apply(lambda x: x.split(", ")))
    .explode("ProteinID")
    .assign(ProteinID = lambda tdf: tdf.ProteinID.apply(lambda x: x.replace("__gene_ENS", "__ENS")))
    .assign(GeneSymbol = lambda tdf: tdf.ProteinID.str.split("__").str[1])
    [["Orthogroup", "OrganismShortName", "GeneSymbol"]]
    .drop_duplicates()
    .pivot_table(
        index="Orthogroup",
        columns="OrganismShortName",
        values="GeneSymbol",
        aggfunc=",".join,
        fill_value="",
    )
    .reset_index()
)
if "human" in others.columns:
    others["Orthogroup"] = np.where(
        others.human != "", others.human, others.Orthogroup
    )

others = (
    others
    .rename(columns = {
        "Orthogroup" : "OrthoSymbol"
    })
    .set_index("OrthoSymbol")
)
print(others)

maf_otx = (
    maf_otx_info
    .assign(OrganismShortName = lambda tdf: tdf.ProteinID.str.split("__").str[0])
    .assign(ProteinID = lambda tdf: tdf.ProteinID.apply(lambda x: x.replace("__gene_ENS", "__ENS")))
    .assign(GeneSymbol = lambda tdf: tdf.ProteinID.str.split("__").str[1])
    [["OrganismShortName", "GeneGroup", "GeneSymbol"]]
    .rename(columns = {"GeneGroup": "OrthoSymbol"})
    .drop_duplicates()
    .pivot_table(
        index="OrthoSymbol",
        columns="OrganismShortName",
        values="GeneSymbol",
        aggfunc=",".join,
        fill_value="",
    )
)
print(maf_otx)

orthotab = (
    pd.concat([maf_otx, others])
)
print(orthotab)
orthotab.to_csv(maf_otx_focused, sep = "\t")
