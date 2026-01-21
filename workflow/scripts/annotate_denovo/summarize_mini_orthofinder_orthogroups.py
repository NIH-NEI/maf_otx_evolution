import pandas as pd
import numpy as np

input = snakemake.input
output = snakemake.output

prot2group = (
    pd.read_csv(input.annot, sep = "\t")
    [["ProteinID", "GeneGroup"]]
)
print(prot2group)

prot2len = (
    pd.read_csv(input.info, sep = "\t")
    [["ProteinID", "ProteinLen"]]
)
print(prot2len)

df = (
    pd.read_csv(input.ogs, sep = "\t")
    .melt(id_vars = "Orthogroup", var_name = "OrganismShortName", value_name = "ProteinID", )
    .dropna()
    .assign(ProteinID = lambda tdf: tdf.ProteinID.apply(lambda x: [y.strip() for y in x.split(",")]))
    .explode("ProteinID")
    .assign(GeneSymbol = lambda tdf: tdf.ProteinID.str.split("__").str[1])
    .dropna(subset = ["ProteinID"])
)
print(df)

tally = (
    df
    .merge(prot2group, how = "left")
    .assign(ngrp = lambda tdf: np.where(pd.isna(tdf.GeneGroup), 0, 1))
)

tot = (
    tally 
    .groupby(["Orthogroup"])
    .agg({"ngrp": sum})
    .rename(columns = {"ngrp": "total"})
    .reset_index()
)

ngrp = (
    tally
    .groupby(["Orthogroup", "GeneGroup"])
    .agg({"ngrp": sum})
    .reset_index()
    .merge(tot, how = "left")
    .assign(percent = lambda tdf: tdf.ngrp * 100 / tdf.total)
)
print(ngrp)

quants = (
    pd.read_csv(input.quant, sep = "\t")
)
print(quants)

THRESHOLD = 1.0

df = (
    df
    .merge(ngrp, how = "left")
    .merge(prot2len, how = "left")
    .query("GeneSymbol == 'Unknown'")
    .merge(quants, how = "left")
    .fillna({"TPM_max": THRESHOLD})
    .query("total > 0")
)
print(df)

df.to_csv(output.long, sep = "\t", index = False)

df = df.query("TPM_max >= @THRESHOLD")
print(df)

df = (
    df
    .query("percent > 25.0")
    [["ProteinID", "GeneGroup", "ProteinLen"]]
    .sort_values(["GeneGroup", "ProteinLen"], ascending=[True, False])
    .drop_duplicates(subset="GeneGroup", keep="first")
)
print(df)
df.to_csv(output.short, sep = "\t", index = False)
