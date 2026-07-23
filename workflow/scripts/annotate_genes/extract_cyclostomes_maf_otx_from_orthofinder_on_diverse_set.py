import pandas as pd
import numpy as np

def get_genes(ogs_path, ogs_dict, gene_group, species_group):
    oglist = ogs_dict[gene_group]
    df = (
        pd.read_csv(ogs_path, sep = "\t")
        .query("Orthogroup in @oglist")
    )
    singletons = [og for og in oglist if og not in df.Orthogroup.tolist()]
    if singletons:
        ogs_txt = ogs_path.replace(".tsv", ".txt")
        with open(ogs_txt, "r") as f:
            for l in f:
                og, rest = l.strip().split(": ", maxsplit=1)
                if og in singletons:
                    row = {"Orthogroup":[og]}
                    proteins = rest.split(" ")
                    for prot in proteins:
                        row[prot.split("__")[0]] = [prot]
                    df = pd.concat([df, pd.DataFrame(row)])
    df = (
        df
        .melt(id_vars = "Orthogroup", var_name = "OrganismShortName", value_name = "ProteinID", )
        .dropna()
        .assign(ProteinID = lambda tdf: tdf.ProteinID.apply(lambda x: [y.strip() for y in x.split(",")]))
        .explode("ProteinID")
        #.assign(ProteinID = lambda tdf: tdf.ProteinID.str.split("__").str[1:3].str.join("__"))
        .dropna(subset = ["ProteinID"])
        .assign(GeneGroup = gene_group)
        .assign(SpeciesGroup = species_group)
        [["OrganismShortName", "SpeciesGroup", "GeneGroup", "ProteinID"]]
    )
    return df

orthogroups = {
    "diverse": {
        "path": "scratch/orthofinder_diverse/results/Results_Feb11/Orthogroups/Orthogroups.tsv",
        "ogs": {
            "MAFL": ["OG0000749"],
            "MAFS": ["OG0001877"],
            "OTX": ["OG0000728"],
        },
    },
}

dfs = []
for species_group, spgrp_dict in orthogroups.items():
    ogs_path = spgrp_dict["path"]
    ogs_dict = spgrp_dict["ogs"]
    for gene_group in ogs_dict.keys():
        dfs.append(get_genes(ogs_path, ogs_dict, gene_group, species_group))

maf_otx = (
    pd.concat(dfs)
)
print(maf_otx)

maf_otx.to_csv("exports/curated_genes/draft_annotated_maf_otx_longformat_cyclostomes.tsv", sep = "\t", index=False)
