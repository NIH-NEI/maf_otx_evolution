import pandas as pd
import numpy as np

orgs = (
    pd.read_csv("configs/species173_table.tsv", sep = "\t")
    [["OrganismShortName", "OrganismID", "OrganismColor"]]
)

print(orgs)

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
    "mammal": {
        "path": "scratch/orthofinder_annotated_mammals/results/Results_Dec01/Orthogroups/Orthogroups.tsv",
        "ogs": {
            "CMAF": ["OG0014168"],
            "MAFA": ["OG0011291"],
            "MAFB": ["OG0008579"],
            "NRL": ["OG0011758", "OG0022597"],
            "MAFF": ["OG0007929"],
            "MAFG": ["OG0014001"],
            "MAFK": ["OG0007930"],
            "OTX1": ["OG0011851"],
            "OTX2": ["OG0002540"],
            "CRX": ["OG0016247"],
        },
    },
    "bird": {
        "path": "scratch/orthofinder_annotated_birds/results/Results_Dec01/Orthogroups/Orthogroups.tsv",
        "ogs": {
            "CMAF": ["OG0007266"],
            "MAFA": ["OG0007267"],
            "MAFB": ["OG0007268"],
            "MAFF": ["OG0011722"],
            "MAFG": ["OG0011723"],
            "MAFK": ["OG0007269"],
            "OTX1": ["OG0011861"],
            "OTX2": ["OG0014511", "OG0016712"],
            "CRX": ["OG0014668"],
        },
    },
    "reptile": {
        "path": "scratch/orthofinder_annotated_nonbird_sauropsids/results/Results_Dec01/Orthogroups/Orthogroups.tsv",
        "ogs": {
            "CMAF": ["OG0006433"],
            "MAFA": ["OG0013301"],
            "MAFB": ["OG0009609"],
            "NRL": ["OG0014298"],
            "MAFF": ["OG0009610"],
            "MAFG": ["OG0005530"],
            "MAFK": ["OG0002641"],
            "OTX1": ["OG0012854"],
            "OTX2": ["OG0002314"],
            "CRX": ["OG0002313"],
        },
    },
    "fish": {
        "path": "scratch/orthofinder_annotated_fish_amphibians/results/Results_Dec01/Orthogroups/Orthogroups.tsv",
        "ogs": {
            "CMAF": ["OG0004727"],
            "MAFA": ["OG0004728"],
            "MAFB": ["OG0004726"],
            "NRL": ["OG0016169", "OG0024312"],
            "MAFF": ["OG0010284", "OG0028082"],
            "MAFG": ["OG0005245"],
            "MAFK": ["OG0014575", "OG0024119"],
            "OTX1": ["OG0010187"],
            "OTX2": ["OG0005123"],
            "CRX": ["OG0005122"],
        },
    },
    "teleost": {
        "path": "scratch/orthofinder_annotated_teleosts/results/Results_Dec03/Orthogroups/Orthogroups.tsv",
        "ogs": {
            "CMAF": ["OG0000593"],
            "MAFA": ["OG0012617"],
            "MAFB": ["OG0000894"],
            "NRL": ["OG0013151"],
            "MAFF": ["OG0005675"],
            "MAFG": ["OG0000560"],
            "MAFK": ["OG0007340"],
            "OTX1": ["OG0007909"],
            "OTX2": ["OG0004178", "OG0007075"],
            "CRX": ["OG0003566"],
        },
    },
    "invertebrate": {
        "path": "scratch/orthofinder_annotated_invertebrates/results/Results_Nov27/Orthogroups/Orthogroups.tsv",
        "ogs": {
            "MAFL": ["OG0010116"],
            "MAFS": ["OG0011228"],
            "OTX": ["OG0000445"],
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

maf_otx.to_csv("exports/curated_genes/draft_annotated_maf_otx_longformat.tsv", sep = "\t", index=False)

df = (
    maf_otx
    .assign(GeneGroup = lambda tdf: np.where(tdf.GeneGroup == "MAFL", "CMAF",
                                    np.where(tdf.GeneGroup == "MAFS", "MAFF", 
                                    np.where(tdf.GeneGroup == "OTX", "OTX1", tdf.GeneGroup))))
    .pivot_table(
        values = "ProteinID",
        index = ["OrganismShortName", "SpeciesGroup"],
        columns = "GeneGroup",
        aggfunc= ','.join,
    )
    .fillna("")
)
df["nMAFL"] = (df[["CMAF", "MAFA", "MAFB", "NRL"]] != "").sum(axis=1)
df["nMAFS"] = (df[["MAFF", "MAFG", "MAFK"]] != "").sum(axis=1)
df["nOTX"] = (df[["OTX1", "OTX2", "CRX"]] != "").sum(axis=1)

print(df)
orgs = orgs.merge(df.reset_index(), how = "left")
print(orgs)
print(orgs.query("OrganismColor == 'Birds'"))
orgs.to_csv("exports/curated_genes/draft_annotated_maf_otx.tsv", sep = "\t")
