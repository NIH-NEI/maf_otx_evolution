import pandas as pd

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
        },
    },
    "bird": {
        "path": "scratch/orthofinder_annotated_birds/results/Results_Dec01/Orthogroups/Orthogroups.tsv",
        "ogs": {
            "CMAF": ["OG0007266"],
            "MAFA": ["OG0007267"],
            "MAFB": ["OG0007268"],
        },
    },
    "reptile": {
        "path": "scratch/orthofinder_annotated_nonbird_sauropsids/results/Results_Dec01/Orthogroups/Orthogroups.tsv",
        "ogs": {
            "CMAF": ["OG0006433"],
            "MAFA": ["OG0013301"],
            "MAFB": ["OG0009609"],
            "NRL": ["OG0014298"],
        },
    },
    "fish": {
        "path": "scratch/orthofinder_annotated_fish_amphibians/results/Results_Dec01/Orthogroups/Orthogroups.tsv",
        "ogs": {
            "CMAF": ["OG0004727"],
            "MAFA": ["OG0004728"],
            "MAFB": ["OG0004726"],
            "NRL": ["OG0016169", "OG0024312"],
        },
    },
    "teleost": {
        "path": "scratch/orthofinder_annotated_teleosts/results/Results_Dec03/Orthogroups/Orthogroups.tsv",
        "ogs": {
            "CMAF": ["OG0000593"],
            "MAFA": ["OG0012617"],
            "MAFB": ["OG0000894"],
            "NRL": ["OG0013151"],
        },
    },
    "invertebrate": {
        "path": "scratch/orthofinder_annotated_invertebrates/results/Results_Nov27/Orthogroups/Orthogroups.tsv",
        "ogs": {
            "CMAF": ["OG0010116"],
        },
    },
}

dfs = []
for species_group, spgrp_dict in orthogroups.items():
    ogs_path = spgrp_dict["path"]
    ogs_dict = spgrp_dict["ogs"]
    for gene_group in ogs_dict.keys():
        dfs.append(get_genes(ogs_path, ogs_dict, gene_group, species_group))


all_genes = (
    pd.concat(dfs)
)
print(all_genes)

all_genes.to_csv("exports/curated_genes/draft_annotated_maf_otx_longformat.tsv", sep = "\t", index=False)
quit()

df = (
    fish_to_mammal_genes
    .pivot_table(
        values = "ProteinID",
        index = "OrganismShortName",
        columns = "GeneGroup",
        aggfunc= ','.join,
    )
    .fillna("")
)
df["nMAF"] = (df != "").sum(axis=1)

print(df)
orgs = orgs.merge(df.reset_index(), how = "left")
print(orgs)
print(orgs.query("OrganismColor == 'Birds'"))
orgs.to_csv("exports/curated_genes/draft_annotated_MAFL.tsv", sep = "\t")
