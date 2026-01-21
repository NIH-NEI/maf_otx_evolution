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

def get_genes(ogs_path, ogs_dict, gene_group):
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
        [["OrganismShortName", "GeneGroup", "ProteinID"]]
    )

    return df


mammal_ogs = {
    "MAFF": ["OG0007929"],
    "MAFG": ["OG0014001"],
    "MAFK": ["OG0007930"],
}
mammal_ogs_path = "scratch/orthofinder_annotated_mammals/results/Results_Dec01/Orthogroups/Orthogroups.tsv"

bird_ogs = {
    "MAFF": ["OG0011722"],
    "MAFG": ["OG0011723"],
    "MAFK": ["OG0007269"],
}
bird_ogs_path = "scratch/orthofinder_annotated_birds/results/Results_Dec01/Orthogroups/Orthogroups.tsv"

reptile_ogs = {
    "MAFF": ["OG0009610"],
    "MAFG": ["OG0005530"],
    "MAFK": ["OG0002641"],
}
reptile_ogs_path = "scratch/orthofinder_annotated_nonbird_sauropsids/results/Results_Dec01/Orthogroups/Orthogroups.tsv"

fish_ogs = {
    "MAFF": ["OG0010284", "OG0028082"],
    "MAFG": ["OG0005245"],
    "MAFK": ["OG0014575", "OG0024119"],
}
fish_ogs_path = "scratch/orthofinder_annotated_fish_amphibians/results/Results_Dec01/Orthogroups/Orthogroups.tsv"

teleost_ogs = {
    "MAFF": ["OG0005675"],
    "MAFG": ["OG0000560"],
    "MAFK": ["OG0007340"],
}
teleost_ogs_path = "scratch/orthofinder_annotated_teleosts/results/Results_Dec03/Orthogroups/Orthogroups.tsv"

invertebrate_ogs = {
    "MAFF": ["OG0011228"],
}
invertebrate_ogs_path = "scratch/orthofinder_annotated_invertebrates/results/Results_Nov27/Orthogroups/Orthogroups.tsv"


mammal_genes = (
    pd.concat([
        get_genes(mammal_ogs_path, mammal_ogs, "MAFF"),
        get_genes(mammal_ogs_path, mammal_ogs, "MAFG"),
        get_genes(mammal_ogs_path, mammal_ogs, "MAFK"),
    ])
)
print(mammal_genes)

bird_genes = (
    pd.concat([
        get_genes(bird_ogs_path, bird_ogs, "MAFF"),
        get_genes(bird_ogs_path, bird_ogs, "MAFG"),
        get_genes(bird_ogs_path, bird_ogs, "MAFK"),
    ])
)
print(bird_genes)

reptile_genes = (
    pd.concat([
        get_genes(reptile_ogs_path, reptile_ogs, "MAFF"),
        get_genes(reptile_ogs_path, reptile_ogs, "MAFG"),
        get_genes(reptile_ogs_path, reptile_ogs, "MAFK"),
    ])
)
print(reptile_genes)

fish_genes = (
    pd.concat([
        get_genes(fish_ogs_path, fish_ogs, "MAFF"),
        get_genes(fish_ogs_path, fish_ogs, "MAFG"),
        get_genes(fish_ogs_path, fish_ogs, "MAFK"),
    ])
)
print(fish_genes)

teleost_genes = (
    pd.concat([
        get_genes(teleost_ogs_path, teleost_ogs, "MAFF"),
        get_genes(teleost_ogs_path, teleost_ogs, "MAFG"),
        get_genes(teleost_ogs_path, teleost_ogs, "MAFK"),
    ])
)
print(teleost_genes)

invertebrate_genes = (
    pd.concat([
        get_genes(invertebrate_ogs_path, invertebrate_ogs, "MAFF"),
    ])
)
print(fish_genes)

fish_to_mammal_genes = (
    pd.concat([
        mammal_genes,
        bird_genes,
        reptile_genes,
        fish_genes,
        teleost_genes,
        invertebrate_genes,
    ])
)

print(fish_to_mammal_genes)
print(fish_to_mammal_genes.query("OrganismShortName == 'barnOwl'"))
print(fish_to_mammal_genes.query("OrganismShortName == 'spottedGar'"))

fish_to_mammal_genes.to_csv("exports/curated_genes/draft_annotated_MAFS_longformat.tsv", sep = "\t", index=False)

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
orgs.to_csv("exports/curated_genes/draft_annotated_MAFS.tsv", sep = "\t", index=False)
