import pandas as pd
from tqdm import tqdm
tqdm.pandas()

input = snakemake.input
output = snakemake.output

orthogroups = pd.read_csv(input.ogs, sep = "\t", dtype = str, low_memory = False)
df = (
    orthogroups
    .melt(id_vars = "Orthogroup", var_name = "OrganismShortName", value_name = "ProteinID")
    .query("ProteinID.notna()")
    .assign(ProteinID = lambda tdf: tdf.ProteinID.apply(lambda x: x.split(", ")))
    .explode("ProteinID")
    .assign(GeneSymbol = lambda tdf: tdf.ProteinID.str.split("__").str[1])
    [["Orthogroup", "OrganismShortName", "GeneSymbol", "ProteinID"]]
    .drop_duplicates()
)
assert(df.loc[df.GeneSymbol.isna(),:].empty)
print(df)

# There could be mutiple proteins for the same gene but assigned to different orthogroups
# Simply ignore them for ortholog table
weird = (
    df
    .groupby(["OrganismShortName", "GeneSymbol"])
    .agg(
        ProteinIDs=("ProteinID", ",".join),
        OrthogroupCount=("Orthogroup", "nunique"),
        Orthogroups=("Orthogroup", ",".join),
    )
    .query("OrthogroupCount > 1")
)
print(weird)
#weird.to_csv("weird.tsv", sep = "\t")
weird_proteins = ",".join(weird.ProteinIDs).split(",")
print(weird_proteins)

df = df.query("ProteinID not in @weird_proteins")

# Now treat multiple proteins with same gene symbol that were assigned the
# same orthogroup as different isoforms
df = (
    df
    [["Orthogroup", "OrganismShortName", "GeneSymbol"]]
    .drop_duplicates()
)
print(df)

organisms = (
    pd.read_csv("configs/species_table.tsv", sep = "\t")
    [["OrganismShortName", "OrganismColor"]]
)
print(organisms)

df = (
    df.merge(organisms, how = "left")
    .assign(OrthoOrg = lambda tdf: tdf.Orthogroup + "__" + tdf.OrganismShortName)
)
print(df)

counts = (
    df
    .groupby(["OrthoOrg", "Orthogroup", "OrganismShortName", "OrganismColor"])
    .agg({"GeneSymbol": "count"})
    .reset_index()
)
print(counts)

unique_mask = counts["GeneSymbol"] <= 1
unique_counts = counts[unique_mask]
print(unique_counts)

multi_orthogroups = df.query("OrthoOrg not in @unique_counts.OrthoOrg")
print(multi_orthogroups)

duplicates_allowed = {
    "Nonchordates" : 1,
    "Protochordates" : 1,
    "Agnathans (hagfish and lampreys)" : 1,
    "Cartilaginous fishes" : 1,
    "Non-teleost ray-finned fishes" : 1,
    "Teleost ray-finned fishes" : 2,
    "Lobe-finned fishes" : 1,
    "Amphibians" : 2,
    "Turtles" : 1,
    "Crocodiles" : 1,
    "Squamata (lizards and snakes)" : 1,
    "Birds" : 1,
    "Mammals" : 1,
}

# get proper count of genes per orthogroup and organism
def count_distinct_genes(grp):
    symbols = grp.GeneSymbol.tolist()
    color = grp.OrganismColor.tolist()[0]
    down_factor = duplicates_allowed[color]
    count = len(symbols)/down_factor
    return count


multi_counts = (
    multi_orthogroups
    .groupby(["OrthoOrg", "Orthogroup", "OrganismShortName", "OrganismColor"])
    .progress_apply(count_distinct_genes)
)
multi_counts.name = "GeneSymbol"
multi_counts = multi_counts.reset_index()
print(multi_counts)

counts = (
    pd.concat([unique_counts, multi_counts])
)
print(counts)

counts = counts.pivot_table(
    index="Orthogroup",
    columns="OrganismShortName",
    values="GeneSymbol",
    aggfunc="max",   # or sum, max, min, etc.
    fill_value=0,
)
print(counts)

# must be present in at least 2/3 rd of the species
must_present_in = counts.shape[1] * 2 / 3
print(counts.shape, must_present_in)

atmost1gene = counts[counts.le(1.0).all(axis=1)]
print(atmost1gene)

atmost1gene_abs_major = atmost1gene[(atmost1gene >= 1).sum(axis = 1) >= must_present_in]
print(atmost1gene_abs_major)

df = orthogroups.query("Orthogroup in @atmost1gene_abs_major.index")
print(df)

df.to_csv(output.atmost1gene, index = None, sep = "\t")
