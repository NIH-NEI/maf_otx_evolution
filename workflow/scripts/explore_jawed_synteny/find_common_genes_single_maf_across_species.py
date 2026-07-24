import pandas as pd
from utils import *
import sys

# Constants
SEARCH_RANGE = 200
LEAST_NMATCH_FOR_OVERLAP = 3

# List of species
OTHER_SPECIES = [
    "opossum",
    "mouse",
    "chicken",
    "frog",
    "gar",
]
SPECIES = ["human"] + OTHER_SPECIES

# Function to read species data
def read_species_data(species):
    """
    Read gene annotation data for a given species.
    
    Args:
    - species: The name of the species.
    
    Returns:
    A pandas DataFrame containing gene annotation data.
    """
    annotation_file = f"scratch/gene_annot_tables/{species}_genes.tsv"
    return (
        pd.read_csv(annotation_file, sep="\t", header=None)
        .rename(columns={0: "chrom", 1: "start", 2: "end", 3: "gene"})
        # Do not uppercase gene names immediately now to
        # preserve species specific symbols
        #.assign(gene = lambda df: df["gene"].str.upper())
        .query("~gene.str.startswith('LOC')")  # Exclude certain genes
        .reset_index(drop=True)
    )

# Function to retrieve genes before and after a given query gene
def get_before_after_genes(species, query_gene, n=SEARCH_RANGE):
    """
    Get genes before and after a given query gene for a specific species.
    
    Args:
    - species: The name of the species.
    - query_gene: The gene to search for.
    - n: The number of genes to retrieve before and after the query gene.
    
    Returns:
    Two lists containing genes before and after the query gene, respectively.
    """
    actual_query = actual_gene[query_gene][species]
    print(species, actual_query)

    df = species_dataframes[species]
    found_genes = df[df.gene == actual_query]
    #print(found_genes)
    
    if found_genes.empty:
        # Query gene not found
        return [], [], pd.DataFrame()
    
    idx = found_genes.index[0]
    return (
        list(reversed(df.iloc[idx - n : idx].gene.to_list())),
        df.iloc[idx + 1 : idx + n + 1, :].gene.to_list(),
        (
            df.iloc[idx - n : idx + n + 1]
            .assign(pos = range(-n, n+1))
            .set_index("pos")
        ),
    )

# Function to assign pseudo order to genes
def assign_pseudo_gene_order(genes_dict, direction="after"):
    """
    Assign pseudo order to genes based on their positions.
    
    Args:
    - genes_dict: A dictionary containing genes for each species.
    - direction: The direction in which to assign order ("before" or "after").
    
    Returns:
    A list of genes in the assigned pseudo order.
    """
    def make_dataframe(species):
        print(species)
        return (
            pd.DataFrame({"gene": [x.upper() for x in genes_dict[species]]})
            .assign(species=species)
            .assign(pos=range(len(genes_dict[species])))
            .assign(pos=lambda df: df.pos + 1)
        )

    tmp = pd.concat([make_dataframe(species) for species in genes_dict.keys()])
    print(tmp)
    species_dataframes = (
        tmp.pivot_table(
            index="gene",
            columns="species",
            values="pos",
            fill_value=float("inf")
        )
        [tmp.species.unique()]
    )
    species_dataframes["spl"] = species_dataframes.apply(lambda x: tuple(x), axis=1)
    species_dataframes["spl2"] = species_dataframes.apply(lambda x: str(x.to_dict()), axis=1)
    species_dataframes = species_dataframes.sort_values("spl")
    return species_dataframes.index.to_list()

# Function to retrieve common genes before and after the query gene across all species
def get_common_genes_before_after(query_gene, directions):
    """
    Get common genes before and after the query gene across all species.
    
    Args:
    - query_gene: The gene to search for.
    
    Returns:
    Two lists containing common genes before and after the query gene, respectively.
    """
    before_genes, after_genes, df = get_before_after_genes("human", query_gene)
    before_gene_lists = {"human": before_genes}
    after_gene_lists = {"human": after_genes}
    annot_dfs = {"human": df}


    for species in OTHER_SPECIES:
        before_genes, after_genes, df = get_before_after_genes(species, query_gene)
        print("after_genes =", after_genes)
        print("before_genes =", before_genes)
        
        if not (before_genes or after_genes):
            # Gene not found
            continue

        if directions[species] is None:
            if (
                has_minimum_gene_overlap(
                    before_genes,
                    before_gene_lists["human"],
                    LEAST_NMATCH_FOR_OVERLAP
                )
                and has_minimum_gene_overlap(
                    after_genes,
                    after_gene_lists["human"],
                    LEAST_NMATCH_FOR_OVERLAP
                )
            ):
                directions[species] = +1
            elif (
                has_minimum_gene_overlap(
                    after_genes,
                    before_gene_lists["human"],
                    LEAST_NMATCH_FOR_OVERLAP
                )
                and has_minimum_gene_overlap(
                    before_genes,
                    after_gene_lists["human"],
                    LEAST_NMATCH_FOR_OVERLAP
                )
            ):
                directions[species] = -1
            else:
                print(before_genes)
                print(after_genes)
                sys.exit("Don't know what to do")

        if (directions[species] == -1):
            before_gene_lists[species] = after_genes
            after_gene_lists[species] = before_genes
            df.index = df.index * -1
            annot_dfs[species] = df
        else:
            before_gene_lists[species] = before_genes
            after_gene_lists[species] = after_genes
            annot_dfs[species] = df


    common_before_genes = get_common_genes(before_gene_lists)
    print_gene_lists(common_before_genes)
    b = assign_pseudo_gene_order(common_before_genes, "before")
    common_after_genes = get_common_genes(after_gene_lists)
    print_gene_lists(common_after_genes)
    a = assign_pseudo_gene_order(common_after_genes, "after")
    annotations = (
        pd.concat(annot_dfs, axis = 1, names = ["species", "gene_info"])
    )
    print(annotations)
    return b, a, annotations

# Main part of the script
if len(sys.argv) < 2:
    sys.exit("Error: specify maf-name")

maf_name = sys.argv[1]
all_maf_names = ["mafa", "mafb", "nrl", "cmaf", "maff", "mafg", "mafk"]
if maf_name not in all_maf_names:
    sys.exit("Error: unknown maf name")

maf_names_data = pd.read_csv("resources/maf_names.csv")
print(maf_names_data)

actual_gene = (
    maf_names_data
    .drop(columns="maftype")
    .melt(id_vars="mafname")
    .groupby("mafname")
    .apply(lambda x: dict(zip(x.variable, x.value)))
    .to_dict()
)
print(actual_gene)

chromosome_directions = (
    pd.read_csv("resources/single_maf_chromosome_direction.csv")
    .set_index("species")
)
print(chromosome_directions)


species_dataframes = {}
for species in SPECIES:
    species_dataframes[species] = read_species_data(species)

for species in SPECIES:
    print(species_dataframes[species])

b, a, df = get_common_genes_before_after(maf_name, chromosome_directions[maf_name].to_dict())
print(b)
print(a)

output_file = f"scratch/aligned_single_maf/{maf_name}_alignment.txt"
print(output_file)

with open(output_file, "w") as f:
    f.write(",".join(b) + "\n")
    f.write(",".join(a) + "\n")


df_file = f"scratch/aligned_single_maf/{maf_name}_alignment.h5"
print(df_file)
df.to_hdf(df_file, key = maf_name)
