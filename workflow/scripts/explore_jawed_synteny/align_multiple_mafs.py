import pandas as pd
from utils import *
import sys

# Function to get a list of genes with a common prefix among all species
def get_genes_with_common_prefix_by_species(capitalized_gene_lists, pos_list):
    # Extract all unique genes from all species
    all_genes = list(set(
        [gene for species in capitalized_gene_lists.keys()
         for gene in capitalized_gene_lists[species]]
    ))

    # Dictionary to store prefix frequencies
    prefix_freq = {}
    for gene in all_genes:
        #prefix = gene[0:3]  # Extract first three characters as prefix
        prefix = remove_trailing_numbers(gene)
        prefix_freq[prefix] = []

    # Populate prefix frequencies for each species
    for species in capitalized_gene_lists.keys():
        for gene in capitalized_gene_lists[species]:
            #prefix = gene[0:3]
            prefix = remove_trailing_numbers(gene)
            # Add the species to the list of species having this prefix
            old = prefix_freq[prefix]
            prefix_freq[prefix] = list(set(old + [species]))

    # Filter out prefixes that occur in only one species or start with "LOC"
    good_prefix_list = [
        prefix
        for prefix in prefix_freq.keys()
        if (len(prefix_freq[prefix]) > 1) and (prefix != "LOC")
    ]

    # Filter genes using good prefixes and return
    return good_prefix_list, {
        species: list(sorted([
            f"{gene}:{pos}"
            for gene,pos in zip(capitalized_gene_lists[species], pos_list[species])
            if remove_trailing_numbers(gene.upper()) in good_prefix_list
        ]))
        for species in capitalized_gene_lists.keys()
    }

# Load maf types from CSV
maf_types = (
    pd.read_csv("configs/maf_names.csv")
    .set_index("maftype")
    .groupby("maftype")
    .agg({"mafname": list})
    ["mafname"]
    .to_dict()
)

# Print available maf types
print("Available MAF types:")
print(maf_types)

# Check if specified maf type is valid
if len(sys.argv) < 2:
    exit("Error: specify maf-type")

maftype = sys.argv[1]

if not(maftype in maf_types.keys()):
   exit("Error: unknown maf-type")

mafs = maf_types[maftype]

print("Selected MAF types:")
print(mafs)

# Load gene lists before and after
before_list = {}
after_list = {}
pos_list = {}

for maf in mafs:
    in_file = f"scratch/aligned_single_maf/{maf}_alignment.txt"
    with open(in_file) as f:
        text = f.readlines()
    before_list[maf] = [x.split(":")[0] for x in text[0].rstrip().split(",")]
    pos_list[maf] = [x.split(":")[1] for x in text[0].rstrip().split(",")]
    #after_list[maf] = text[1].rstrip().split(",")

# Print loaded gene lists
print("Before Alignment:")
print(before_list)
print("After Alignment:")
print(after_list)

# Select the first maf as reference and compare with other mafs
reference_maf = mafs[0]
other_mafs = mafs[1:]

#for maf in other_mafs:
#    b = before_list[maf]
#    a = after_list[maf]
#    if (has_minimum_trimmed_gene_overlap(b, before_list[reference_maf]) and
#            has_minimum_trimmed_gene_overlap(a, after_list[reference_maf])):
#        before_list[maf] = b
#        after_list[maf] = a
#    elif (has_minimum_trimmed_gene_overlap(a, before_list[reference_maf]) and
#            has_minimum_trimmed_gene_overlap(b, after_list[reference_maf])):
#        before_list[maf] = a
#        after_list[maf] =  b
#    else:
#        print("Gene lists before and after:")
#        print(b)
#        print(a)
#        exit("Error: Cannot find overlap")

# Get gene lists with common prefix
prefix_before, res_before = get_genes_with_common_prefix_by_species(before_list, pos_list)
#prefix_after, res_after = get_genes_with_common_prefix_by_species(after_list)

# Print final gene lists
print("=======================  Final =================")
print("Leftside of Alignment:")
print_gene_lists(res_before)
#print("Rightside of Alignment:")
#print_gene_lists(res_after)

output_file = f"scratch/aligned_maf_groups/{maftype}_prefixes.txt"
print(output_file)

with open(output_file, "w") as f:
    f.write(",".join(prefix_before) + "\n")
#    f.write(",".join(prefix_after) + "\n")

output_file = f"scratch/aligned_maf_groups/{maftype}_alignment.txt"
print(output_file)

with open(output_file, "w") as f:
    for maf, genes in res_before.items():
        f.write(f"before,{maf}|" + ",".join(genes) + "\n")
#    for maf, genes in res_after.items():
#        f.write(f"after,{maf}:" + ",".join(genes) + "\n")

