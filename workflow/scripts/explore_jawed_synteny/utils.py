import re

def capitalize_gene_lists(list_of_gene_lists):
    """
    Capitalizes all gene names in each species' gene list.

    Args:
        list_of_gene_lists (dict): A dictionary where keys are species and values are lists of genes.

    Returns:
        dict: A dictionary with species as keys and capitalized gene lists as values.
    """
    capitalized_gene_lists = {
        species: [gene.upper() for gene in gene_list]
        for species, gene_list in list_of_gene_lists.items()
    }
    return capitalized_gene_lists


def print_gene_lists(gene_lists):
    """
    Prints the species and their respective gene lists.

    Args:
        gene_lists (dict): A dictionary where keys are species and values are lists of genes.
    """
    for species, genes in gene_lists.items():
        print(f"{species}:\t{genes}")


def remove_trailing_numbers(input_string):
    """
    Removes trailing numbers from a string.

    Args:
        input_string (str): The string to process.

    Returns:
        str: The string without trailing numbers.
    """
    return re.sub(r'\d+$', '', input_string)


def get_common_genes_one_side(list_of_gene_lists):
    """
    Identifies genes common among multiple species.

    Args:
        list_of_gene_lists (dict): A dictionary where keys are species and values are lists of genes.

    Returns:
        dict: A dictionary with species as keys and lists of common genes as values.
    """
    capitalized_gene_lists = capitalize_gene_lists(list_of_gene_lists)

    # Get all unique genes across all species
    all_genes = list(set([
        gene
        for species in capitalized_gene_lists.keys()
        for gene in capitalized_gene_lists[species]
    ]))

    # Determine presence of each gene in different species
    presence = {gene: [] for gene in all_genes}
    for species in capitalized_gene_lists.keys():
        for gene in capitalized_gene_lists[species]:
            old_species = presence[gene]
            presence[gene] = list(set(old_species + [species]))

    # Filter genes present in at least two species
    common_genes = [
        gene
        for gene, species_list in presence.items()
        if len(species_list) > 1
    ]

    # Filter gene lists based on common genes
    common_gene_lists = {
        species: [gene for gene in gene_list if gene.upper() in common_genes]
        for species, gene_list in list_of_gene_lists.items()
    }

    return common_gene_lists

def get_common_genes(list_of_gene_lists):
    """
    Identifies genes common among multiple species.

    Args:
        list_of_gene_lists (dict): A dictionary where keys are species and values are lists of genes.

    Returns:
        dict: A dictionary with species as keys and lists of common genes as values.
    """
    capitalized_gene_lists = capitalize_gene_lists(list_of_gene_lists)

    # Get all unique genes across all species
    all_genes = list(set([
        gene
        for species in capitalized_gene_lists.keys()
        for gene in capitalized_gene_lists[species]
    ]))

    # Determine presence of each gene in different species
    presence = {gene: [] for gene in all_genes}
    for species in capitalized_gene_lists.keys():
        for gene in capitalized_gene_lists[species]:
            old_species = presence[gene]
            presence[gene] = list(set(old_species + [species]))

    # Filter genes present in at least two species
    common_genes = [
        gene
        for gene, species_list in presence.items()
        if len(species_list) > 1
    ]
    print(common_genes)
    print(len(common_genes))
    import pandas as pd
    df = (
        pd.DataFrame({
            "gene": presence.keys(),
            "presence": [len(species_list) for gene, species_list in
                presence.items()]
        })
        .query("presence > 1")
        .sort_values(by = "presence", ascending=False)
    )
    print(df)

    # Filter gene lists based on common genes
    common_gene_lists = {
        species: [gene for gene in gene_list if gene.upper() in common_genes]
        for species, gene_list in list_of_gene_lists.items()
    }

    return df.gene.to_list(), common_gene_lists


def has_minimum_gene_overlap(list_a, list_b, min_overlap):
    """
    Checks if two gene lists have a minimum overlap.

    Args:
        list_a (list): First list of genes.
        list_b (list): Second list of genes.
        min_overlap (int): Minimum number of common genes required for overlap.

    Returns:
        bool: True if overlap meets the minimum requirement, False otherwise.
    """
    n_found = 0
    for gene_a in sorted(list_a):
        for gene_b in sorted(list_b):
            if gene_a.upper().startswith(gene_b.upper()) or gene_b.upper().startswith(gene_a.upper()):
                n_found += 1
    return n_found >= min_overlap


def has_minimum_trimmed_gene_overlap(list_a, list_b):
    """
    Checks if two lists of genes, after removing trailing numbers, have at least one common gene.

    Args:
        list_a (list): First list of genes.
        list_b (list): Second list of genes.

    Returns:
        bool: True if at least one common gene is found after trimming, False otherwise.
    """
    # Trim trailing numbers from gene names in both lists
    trimmed_list_a = [remove_trailing_numbers(gene) for gene in list_a]
    trimmed_list_b = [remove_trailing_numbers(gene) for gene in list_b]

    # Initialize counter for common genes
    num_common_genes = 0

    # Loop through sorted and capitalized genes in trimmed lists
    for gene_a in sorted(trimmed_list_a):
        for gene_b in sorted(trimmed_list_b):
            # Capitalize gene names for comparison
            gene_a_upper = gene_a.upper()
            gene_b_upper = gene_b.upper()

            # Check if one gene starts with the other
            if gene_a_upper.startswith(gene_b_upper) or gene_b_upper.startswith(gene_a_upper):
                # Increment counter if a common gene is found
                num_common_genes += 1

    # Return True if at least one common gene is found, otherwise False
    return num_common_genes >= 1
