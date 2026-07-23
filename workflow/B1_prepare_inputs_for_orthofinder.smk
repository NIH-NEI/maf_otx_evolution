import pandas as pd
from os.path import basename

use_annotations = [
    "GCF",
    "Ensembl",
    "Custom"
]

annotated = (
    pd.read_csv("configs/species181_table.tsv", sep = "\t")
    .query("AnnotationSource in @use_annotations")
)
print(annotated)
print(annotated.head(20))

invertebrates = (
    annotated
    .query("OrganismColor in ['Nonchordates', 'Protochordates']")
)
print(invertebrates)

mammals = (
    annotated
    .query("OrganismColor in ['Mammals']")
)
print(mammals)

birds = (
    annotated
    .query("OrganismColor in ['Birds']")
)
print(birds)

nonbird_sauropsids = (
    annotated
    .query("OrganismColor in ['Squamata (lizards and snakes)', 'Crocodiles', 'Turtles']")
)
print(nonbird_sauropsids)

fish_amphibians = (
    annotated
    .query("OrganismColor in ['Amphibians', 'Cartilaginous fishes', 'Lobe-finned fishes', 'Non-teleost ray-finned fishes']")
)
print(fish_amphibians)

teleosts = (
    annotated
    .query("OrganismColor in ['Teleost ray-finned fishes']")
)
print(teleosts)

OF_groups = {
    "annotated_mammal": mammals,
    "annotated_bird": birds,
    "annotated_reptile": nonbird_sauropsids,
    "annotated_fish": fish_amphibians,
    "annotated_teleost": teleosts,
    "annotated_invertebrates": invertebrates,
}

group_org_pairs = []
for grp, df in OF_groups.items():
    for org in df.OrganismShortName.unique():
        group_org_pairs.append((grp, org))
group_org_pairs = pd.DataFrame(group_org_pairs, columns = ["grp", "org"])
print(group_org_pairs)


rule all:
    input:
        expand("scratch/orthofinder_{group}/peptides/{org}.faa", zip, group=group_org_pairs.grp, org=group_org_pairs.org),
        expand("scratch/orthofinder_{group}/topology_{group}.nwk", group = OF_groups.keys()),

rule create_pep_links_for_orthofinder:
    input:
        "scratch/canonical_peptides/{org}.faa"
    output:
        "scratch/orthofinder_{group}/peptides/{org}.faa"
    shell:
        """
        ln -rs {input} {output}
        """

rule create_input_topology_for_orthofinder:
    input:
        intree = "imports/trees/species173_topology_with_scientific_names.nwk",
        rscrpt = "workflow/scripts/prepare_orthofinder/build_input_topology_for_orthofinder.R",
    output:
        outree = "scratch/orthofinder_{group}/topology_{group}.nwk",
        orglist = "scratch/orthofinder_{group}/orglist_{group}.txt",
    run:
        print(wildcards.group)
        df = OF_groups.get(wildcards.group, None)
        print(df)
        df["OrganismShortName"].to_csv(output.orglist, header = None, index = None)
        shell("Rscript {input.rscrpt} {output.orglist} {input.intree} {output.outree}")


