import pandas as pd
from os.path import basename

use_annotations = [
    "GCF",
    "Ensembl",
    "Custom"
]

annotated = (
    pd.read_csv("configs/species_table.tsv", sep = "\t")
    .query("SurrogateAnnotationSource in @use_annotations")
    .query("OrganismOrder177 > 27")
)
print(annotated)

bulkrna = (
    annotated
    .query("bulkrna107 == 'Yes'")
)
print(bulkrna)

scrna = (
    annotated
    .query("scrna35 == 'Yes'")
)
print(scrna)

rnaseq = (
    annotated
    .query("(scrna35 == 'Yes') or (bulkrna107 == 'Yes')")
)
print(rnaseq)



OF_groups = {
    "surrogate_bulkrna": bulkrna,
    "surrogate_scrna": scrna,
    "surrogate_rnaseq": rnaseq,
}

group_org_pairs = []
for grp, df in OF_groups.items():
    for org in df.SurrogateShortName.unique():
        group_org_pairs.append((grp, org))
group_org_pairs = pd.DataFrame(group_org_pairs, columns = ["grp", "org"])
print(group_org_pairs)



rule all:
    input:
        expand("scratch/orthofinder_{group}/peptides/{org}.faa", zip, group=group_org_pairs.grp, org=group_org_pairs.org),
        expand("scratch/orthofinder_{group}/topology_{group}.nwk", group = OF_groups.keys()),
        expand("scratch/orthofinder_{group}/run_orthofinder_{group}.sh", group = OF_groups.keys()),

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
        intree = "imports/trees/species177_topology_with_scientific_names.nwk",
        rscrpt = "workflow/scripts/prepare_orthofinder/build_input_topology_for_orthofinder.R",
    output:
        outree = "scratch/orthofinder_{group}/topology_{group}.nwk",
        orglist = "scratch/orthofinder_{group}/orglist_{group}.txt",
    run:
        print(wildcards.group)
        df = OF_groups.get(wildcards.group, None)
        print(df)
        df["SurrogateShortName"].drop_duplicates().to_csv(output.orglist, header = None, index = None)
        shell("Rscript {input.rscrpt} {output.orglist} {input.intree} {output.outree}")


rule create_orthofinder_sbatch:
    input:
        outree = "scratch/orthofinder_{group}/topology_{group}.nwk",
        orglist = "scratch/orthofinder_{group}/orglist_{group}.txt",
    output:
        sbatch = "scratch/orthofinder_{group}/run_orthofinder_{group}.sh",
    params:
        wd = "scratch/orthofinder_{group}",
        threads = 32,
        alg_threads = 16,
    shell:
        """
        ABS_WD=$(readlink -f {params.wd})
        ABS_PEPTIDES=$(readlink -f {params.wd}/peptides)
        ABS_TREE=$(readlink -f {input.outree})

        cat << EOF > {output.sbatch}
#!/bin/bash
#SBATCH --mem=256g
#SBATCH --cpus-per-task={params.threads}
#SBATCH --time=96:00:00

module load diamond OrthoFinder

cd $ABS_WD

orthofinder \\
    -f $ABS_PEPTIDES \\
    -s $ABS_TREE \\
    -t {params.threads} \\
    -a {params.alg_threads} \\
    -o results
EOF
        chmod +x {output.sbatch}
        """
