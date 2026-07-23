import pandas as pd
from os.path import basename

orthogroups_path = dict(
    invertebrate = "scratch/orthofinder_surrogate_bulk_invertebrate/OrthoFinder_results/Results_Feb25/Orthogroups/Orthogroups.tsv",
    jawless = "scratch/orthofinder_surrogate_bulk_jawless/OrthoFinder_results/Results_Feb25/Orthogroups/Orthogroups.tsv",
    jawed = "scratch/orthofinder_surrogate_bulk_jawed/results/Results_Dec31/Orthogroups/Orthogroups.tsv",
)
print(orthogroups_path)

rule all:
    input:
        expand("scratch/orthofinder_surrogate_bulk_{spgrp}/ortholog_tables/orthotab_atmost1gene.txt",
                spgrp = orthogroups_path.keys()),
        expand("exports/ortholog_tables/{spgrp}_orthotab_maf_focused.tsv",
                spgrp = orthogroups_path.keys()),

rule prepare_ortho_table:
    input:
        ogs = lambda wc: orthogroups_path[wc.spgrp]
    output:
        atmost1gene = "scratch/orthofinder_surrogate_bulk_{spgrp}/ortholog_tables/orthotab_atmost1gene.txt",
    script:
        "scripts/quantify_rnaseq/prepare_orthotable_bulkrna.py"


rule augment_maf_otx_to_ortho_table:
    input:
        atmost1gene = "scratch/orthofinder_surrogate_bulk_{spgrp}/ortholog_tables/orthotab_atmost1gene.txt",
        grand_list = "scratch/grand_list/grand_list.tsv",
    output:
        maf_otx_focused = "exports/ortholog_tables/{spgrp}_orthotab_maf_focused.tsv",
    script:
        "scripts/quantify_rnaseq/augment_maf_otx.py"

