import pandas as pd

ALL_SAMPLES_FILE = "configs/species107_rnaseq_samples.tsv"
ALL_SPECIES_FILE = "configs/species181_table.tsv"


organisms = pd.read_csv(ALL_SPECIES_FILE, sep = "\t")

spgroups = dict(
    invertebrate = [
        "Nonchordates",
        "Protochordates",
    ],
    jawless = [
        "Agnathans (hagfish and lampreys)",
    ],
    jawed = [
        "Cartilaginous fishes",
        "Non-teleost ray-finned fishes",
        "Teleost ray-finned fishes",
        "Lobe-finned fishes",
        "Amphibians",
        "Turtles",
        "Crocodiles",
        "Squamata (lizards and snakes)",
        "Birds",
        "Mammals",
    ]
)

samples = (
    pd.read_csv(ALL_SAMPLES_FILE, sep = "\t")
    [["SampleID", "OrganismID"]]
    .merge(organisms, how = "left")
    .set_index("SampleID", drop = True)
)
#print(samples)


orthotab_path = dict(
    invertebrate = "exports/ortholog_tables/invertebrate_orthotab_maf_focused.tsv",
    jawless = "exports/ortholog_tables/jawless_orthotab_maf_focused.tsv",
    jawed = "exports/ortholog_tables/jawed_orthotab_maf_focused.tsv",
)
print(orthotab_path)

rule all:
    input:
        expand(
            "exports/kallisto_normalized/expr_matrix__{spgrp}.tsv",
            spgrp = ["jawed", "jawless", "invertebrate"],
        ),


def dict_expand_interspecies(wildcards):
    fpath = "scratch/kallisto_gene_quants/gene_expr__{sample}.tsv"
    suborgs = spgroups[wildcards.spgrp] 
    subsamples = samples.query("OrganismColor in @suborgs")
    return dict(zip(subsamples.index, expand(fpath, sample = subsamples.index)))

rule normalize_expression_interspecies:
    input:
        unpack(dict_expand_interspecies),
        ortho = lambda wc: orthotab_path[wc.spgrp],
        samples = ALL_SAMPLES_FILE,
        orgs = ALL_SPECIES_FILE,
    output:
        norm = "exports/kallisto_normalized/expr_matrix__{spgrp}.tsv",
        dge = "scratch/kallisto_normalized_{spgrp}/expr_matrix__{spgrp}.rds",
        factor = "scratch/kallisto_normalized_{spgrp}/norm_factors__{spgrp}.tsv",
    script:
        "scripts/quantify_rnaseq/normalize_interspecies_tmm_cpm.R"


