import pandas as pd

configfile: "configs/config.yaml"
#ALL_SAMPLES_FILE = config["samples_file"]
DENOVO_ASSEMBLY_DIR = config["denovo_dir"]
TRIMMED_FASTQ_DIR = config["trimmed_dir"]
NCBI_DOWNLOADS_DIR = config["ncbi_downloads_dir"]

ALL_SAMPLES_FILE = "configs/species107_rnaseq_samples.tsv"
ALL_SPECIES_FILE = "configs/species_table.tsv"


organisms = pd.read_csv(ALL_SPECIES_FILE, sep = "\t")

samples = (
    pd.read_csv(ALL_SAMPLES_FILE, sep = "\t")
    [["SampleID", "OrganismID"]]
    .merge(organisms, how = "left")
    .query("OrganismOrder177 > 27")
    .set_index("SampleID", drop = True)
)
#print(samples)

organisms = organisms.set_index("OrganismShortName", drop = False)

multi_sample_species = (
    samples.reset_index().groupby("OrganismID").agg({"SampleID": "count"})
    .query("SampleID > 1")
    .index.tolist()
)
print(multi_sample_species)


ortho_files = dict(
    single_copy = "",
    atmost1gene = "",
    atmost1gene_atleast2org = "",
    atmost1gene_atleast2org_eachgrp = "",
    maf_focused = "scratch/orthofinder_surrogate_bulkrna/ortholog_tables/orthotab_maf_focused.tsv",
)

rule all:
    input:
        expand("scratch/kallisto_gene_quants/gene_expr__{sample}.tsv",
                sample = samples.index),
        #expand("scratch/kallisto_normalized_intraspecies/{org}/norm_expr_{org}.tsv",
        #        org = multi_sample_species),
        #"scratch/OrthoFinder/one_to_one_ortholog.tsv",
        #"scratch/kallisto_quants/gene_expr__Columba_livia_1.tsv",
        #"scratch/kallisto_normalized_quants/ortho~maf_focused/norm~tmm_tpm_edger/expr_samples__tmm_tpm_edger__maf_focused.tsv",
        #"scratch/kallisto_normalized_quants/ortho~maf_focused/norm~tmm_abundance_edger/expr_samples__tmm_abundance_edger__maf_focused.tsv",
        #"scratch/kallisto_normalized_quants/ortho~atmost1gene/norm~tmm_abundance_edger/expr_samples__tmm_abundance_edger__atmost1gene.tsv",
        "scratch/kallisto_normalized_interspecies/ortho~maf_focused/norm~tmm_abundance_edger/expr_matrix__tmm_abundance_edger__maf_focused.tsv",

def get_seq2gene(wc):
    sorg = samples['SurrogateShortName'][wc.sample]
    return f"scratch/transcripts_annotation/{sorg}/{sorg}_cds_pep_gene_info.tsv"

rule quantify_genes_single_sample:
    input:
        quants = "scratch/kallisto_quants/{sample}/abundance.tsv",
        seq2gene = get_seq2gene,
    output:
        gene = "scratch/kallisto_gene_quants/gene_expr__{sample}.tsv",
    script:
        "scripts/quantify_rnaseq/quantify_genes_tximport.R"


def dict_expand_intraspecies(wildcards):
    mysamples = samples.query("OrganismID == @wildcards.org")
    fpath = "scratch/kallisto_quants_mafs/gene_expr__{sample}.tsv"
    return dict(zip(mysamples.index, expand(fpath, sample = mysamples.index)))

rule normalize_expression_intraspecies:
    input:
        unpack(dict_expand_intraspecies),
        samples = ALL_SAMPLES_FILE,
    output:
        norm = "scratch/kallisto_normalized_intraspecies/{org}/norm_expr_{org}.tsv",
        dge = "scratch/kallisto_normalized_intraspecies/{org}/dge_{org}.rds",
        factor = "scratch/kallisto_normalized_intraspecies/{org}/norm_factors_{org}.tsv",
    script:
        "scripts/quantify_rnaseq/normalize_intraspecies_tmm.R"

def dict_expand_interspecies(wildcards):
    fpath = "scratch/kallisto_gene_quants/gene_expr__{sample}.tsv"
    return dict(zip(samples.index, expand(fpath, sample = samples.index)))

rule normalize_expression_interspecies:
    input:
        unpack(dict_expand_interspecies),
        ortho = lambda wc: ortho_files[wc.scportho],
        samples = ALL_SAMPLES_FILE,
        orgs = ALL_SPECIES_FILE,
    output:
        norm = "scratch/kallisto_normalized_interspecies/ortho~{scportho}/norm~{norm}/expr_matrix__{norm}__{scportho}.tsv",
        dge = "scratch/kallisto_normalized_interspecies/ortho~{scportho}/norm~{norm}/expr_matrix__{norm}__{scportho}.rds",
        factor = "scratch/kallisto_normalized_interspecies/ortho~{scportho}/norm~{norm}/norm_factors__{norm}__{scportho}.tsv",
    script:
        "scripts/quantify_rnaseq/normalize_interspecies_tmm_cpm.R"


