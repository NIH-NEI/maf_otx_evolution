import pandas as pd
import numpy as np

configfile: "configs/config.yaml"
TRIMMED_FASTQ_DIR = config["trimmed_dir"]


ALL_SAMPLES_FILE = "configs/species107_rnaseq_samples.tsv"
ALL_SPECIES_FILE = "configs/species181_table.tsv"


organisms = pd.read_csv(ALL_SPECIES_FILE, sep = "\t")

samples = (
    pd.read_csv(ALL_SAMPLES_FILE, sep = "\t")
    [["SampleID", "OrganismID", "RunType"]]
    .merge(organisms, how = "left")
    .set_index("SampleID", drop = True)
)
#print(samples)

organisms = organisms.set_index("OrganismShortName", drop = False)

rule all:
    input:
        expand("scratch/kallisto_index/{org}_kallisto.idx",
                org = samples.OrganismShortName.unique()),
        expand("scratch/kallisto_quants/{sample}/abundance.tsv",
                sample = samples.index),
        expand("scratch/kallisto_gene_quants/gene_expr__{sample}.tsv",
                sample = samples.index),
        "scratch/kallisto_quants/all_samples_quality.tsv",


def get_cds_file(wc):
    annot_src = organisms["SurrogateAnnotationSource"][wc.org]
    annot_rel =  organisms["SurrogateAnnotationRelease"][wc.org]
    annot_name =  organisms["SurrogateAnnotationName"][wc.org]
    annot_org =  organisms["SurrogateShortName"][wc.org]
    annot_id =  organisms["SurrogateID"][wc.org]
    if annot_src == "GCF":
        return f"imports/genomes_annotations/{annot_org}/{annot_rel}_{annot_name}_cds_from_genomic.fna.gz"
    if annot_src == "Ensembl":
        return f"imports/genomes_annotations/{annot_org}/{annot_id}.{annot_name}.cds.all.fa.gz"
    if annot_src == "Custom":
        return {
            "haha": "hoho",
        }[annot_org]
    # It should not come here at all
    print("Error! Unknown annotation source")
    quit(-1)

rule kallisto_index:
    input:
        get_cds_file,
    output:
        "scratch/kallisto_index/{org}_kallisto.idx",
    log:
        out = 'scratch/kallisto_index/{org}_kallisto_stdout.log',
        err = 'scratch/kallisto_index/{org}_kallisto_stderr.log'
    shell:
        """
        module load kallisto/0.48.0
        kallisto index -i {output} {input} 2>{log.err} >{log.out}
        """

def get_input_for_kallisto(wc):
    run_type = samples["RunType"][wc.sample]
    my_species = samples["OrganismShortName"][wc.sample]
    infiles = {
        "idx": f"scratch/kallisto_index/{my_species}_kallisto.idx"
    }
    if run_type == "PE":
        infiles["fq1"] = f"{TRIMMED_FASTQ_DIR}/{wc.sample}_R1-trimmed.fastq.gz",
        infiles["fq2"] = f"{TRIMMED_FASTQ_DIR}/{wc.sample}_R2-trimmed.fastq.gz",
    else:
        infiles["fq1"] = f"{TRIMMED_FASTQ_DIR}/{wc.sample}-trimmed.fastq.gz",

    return infiles


rule kallisto_quant:
    input:
        unpack(get_input_for_kallisto)
    output:
        "scratch/kallisto_quants/{sample}/abundance.tsv",
    threads: 8
    log:
        out = 'scratch/kallisto_quants/{sample}/kallisto_stdout.log',
        err = 'scratch/kallisto_quants/{sample}/kallisto_stderr.log',
    params:
        extra = "-b 30",
        run_type = lambda wc: samples["RunType"][wc.sample]
    run:
        import os
        options = params.extra + f" --threads={threads}"
        infiles = f"{input.fq1}"
        options += " --single --single-overhang -l 200 -s 20" if params.run_type == "SE" else ""
        infiles += f" {input.fq2}" if params.run_type == "PE" else ""
        outdir = os.path.dirname(os.path.realpath(f"{output}"))
        cmd = "module load kallisto/0.48.0; "
        cmd += f"kallisto quant -i {input.idx} -o {outdir} {options} {infiles} "
        cmd += f"2>{log.err} >{log.out}"
        shell(cmd)


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


rule get_alignment_quality:
    input:
        expand("scratch/kallisto_quants/{sample}/kallisto_stderr.log",
                sample = samples.index)
    output:
        "scratch/kallisto_quants/all_samples_quality.tsv",
    script:
        "scripts/quantify_rnaseq/gather_kallisto_alignment_percentage.py"

