import pandas as pd
import numpy as np

configfile: "configs/config.yaml"
TRIMMED_FASTQ_DIR = config["trimmed_dir"]


ALL_SAMPLES_FILE = "configs/species107_rnaseq_samples.tsv"

samples = (
    pd.read_csv(ALL_SAMPLES_FILE, sep = "\t")
    .query("Source == 'Mine'")
    [["SampleID", "OrganismID", "RunType"]]
    .set_index("SampleID", drop = True)
)
print(samples)

rule all:
    input:
        expand("scratch/geo_submission/{org}_isoform_counts.csv",
                org = samples.OrganismID.unique()),
        expand("scratch/geo_submission/{org}_gene_abundance.csv",
                org = samples.OrganismID.unique()),


def get_isoform_quants_file(wc):
    organism = wc.org
    org_samples = samples.query("OrganismID == @organism")
    infiles = {
        sid: f"scratch/kallisto_quants/{sid}/abundance.tsv"
        for sid in org_samples.index
    }
    return infiles

def get_isoform_quants(sid, infile):
    df = (
        pd.read_csv(infile, sep = "\t")
        .rename(columns = {
            "est_counts": sid,
            "target_id": "transript_id",
        })
        .set_index("transript_id")
        [[sid]]
    )
    return(df)

rule combine_isoform_quants_per_org:
    input:
       unpack(get_isoform_quants_file)
    output:
        "scratch/geo_submission/{org}_kallisto_estimated_counts.csv",
    run:
        import pandas as pd
        df = (
            pd.concat([
                get_isoform_quants(sid, infile)
                for sid, infile in  input.items()
            ], axis = 1)
        )
        print(df)
        df.to_csv(output[0])


def get_gene_quants_file(wc):
    organism = wc.org
    org_samples = samples.query("OrganismID == @organism")
    infiles = {
        sid: f"scratch/kallisto_gene_quants/gene_expr__{sid}.tsv"
        for sid in org_samples.index
    }
    return infiles

def get_gene_quants(sid, infile):
    df = (
        pd.read_csv(infile, sep = "\t")
        .rename(columns = {
            "abndance": sid,
        })
        .set_index("gene_id")
        [[sid]]
    )
    return(df)

rule combine_gene_quants_per_org:
    input:
       unpack(get_gene_quants_file)
    output:
        "scratch/geo_submission/{org}_gene_abundance.csv",
    run:
        import pandas as pd
        df = (
            pd.concat([
                get_gene_quants(sid, infile)
                for sid, infile in  input.items()
            ], axis = 1)
        )
        print(df)
        df.to_csv(output[0])


