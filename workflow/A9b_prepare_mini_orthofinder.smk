print("Snakemake Python:", sys.executable)
import pandas as pd
import glob
from io import StringIO

TRANSDECODER_RES_TOP = "/data/VisionEvo/TransDecoderOutput/OrganismLevel"
INTERPROSCAN_TOP_DIR = "/data/VisionEvo/Interproscan"

rnaseq = (
    pd.read_csv("configs/species173_table.tsv", sep = "\t")
    .query("bulkrna107 == 'Yes'")
#    .iloc[0:1,:]
)
print(rnaseq)

vertebrates = rnaseq.query("IsVertebrate == 'Yes'")
invertebrates = rnaseq.query("IsVertebrate == 'No'")

short2id = rnaseq.set_index("OrganismShortName")["OrganismID"].to_dict()

rule all:
    input:
        "scratch/mini_orthofinder/annotated_vertebrate/fasta_list.txt",
        "scratch/mini_orthofinder/annotated_invertebrate/fasta_list.txt",
        expand("scratch/mini_orthofinder/denovo_vertebrate/{org}/fasta_list.txt", org = vertebrates.OrganismShortName),
        expand("scratch/mini_orthofinder/denovo_invertebrate/{org}/fasta_list.txt", org = invertebrates.OrganismShortName),


rule prepare_annotated_fasta_for_mini_orthofinder:
    input:
        tsv="scratch/annotated_proteins/annotated_maf_otx.tsv",
        ps="workflow/scripts/annotate_denovo/prepare_mini_orthofinder.py",
    output:
        "scratch/mini_orthofinder/annotated_{spgrp}/fasta_list.txt",
    shell:
        """
        set -euo pipefail
        txt={output}
        topdir=$(dirname $txt)
        fastadir="$topdir/peptides"
        python {input.ps} {input.tsv} $fastadir {wildcards.spgrp}
        ls $fastadir > {output}
        """

rule prepare_denovo_fasta_for_mini_orthofinder:
    input:
        tsv="scratch/mini_orthofinder/denovo_blast_interpro/{org}_denovo.tsv",
        ps="workflow/scripts/annotate_denovo/prepare_mini_orthofinder.py",
    output:
        "scratch/mini_orthofinder/denovo_{spgrp}/{org}/fasta_list.txt",
    shell:
        """
        set -euo pipefail
        txt={output}
        topdir=$(dirname $txt)
        fastadir="$topdir/peptides"
        python {input.ps} {input.tsv} $fastadir 
        ls $fastadir > {output}
        """
