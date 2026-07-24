from io import StringIO
import pandas as pd
import glob


rule all:
    input:
        "scratch/grand_list_interpro/grand_list_interproscan.tsv",

rule extract_mafs_otx_from_interproscan:
    input:
        sig = "configs/interpro_maf_otx_signatures.txt"
    output:
        ips = "scratch/grand_list_interpro/grand_list_interproscan.tsv",
    params:
        indir = "scratch/grand_list_interpro/interpro_run_dir",
    shell:
        """
        (grep -H -w -f {input.sig} {params.indir}/*.tsv > {output}) || true
        """

