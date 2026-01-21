from io import StringIO
import pandas as pd
import glob

INTERPROSCAN_TOP_DIR = "/data/VisionEvo/Interproscan"

rnaseq = (
    pd.read_csv("configs/species173_table.tsv", sep = "\t")
    .query("bulkrna107 == 'Yes'")
)
print(rnaseq)

rule all:
    input:
        expand(INTERPROSCAN_TOP_DIR + "/{organism}/{organism}_denovo.fasta.transdecoder_nostar.tsv",
                organism = rnaseq.OrganismID)

rule fish_mafs_otx_from_interproscan:
    input:
        sig = "configs/interpro_maf_otx_signatures.txt"
    output:
        INTERPROSCAN_TOP_DIR + "/{organism}/{organism}_denovo.fasta.transdecoder_nostar.tsv",
    params:
        indir = INTERPROSCAN_TOP_DIR + "/{organism}/{organism}_denovo.fasta.transdecoder_nostar",
    shell:
        """
        (grep -H -w -f {input.sig} {params.indir}/*.tsv > {output}) || true
        """

