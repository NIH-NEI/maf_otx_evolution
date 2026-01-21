from io import StringIO
import pandas as pd
import glob

TRANSCODER_OUTPUT_DIR = "/data/VisionEvo/TransDecoderOutput/OrganismLevel"
INTERPROSCAN_TOP_DIR = "/data/VisionEvo/Interproscan"

category_list = "Amphioxus|Arthropod|Hagfish|Jellyfish|Lamprey|Mollusc|Tunicate".split("|")

rnaseq = (
    pd.read_csv("configs/species173_table.tsv", sep = "\t")
#    .merge((
#        pd.read_csv("configs/species173_annotations.tsv", sep = "\t")
#        [["OrganismID", "AnnotationSource", "AnnotationName", "AnnotationRelease"]]
#    ), how = "left")
    .query("bulkrna107 == 'Yes'")
)
print(rnaseq)

rule all:
    input:
        expand(
            INTERPROSCAN_TOP_DIR + "/{organism}/{organism}_denovo.fasta.transdecoder_nostar.txt",
            zip,
            organism = rnaseq.OrganismID,
        )


rule run_interproscan:
    input:
        TRANSCODER_OUTPUT_DIR + "/{sample}/{sample}_denovo.fasta.transdecoder.pep",
    output:
        tmp = TRANSCODER_OUTPUT_DIR + "/{organism}/{sample}_denovo.fasta.transdecoder_nostar.pep",
        start = INTERPROSCAN_TOP_DIR + "/{organism}/{sample}_denovo.fasta.transdecoder_nostar.txt",
    params:
        outdir = INTERPROSCAN_TOP_DIR + "/{organism}",
        prefix = "{sample}_denovo.fasta.transdecoder_nostar",
    shell:
        """
        module load interproscan
        cat {input} | sed 's/*//' > {output.tmp}
        cd {params.outdir}
        interproscan {output.tmp} {params.prefix} 1000
        cd {params.prefix}
        sbatch interproscan.batch
        touch {output.start}
        """
