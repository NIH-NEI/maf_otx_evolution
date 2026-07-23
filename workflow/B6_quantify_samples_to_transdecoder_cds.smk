print("Snakemake Python:", sys.executable)
import pandas as pd
import glob
from io import StringIO

TRANSDECODER_RES_TOP = "/data/VisionEvo/TransDecoderOutput/OrganismLevel"
TRANSDECODER_QUANT_TOP = "/data/VisionEvo/TransDecoderOutput/QuantSalmon"
TRIMMED_FASTQ_DIR = "/data/VisionEvo/NNRL_sequenced_retinal_transcriptomes/trimmed_fastqs"
RNASEQ_SAMPLES_INFO = "/data/VisionEvo/RNAseq_samples_info"

rnaseq = (
    pd.read_csv("configs/species181_table.tsv", sep = "\t")
    .query("bulkrna107 == 'Yes'")
)
print(rnaseq)

short2id = rnaseq.set_index("OrganismShortName")["OrganismID"].to_dict()
vertebrates = rnaseq.query("IsVertebrate == 'Yes'").iloc[0:1,:]
invertebrates = rnaseq.query("IsVertebrate == 'No'")
print(vertebrates.shape)

spgrp = {}
for idx, row in rnaseq.iterrows():
    spgrp[row.OrganismShortName] = "vertebrate" if row.IsVertebrate == "Yes" else "invertebrate"
print(spgrp)


rule all:
    input:
        expand(TRANSDECODER_RES_TOP + "/{org}/{org}_denovo.fasta.transdecoder.cds.gene_cds_map",
                org = rnaseq.OrganismID),
        expand(TRANSDECODER_QUANT_TOP + "/{org}/quants.done", org = rnaseq.OrganismID),

rule map_cds_to_gene:
    input:
        cds = TRANSDECODER_RES_TOP + "/{org}/{org}_denovo.fasta.transdecoder.cds",
    output:
        c2g = TRANSDECODER_RES_TOP + "/{org}/{org}_denovo.fasta.transdecoder.cds.gene_cds_map",
    shell:
        """
        cat {input.cds} \
            | grep ">" | sed 's/>//' \
            | sed 's/~~/\\t/' | sed 's/ /\\t/' \
            | awk '{{print $2"\\t"$1}}' \
            > {output.c2g}
        """

rule generate_samples_file:
    input:
        samps = "configs/species107_rnaseq_samples.tsv",
        orgs = "configs/species181_table.tsv",
    output:
        RNASEQ_SAMPLES_INFO + "/{org}_samples.tsv",
    run:
        import pandas as pd
        df = (
            pd.read_csv(input.samps, sep = "\t")
            .merge(
                (
                    pd.read_csv(input.orgs, sep = "\t")
                    [["OrganismID", "OrganismShortName"]]
                ),
                how = "left"
            )
            .query("OrganismID in @wildcards.org")
            .query("RunType == 'PE'")
            .assign(fq1 = lambda tdf: TRIMMED_FASTQ_DIR + "/" + tdf.SampleID + "_R1-trimmed.fastq.gz")
            .assign(fq2 = lambda tdf: TRIMMED_FASTQ_DIR + "/" + tdf.SampleID + "_R2-trimmed.fastq.gz")
            [["OrganismShortName", "SampleID", "fq1", "fq2"]]
        )
        print(df)
        df.to_csv(output[0], sep = "\t", index = False, header = False)

rule quantify_using_trinity_utility_salmon:
    input:
        cds = TRANSDECODER_RES_TOP + "/{org}/{org}_denovo.fasta.transdecoder.cds",
        c2g = TRANSDECODER_RES_TOP + "/{org}/{org}_denovo.fasta.transdecoder.cds.gene_cds_map",
        samples = RNASEQ_SAMPLES_INFO + "/{org}_samples.tsv",
    output:
        TRANSDECODER_QUANT_TOP + "/{org}/quants.done",
    params:
        wd = TRANSDECODER_QUANT_TOP + "/{org}",
    shell:
        """
        set -euo pipefail
        module load trinity
        cd {params.wd}
        $TRINITY_HOME/util/align_and_estimate_abundance.pl \
                --seqType fq  \
                --samples_file {input.samples} \
                --transcripts {input.cds} \
                --est_method salmon \
                --gene_trans_map {input.c2g} \
                --prep_reference
        touch {output}
        """

