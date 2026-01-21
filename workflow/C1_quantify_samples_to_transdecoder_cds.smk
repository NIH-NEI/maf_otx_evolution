print("Snakemake Python:", sys.executable)
import pandas as pd
import glob
from io import StringIO

TRANSDECODER_RES_TOP = "/data/VisionEvo/TransDecoderOutput/OrganismLevel"
TRANSDECODER_QUANT_TOP = "/data/VisionEvo/TransDecoderOutput/QuantSalmon"
TRIMMED_FASTQ_DIR = "/data/VisionEvo/NNRL_sequenced_retinal_transcriptomes/trimmed_fastqs"
RNASEQ_SAMPLES_INFO = "/data/VisionEvo/RNAseq_samples_info"

rnaseq = (
    pd.read_csv("configs/species173_table.tsv", sep = "\t")
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
        orgs = "configs/species173_table.tsv",
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

#rule summarize_mini_orthofinder_denovo_vertebrate:
#    input:
#        info = "scratch/mini_orthofinder/denovo_blast_interpro/{org}_denovo.tsv",
#        ogs = lambda wc: f"scratch/mini_orthofinder/denovo_{spgrp[wc.org]}/{wc.org}/Orthogroups.tsv",
#        annot = "exports/curated_genes/draft_annotated_maf_otx_longformat.tsv",
#    output:
#        short = "scratch/mini_orthofinder/denovo_summary/{org}_mini_orthofinder_conclusion.tsv",
#        long = "scratch/mini_orthofinder/denovo_summary/{org}_mini_orthofinder_details.tsv",
#    script:
#        "scripts/annotate_denovo/summarize_mini_orthofinder_orthogroups.py"
#
#rule combine_denovo_finds:
#    input:
#        expand("scratch/mini_orthofinder/denovo_summary/{org}_mini_orthofinder_conclusion.tsv",
#                org = rnaseq.OrganismShortName),
#    output:
#        "exports/curated_genes/draft_denovo_maf_otx_longformat.tsv",
#    run:
#        import pandas as pd
#        df = (
#            pd.concat([
#                pd.read_csv(infile, sep = "\t")
#                for infile in input
#            ])
#        )
#        print(df)
#        df.to_csv(output[0], sep = "\t", index = False)
#
#rule merge_annotated_denovo:
#    input:
#        annotated = "exports/curated_genes/draft_annotated_maf_otx_longformat.tsv",
#        denovo = "exports/curated_genes/draft_denovo_maf_otx_longformat.tsv",
#    output:
#        long = "exports/curated_genes/draft_combined_maf_otx_longformat.tsv",
#        wide = "exports/curated_genes/draft_combined_maf_otx.tsv",
#    run:
#        import pandas as pd
#        import numpy as np
#        df = (
#            pd.concat(
#                {
#                    key: pd.read_csv(infile, sep = "\t").set_index("ProteinID")
#                    for key, infile in input.items()
#                },
#                names = ["AnnotationType", "ProteinID"],
#            )
#            .reset_index()
#            .assign(OrganismShortName = lambda tdf: tdf.ProteinID.str.split("__").str[0].str.replace("Denovo", ""))
#        )
#        print(df)
#        df.to_csv(output.long, sep = "\t", index = False)
#
#        large_mafs = ["CMAF", "MAFA", "MAFB", "NRL"]
#        small_mafs = ["MAFF", "MAFG", "MAFK"]
#        otxs = ["OTX1", "OTX2", "CRX"]
#
#        df = (
#            df
#            .assign(GeneGroup = lambda tdf: np.where(tdf.GeneGroup == "MAFL", "CMAF",
#                                            np.where(tdf.GeneGroup == "MAFS", "MAFF",
#                                            np.where(tdf.GeneGroup == "OTX", "OTX1", tdf.GeneGroup))))
#            .pivot_table(
#                values = "ProteinID",
#                index = ["OrganismShortName", "AnnotationType"],
#                columns = "GeneGroup",
#                aggfunc= ','.join,
#            )
#            .fillna("")
#            [large_mafs + small_mafs + otxs]
#        )
#        df["nMAFL"] = (df[large_mafs] != "").sum(axis=1)
#        df["nMAFS"] = (df[small_mafs] != "").sum(axis=1)
#        df["nOTX"] = (df[otxs] != "").sum(axis=1)
#        print(df)
#
#        annotated_sources = [
#            "GCF", "Ensembl", "Custom"
#        ]
#
#        orgs = (
#            pd.read_csv("configs/species173_table.tsv", sep = "\t")
#            [["OrganismShortName", "OrganismID", "OrganismColor", "bulkrna107"]]
#            .merge(
#                (
#                    pd.read_csv("configs/species173_annotations.tsv", sep = "\t")
#                    [["OrganismShortName", "AnnotationSource"]]
#                ),
#                how  ="left",
#            )
#            .query("((AnnotationSource in @annotated_sources) or (bulkrna107 == 'Yes'))")
#        )
#        print(orgs)
#
#        orgs = orgs.merge(df.reset_index(), how = "left")
#        print(orgs)
#        orgs.to_csv(output.wide, sep = "\t", index = False)
#
