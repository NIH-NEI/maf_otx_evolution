print("Snakemake Python:", sys.executable)
import pandas as pd
import glob
from io import StringIO
import os

TRANSDECODER_RES_TOP = "/data/VisionEvo/TransDecoderOutput/OrganismLevel"
INTERPROSCAN_TOP_DIR = "/data/VisionEvo/Interproscan"
TRANSDECODER_QUANT_TOP = "/data/VisionEvo/TransDecoderOutput/QuantSalmon"

rnaseq = (
    pd.read_csv("configs/species173_table.tsv", sep = "\t")
    .query("bulkrna107 == 'Yes'")
#    .query("OrganismShortName in ['chicken', 'quail']")
)
print(rnaseq)

short2id = rnaseq.set_index("OrganismShortName")["OrganismID"].to_dict()
id2short = rnaseq.set_index("OrganismID")["OrganismShortName"].to_dict()

vertebrates = rnaseq.query("IsVertebrate == 'Yes'").iloc[0:1,:]
invertebrates = rnaseq.query("IsVertebrate == 'No'")
print(vertebrates.shape)

print(glob.glob(TRANSDECODER_QUANT_TOP + "/Coturnix_japonica/*/quant.sf"))

spgrp = {}
for idx, row in rnaseq.iterrows():
    spgrp[row.OrganismShortName] = "vertebrate" if row.IsVertebrate == "Yes" else "invertebrate"
print(spgrp)

rule all:
    input:
        "exports/curated_genes/draft_combined_maf_otx.tsv",
        "scratch/denovo_summary_combined/denovo_maf_otx_details.tsv",
        TRANSDECODER_QUANT_TOP + "/Coturnix_japonica/quants_summary.tsv",

rule consolidate_denovo_quantification:
    input:
        lambda wc: glob.glob(TRANSDECODER_QUANT_TOP + f"/{wc.org}/*/quant.sf"),
    output:
        TRANSDECODER_QUANT_TOP + "/{org}/quants_summary.tsv",
    run:
        import pandas as pd
        import os
        protein_prefix = f"{id2short[wildcards.org]}Denovo__Unknown__"
        print(protein_prefix)

        dfs = []
        for infile in input:
            tdf = (
                pd.read_csv(infile, sep = "\t")
                .assign(SampleID = os.path.basename(os.path.dirname(infile)))
            )
            if not tdf.empty:
                dfs.append(tdf)

        print(dfs)

        if dfs:
            df = (
                pd.concat(dfs)
                .assign(ProteinID = lambda tdf: protein_prefix + tdf.Name)
                .groupby(["ProteinID"])
                .agg({"TPM": ["min", "median", "max"],
                    "NumReads": ["min", "median", "max"]})
            )
            df.columns = ["_".join(x) for x in df.columns]
        else:
            df = pd.DataFrame([], columns = [
                "ProteinID", "TPM_min", "TPM_median", "TPM_max",
                "NumReads_min", "NumReads_median", "NumReads_max"
            ])
        print(df)

        df.to_csv(output[0], sep = "\t")

rule summarize_mini_orthofinder_denovo_vertebrate:
    input:
        info = "scratch/mini_orthofinder/denovo_blast_interpro/{org}_denovo.tsv",
        ogs = lambda wc: f"scratch/mini_orthofinder/denovo_{spgrp[wc.org]}/{wc.org}/Orthogroups.tsv",
        annot = "exports/curated_genes/draft_annotated_maf_otx_longformat.tsv",
        quant = lambda wc: TRANSDECODER_QUANT_TOP + f"/{short2id[wc.org]}/quants_summary.tsv",
    output:
        short = "scratch/mini_orthofinder/denovo_summary/{org}_mini_orthofinder_conclusion.tsv",
        long = "scratch/mini_orthofinder/denovo_summary/{org}_mini_orthofinder_details.tsv",
    script:
        "scripts/annotate_denovo/summarize_mini_orthofinder_orthogroups.py"

rule combine_denovo_finds_details:
    input:
        expand("scratch/mini_orthofinder/denovo_summary/{org}_mini_orthofinder_details.tsv",
                org = rnaseq.OrganismShortName),
    output:
        "scratch/denovo_summary_combined/denovo_maf_otx_details.tsv",
    run:
        import pandas as pd
        df = (
            pd.concat([
                pd.read_csv(infile, sep = "\t")
                for infile in input
            ])
        )
        print(df)
        df.to_csv(output[0], sep = "\t", index = False)

rule combine_denovo_finds:
    input:
        expand("scratch/mini_orthofinder/denovo_summary/{org}_mini_orthofinder_conclusion.tsv",
                org = rnaseq.OrganismShortName),
    output:
        "exports/curated_genes/draft_denovo_maf_otx_longformat.tsv",
    run:
        import pandas as pd
        df = (
            pd.concat([
                pd.read_csv(infile, sep = "\t")
                for infile in input
            ])
        )
        print(df)
        df.to_csv(output[0], sep = "\t", index = False)

rule merge_annotated_denovo:
    input:
        annotated = "exports/curated_genes/draft_annotated_maf_otx_longformat.tsv",
        denovo = "exports/curated_genes/draft_denovo_maf_otx_longformat.tsv",
    output:
        long = "exports/curated_genes/draft_combined_maf_otx_longformat.tsv",
        wide = "exports/curated_genes/draft_combined_maf_otx.tsv",
    run:
        import pandas as pd
        import numpy as np
        df = (
            pd.concat(
                {
                    key: pd.read_csv(infile, sep = "\t").set_index("ProteinID")
                    for key, infile in input.items()
                },
                names = ["AnnotationType", "ProteinID"],
            )
            .reset_index()
            .assign(OrganismShortName = lambda tdf: tdf.ProteinID.str.split("__").str[0].str.replace("Denovo", ""))
        )
        print(df)
        df.to_csv(output.long, sep = "\t", index = False)

        large_mafs = ["CMAF", "MAFA", "MAFB", "NRL"]
        small_mafs = ["MAFF", "MAFG", "MAFK"]
        otxs = ["OTX1", "OTX2", "CRX"]

        df = (
            df
            .assign(GeneGroup = lambda tdf: np.where(tdf.GeneGroup == "MAFL", "CMAF",
                                            np.where(tdf.GeneGroup == "MAFS", "MAFF",
                                            np.where(tdf.GeneGroup == "OTX", "OTX1", tdf.GeneGroup))))
            .pivot_table(
                values = "ProteinID",
                index = ["OrganismShortName", "AnnotationType"],
                columns = "GeneGroup",
                aggfunc= ','.join,
            )
            .fillna("")
            [large_mafs + small_mafs + otxs]
        )
        df["nMAFL"] = (df[large_mafs] != "").sum(axis=1)
        df["nMAFS"] = (df[small_mafs] != "").sum(axis=1)
        df["nOTX"] = (df[otxs] != "").sum(axis=1)
        print(df)

        annotated_sources = [
            "GCF", "Ensembl", "Custom"
        ]

        orgs = (
            pd.read_csv("configs/species173_table.tsv", sep = "\t")
            [["OrganismShortName", "OrganismID", "OrganismColor", "bulkrna107"]]
            .merge(
                (
                    pd.read_csv("configs/species173_annotations.tsv", sep = "\t")
                    [["OrganismShortName", "AnnotationSource"]]
                ),
                how  ="left",
            )
            .query("((AnnotationSource in @annotated_sources) or (bulkrna107 == 'Yes'))")
        )
        print(orgs)

        orgs = orgs.merge(df.reset_index(), how = "left")
        print(orgs)
        orgs.to_csv(output.wide, sep = "\t", index = False)







