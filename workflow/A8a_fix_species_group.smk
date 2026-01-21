from io import StringIO
import pandas as pd
import glob

TRANSCODER_OUTPUT_DIR = "/data/VisionEvo/TransDecoderOutput/OrganismLevel"
INTERPROSCAN_TOP_DIR = "/data/VisionEvo/Interproscan"

organisms = (
    pd.read_csv("configs/species173_table.tsv", sep = "\t")
    .query("OrganismColor in ['Cartilaginous fishes', 'Non-teleost ray-finned fishes']")
)
print(organisms)
rnaseq = organisms.query("bulkrna107 == 'Yes'")
print(rnaseq)

annotated = (
    pd.read_csv("exports/curated_genes/final_annotated_MAFL.tsv", sep = "\t")
    [["OrganismShortName", "CMAF", "MAFA", "MAFB", "NRL"]]
    .melt(id_vars = "OrganismShortName", value_name = "ProteinID", var_name = "GeneGroup") 
    .dropna(subset="ProteinID")
    .assign(ProteinID = lambda tdf: tdf.ProteinID.apply(lambda x: x.split(",")))
    .explode("ProteinID")
    [["OrganismShortName", "ProteinID"]]
)
print(annotated)

def get_denovo_proteins(org):
    interpro_res_file = f"/data/VisionEvo/Interproscan/{org}/{org}_denovo.fasta.transdecoder_nostar.tsv"
    df = (
        pd.read_csv(interpro_res_file, header = None, usecols = [0,5], names=["ProteinID", "signature"])
        .fillna("")
#        .query('signature.str.contains("bZIP_Maf")')
    )
    print(df)
    print(df.signature.unique())
    quit()

denovo = (
    pd.concat([
        get_denovo_proteins(org_id)
        for org_id in rnaseq.OrganismID
    ])
)
print(denovo)
quit()

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
        cd {params.prefix}
        sbatch interproscan.batch
        interproscan {output.tmp} {params.prefix} 1000
        touch {output.start}
        """
