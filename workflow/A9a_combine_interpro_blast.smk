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

short2id = rnaseq.set_index("OrganismShortName")["OrganismID"].to_dict()

rule all:
    input:
        expand("scratch/mini_orthofinder/denovo_interpro/{org}_interpro.tsv", org = rnaseq.OrganismShortName),
        expand("scratch/mini_orthofinder/denovo_blast/{org}_blast.tsv", org = rnaseq.OrganismShortName),
        expand("scratch/mini_orthofinder/denovo_blast_interpro/{org}_denovo.tsv", org = rnaseq.OrganismShortName),
        "scratch/mini_orthofinder/annotated_vertebrate/fasta_list.txt",
        "scratch/mini_orthofinder/annotated_invertebrate/fasta_list.txt",


rule finalize_interpro_result:
    input:
        ips = lambda wc: f"{INTERPROSCAN_TOP_DIR}/{short2id[wc.org]}/{short2id[wc.org]}_denovo.fasta.transdecoder_nostar.tsv",
    output:
        "scratch/mini_orthofinder/denovo_interpro/{org}_interpro.tsv",
    run:
        import pandas as pd
        df = (
            pd.read_csv(input.ips, sep = "\t", header=None,
                usecols = [0, 2, 5],
                names = ["C0", "ProteinLen", "Remark"],
            )
            .assign(pep = lambda tdf: tdf.C0.str.split(".tsv:").str[0])
            .assign(ProteinID = lambda tdf: tdf.C0.str.split(".tsv:").str[1])
            [["ProteinID", "ProteinLen", "Remark", "pep"]]
        )
        print(df)
        df.to_csv(output[0], index = False, sep = "\t")


rule finalize_blast_result:
    input:
        "scratch/blast_query_proteins/{org}/best_db_{org}__query_maf_otx.tsv",
    output:
        "scratch/mini_orthofinder/denovo_blast/{org}_blast.tsv",
    run:
        import pandas as pd
        organism = short2id[wildcards.org]
        df = (
            pd.read_csv(input[0], sep = "\t")
            .rename(columns = {
                "saccver": "ProteinID",
                "slen": "ProteinLen",
            })
            .assign(Remark = "BlastFind")
            .assign(pep = f"{TRANSDECODER_RES_TOP}/{organism}/{organism}_denovo.fasta.transdecoder.pep")
            .assign(ProteinID = lambda tdf: tdf.ProteinID.str.split("__").str[1])
            .reset_index()
            [["ProteinID", "ProteinLen", "Remark", "pep"]]
        )
        print(df)
        df.to_csv(output[0], sep="\t", index=False)


rule consolidate_all_denovo:
    input:
        "scratch/mini_orthofinder/denovo_interpro/{org}_interpro.tsv",
        "scratch/mini_orthofinder/denovo_blast/{org}_blast.tsv",
    output:
        "scratch/mini_orthofinder/denovo_blast_interpro/{org}_denovo.tsv",
    run:
        import pandas as pd
        from Bio import SeqIO
        from Bio.Seq import Seq
        from Bio.SeqRecord import SeqRecord

        def add_protein_seqs(group: pd.DataFrame) -> pd.DataFrame:
            # group.name is the value of the "pep" column for this group
            pep_file = group.name

            print(f"Processing {pep_file} ...")

            required = group.ProteinID.tolist()

            # Read all sequences from this fasta once
            seq_dict = {}
            with open(pep_file, "r") as spe:
                for record in SeqIO.parse(spe, "fasta"):
                    if record.id in required:
                        seq_dict[record.id] = str(record.seq)

            if len(required) != len(seq_dict):
                print("Error: Not all sequences are found!")
                print("Needed:", required)
                print("Found:", seq_dict.keys())
                exit(-1)

            # Map ProteinID to sequences
            group["ProteinSeq"] = group["ProteinID"].map(seq_dict.get)

            return group

        df = (
            pd.concat([
                pd.read_csv(infile, sep = "\t")
                for infile in input
            ])
        )
        # keep best remark
        remarks = (
            df
            .groupby("ProteinID")
            .agg({"Remark": lambda x: ";".join(x.unique())})
            .reset_index()
        )
        print(remarks)

        df = (
            df[["ProteinID", "ProteinLen", "pep"]]
            .merge(remarks, how = "left")
            .assign(src_gene = lambda xdf: xdf.ProteinID.str.split("_").str[:-1].str.join("_"))
            .sort_values(["src_gene", "ProteinLen"], ascending=[True, False])
            .drop_duplicates("src_gene", keep="first")
            .groupby("pep", group_keys=False)
            .apply(add_protein_seqs)
            .assign(ProteinID = lambda tdf: f"{wildcards.org}Denovo__Unknown__" + tdf.ProteinID)
            [['ProteinID', 'ProteinLen', 'Remark', 'ProteinSeq']]
        )
        print(df)
        df.to_csv(output[0], sep ="\t", index = False)

