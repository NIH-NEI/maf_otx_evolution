
import pandas as pd
import glob
from io import StringIO

TRANSDECODER_RES_TOP = "/data/VisionEvo/TransDecoderOutput/OrganismLevel"
INTERPROSCAN_TOP_DIR = "/data/VisionEvo/Interproscan"

rnaseq = (
    pd.read_csv("configs/species173_table.tsv", sep = "\t")
    .query("bulkrna107 == 'Yes'")
)
print(rnaseq)

short2id = rnaseq.set_index("OrganismShortName")["OrganismID"].to_dict()
vertebrates = rnaseq.query("IsVertebrate == 'Yes'").iloc[0:1,:]
invertebrates = rnaseq.query("IsVertebrate == 'No'")
print(vertebrates.shape)

rule all:
    input:
        "scratch/grand_list/grand_list_with_tpm.tsv",
        #"scratch/grand_list/grand_list.faa",


rule consolidate_mini_ortho_probables:
    input:
        finds = "scratch/mini_orthofinder/denovo_summary/{org}_mini_orthofinder_details.tsv",
        info = "scratch/mini_orthofinder/denovo_blast_interpro/{org}_denovo.tsv",
    output:
        "scratch/mini_orthofinder/denovo_summary/{org}_mini_orthofinder_with_seq.tsv",
    run:
        import pandas as pd
        df = (
            pd.read_csv(input.finds, sep = "\t")
            .groupby("ProteinID")
            .agg({"GeneGroup": " or ".join, "TPM_max": "max"})
            .rename(columns = {"GeneGroup": "Remark"})
            .reset_index()
        )
        print(df)
        info = (
            pd.read_csv(input.info, sep = "\t")
            [["ProteinID", "ProteinLen", "ProteinSeq"]]
        )
        print(info)

        df = (
            df.merge(info, how = "left")
            [["ProteinID", "ProteinLen", "Remark", "TPM_max", "ProteinSeq"]]
        )
        print(df)

        df.to_csv(output[0], sep = "\t", index = False)


rule make_grand_list:
    input:
        "scratch/annotated_proteins/annotated_maf_otx.tsv",
        expand("scratch/mini_orthofinder/denovo_summary/{org}_mini_orthofinder_with_seq.tsv", org = rnaseq.OrganismShortName),
    output:
        "scratch/grand_list/grand_list_with_tpm.tsv"
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
        print(df)
        print(df.columns)

        outgroups = (
            pd.read_csv("configs/outgroups.tsv", sep = "\t")
        )
        print(outgroups)

        outgroups_to_add = list(set(outgroups.Outgroup) - set(df.ProteinID))
        print(outgroups_to_add)

        outgroup_orgs = [x.split("__")[0] for x in outgroups_to_add]
        print(outgroup_orgs)


        outgroups_dict = {
            "ProteinID": outgroups_to_add,
            "pep": [
                f"scratch/canonical_peptides/{org}.faa"
                for org in outgroup_orgs
            ]
        }

        outgroups = (
            pd.DataFrame(outgroups_dict)
            .groupby("pep", group_keys=False)
            .apply(add_protein_seqs)
            .assign(ProteinLen = lambda tdf: tdf.ProteinSeq.apply(lambda x: len(x)))
            .assign(Remark = "Outgroup")
            [['ProteinID', 'ProteinLen', 'Remark', 'ProteinSeq']]
        )
        print(outgroups)

        noor_amey_file = "/data/pals2/old_maf_zenodo/amphi_tuni_hag_lamp/noor_amy_combined.faa"
        noor_amey_headers  = [record.id for record in SeqIO.parse(noor_amey_file, "fasta")]
        noor_amey_dict = {
            "ProteinID": noor_amey_headers,
            "pep": [
                noor_amey_file
                for org in noor_amey_headers
            ]
        }
        noor_amey = (
            pd.DataFrame(noor_amey_dict)
            .groupby("pep", group_keys=False)
            .apply(add_protein_seqs)
            .assign(ProteinLen = lambda tdf: tdf.ProteinSeq.apply(lambda x: len(x)))
            .assign(Remark = "NoorAmey")
            [['ProteinID', 'ProteinLen', 'Remark', 'ProteinSeq']]
        )
        print(noor_amey)

        df = (
            pd.concat([df, outgroups, noor_amey])
            .assign(ProteinID = lambda tdf: tdf.ProteinID.apply(lambda x: x.replace("Denovo__Unknown__", "__Unknown__")))
        )
        print(df)
        df.to_csv(output[0], sep = "\t", index = False)


rule tsv2fasta:
    input:
        "scratch/grand_list/grand_list.tsv",
    output:
        "scratch/grand_list/grand_list.faa",
    run:
        import pandas as pd
        from Bio import SeqIO
        from Bio.Seq import Seq
        from Bio.SeqRecord import SeqRecord
        df = (
            pd.read_csv(input[0], sep = "\t")
        )
        print(df)

        outfile = output[0]
        records = []
        for _, row in df.iterrows():
            if pd.isna(row["ProteinSeq"]):
                continue  # skip missing sequences

            rec = SeqRecord(
                Seq(row["ProteinSeq"]),
                id=row["ProteinID"],
                description=""   # no trailing description
            )
            records.append(rec)

        SeqIO.write(records, outfile, "fasta")
        print(f"Wrote {len(records)} sequences to {outfile}")

