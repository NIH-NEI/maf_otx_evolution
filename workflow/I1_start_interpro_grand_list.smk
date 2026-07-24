
import pandas as pd
import glob
from io import StringIO

rule all:
    input:
        "scratch/grand_list_interpro/interpro_run_dir/interproscan.batch",

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

rule start_interproscan:
    input:
        "scratch/grand_list/grand_list.faa",
    output:
        scrpt = "scratch/grand_list_interpro/interpro_run_dir/interproscan.batch",
        faa = "scratch/grand_list/grand_list_nostar.faa",
    params:
        outdir = "scratch/grand_list_interpro/interpro_run_dir",
    shell:
        """
        module load interproscan
        sed 's/*//g' {input} > {output.faa}
        rm -rf {params.outdir}
        interproscan --goterms --pathways -f tsv,json {output.faa} {params.outdir} 600
        """

