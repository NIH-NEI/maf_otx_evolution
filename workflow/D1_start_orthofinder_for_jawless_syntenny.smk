
import pandas as pd
import glob
from io import StringIO

orgs = [
    "acornWorm",
    "seaurchin",
    "echinoderm",
    "amphioxus",
    "vaseTunicate",
    "seaSquirt",
    "brownHagfish",
    "atlanticHagfish",
    "inshoreHagfish",
    "seaLamprey",
    "europeanRiverLamprey",
    "europeanBrookLamprey",
    "fareasternBrookLamprey",
    "bambooShark",
    "spottedGar",
    "zebrafish",
    "northernPike",
    "chicken",
    "human",
]

rule all:
    input:
        expand("scratch/canonical_peptides_unique/{org}.faa", org = orgs),
        expand("scratch/orthofinder_diverse/peptides/{org}.faa", org = orgs),
        "scratch/orthofinder_diverse/run_orthofinder_diverse.sh",

rule keep_unique_proteins:
    input:
        faa = "scratch/canonical_peptides/{org}.faa",
        tsv = "scratch/canonical_peptides/{org}_info.tsv",
    output:
        faa = "scratch/canonical_peptides_unique/{org}.faa",
        tsv = "scratch/canonical_peptides_unique/{org}_info.tsv",
        dup = "scratch/canonical_peptides_unique/{org}_duplicates_per_symbol.tsv",
    run:
        import pandas as pd
        from Bio import SeqIO
        from Bio.Seq import Seq
        from Bio.SeqRecord import SeqRecord
        import re

        organism = wildcards.org

        df = (
            pd.read_csv(input.tsv, sep = "\t")
            .assign(ShortPID = lambda tdf: tdf.ProteinID)
            .assign(ProteinID = lambda tdf: organism + "__" + tdf.GeneSymbol + "__" + tdf.ProteinID)
            [["GeneSymbol", "ProteinID", "ShortPID"]]
        )
        print(df)
        print(df.query("GeneSymbol == 'ZYX'"))
        df = df.drop_duplicates()
        print(df)
        print(df.query("GeneSymbol == 'ZYX'"))

        df = (
            df
            .sort_values(by = ["GeneSymbol", "ProteinID"])
            .groupby("GeneSymbol")
            .agg({
                "ProteinID": ["count", "first"],
                "ShortPID": ",".join,
            })
        )
        df.columns = ["nProtein", "ProteinID", "ShortPIDs"]
        dup = df.query("nProtein > 1")
        if not df.query("nProtein > 1").empty:
            print("!!!! WARNING !!! There are multiple proteins with same gene symbol!")
            print(dup)
        dup.to_csv(output.dup, sep = "\t", index = True)

        df = df.reset_index()
        print(df.query("GeneSymbol == 'NRL'"))


        results = []
        for record in SeqIO.parse(input.faa, "fasta"):
            #print(record)
            seq = str(record.seq)
            cleaned = re.sub(r'[^ACDEFGHIKLMNPQRSTVWYXBZJU]', '', seq.upper())
            #print("orignal:", seq)
            #print("cleaned:", cleaned)
            #results.append([record.id, str(record.seq).replace("*", "")])
            results.append([record.id, int(seq==cleaned), cleaned])

        sequences = pd.DataFrame(
            results, columns = ["ProteinID", "IsCleaned", "ProteinSeq"]
        )
        print(sequences)

        df = df.merge(sequences, how = "left")
        print(df)
        print(df.query("GeneSymbol == 'NRL'"))

        def avoid_generic_symbol(pids):
            x = [s.split("__", 1)[1] for s in pids]
            generic = [s for s in x if s.startswith("LOC")]
            given = [s for s in x if not s.startswith("LOC")]
            #if len(pids) > 1:
            #    print(pids)
            #    print(generic, given)
            return ",".join(given + generic)

        df = (
            df.groupby("ProteinSeq")
            .agg({
                "GeneSymbol": "count",
                "ProteinID": avoid_generic_symbol,
            })
            .reset_index()
            .assign(AlternatePIDs = lambda tdf: tdf.ProteinID.str.split(",", n=1).str[1])
            .assign(ProteinID = lambda tdf: organism + "__" + tdf.ProteinID.str.split(",", n=1).str[0])
            [["ProteinID", "AlternatePIDs", "GeneSymbol", "ProteinSeq"]]
            .rename(columns = {"GeneSymbol": "nPIDs"})
        )
        print(df.query("nPIDs > 1"))
        print(df.nPIDs.max())
        df.to_csv(output.tsv, index = None, sep = "\t")

        records = []
        outfile = output.faa
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


rule create_pep_links_for_orthofinder:
    input:
        "scratch/canonical_peptides_unique/{org}.faa"
    output:
        "scratch/orthofinder_diverse/peptides/{org}.faa"
    shell:
        """
        ln -rs {input} {output}
        """


rule create_orthofinder_sbatch:
    input:
        expand("scratch/orthofinder_diverse/peptides/{org}.faa", org = orgs),
    output:
        sbatch = "scratch/orthofinder_diverse/run_orthofinder_diverse.sh",
    params:
        wd = "scratch/orthofinder_diverse",
        threads = 32,
        alg_threads = 16,
    shell:
        """
        ABS_WD=$(readlink -f {params.wd})
        ABS_PEPTIDES=$(readlink -f {params.wd}/peptides)

        cat << EOF > {output.sbatch}
#!/bin/bash
#SBATCH --mem=256g
#SBATCH --cpus-per-task={params.threads}
#SBATCH --time=96:00:00

module load diamond OrthoFinder

cd $ABS_WD

orthofinder \\
    -f $ABS_PEPTIDES \\
    -t {params.threads} \\
    -a {params.alg_threads} \\
    -o results
EOF
        chmod +x {output.sbatch}
        """
