import pandas as pd
from os.path import basename
from Bio import SeqIO


TOP_DIR = "scratch/hyphy_analysis"

genes = [
    "NRL",
    "MAFA",
    "MAFB",
    "CMAF",
]

rule all:
    input:
        expand(
            TOP_DIR + "/{basename}/{basename}__{part}__{aligner}__{postalign}.fna",
            basename = [f"jawed_{gene}" for gene in genes],
            part = ["full"],
            aligner = ["clustalo"], #, "mafft"],
            postalign = ["fullaln"],
        ),
        expand(
            TOP_DIR + "/{basename}/{basename}__{part}__{aligner}__{postalign}__inframe.fna",
            basename = [f"jawed_{gene}" for gene in genes],
            part = ["full"],
            aligner = ["clustalo"], #, "mafft"],
            postalign = ["fullaln"],
        ),

rule create_fasta:
    input:
        tsv = "scratch/final_list/{basename}.tsv",
    output:
        fasta = TOP_DIR + "/{basename}/{basename}__full.faa",
    run:
        import pandas as pd
        from Bio.SeqRecord import SeqRecord
        from Bio.Seq import Seq
        df = pd.read_csv(input.tsv, sep = "\t")
        outfile = output.fasta
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


rule align_by_clustal:
    input:
        fasta = TOP_DIR + "/{basename}__{part}.faa",
    output:
        fasta = TOP_DIR + "/{basename}__{part}__clustalo__fullaln.fasta",
    threads: 16
    shell:
        """
        module load clustalo
        clustalo -i {input.fasta} -o {output.fasta} --threads={threads} --outfmt=fasta
        """



rule subset_fna:
    input:
        pep_aln = TOP_DIR + "/{basename}__{part}__{aligner}__{postalign}.fasta",
        grnd_fna = "scratch/grand_list/grand_list.fna",
    output:
        cds_fna = TOP_DIR + "/{basename}__{part}__{aligner}__{postalign}.fna",
    shell:
        """
        module load seqkit

        # 1) Extract sequence IDs from peptide alignment headers
        #    (takes first token after '>' which matches Biopython rec.id / seqkit default name)
        awk '/^>/{{print substr($1,2)}}' {input.pep_aln} > {output.cds_fna}.ids

        # 2) Subset grand CDS fasta to only those IDs
        #    -n: match by name (header id), -f: file of ids, -o: output fasta
        seqkit grep -n -f {output.cds_fna}.ids {input.grnd_fna} -o {output.cds_fna}

        # Optional sanity check: ensure we got everything
        npep=$(grep -c '^>' {input.pep_aln} || true)
        ncds=$(grep -c '^>' {output.cds_fna} || true)
        echo "pep_aln seqs: $npep  cds subset seqs: $ncds"

        if [[ "$npep" -ne "$ncds" ]]; then
            echo "ERROR: count mismatch (npep=$npep, ncds=$ncds) for {wildcards.basename} {wildcards.part} {wildcards.aligner} {wildcards.postalign}" >&2
            sdiff <(grep '^>' {input.pep_aln} | sort) <(grep '^>' {output.cds_fna} | sort)
            exit 1
        fi
        """


rule align_by_pal2nal:
    input:
        pep_aln = TOP_DIR + "/{basename}__{part}__{aligner}__{postalign}.fasta",
        cds_fna = TOP_DIR + "/{basename}__{part}__{aligner}__{postalign}.fna",
    output:
        cds_aln = TOP_DIR + "/{basename}__{part}__{aligner}__{postalign}__inframe.fna",
    shell:
        r"""
        # 3) Run pal2nal on the subset CDS
        echo "Running pal2nl"
        echo "perl imports/pal2nal.v14/pal2nal.pl {input.pep_aln} {input.cds_fna} -output fasta > {output.cds_aln}"
        perl imports/pal2nal.v14/pal2nal.pl {input.pep_aln} {input.cds_fna} -output fasta > {output.cds_aln}
        """

rule check_cds_peptides_nucleotides:
    input:
        faa = "imports/peptides/{org}.fa",
        fna = "imports/nucleotides/{org}.fna",
    output:
        "imports/match_nuc_pep/{org}.tsv"
    run:
        from Bio.Seq import Seq
        from Bio import SeqIO
        for fna_rec, faa_rec in zip(SeqIO.parse(input.fna, "fasta"), SeqIO.parse(input.faa, "fasta")):
            if fna_rec.id != faa_rec.id:
                print("Names do not match: nuc", fna_rec.id, "vs pep", faa_rec.id)
            protein = fna_rec.seq.translate(to_stop=False)
            if "*" in protein[:-1]:  # Ignore final codon
                print(f"Internal stop codon found in {fna_rec.id}")
            elif protein[-1] != "*":
                print(f"No stop codon found in {fna_rec.id}")
            else:
                protein = protein[:-1]

            if str(faa_rec.seq) != protein:
                print(">", fna_rec.id, "translation do not match: nuc")
                print(protein)
                print("vs pep")
                print(faa_rec.seq)
