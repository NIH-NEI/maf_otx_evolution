import pandas as pd
import glob

TRANSDECODER_RES_TOP = "/data/VisionEvo/TransDecoderOutput/OrganismLevel"
DENOVO_ASSEMBLY_DIR = "/data/VisionEvo/NNRL_sequenced_retinal_transcriptomes/Analysis/Denovo_assemblies-incomplete/Nucleotide_All"

rnaseq = (
    pd.read_csv("configs/species181_table.tsv", sep = "\t")
    .query("bulkrna107 == 'Yes'")
)
print(rnaseq)

rule all:
    input:
        expand("scratch/blastdb_per_organism/{org}/{org}", org  = rnaseq.OrganismID.unique()),

#
rule make_blast_db:
    input:
        TRANSDECODER_RES_TOP + "/{org}/{org}_denovo.fasta.transdecoder.pep",
    output:
        touch("scratch/blastdb_per_organism/{org}/{org}")
    shell:
        """
        module load blast
        organism={wildcards.org}
        cat {input} \
            | awk -v prefix=$organism '/>/{{sub(">","&"prefix"__")}}1' \
            | makeblastdb -blastdb_version 5 \
                -title "Denovo proteins of $organism" \
                -input_type "fasta" \
                -dbtype prot \
                -out {output}
        """

rule enlist_sequences_in_blast_db:
	input:
        "scratch/blastdb_per_organism/{org}/{org}"
	output:
        "scratch/blastdb_per_organism/{org}/{org}_seqnames.txt"
	shell:
		"""
		module load blast
		blastdbcmd \
			-db {input} \
			-entry all -outfmt "%o %t" \
			> {output}
		"""
