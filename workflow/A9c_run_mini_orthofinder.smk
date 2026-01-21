print("Snakemake Python:", sys.executable)
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
vertebrates = rnaseq.query("IsVertebrate == 'Yes'")
invertebrates = rnaseq.query("IsVertebrate == 'No'")
print(vertebrates.shape)

rule all:
    input:
        #"scratch/mini_orthofinder/annotated_invertebrate/Orthogroups.txt",
        #"scratch/mini_orthofinder/annotated_vertebrate/Orthogroups.txt",
        expand("scratch/mini_orthofinder/denovo_vertebrate/{org}/Orthogroups.tsv", org = vertebrates.OrganismShortName),
        expand("scratch/mini_orthofinder/denovo_invertebrate/{org}/Orthogroups.tsv", org = invertebrates.OrganismShortName),

rule run_mini_orthofinder_annotated_invertebrate:
    input:
        peptides = "scratch/mini_orthofinder/annotated_invertebrate/peptides",
    output:
        "scratch/mini_orthofinder/annotated_invertebrate/done.txt",
    threads: 32
    shell:
        """
        module load diamond OrthoFinder
        orthofinder \
            -f {input.peptides} \
            -t 32 -a 16 \
            -I 1.5 \
            -o "scratch/mini_orthofinder/annotated_invertebrate/OrthoFinder_I1.5" -n "Fixed"
        touch {output}
        """

rule run_mini_orthofinder_annotated_vertebrate:
    input:
        peptides = "scratch/mini_orthofinder/annotated_vertebrate/peptides",
    output:
        "scratch/mini_orthofinder/annotated_vertebrate/done.txt",
    threads: 32
    shell:
        """
        module load diamond OrthoFinder
        orthofinder \
            -f {input.peptides} \
            -t 32 -a 16 \
            -I 1.7 \
            -o "scratch/mini_orthofinder/annotated_vertebrate/OrthoFinder" -n "Fixed_I1.7"
        touch {output}
        """

rule run_mini_orthofinder_denovo_vertebrate:
    input:
        peptides = "scratch/mini_orthofinder/denovo_vertebrate/{org}/peptides",
        prev = "scratch/mini_orthofinder/annotated_vertebrate/OrthoFinder/Results_Fixed_I1.7",
    output:
        "scratch/mini_orthofinder/denovo_vertebrate/{org}/Orthogroups.tsv"
    params:
        resdir="scratch/mini_orthofinder/denovo_vertebrate/{org}"
    threads: 32
    shell:
        """
        set -euo pipefail
        module load diamond OrthoFinder
        orthofinder \
            -b {input.prev} \
            -f {input.peptides} \
            -t 32 -a 16 \
            -I 5.0
        derived=$(ls -1td $(readlink -f {input.prev}/../Results_*) | grep -v 'Results_Fixed' | head -1)
        echo $derived
        cp $derived/Orthogroups/Orthogroups.t* {params.resdir}/
        rm -rf $derived
        """

rule run_mini_orthofinder_denovo_invertebrate:
    input:
        peptides = "scratch/mini_orthofinder/denovo_invertebrate/{org}/peptides",
        prev = "scratch/mini_orthofinder/annotated_invertebrate/OrthoFinder_I1.5/Results_Fixed",
    output:
        "scratch/mini_orthofinder/denovo_invertebrate/{org}/Orthogroups.tsv"
    params:
        resdir="scratch/mini_orthofinder/denovo_invertebrate/{org}"
    threads: 32
    shell:
        """
        set -euo pipefail
        module load diamond OrthoFinder
        orthofinder \
            -b {input.prev} \
            -f {input.peptides} \
            -t 32 -a 16 \
            -I 5.0
        derived=$(ls -1td $(readlink -f {input.prev}/../Results_*) | grep -v 'Results_Fixed' | head -1)
        echo $derived
        cp $derived/Orthogroups/Orthogroups.t* {params.resdir}/
        rm -rf $derived
        """

