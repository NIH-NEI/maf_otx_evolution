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

localrules: all, filter_vertebrate_peptides, subset_fna, align_by_clustal, align_by_pal2nal, extract_subtree

rule all:
    input:
        expand(
            TOP_DIR + "/{basename}/{basename}__{part}__{aligner}__{postalign}__{constr}__{treetype}.fna.ABSREL.json",
            basename = [f"jawed_{gene}" for gene in genes],
            part = ["full"],
            aligner = ["clustalo"],
            postalign = ["fullaln"],
            constr = ["species_topology"],
            treetype = ["given"],
        ),


rule extract_subtree:
    input:
        aln = TOP_DIR + "/{basename}__{part}__{aligner}__{postalign}__inframe.fna",
        tree = "imports/trees/species181_topology_shortnames.nwk",
    output:
        tree = TOP_DIR + "/{basename}__{part}__{aligner}__{postalign}__species_topology__given.newick",
    script:
        "scripts/draw_tree_alignment/create_protein_topology_from_species.R"

rule run_hyphy_absrel:
    input:
        aln = TOP_DIR + "/{basename}__{part}__{aligner}__{postalign}__inframe.fna",
        tree = TOP_DIR + "/{basename}__{part}__{aligner}__{postalign}__{constr}__{treetype}.newick",
    output:
        json = TOP_DIR + "/{basename}__{part}__{aligner}__{postalign}__{constr}__{treetype}.fna.ABSREL.json",
        aln = TOP_DIR + "/{basename}__{part}__{aligner}__{postalign}__{constr}__{treetype}.fna",
    threads: 32
    shell:
        """
        echo "module load hyphy; cp {input.aln} {output.aln}; hyphy CPU={threads} absrel --alignment {output.aln} --tree {input.tree}"
        """
