import pandas as pd
from Bio import SeqIO

constraints = dict(
    LargeMAF = [
        "noconstr",
#        "constr_largemaf",
#        "constr_lamprey",
        #"constr_nrl_mafa_cmaf_mafb",
        #"constr_nrl_mafb_cmaf_mafa",
 #       "constr_nrl_cmaf_mafb_mafa",
 #       "constr_nrl_mafa_cmaf_mafb_tj",
 #       "constr_nrl_mafb_cmaf_mafa_tj",
 #       "constr_nrl_cmaf_mafb_mafa_tj",
 #       "constr_humanC",
 #       "constr_humanC_lamprey1",
 #       "constr_nrl_subtree",
#        "constr_humanC_lamprey2",
#        "constr_humanC_lamprey3",
#        "constr_humanC_lamprey4",
#        "constr_humanC_lamprey5",
#        "constr_humanC_lamprey6",
    ],
    SmallMAF = [
#        "noconstr",
#        "constr_smallmaf",
        #"constr_maff_mafg_mafk",
        #"constr_maff_mafg_mafk_mafs",
#        "constr_humanF",
#        "constr_humanF_mixed1",
#        "constr_humanF_mixed2",
    ],
    AllMAF = [
#        "noconstr",
        "constr_allmaf",
#
    ],
    AllOTX = [
#        "noconstr",
#        "constr_allotx",
 #       "constr_otx1_oc",
        "constr_otx_jawless",
 #       "constr_otx1_otx2_crx_oc",
        #"constr_otx1_otx2_crx_oc_hydra",
    ],
    Unknown = [
        "noconstr",
    ],
)

from pathlib import Path


infile = config["infile"]
print(infile)

basename =  Path(infile).stem
print(basename)

TOP_DIR = config.get("topdir", f"scratch/temp_tree_alignment/{basename}")
print(TOP_DIR)

gene_family = config.get("family", "Unknown")
print(gene_family)

constraint_list = constraints[gene_family]
print(constraint_list)


rule all:
    input:
        expand(
            TOP_DIR + "/{basename}__{part}__{aligner}__{postalign}__{constr}__{treetype}.pdf",
            basename = basename,
            part = ["full"],
            aligner = ["clustalo"], #["clustalo"], #, "mafft"],
            postalign = ["fullaln"], #["trimaln"], #["fullaln"],
            constr = constraint_list, #["noconstr"],
            #treetype = ["iq"], #["fast", "iq", "raxml"],
            treetype = ["fast"], #["fast", "iq", "raxml"],
        ),


rule create_fasta:
    input:
        tsv = infile,
    output:
        fasta = TOP_DIR + "/{basename}__full.faa",
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

rule align_by_mafft:
    input:
        fasta = TOP_DIR + "/{basename}__{part}.faa",
    output:
        fasta = TOP_DIR + "/{basename}__{part}__mafft__fullaln.fasta",
    threads: 16
    shell:
        """
        module load mafft
        mafft --maxiterate 1000 --localpair --thread {threads} {input.fasta} > {output.fasta}
        """
        #mafft --auto --thread {threads} {input.fasta} > {output.fasta}

rule trim_aligned_sequences:
    input:
        fasta = TOP_DIR + "/{basename}__{part}__{aligner}__fullaln.fasta",
    output:
        tmp = TOP_DIR + "/{basename}__{part}__{aligner}__trimtmp.fasta",
        fasta = TOP_DIR + "/{basename}__{part}__{aligner}__trimaln.fasta",
    shell:
        """
        module load trimal
        cat {input.fasta} | tr "*" "-" > {output.tmp}
        trimal -in {output.tmp} -out {output.fasta} -automated1
        """
        #trimal -in {output.tmp} -out {output.fasta} -gt 0.5 -cons 50

# costraints in fasttree
# http://www.microbesonline.org/fasttree/constrained.html#Specifying
# Use the perl file /data/pals2/eye_evolution/code/TreeToConstraints.pl
rule reformat_constraint:
    input:
        "configs/constraints/{constr}.newick",
    output:
        "configs/constraints/{constr}.fasta",
    shell:
        """
        perl workflow/scripts/draw_tree_alignment/TreeToConstraints.pl  < {input} > {output}
        """

rule create_unconstrained_fasttree:
    input:
        fasta = TOP_DIR + "/{basename}__{part}__{aligner}__{postalign}.fasta",
    output:
        tree = TOP_DIR + "/{basename}__{part}__{aligner}__{postalign}__noconstr__fast.newick",
    threads: 8
    shell:
        """
        module load FastTree
        export OMP_NUM_THREADS={threads}
        FastTreeMP < {input.fasta} > {output}
        """

rule create_constrained_fasttree:
    input:
        fasta = TOP_DIR + "/{basename}__{part}__{aligner}__{postalign}.fasta",
        constr = "configs/constraints/{constr}.fasta",
    output:
        tree = TOP_DIR + "/{basename}__{part}__{aligner}__{postalign}__{constr}__fast.newick",
    threads: 1
    shell:
        """
        module load FastTree
        export OMP_NUM_THREADS={threads}
        FastTreeMP -constraints {input.constr} < {input.fasta} > {output}
        """
        #FastTree -pseudo -constraints {input.constr} < {input.fasta} > {output}
        #FastTreeMP -constraints {input.constr} < {input.fasta} > {output}

#rule create_raxmltree:
#    input:
#        fasta = TOP_DIR + "/{basename}__{part}__{aligner}__{postalign}.fasta",
#        constr = "imports/maf_evolution/constraints/{constr}.newick",
#    output:
#        tree =
#        TOP_DIR + "/{basename}__{part}__{aligner}__{postalign}__{constr}__raxml.newick",
#    shell:
#        """
#        module load raxml-ng
#        raxml-ng --msa {input.fasta} --msa-format FASTA --data-type AA --model PROTGTR --tree-constraint {input.constr}
#        false
#        """
#
rule create_iqtree_unconstrained:
    input:
        fasta = TOP_DIR + "/{basename}__{part}__{aligner}__{postalign}.fasta",
    output:
        tree = TOP_DIR + "/{basename}__{part}__{aligner}__{postalign}__noconstr__iq.newick",
    threads: 16
    shell:
        """
        module load iqtree
        iqtree2 -s {input.fasta} -m 'WAG' -T {threads} -B 1000 --redo-tree
        cp {input.fasta}.treefile {output.tree}
        """
        #iqtree2 -s {input} -m 'WAG+R4' -T {threads} -B 1000 --redo-tree

rule create_iqtree_constrained:
    input:
        fasta = TOP_DIR + "/{basename}__{part}__{aligner}__{postalign}.fasta",
        constr = "configs/constraints/{constr}.newick",
        init = TOP_DIR + "/{basename}__{part}__{aligner}__{postalign}__{constr}__fast.newick",
    output:
        tree = TOP_DIR + "/{basename}__{part}__{aligner}__{postalign}__{constr}__iq.newick",
    threads: 32
    shell:
        """
        module load iqtree
        iqtree2 -s {input.fasta} -g {input.constr} -m LG+C20+F+G -T {threads} -B 1000 -alrt 1000 --redo-tree
        cp {input.fasta}.treefile {output.tree}
        """
        #iqtree2 -s {input.fasta} -g {input.constr} -m LG+C20+F+G -T {threads} -B 1000 -alrt 1000 --redo-tree -ft {input.init}
        #iqtree2 -s {input.fasta} -g {input.constr} -m 'WAG' -T {threads} -B 1000 --redo-tree -nt {threads}
        #iqtree2 -s {input} -m 'WAG+R4' -T {threads} -B 1000 --redo-tree

#rule root_phylogeny_ape:
#    input:
#        tree = TOP_DIR + "/{basename}__{part}__{aligner}__{postalign}__{constr}__{treetype}.newick",
#        outgroup = "configs/outgroup_{basename}.tsv",
#    output:
#        tree = TOP_DIR + "/{basename}__{part}__{aligner}__{postalign}__{constr}__{treetype}_rooted.newick",
#    script:
#        "scripts/create_tree_alignment/root_tree.R"


rule plot_phylogeny:
    input:
        aln = TOP_DIR + "/{basename}__full__{aligner}__{postalign}.fasta",
        tree = TOP_DIR + "/{basename}__{part}__{aligner}__{postalign}__{constr}__{treetype}.newick",
        info = infile,
    output:
        plot = TOP_DIR + "/{basename}__{part}__{aligner}__{postalign}__{constr}__{treetype}.pdf",
#        tree = TOP_DIR + "/{basename}__{part}__{aligner}__{postalign}__{constr}__{treetype}.itol.newick",
#        itol = TOP_DIR + "/{basename}__{part}__{aligner}__{postalign}__{constr}__{treetype}.itol.csv",
#        fasta = TOP_DIR + "/{basename}__{part}__{aligner}__{postalign}__{constr}__{treetype}.fasta",
#        tip = TOP_DIR + "/{basename}__{part}__{aligner}__{postalign}__{constr}__{treetype}.tip.txt",
    script:
        "scripts/draw_tree_alignment/draw_tree.R"

