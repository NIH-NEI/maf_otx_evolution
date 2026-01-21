import pandas as pd
from Bio import SeqIO

constraints = dict(
    LargeMAF = [
        "noconstr",
        #"constr_nrl_mafa_cmaf_mafb",
        #"constr_nrl_mafb_cmaf_mafa",
        "constr_nrl_cmaf_mafb_mafa",
        "constr_nrl_mafa_cmaf_mafb_tj",
        "constr_nrl_mafb_cmaf_mafa_tj",
        "constr_nrl_cmaf_mafb_mafa_tj",
    ],
    SmallMAF = [
        "noconstr",
        "constr_maff_mafg_mafk",
        "constr_maff_mafg_mafk_mafs",
    ],
    AllOTX = [
        "noconstr",
        "constr_otx1_otx2_crx_oc",
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


rule all:
    input:
        expand(
            TOP_DIR + "/{basename}__{part}__{aligner}__{postalign}__{constr}__{treetype}.pdf",
            basename = basename,
            part = ["full"],
            aligner = ["clustalo"], #, "mafft"],
            postalign = ["fullaln"],
            constr = constraint_list, #["noconstr"],
            treetype = ["fast"], #, "raxml"], #"iq"],
        ),

        #        "scratch/protein_tree_aln/species50/ALLMAFS__full.faa",
#        expand(
#            TOP_DIR + "/{basename}__{part}__{aligner}__{postalign}__{constr}__{treetype}.{ext}",
#            topdir = "scratch/protein_tree_aln",
#            collection = "species50",
#            basename = ["MAFL"],
#            part = ["full"], #, "bzip"],
#            aligner = ["clustalo"], #, "mafft"],
#            postalign = ["fullaln"],
#            constr = constraints["MAFL"],
#            treetype = ["fast"], #, "raxml"], #"iq"],
#            ext = ["newick", "pdf"],
#        ),
#        expand(
#            TOP_DIR + "/comparison/{basename}__{part}__{aligner}__{postalign}__{treetype}__ztable.tsv",
#            topdir = "scratch/protein_tree_aln",
#            collection = "species50",
#            basename = ["MAFL"], #, "MAFS"],
#            part = ["full"], #, "bzip"],
#            aligner = ["clustalo"], #, "mafft"],
#            postalign = ["fullaln"],
#            treetype = ["fast"], #, "raxml"], #"iq"],
#        ),
#        expand(
#            TOP_DIR + "/{basename}__{part}__{aligner}__{postalign}__{constr}__{treetype}.pdf",
#            topdir = "scratch/protein_tree_aln",
#            collection = "species50",
#            basename = ["MAFL"],
#            part = ["full"], #, "bzip"],
#            aligner = ["clustalo"], #, "mafft"],
#            postalign = ["fullaln"], #, "trimaln"],
#            constr = constraints["MAFL"],
#            treetype = ["fast"], #, "raxml"], #"iq"],
#        ),
#        expand(
#            TOP_DIR + "/comparison/{basename}__{part}__{aligner}__{postalign}__{treetype}__ztable.tsv",
#            topdir = "scratch/protein_tree_aln",
#            collection = "species50",
#            basename = ["MAFS"], #, "MAFS"],
#            part = ["full"], #, "bzip"],
#            aligner = ["clustalo"], #, "mafft"],
#            postalign = ["fullaln"],
#            treetype = ["fast"], #, "raxml"], #"iq"],
#        ),
#        expand(
#            TOP_DIR + "/{basename}__{part}__{aligner}__{postalign}__{constr}__{treetype}.{ext}",
#            topdir = "scratch/protein_tree_aln",
#            collection = "species50",
#            basename = ["OTX"],
#            part = ["full"], #, "bzip"],
#            aligner = ["clustalo"], #, "mafft"],
#            postalign = ["fullaln"],
#            constr = constraints["OTX"],
#            treetype = ["fast"], #, "raxml"], #"iq"],
#            ext = ["newick", "pdf"],
#        ),
#        expand(
#            TOP_DIR + "/comparison/{basename}__{part}__{aligner}__{postalign}__{treetype}__ztable.tsv",
#            topdir = "scratch/protein_tree_aln",
#            collection = "species50",
#            basename = ["OTX"], #, "MAFS"],
#            part = ["full"], #, "bzip"],
#            aligner = ["clustalo"], #, "mafft"],
#            postalign = ["fullaln"],
#            treetype = ["fast"], #, "raxml"], #"iq"],
#        ),
#        expand(
#            TOP_DIR + "/{basename}__{part}__{aligner}__{postalign}__{constr}__{treetype}.pdf",
#            topdir = "scratch/protein_tree_aln",
#            collection = "species50",
#            basename = ["GNB"],
#            part = ["full"], #, "bzip"],
#            aligner = ["clustalo"], #, "mafft"],
#            postalign = ["fullaln"], #, "trimaln"],
#            constr = [
#                "constr_none",
#                #"constr_maff_mafg_mafk",
#                #"constr_maff_mafg_mafk_mafs",
#            ],
#            treetype = ["fast"], #, "raxml"], #"iq"],
#        ),

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

#rule align_by_mafft:
#    input:
#        fasta = TOP_DIR + "/{basename}__{part}.faa",
#    output:
#        fasta = TOP_DIR + "/{basename}__{part}__mafft__fullaln.fasta",
#    threads: 16
#    shell:
#        """
#        module load mafft
#        mafft --maxiterate 1000 --localpair --thread {threads} {input.fasta} > {output.fasta}
#        """
#        #mafft --auto --thread {threads} {input.fasta} > {output.fasta}
#
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
        trimal -in {output.tmp} -out {output.fasta} -gt 0.9 -cons 1
        """
        #trimal -in {output.tmp} -out {output.fasta} -automated1

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
        FastTree -pseudo -constraints {input.constr} < {input.fasta} > {output}
        """
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
#rule create_iqtree:
#    input:
#        fasta = TOP_DIR + "/{basename}__{part}__{aligner}__{postalign}.fasta",
#        constr = "imports/maf_evolution/constraints/{constr}.newick",
#    output:
#        tree = TOP_DIR + "/{basename}__{part}__{aligner}__{postalign}__{constr}__iq.newick",
#    threads: 16
#    shell:
#        """
#        module load iqtree
#        iqtree2 -s {input.fasta} -g {input.constr} -m 'WAG' -T {threads} -B 1000 --redo-tree
#        cp {input.fasta}.treefile {output.tree}
#        """
#        #iqtree2 -s {input} -m 'WAG+R4' -T {threads} -B 1000 --redo-tree
#
#def get_input_files_for_iq_minus_z(wc):
#    infiles = dict(
#       fasta = TOP_DIR + "/{basename}__{part}__{aligner}__{postalign}.fasta",
#    )
#    for i, constraint in enumerate(constraints[wc.basename]):
#        print(i, constraint)
#        key = constraint
#        value = (
#            TOP_DIR + "/{basename}__{part}__{aligner}__{postalign}__"
#            + constraint
#            + "__{treetype}.newick"
#        )
#        infiles[key] = value
#    return infiles
#
#rule justify_constrained_tree:
#    input:
#        unpack(get_input_files_for_iq_minus_z),
#    output:
#        zinfo = TOP_DIR + "/comparison/{basename}__{part}__{aligner}__{postalign}__{treetype}__zinfo.tsv",
#        ztree = TOP_DIR + "/comparison/{basename}__{part}__{aligner}__{postalign}__{treetype}__zinfo.trees",
#        result = TOP_DIR + "/comparison/{basename}__{part}__{aligner}__{postalign}__{treetype}.iqtree",
#    params:
#        prefix = TOP_DIR + "/comparison/{basename}__{part}__{aligner}__{postalign}__{treetype}",
#    run:
#        import pandas as pd
#        import csv
#        def get_single_line_newick(file_path):
#            with open(file_path) as f:
#                lines = f.readlines()
#            return ''.join([l.strip().replace("\s+", "") for l in lines])
#
#        df = (
#            pd.DataFrame(input.items(), columns = ["Constraint", "NewickFile"])
#            .query("Constraint != 'fasta'")
#            .assign(Tree = lambda tdf: range(1, len(tdf) + 1))
#            .assign(TreeString = lambda tdf: tdf.NewickFile.apply(get_single_line_newick))
#            [["Tree", "Constraint", "NewickFile", "TreeString"]]
#        )
#        print(df)
#        df.to_csv(output.zinfo, sep = "\t", index=False)
#        df["TreeString"].to_csv(
#            output.ztree,
#            index=False,
#            header=False,
#            sep="\t",
#            quoting=csv.QUOTE_NONE,
#            escapechar='\\',
#        )
#        cmd = """
#        module load iqtree;
#        iqtree2 -redo -nt 16 -s {input.fasta} -m 'WAG' -z {output.ztree} -n 0 -zb 10000 -au -pre {params.prefix};
#        """
#        print(cmd)
#        shell(cmd)
#
#rule reformat_comparison:
#    input:
#        result = TOP_DIR + "/comparison/{basename}__{part}__{aligner}__{postalign}__{treetype}.iqtree",
#        zinfo = TOP_DIR + "/comparison/{basename}__{part}__{aligner}__{postalign}__{treetype}__zinfo.tsv",
#    output:
#        table = TOP_DIR + "/comparison/{basename}__{part}__{aligner}__{postalign}__{treetype}__ztable.tsv",
#    script:
#        "scripts/create_tree_alignment/reformat_iqtree_zoption_results.py"
#
#rule root_phylogeny_ete3:
#    input:
#        tree = TOP_DIR + "/{basename}__{part}__{aligner}__{postalign}__{constr}__{treetype}.newick",
#        outgroup = "configs/outgroup_{basename}.tsv",
#    output:
#        tree =
#        TOP_DIR + "/{basename}__{part}__{aligner}__{postalign}__{constr}__{treetype}_ete3rooted.newick",
#    run:
#        import pandas as pd
#        outgroups = (
#            pd.read_csv(input.outgroup, sep = "\t")
#            .assign(ProteinName = lambda tdf: tdf.ShortName + "__" + tdf.GeneSymbol)
#            ["ProteinName"].tolist()
#        )
#        print(outgroups)
#        from ete3 import Tree
#        t = Tree(input.tree, format=1)  # format=1 for internal node support
#        print("Is the tree rooted?", t.get_tree_root().is_root())
#        ancestor = t.get_common_ancestor(outgroups)
#        t.set_outgroup(ancestor)
#        t.write(outfile = output.tree)
#
#
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

