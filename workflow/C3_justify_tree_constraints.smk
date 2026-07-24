import pandas as pd

constraints = dict(
    LargeMAF = [
        "noconstr",
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
            TOP_DIR + "/comparison/{basename}__{part}__{aligner}__{postalign}__{treetype}__ztable.tsv",
            basename = basename,
            part = ["full"], #, "bzip"],
            aligner = ["clustalo"], #, "mafft"],
            postalign = ["fullaln"],
            treetype = ["fast"], #, "raxml"], #"iq"],
        ),

def get_input_files_for_iq_minus_z(wc):
    infiles = dict(
       fasta = TOP_DIR + "/{basename}__{part}__{aligner}__{postalign}.fasta",
    )
    for i, constraint in enumerate(constraint_list):
        print(i, constraint)
        key = constraint
        value = (
            TOP_DIR + "/{basename}__{part}__{aligner}__{postalign}__"
            + constraint
            + "__{treetype}.newick"
        )
        infiles[key] = value
    return infiles


rule justify_constrained_tree:
    input:
        unpack(get_input_files_for_iq_minus_z),
    output:
        zinfo = TOP_DIR + "/comparison/{basename}__{part}__{aligner}__{postalign}__{treetype}__zinfo.tsv",
        ztree = TOP_DIR + "/comparison/{basename}__{part}__{aligner}__{postalign}__{treetype}__allconstraints.newick",
        result = TOP_DIR + "/comparison/{basename}__{part}__{aligner}__{postalign}__{treetype}.iqtree",
    params:
        prefix = TOP_DIR + "/comparison/{basename}__{part}__{aligner}__{postalign}__{treetype}",
    run:
        import pandas as pd
        import csv
        def get_single_line_newick(file_path):
            with open(file_path) as f:
                lines = f.readlines()
            return ''.join([l.strip().replace("\s+", "") for l in lines])

        df = (
            pd.DataFrame(input.items(), columns = ["Constraint", "NewickFile"])
            .query("Constraint != 'fasta'")
            .assign(Tree = lambda tdf: range(1, len(tdf) + 1))
            .assign(TreeString = lambda tdf: tdf.NewickFile.apply(get_single_line_newick))
            [["Tree", "Constraint", "NewickFile", "TreeString"]]
        )
        print(df)
        df.to_csv(output.zinfo, sep = "\t", index=False)
        df["TreeString"].to_csv(
            output.ztree,
            index=False,
            header=False,
            sep="\t",
            quoting=csv.QUOTE_NONE,
            escapechar='\\',
        )
        cmd = """
        module load iqtree;
        iqtree2 -redo -nt 16 -s {input.fasta} -m 'WAG' -z {output.ztree} -n 0 -zb 10000 -au -pre {params.prefix};
        """
        print(cmd)
        shell(cmd)


rule reformat_comparison:
    input:
        result = TOP_DIR + "/comparison/{basename}__{part}__{aligner}__{postalign}__{treetype}.iqtree",
        zinfo = TOP_DIR + "/comparison/{basename}__{part}__{aligner}__{postalign}__{treetype}__zinfo.tsv",
    output:
        table = TOP_DIR + "/comparison/{basename}__{part}__{aligner}__{postalign}__{treetype}__ztable.tsv",
    script:
        "scripts/justify_constraint_trees/reformat_iqtree_zoption_results.py"
