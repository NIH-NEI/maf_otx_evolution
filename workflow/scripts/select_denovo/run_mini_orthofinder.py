import pandas as pd

#protein_info_file = snakemake.input["protein_info"]
blast_res_file = snakemake.input["blast_res"]
output_file = snakemake.output[0]

#protein_info = (
#    pd.read_table(protein_info_file)
#    [["symbol", "product_accession"]]
#)

col_names = "qseqid qgi qacc qaccver qlen sseqid sallseqid sgi sallgi sacc saccver sallacc slen qstart qend sstart send qseq sseq evalue bitscore score length pident nident mismatch positive gapopen gaps ppos frames qframe sframe btop staxid ssciname scomname sblastname sskingdom staxids sscinames scomnames sblastnames sskingdoms stitle salltitles sstrand qcovs qcovhsp qcovus".split()

def get_sample_results(path):
    df = (
        pd.read_table(path, header = None, comment="#", names = col_names)
        .assign(subject = lambda xdf: xdf.saccver.str.split("__").str[0])
        .assign(contig = lambda xdf: xdf.saccver.str.split("__").str[1])
        .assign(qalignlen = lambda xdf: xdf.qseq.apply(len))
        .assign(coverage = lambda xdf:  (xdf.qalignlen / xdf.qlen))
        .assign(normident = lambda xdf: xdf.pident * xdf.coverage)
        [["qaccver", "saccver", "subject", "contig", "normident", "pident", "evalue", "bitscore", "coverage", "slen"]]
        .query("evalue < 1e-6")
        .query("pident > 30.0")
#        .assign(src_gene = lambda xdf: xdf.saccver.str.split("_").str[:-1].str.join("_"))
    )
    print(df)
#    print(df.src_gene.unique())
#
#    best = (
#        df.groupby(["src_gene"])
#        .agg({
#            "slen": ["max", "idxmax"],
#        })
#    )
#    best.columns = ["_".join(x) for x in best.columns]
#
#    keep = df.iloc[best["slen_idxmax"],:]
#    print(keep)
#    quit()
#
#    best["normident_best_contig"] = df.iloc[best["normident_idxmax"],:]["contig"].values
##    best["pident_best_contig"] = df.iloc[best["pident_idxmax"],:]["contig"].values
#    best["evalue_best_contig"] = df.iloc[best["evalue_idxmin"],:]["contig"].values
#    best["bitscore_best_contig"] = df.iloc[best["bitscore_idxmax"],:]["contig"].values
#
#    best = (
#        best.reset_index()
#        .rename(columns = {
#            "qaccver": "product_accession",
#            "normident_max": "normident",
#            "pident_max": "pident",
#            "evalue_min": "evalue",
#            "bitscore_max": "bitscore",
#        })
#        .drop(columns = ["normident_idxmax", "pident_idxmax", "evalue_idxmin", "bitscore_idxmax"])
#    )

    #return protein_info.merge(best, how = "inner")
    return df

df = get_sample_results(blast_res_file)
print(df)

df.to_csv(output_file, index=False, sep="\t")

