import pandas as pd

blast_res_file = snakemake.input["blast_res"]
output_file = snakemake.output[0]

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
    return df

df = get_sample_results(blast_res_file)
print(df)

df.to_csv(output_file, index=False, sep="\t")

