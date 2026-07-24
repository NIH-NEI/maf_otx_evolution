import pandas as pd
from os.path import basename

annotated_sources = [
    "GCF",
    "Ensembl",
    "Custom"
]

organisms = (
    pd.read_csv("configs/species181_annotations.tsv", sep = "\t")
    .query("AnnotationSource in @annotated_sources")
    .set_index("OrganismShortName", drop = False)
)
print(organisms)

rule all:
    input:
        "scratch/grand_list/grand_list_nt.tsv",

def get_cds_file(org):
    annot_src = organisms["AnnotationSource"][org]
    annot_rel =  organisms["AnnotationRelease"][org]
    annot_name =  organisms["AnnotationName"][org]
    annot_org =  organisms["OrganismShortName"][org]
    annot_id =  organisms["OrganismID"][org]
    #print("Current org:", wc.org)
    #print("Annotation org:", annot_org)
    #print("Annotation org id:", annot_id)
    #print("Annotation source:", annot_src)
    #print("Annotation release:", annot_rel)
    #print("Annotation name:", annot_name)
    if annot_src == "GCF":
        return f"imports/genomes_annotations/{annot_org}/{annot_rel}_{annot_name}_cds_from_genomic.fna.gz"
    if annot_src == "Ensembl":
        return f"imports/genomes_annotations/{annot_org}/{annot_id}.{annot_name}.cds.all.fa.gz"
    if annot_src == "Custom":
        return {
            "haha": "hoho",
        }[annot_org]
    # It should not come here at all
    print("Error! Unknown annotation source")
    quit(-1)

def get_handle(fpath):
    if fpath[-3:] == ".gz":
        import gzip
        return gzip.open(fpath, "rt")
    return open(fpath, "r")

def get_cds_sequences(grp):
    from Bio import SeqIO
    import re
    organism = grp.name
    fasta_file = get_cds_file(organism)
    grp = grp.copy()
    grp["tgt"] = grp["ProteinID"].str.split("__").str[2]
    print(fasta_file)
    target_names = grp.tgt.tolist()
    print(target_names)

    annot_src = organisms["AnnotationSource"][organism]
    print(annot_src)
    if annot_src  == "Ensembl":
        p2t_file = f"scratch/canonical_peptides/{organism}.p2t"
        p2t = (
            pd.read_csv(p2t_file, header = None, sep = "\t", names = ["pid", "tid"])
            .assign(pid = lambda tdf: tdf["pid"].str.split(".").str[0])
            .set_index("pid")["tid"].to_dict()
        )
    else:
        p2t = {tgt:tgt for tgt in target_names}

    results = []
    for record in SeqIO.parse(get_handle(fasta_file), "fasta"):
        #print(record.id)
        for tgt in target_names:
            if p2t[tgt] in record.id:
                results.append([tgt, record.id, str(record.seq)])
    print(results)

    fasta_df = pd.DataFrame(results, columns = ["tgt", "recordid", "CDS"])
    print(fasta_df)

    fasta_df = fasta_df[["tgt", "CDS"]].drop_duplicates()

    # 3) Validate FASTA has exactly one CDS per tgt (for this organism)
    counts = fasta_df.groupby("tgt").size()
    dup = counts[counts > 1]

    problem = fasta_df.loc[fasta_df.tgt.isin(dup.index),:]
    for x in problem["CDS"].to_list():
        print(x)
    print(problem)
    if not dup.empty:
        print(f"!!!!!!WARNING!!!!! {organism}: multiple CDS records for targets: {dup.index.tolist()[:10]}")

    fasta_df = fasta_df.groupby("tgt").agg("first").reset_index()

    # 4) Merge: keep original ProteinID, attach CDS via tgt
    out = grp.merge(fasta_df, on="tgt", how="left")

    # 5) Validate all targets found
    missing = out.loc[out["CDS"].isna(), "tgt"].unique()
    if len(missing):
        print(f"!!!!!!WARNING!!!!! {organism}: missing CDS for targets: {missing[:10].tolist()}")

    return out.drop(columns=["tgt"])

rule get_nucleotides_for_peptides:
    input:
        tsv = "scratch/grand_list/grand_list.tsv",
    output:
        tsv = "scratch/grand_list/grand_list_nt.tsv",
        fna = "scratch/grand_list/grand_list.fna",
    run:
        from Bio import SeqIO
        import pandas as pd
        import gzip
        toskip = [
            "fareasternBrookLamprey",
            "grassMouse",
            "treeShrew",
        ]

        selected = [
        ]

        df = (
            pd.read_csv(input.tsv, sep = "\t")
            .query("ProteinSource != 'Denovo'")
            .query("OrganismShortName not in @toskip")
            #.query("OrganismShortName in @selected")
            [["ProteinID", "OrganismShortName"]]
            .groupby("OrganismShortName")
            .apply(get_cds_sequences)
            .reset_index(drop = True)
        )
        print(df)
        df.to_csv(output.tsv, index = False, sep = "\t")
        with open(output.fna, "w") as f:
            for i, row in df.iterrows():
                f.write(f">{row.ProteinID}\n{row.CDS}\n")

