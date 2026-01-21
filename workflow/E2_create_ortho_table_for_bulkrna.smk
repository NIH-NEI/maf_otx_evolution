import pandas as pd
from os.path import basename

print("Snakemake Python:", sys.executable)

use_annotations = [
    "GCF",
    "Ensembl",
    "Custom"
]

annotated = (
    pd.read_csv("configs/species_table.tsv", sep = "\t")
    .query("SurrogateAnnotationSource in @use_annotations")
    .query("OrganismOrder177 > 27")
)
print(annotated)

bulkrna = (
    annotated
    .query("bulkrna107 == 'Yes'")
)
print(bulkrna)

scrna = (
    annotated
    .query("scrna35 == 'Yes'")
)
print(scrna)

rnaseq = (
    annotated
    .query("(scrna35 == 'Yes') or (bulkrna107 == 'Yes')")
)
print(rnaseq)



OF_groups = {
    "surrogate_bulkrna": bulkrna,
    "surrogate_scrna": scrna,
    "surrogate_rnaseq": rnaseq,
}

group_org_pairs = []
for grp, df in OF_groups.items():
    for org in df.SurrogateShortName.unique():
        group_org_pairs.append((grp, org))
group_org_pairs = pd.DataFrame(group_org_pairs, columns = ["grp", "org"])
print(group_org_pairs)



rule all:
    input:
        #expand("scratch/orthofinder_{group}/peptides/{org}.faa", zip, group=group_org_pairs.grp, org=group_org_pairs.org),
        #expand("scratch/orthofinder_{group}/topology_{group}.nwk", group = OF_groups.keys()),
        #expand("scratch/orthofinder_{group}/run_orthofinder_{group}.sh", group = OF_groups.keys()),
        atmost1gene = "scratch/orthofinder_surrogate_bulkrna/ortholog_tables/orthotab_atmost1gene.txt",
        maf_focused = "scratch/orthofinder_surrogate_bulkrna/ortholog_tables/orthotab_maf_focused.tsv",
        glist = "scratch/orthofinder_surrogate_bulkrna/ortholog_tables/patch_grand.tsv",

rule prepare_ortho_table:
    input:
        ogs = "scratch/orthofinder_surrogate_{spgrp}/results/Results_Dec31/Orthogroups/Orthogroups.tsv",
    output:
        atmost1gene = "scratch/orthofinder_surrogate_{spgrp}/ortholog_tables/orthotab_atmost1gene.txt",
    run:
        import pandas as pd
        from tqdm import tqdm
        tqdm.pandas()

        orthogroups = pd.read_csv(input.ogs, sep = "\t", dtype = str, low_memory = False)
        df = (
            orthogroups
            .melt(id_vars = "Orthogroup", var_name = "OrganismShortName", value_name = "ProteinID")
            .query("ProteinID.notna()")
            .assign(ProteinID = lambda tdf: tdf.ProteinID.apply(lambda x: x.split(", ")))
            .explode("ProteinID")
            .assign(GeneSymbol = lambda tdf: tdf.ProteinID.str.split("__").str[1])
            [["Orthogroup", "OrganismShortName", "GeneSymbol", "ProteinID"]]
            .drop_duplicates()
        )
        assert(df.loc[df.GeneSymbol.isna(),:].empty)
        print(df)

        # There could be mutiple proteins for the same gene but assigned to different orthogroups
        # Simply ignore them for ortholog table
        weird = (
            df
            .groupby(["OrganismShortName", "GeneSymbol"])
            .agg(
                ProteinIDs=("ProteinID", ",".join),
                OrthogroupCount=("Orthogroup", "nunique"),
                Orthogroups=("Orthogroup", ",".join),
            )
            .query("OrthogroupCount > 1")
        )
        print(weird)
        #weird.to_csv("weird.tsv", sep = "\t")
        weird_proteins = ",".join(weird.ProteinIDs).split(",")
        print(weird_proteins)

        df = df.query("ProteinID not in @weird_proteins")

        # Now treat multiple proteins with same gene symbol that were assigned the
        # same orthogroup as different isoforms
        df = (
            df
            [["Orthogroup", "OrganismShortName", "GeneSymbol"]]
            .drop_duplicates()
        )
        print(df)

        # african frog has proteins with suffix .L and .S treat them as same gene
        df["GeneSymbol"] = df["GeneSymbol"].str.replace(
            r"\.(L|S)$", "", regex=True
        )

        organisms = (
            pd.read_csv("configs/species_table.tsv", sep = "\t")
            [["OrganismShortName", "OrganismColor"]]
        )
        print(organisms)

        df = (
            df.merge(organisms, how = "left")
            .assign(OrthoOrg = lambda tdf: tdf.Orthogroup + "__" + tdf.OrganismShortName)
        )
        print(df)

        counts = (
            df
            .groupby(["OrthoOrg", "Orthogroup", "OrganismShortName", "OrganismColor"])
            .agg({"GeneSymbol": "count"})
            .reset_index()
        )
        print(counts)

        unique_mask = counts["GeneSymbol"] <= 1
        unique_counts = counts[unique_mask]
        print(unique_counts)

        multi_orthogroups = df.query("OrthoOrg not in @unique_counts.OrthoOrg")
        print(multi_orthogroups)

        duplicates_allowed = [
            "Nonchordates",
            "Protochordates",
            "Agnathans (hagfish and lampreys)",
            "Cartilaginous fishes",
            "Non-teleost ray-finned fishes",
            "Teleost ray-finned fishes",
            "Lobe-finned fishes",
            "Amphibians",
        ]

        # get proper count of genes per orthogroup and organism
        def count_distinct_genes(grp):
            symbols = grp.GeneSymbol.tolist()
            org = grp.OrganismShortName.tolist()[0]
            print(org)

            if org == "human":
                return len(symbols)

            color = grp.OrganismColor.tolist()[0]
            ortho_org = grp.OrthoOrg.tolist()[0]
            #print(symbols)

            #symbols = [x for x in symbols if not x.startswith("gene_ENS")]
            #print(symbols)

            locs = [x for x in symbols if x.startswith("LOC")]
            #print(locs)

            rest = [x for x in symbols if not x.startswith("LOC")]
            #print(rest)

            down_factor = 5 if color in duplicates_allowed else 1

            count = 1 if len(rest) <= 1 else len(symbols)/down_factor
            #print(ortho_org, locs, rest, count)

            return count


        multi_counts = (
            multi_orthogroups
            .groupby(["OrthoOrg", "Orthogroup", "OrganismShortName", "OrganismColor"])
            .progress_apply(count_distinct_genes)
        )
        multi_counts.name = "GeneSymbol"
        multi_counts = multi_counts.reset_index()
        print(multi_counts)

        counts = (
            pd.concat([unique_counts, multi_counts])
        )
        print(counts)

        counts = counts.pivot_table(
            index="Orthogroup",
            columns="OrganismShortName",
            values="GeneSymbol",
            aggfunc="max",   # or sum, max, min, etc.
            fill_value=0,
        )
        print(counts)

        atmost1gene = counts[(counts <= 1).all(axis=1)]
        print(atmost1gene)

        atmost1gene_human = atmost1gene[atmost1gene.human == 1]
        print(atmost1gene_human)

        df = (
            atmost1gene_human
            .reset_index()
        )
        df[["Orthogroup"]].to_csv(output.atmost1gene, index = None, header = None)

rule create_maf_focussed_ortho_table:
    input:
        ogs = "scratch/orthofinder_surrogate_{spgrp}/results/Results_Dec31/Orthogroups/Orthogroups.tsv",
        atmost1gene = "scratch/orthofinder_surrogate_{spgrp}/ortholog_tables/orthotab_atmost1gene.txt",
        maf_otx_ogs = "scratch/orthofinder_surrogate_{spgrp}/ortholog_tables/maf_otx_ogs.tsv",
        grand_list = "scratch/grand_list/grand_list.tsv",
    output:
        maf_otx_focused = "scratch/orthofinder_surrogate_{spgrp}/ortholog_tables/orthotab_maf_focused.tsv",
    run:
        selected = (
            pd.read_csv(input.atmost1gene, sep = "\t", header = None, names = ["Orthogroup"])
        )
        print(selected)

        # this list comes from manual inspection, cannot be automated for the time
        maf_otx_ogs = (
            pd.read_csv(input.maf_otx_ogs, sep = "\t")
        )
        print(maf_otx_ogs)

        ogs = (
            pd.read_csv(input.ogs, sep = "\t", dtype = str, low_memory = False)
            .query("Orthogroup in @selected.Orthogroup")
            .query("Orthogroup not in @maf_otx_ogs.Orthogroup")
            .set_index("Orthogroup")
        )
        print(ogs)

        organisms = ogs.columns
        print(organisms)

        others = (
            ogs.reset_index()
            .melt(id_vars = "Orthogroup", var_name = "OrganismShortName", value_name = "ProteinID")
            .query("ProteinID.notna()")
            .assign(ProteinID = lambda tdf: tdf.ProteinID.apply(lambda x: x.split(", ")))
            .explode("ProteinID")
            .assign(ProteinID = lambda tdf: tdf.ProteinID.apply(lambda x: x.replace("__gene_ENS", "__ENS")))
            .assign(GeneSymbol = lambda tdf: tdf.ProteinID.str.split("__").str[1])
            [["Orthogroup", "OrganismShortName", "GeneSymbol"]]
            .drop_duplicates()
            .pivot_table(
                index="Orthogroup",
                columns="OrganismShortName",
                values="GeneSymbol",
                aggfunc=",".join,
                fill_value="",
            )
            .assign(OrthoSymbol = lambda tdf: tdf.human)
            .reset_index(drop = True)
            .set_index("OrthoSymbol")
        )
        print(others)

        maf_otx = (
            pd.read_csv(input.grand_list, sep = "\t")
            .query("GeneFamily in ['AllOTX', 'LargeMAF', 'SmallMAF']")
            .query("GeneGroup in ['CMAF', 'NRL', 'MAFA', 'MAFB', 'MAFF', 'MAFK', 'MAFG', 'OTX1', 'OTX2', 'CRX']")
            .query("OrganismShortName in @organisms")
            .query("Status == 'Final'")
            .query("ProteinSource == 'Annotated'")
            .assign(ProteinID = lambda tdf: tdf.ProteinID.apply(lambda x: x.replace("__gene_ENS", "__ENS")))
            .assign(GeneSymbol = lambda tdf: tdf.ProteinID.str.split("__").str[1])
            [["OrganismShortName", "GeneGroup", "GeneSymbol"]]
            .rename(columns = {"GeneGroup": "OrthoSymbol"})
            .drop_duplicates()
            .pivot_table(
                index="OrthoSymbol",
                columns="OrganismShortName",
                values="GeneSymbol",
                aggfunc=",".join,
                fill_value="",
            )
        )
        print(maf_otx)

        orthotab = (
            pd.concat([maf_otx, others])
        )
        print(orthotab)
        orthotab.to_csv(output.maf_otx_focused, sep = "\t")

rule grand_list_for_three_extra_species:
    input:
        ogs = "scratch/orthofinder_surrogate_{spgrp}/results/Results_Dec31/Orthogroups/Orthogroups.tsv",
        maf_otx_ogs = "scratch/orthofinder_surrogate_{spgrp}/ortholog_tables/maf_otx_ogs.tsv",
    output:
        glist = "scratch/orthofinder_surrogate_{spgrp}/ortholog_tables/patch_grand.tsv",
    run:
        import pandas as pd
        from Bio import SeqIO
        from Bio.Seq import Seq
        from Bio.SeqRecord import SeqRecord

        def add_protein_seqs(group: pd.DataFrame) -> pd.DataFrame:
            # group.name is the value of the "pep" column for this group
            pep_file = group.name

            print(f"Processing {pep_file} ...")

            required = group.ProteinID.tolist()

            # Read all sequences from this fasta once
            seq_dict = {}
            with open(pep_file, "r") as spe:
                for record in SeqIO.parse(spe, "fasta"):
                    record.id = record.id.replace("__gene:ENS", "__ENS")
                    if record.id in required:
                        seq_dict[record.id] = str(record.seq)

            if len(required) != len(seq_dict):
                print("Error: Not all sequences are found!")
                print("Needed:", required)
                print("Found:", seq_dict.keys())
                exit(-1)

            # Map ProteinID to sequences
            group["ProteinSeq"] = group["ProteinID"].map(seq_dict.get)

            return group

        maf_otx_ogs = (
            pd.read_csv(input.maf_otx_ogs, sep = "\t")
        )
        print(maf_otx_ogs)

        ogs = (
            pd.read_csv(input.ogs, sep = "\t", dtype = str, low_memory = False)
            .query("Orthogroup in @maf_otx_ogs.Orthogroup")
            .set_index("Orthogroup")
            [["americanEel", "sunangel", "orientalScopsOwl"]]
        )
        print(ogs)

        organisms = ogs.columns
        print(organisms)

        others = (
            ogs.reset_index()
            .melt(id_vars = "Orthogroup", var_name = "OrganismShortName", value_name = "ProteinID")
            .query("ProteinID.notna()")
            .assign(ProteinID = lambda tdf: tdf.ProteinID.apply(lambda x: x.split(", ")))
            .explode("ProteinID")
            .assign(ProteinID = lambda tdf: tdf.ProteinID.apply(lambda x: x.replace("__gene_ENS", "__ENS")))
            .assign(pep = lambda tdf: "scratch/canonical_peptides/" + tdf.OrganismShortName + ".faa")
            .groupby("pep", group_keys=False)
            .apply(add_protein_seqs)
            .assign(ProteinLen = lambda tdf: tdf.ProteinSeq.apply(lambda x: len(x)))
            .assign(ProteinSource = "Annotated")
            [['ProteinID', 'ProteinLen', 'ProteinSource', 'ProteinSeq']]

        )
        print(others)
        others.to_csv(output.glist, sep = "\t", index = False)

