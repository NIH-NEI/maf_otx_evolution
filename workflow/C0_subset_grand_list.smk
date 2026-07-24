import pandas as pd

representative_species = (
    pd.read_csv("configs/species181_table.tsv", sep = "\t")
    .query("Representative == 'Yes'")
    ["OrganismShortName"]
    .tolist()
)
print(representative_species)

#skip_species = [
#    "treeShrew",
#    "grassMouse",
#]

all_valid_groups = dict(
    CMAF = ["CMAF"],
    NRL = ["NRL"],
    MAFA = ["MAFA"],
    MAFB = ["MAFB"],
    LargeMAF = ["CMAF", "NRL", "MAFA", "MAFB", "MAFL","MAFL1A","MAFL1B","MAFL2A","MAFL2B"],
    MAFF = ["MAFF"],
    MAFG = ["MAFG"],
    MAFK = ["MAFK"],
    SmallMAF = ["MAFF", "MAFG", "MAFK", "MAFS","MAFS1A"],
    AllMAF = [
        "CMAF", "NRL", "MAFA", "MAFB", "MAFL","MAFL1A","MAFL1B","MAFL2A","MAFL2B",
        "MAFF", "MAFG", "MAFK", "MAFS", "MAFS1A"
    ],
    OTX1 = ["OTX1"],
    OTX2 = ["OTX2"],
    CRX = ["CRX"],
    AllOTX =["OTX1", "OTX2", "CRX", "OTX","OTX1C","OTX2A","OTX2B","OTX2C"],
)

rule all:
    input:
        expand("scratch/final_list/full_{genegrp}.tsv", genegrp = [
            "CMAF", "NRL", "MAFA", "MAFB",
            "OTX1", "OTX2", "CRX",
            "LargeMAF",
            "SmallMAF",
            "AllMAF",
            "AllOTX",
        ]),
        expand("scratch/final_list/annotated_{genegrp}.tsv", genegrp = [
            "CMAF", "NRL", "MAFA", "MAFB",
            "OTX1", "OTX2", "CRX",
            "LargeMAF",
            "SmallMAF",
            "AllMAF",
            "AllOTX",
        ]),
        expand("scratch/final_list/representative_{genegrp}.tsv", genegrp = [
            "LargeMAF",
            "SmallMAF",
            "AllMAF",
            "AllOTX",
        ]),
        "scratch/final_list/final_exported_list.tsv",

rule subset_full:
    input:
        grand = "scratch/grand_list/grand_list.tsv",
        outgrp = "configs/outgroups.tsv",
    output:
        "scratch/final_list/full_{selection}.tsv",
    run:
        import pandas as pd
        selection = wildcards.selection
        valid_groups = all_valid_groups[selection]
        print(valid_groups)

        outgroups = (
            pd.read_csv(input.outgrp, sep = "\t")
            .query("GeneGroup == @selection")
            .Outgroup.tolist()
        )
        print(outgroups)

        df = pd.read_csv(input.grand, sep = "\t")
        print(df.query("ProteinID == 'drosophila__oc__NP_001356934.1'"))

        df_main = (
            df
            .query("GeneGroup in @valid_groups")
            #.query("ProteinSource != 'Denovo'")
            .query("Status != 'Undecided'")
        )
        print(df_main.query("ProteinID == 'drosophila__oc__NP_001356934.1'"))

        df_outgrp = (
            df.
            query("ProteinID in @outgroups")
            .assign(GeneGroup = "Outgroup")
        )

        df = pd.concat([df_main, df_outgrp])
        print(df)
        print(df.query("ProteinID == 'drosophila__oc__NP_001356934.1'"))

        df.to_csv(output[0], index = False, sep = "\t")


rule subset_annotated:
    input:
        full = "scratch/final_list/full_{selection}.tsv",
    output:
        "scratch/final_list/annotated_{selection}.tsv",
    run:
        import pandas as pd

        df = (
            pd.read_csv(input.full, sep = "\t")
            .query("(ProteinSource != 'Denovo') or (GeneGroup == 'Outgroup')")
        )

        df.to_csv(output[0], index = False, sep = "\t")



rule subset_representative:
    input:
        annot = "scratch/final_list/annotated_{selection}.tsv",
    output:
        "scratch/final_list/representative_{selection}.tsv",
    run:
        import pandas as pd

        df = (
            pd.read_csv(input.annot, sep = "\t")
            .query("(OrganismShortName in @representative_species) or (GeneGroup == 'Outgroup')")
        )


        df.to_csv(output[0], index = False, sep = "\t")

rule export_table:
    input:
        grand = "scratch/grand_list/grand_list.tsv",
    output:
        "scratch/final_list/final_exported_list.tsv",
    run:
        import pandas as pd

        df = (
            pd.read_csv(input.grand, sep = "\t")
            .query("GeneGroup not in ['Outgroup', 'Ignore']")
            .query("Status != 'Undecided'")
            .query("LengthCriteria not in ['Truncated', 'VeryLong']")
            .assign(ProteinLen = lambda tdf: tdf.ProteinSeq.apply(lambda x: len(x)))
        )
        print(df)
        print(df[df.ProteinSeq.isna()])
        print(df[df.ProteinLen.isna()])

        df.to_csv(output[0], index = False, sep = "\t")


