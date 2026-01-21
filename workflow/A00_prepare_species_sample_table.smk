
rule all:
    input:
        "configs/species_table.tsv"

rule prepare_species_table:
    input:
        names = "configs/species177_names_orders.tsv",
        datasets = "configs/species177_datasets.tsv",
        annotations = "configs/species177_annotations.tsv",
    output:
        table = "configs/species_table.tsv"
    run:
        import pandas as pd
        import numpy as np
        annot = pd.read_csv(input.annotations, sep = "\t")
        print(annot)
        annot["SurrogateID"] = annot["SurrogateID"].fillna(annot["OrganismID"])
        print(annot)
        annot2 = (
            annot[[
                "OrganismID",
                "OrganismShortName",
                "AnnotationSource",
                "AnnotationRelease",
                "AnnotationName",
            ]].rename(columns = {
                "OrganismID": "SurrogateID",
                "OrganismShortName": "SurrogateShortName",
                "AnnotationSource": "SurrogateAnnotationSource",
                "AnnotationRelease": "SurrogateAnnotationRelease",
                "AnnotationName": "SurrogateAnnotationName",
            })
        )
        print(annot2)
        annot = annot.merge(annot2, how = "left")
        print(annot)
        df = (
            pd.read_csv(input.names, sep = "\t")
            .merge(pd.read_csv(input.datasets, sep = "\t"), how = "left")
            .merge(annot, how = "left")
            .assign(scrna35 = lambda tdf: np.where(
                (tdf.hahnscrna17 == "Yes") | (tdf.liscrna24 == "Yes"),
                "Yes",
                "No"
            ))
            [[
                "OrganismOrder177",
                "OrganismShortName",
                "OrganismID",
                "OrganismCommonName",
                "OrganismScientificName",
                "OrganismColor",
                "IsVertebrate",
                "Datasets",
                "OrganismOrder50",
                "annotated50",
                "bulkrna107",
                "scrna35",
                "hahnscrna17",
                "liscrna24",
                "near29",
                "AnnotationSource",
                "AnnotationRelease",
                "AnnotationName",
                "AnnotationComment",
                "SurrogateID",
                "SurrogateShortName",
                "SurrogateAnnotationSource",
                "SurrogateAnnotationRelease",
                "SurrogateAnnotationName",
            ]]
        )
        print(df)
        df.to_csv(output[0], sep = "\t", index = False)



