
rule all:
    input:
        "configs/species181_table.tsv"

rule prepare_species_table:
    input:
        names = "configs/species181_names_orders.tsv",
        datasets = "configs/species181_datasets.tsv",
        annotations = "configs/species181_annotations.tsv",
        represent = "configs/species181_representatives.tsv",
    output:
        table = "configs/species_table.tsv"
    run:
        import pandas as pd
        import numpy as np
        annot = pd.read_csv(input.annotations, sep = "\t")
        annot["SurrogateID"] = annot["SurrogateID"].fillna(annot["OrganismID"])
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
        annot = annot.merge(annot2, how = "left")
        df = (
            pd.read_csv(input.names, sep = "\t")
            .merge(pd.read_csv(input.represent, sep = "\t"), how = "left")
            .merge(pd.read_csv(input.datasets, sep = "\t"), how = "left")
            .merge(annot, how = "left")
            .assign(scrna35 = lambda tdf: np.where(
                (tdf.hahnscrna17 == "Yes") | (tdf.liscrna24 == "Yes"),
                "Yes",
                "No"
            ))
            [[
                "OrganismOrder181",
                "OrganismShortName",
                "OrganismID",
                "OrganismCommonName",
                "OrganismScientificName",
                "OrganismColor",
                "IsVertebrate",
                "Representative",
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
        df.to_csv(output[0], sep = "\t", index = False)



