

rule all:
    input:
        "exports/tree_alignment_bundle/species_maf_otx_trees_alignments.nex"


rule create_nexus:
    input:
        index = "configs/nexus_bundle_index.csv",
        scrpt = "workflow/scripts/draw_tree_alignment/make_combined_nexus.py",
    output:
        "exports/tree_alignment_bundle/species_maf_otx_trees_alignments.nex"
    shell:
        """
        python {input.scrpt} \
            --index {input.index} \
            --output {output} \
            --datatype PROTEIN
        """
