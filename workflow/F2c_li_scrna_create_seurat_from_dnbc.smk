import pandas as pd

LISCRNA_DNBC_COUNTS_TOP = "/data/VisionDevo/LiRetinaScRNA/dnbc4tools/counts" #/GSM8246586_SRR28879684/SRR28879684_1.fastq.gz"

organisms = (
    pd.read_csv("configs/species_table.tsv", sep = "\t")
    #[["OrganismShortName", "OrganismID", "SurrogateShortName", "SurrogateID", "SurrogateAnnotationSource"]]
)
print(organisms)

samples = (
    pd.read_csv("configs/scrna_data.tsv", sep = "\t")
    .query("Dataset == 'liscrna24'")
    #.query("OrganismShortName not in ['human']")
    [["OrganismShortName", "SampleID", "GSM"]]
    .merge(organisms, how = "left")
    .set_index("SampleID", drop = False)
    .query("OrganismShortName == SurrogateShortName")
    .query("SurrogateAnnotationSource == 'GCF'")
    .query("OrganismShortName != 'human'")
)
print(samples)

rule all:
    input:
        #"scratch/scrna_dotplots/allorgs_dotplot.pdf"
        expand("scratch/liscrna_dnbc_seurat/{org}/{sid}__basic_seurat.rds", zip, org = samples.OrganismShortName, sid = samples.SampleID),


rule create_basic_seurat_objects:
    input:
        barcodes = LISCRNA_DNBC_COUNTS_TOP + "/{sid}/outs/filter_matrix/barcodes.tsv.gz",
        features = LISCRNA_DNBC_COUNTS_TOP + "/{sid}/outs/filter_matrix/features.tsv.gz",
        matrix = LISCRNA_DNBC_COUNTS_TOP + "/{sid}/outs/filter_matrix/matrix.mtx.gz",
    params:
        data_dir = LISCRNA_DNBC_COUNTS_TOP + "/{sid}/outs/filter_matrix",
    output:
        seu = "scratch/liscrna_dnbc_seurat/{org}/{sid}__basic_seurat.rds",
        pdf = "scratch/liscrna_dnbc_seurat/{org}/{sid}__basic_seurat.pdf",
    script:
        "scripts/process_liscrna/create_basic_seurat_object.R"

#rule draw_dot_plots:
#    input:
#        seu = lambda wc: SHEKHAR_TOP + f"/SeuratObjects/{my2hahn[wc.org]}_initial.rds",
#    output:
#        "scratch/scrna_dotplots/{org}_dotplot.pdf"
#    script:
#        "scripts/visualize_scrna/draw_hahn_dotplots_maf_otx.R"
#
#rule extract_dot_plot_data:
#    input:
#        seu = lambda wc: SHEKHAR_TOP + f"/SeuratObjects/{my2hahn[wc.org]}_initial.rds",
#    output:
#        "scratch/scrna_dotplots/{org}_dotplot.tsv"
#    script:
#        "scripts/visualize_scrna/extract_dotplot_data_maf_otx.R"
#
#rule combine_dotplots_data:
#    input:
#        expand("scratch/scrna_dotplots/{org}_dotplot.tsv", org = samples.OrganismShortName),
#    output:
#        "scratch/scrna_dotplots/allorgs_dotplot.tsv"
#    run:
#        import pandas as pd
#        df = pd.concat([
#            pd.read_csv(infile, sep = "\t")
#            for infile in input
#        ])
#        df.to_csv(output[0], sep = "\t")
#
#rule create_combine_dotplot:
#    input:
#        "scratch/scrna_dotplots/allorgs_dotplot.tsv"
#    output:
#        "scratch/scrna_dotplots/allorgs_dotplot.pdf"
#    script:
#        "scripts/visualize_scrna/draw_hahn_combined_orgs_dotplots_maf_otx.R"
