import pandas as pd

#    "https://www.ncbi.nlm.nih.gov/geo/download/?acc={GSM}&format=file&file={GSM}_{SampleCode}_{ftype}.gz"

LI_URL_FMT = (
    "https://ftp.ncbi.nlm.nih.gov/geo/samples/"
    "{gsm_nnn}/{GSM}/suppl/{GSM}_{SampleCode}_{ftype}.gz"
)

samples = (
    pd.read_csv("configs/scrna_data.tsv", sep = "\t")
    .query("Dataset == 'liscrna24'")
    .set_index("SampleID", drop = False)
#    .query("OrganismShortName == 'crabMacaque'")
)

print(samples)

rule all:
    input:
        #"scratch/scrna_dotplots/allorgs_dotplot.pdf"
        expand("scratch/liscrna_seurat/{org}/{sid}__basic_seurat.rds", zip, org = samples.OrganismShortName, sid = samples.SampleID),


rule download_geo_files_liscrna:
    output:
        "imports/liscrna/{org}/{sid}/{ftype}.gz",
    run:
        GSM = samples["GSM"][wildcards.sid]
        SampleCode = samples["SampleCode"][wildcards.sid]
        ftype = wildcards.ftype  # e.g. "matrix.mtx" or "barcodes.tsv" etc, depending on your naming
        # GSM8246593 -> GSM8246nnn
        gsm_nnn = GSM[:-3] + "nnn"
        url = LI_URL_FMT.format(
            GSM = GSM,
            SampleCode = SampleCode,
            ftype = ftype,
            gsm_nnn = gsm_nnn
        )
        cmd = "wget -O " + output[0] + " " + url
        print(cmd)
        shell(cmd)

rule create_basic_seurat_objects:
    input:
        barcodes = "imports/liscrna/{org}/{sid}/barcodes.tsv.gz",
        features = "imports/liscrna/{org}/{sid}/features.tsv.gz",
        matrix = "imports/liscrna/{org}/{sid}/matrix.mtx.gz",
    params:
        data_dir = "imports/liscrna/{org}/{sid}",
    output:
        seu = "scratch/liscrna_seurat/{org}/{sid}__basic_seurat.rds",
        pdf = "scratch/liscrna_seurat/{org}/{sid}__basic_seurat.pdf",
    script:
        "scripts/process_liscrna/create_basic_seurat_object.R"

rule draw_dot_plots:
    input:
        seu = lambda wc: SHEKHAR_TOP + f"/SeuratObjects/{my2hahn[wc.org]}_initial.rds",
    output:
        "scratch/scrna_dotplots/{org}_dotplot.pdf"
    script:
        "scripts/visualize_scrna/draw_hahn_dotplots_maf_otx.R"

rule extract_dot_plot_data:
    input:
        seu = lambda wc: SHEKHAR_TOP + f"/SeuratObjects/{my2hahn[wc.org]}_initial.rds",
    output:
        "scratch/scrna_dotplots/{org}_dotplot.tsv"
    script:
        "scripts/visualize_scrna/extract_dotplot_data_maf_otx.R"

rule combine_dotplots_data:
    input:
        expand("scratch/scrna_dotplots/{org}_dotplot.tsv", org = samples.OrganismShortName),
    output:
        "scratch/scrna_dotplots/allorgs_dotplot.tsv"
    run:
        import pandas as pd
        df = pd.concat([
            pd.read_csv(infile, sep = "\t")
            for infile in input
        ])
        df.to_csv(output[0], sep = "\t")

rule create_combine_dotplot:
    input:
        "scratch/scrna_dotplots/allorgs_dotplot.tsv"
    output:
        "scratch/scrna_dotplots/allorgs_dotplot.pdf"
    script:
        "scripts/visualize_scrna/draw_hahn_combined_orgs_dotplots_maf_otx.R"
