import pandas as pd
import os

LISCRNA_SRA_TOP = "/data/VisionDevo/LiRetinaScRNA" #/GSM8246586_SRR28879684/SRR28879684_1.fastq.gz"

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
)
print(samples)

organisms = (
    organisms
    .set_index("OrganismShortName", drop = False)
)

rule all:
    input:
        #expand("results/seurat/{species}_seurat.rds", species=SPECIES),
        #"results/quality_control/multiqc_report.html"
        expand(LISCRNA_SRA_TOP+"/BySample/{sample}", sample = samples.SampleID),
        expand(LISCRNA_SRA_TOP+"/cellranger/filtered_gtf/{org}_filtered.gtf", org = samples.OrganismShortName.unique()),
        expand(LISCRNA_SRA_TOP+"/cellranger/reference_genome/{org}_genome/reference.json", org = samples.OrganismShortName.unique()),
        expand(LISCRNA_SRA_TOP+"/cellranger/counts/{sample}/outs/filtered_feature_bc_matrix/matrix.mtx.gz", sample = samples.SampleID),

rule merge_fastqs:
    """Merge multiple FASTQ files per species if needed"""
    input:
        r1=lambda wildcards: expand("{prefix}_1.fastq.gz", prefix=samples["prefix"][wildcards.sample]),
        r2=lambda wildcards: expand("{prefix}_2.fastq.gz", prefix=samples["prefix"][wildcards.sample]),
    output:
        r1=LISCRNA_SRA_TOP+"/merged_fastqs_sample_wise/{sample}_R1.fastq.gz",
        r2=LISCRNA_SRA_TOP+"/merged_fastqs_sample_wise/{sample}_R2.fastq.gz",
    shell:
        """
        cat {input.r1} > {output.r1}
        cat {input.r2} > {output.r2}
        """

rule rename_fastq_for_cell_ranger:
    input:
        root = lambda wc: LISCRNA_SRA_TOP+"/ByGSM/"+samples["GSM"][wc.sample],
        info = "configs/scrna_data.tsv",
    output:
        LISCRNA_SRA_TOP+"/BySample/{sample}/cellranger_fastq_links.summary.tsv",
    shell:
        """
        python workflow/scripts/process_liscrna/rename_fastq_for_cellranger.py \
            --root {input.root} \
            --map_tsv {input.info} \
            --out {output} \
            --mode symlink \
            #--dry_run \
            --reads_for_lane 200
        """
def get_gtf(wc):
    surrogate_org = organisms["SurrogateShortName"][wc.org]
    annotation_name = organisms["SurrogateAnnotationName"][wc.org]
    annotation_rel = organisms["SurrogateAnnotationRelease"][wc.org]
    gtf = f"imports/genomes_annotations/{surrogate_org}/{annotation_rel}_{annotation_name}_genomic.gtf"
    gtfgz = gtf + ".gz"
    return (gtf, gtfgz)

rule make_cellranger_gtf:
    input:
        gtfgz=lambda wc: get_gtf(wc)[1],
    output:
        filtered=LISCRNA_SRA_TOP + "/cellranger/filtered_gtf/{org}_filtered.gtf",
    params:
        tmpgtf=lambda wc: get_gtf(wc)[0],
    shell:
        """
        module load cellranger
        gunzip -c {input.gtfgz} > {params.tmpgtf}
        cellranger mkgtf \
            {params.tmpgtf} \
            {output.filtered} \
            --attribute=gene_biotype:protein_coding
        rm -f {params.tmpgtf}
        """

def get_fasta(wc):
    surrogate_org = organisms["SurrogateShortName"][wc.org]
    annotation_name = organisms["SurrogateAnnotationName"][wc.org]
    annotation_rel = organisms["SurrogateAnnotationRelease"][wc.org]
    return f"imports/genomes_annotations/{surrogate_org}/{annotation_rel}_{annotation_name}_genomic.fna.gz"

rule make_cellranger_reference:
    input:
        fasta_gz=lambda wc: get_fasta(wc),
        filtered_gtf=LISCRNA_SRA_TOP + "/cellranger/filtered_gtf/{org}_filtered.gtf",
    output:
        refjson=LISCRNA_SRA_TOP + "/cellranger/reference_genome/{org}_genome/reference.json",
    params:
        ref_parent=LISCRNA_SRA_TOP + "/cellranger/reference_genome",
        genome_name=lambda wc: f"{wc.org}_genome",
    shell:
        r"""
        module load cellranger

        tmpdir=/lscratch/$SLURM_JOB_ID
        mkdir -p $tmpdir

        tmp_fasta=$(mktemp $tmpdir/{wildcards.org}.XXXXXX.fa)
        gunzip -c {input.fasta_gz} > $tmp_fasta

        cd {params.ref_parent}

        rm -rf {params.genome_name}

        cellranger mkref \
            --genome={params.genome_name} \
            --fasta=$tmp_fasta \
            --genes={input.filtered_gtf} \
            --memgb=100
        """

#rule fastqc:
#    """Quality control of raw reads"""
#    input:
#        "data/fastq/{species}/{srr}_{read}.fastq.gz"
#    output:
#        html="results/qc/fastqc/{species}/{srr}_{read}_fastqc.html",
#        zip="results/qc/fastqc/{species}/{srr}_{read}_fastqc.zip"
#    params:
#        outdir="results/qc/fastqc/{species}"
#    threads: 2
#    conda:
#        "envs/qc.yaml"
#    shell:
#        """
#        mkdir -p {params.outdir}
#        fastqc -t {threads} -o {params.outdir} {input}
#        """

rule cellranger_count:
    """Run Cell Ranger count for each sample"""
    input:
        fastq_info=LISCRNA_SRA_TOP + "/BySample/{sample}/cellranger_fastq_links.summary.tsv",
        refjson=lambda wc: (
            LISCRNA_SRA_TOP
            + "/cellranger/reference_genome/"
            + samples["OrganismShortName"][wc.sample]
            + "_genome/reference.json"
        ),
    output:
        bam=LISCRNA_SRA_TOP + "/cellranger/counts/{sample}/outs/possorted_genome_bam.bam",
        matrix=LISCRNA_SRA_TOP + "/cellranger/counts/{sample}/outs/filtered_feature_bc_matrix/matrix.mtx.gz",
        features=LISCRNA_SRA_TOP + "/cellranger/counts/{sample}/outs/filtered_feature_bc_matrix/features.tsv.gz",
        barcodes=LISCRNA_SRA_TOP + "/cellranger/counts/{sample}/outs/filtered_feature_bc_matrix/barcodes.tsv.gz",
        metrics=LISCRNA_SRA_TOP + "/cellranger/counts/{sample}/outs/metrics_summary.csv",
    params:
        run_id=lambda wc: wc.sample,
        transcriptome=lambda wc: (
            LISCRNA_SRA_TOP
            + "/cellranger/reference_genome/"
            + samples["OrganismShortName"][wc.sample]
            + "_genome"
        ),
        fastq_dir=LISCRNA_SRA_TOP + "/BySample/{sample}",
        outdir=LISCRNA_SRA_TOP + "/cellranger/counts",
    threads: 16
    resources:
        mem_mb=64000
    shell:
        r"""
        module load cellranger

        cd {params.outdir}

        # Cell Ranger CANNOT overwrite — clean stale dirs
        if [ -d "{params.run_id}" ]; then
            echo "Removing existing Cell Ranger output {params.run_id}"
            rm -rf "{params.run_id}"
        fi

        cellranger count \
            --create-bam true \
            --id={params.run_id} \
            --transcriptome={params.transcriptome} \
            --fastqs={params.fastq_dir} \
            --sample={params.run_id} \
            --localcores={threads} \
            --localmem=60 \
            --chemistry=auto
        """

#rule create_seurat_object:
#    """Create Seurat object from Cell Ranger output"""
#    input:
#        matrix="results/cellranger/{species}/outs/filtered_feature_bc_matrix/matrix.mtx.gz",
#        features="results/cellranger/{species}/outs/filtered_feature_bc_matrix/features.tsv.gz",
#        barcodes="results/cellranger/{species}/outs/filtered_feature_bc_matrix/barcodes.tsv.gz",
#        metrics="results/cellranger/{species}/outs/metrics_summary.csv"
#    output:
#        seurat="results/seurat/{species}_seurat.rds",
#        qc_plots="results/seurat/{species}_qc_plots.pdf"
#    params:
#        species="{species}",
#        min_features=200,
#        max_features=7000,
#        max_mt_percent=20
#    threads: 4
#    resources:
#        mem_mb=32000
#    conda:
#        "envs/seurat.yaml"
#    script:
#        "scripts/create_seurat_object.R"
#
#rule multiqc:
#    """Aggregate QC reports"""
#    input:
#        fastqc=expand("results/qc/fastqc/{species}/{srr}_{read}_fastqc.zip",
#                     species=SPECIES, srr=lambda w: SAMPLES[SPECIES[0]], read=["1", "2"]),
#        cellranger=expand("results/cellranger/{species}/outs/metrics_summary.csv",
#                         species=SPECIES)
#    output:
#        "results/quality_control/multiqc_report.html"
#    params:
#        search_dir="results"
#    conda:
#        "envs/qc.yaml"
#    shell:
#        """
#        multiqc {params.search_dir} -o results/quality_control -f
#        """
#
#rule create_integrated_atlas:
#    """Create integrated cross-species atlas"""
#    input:
#        seurat_objects=expand("results/seurat/{species}_seurat.rds", species=SPECIES)
#    output:
#        integrated="results/integrated/cross_species_atlas.rds",
#        plots="results/integrated/integration_plots.pdf"
#    threads: 8
#    resources:
#        mem_mb=128000
#    conda:
#        "envs/seurat.yaml"
#    script:
#        "scripts/integrate_species.R"
