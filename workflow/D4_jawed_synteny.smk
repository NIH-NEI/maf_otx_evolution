gtf_files = dict(
    human = "imports/genomes_annotations/human/GCF_000001405.40_GRCh38.p14_genomic.gtf.gz",
    mouse = "imports/genomes_annotations/labMouse/GCF_000001635.27_GRCm39_genomic.gtf.gz",
    opossum = "imports/genomes_annotations/opossum/GCF_027887165.1_mMonDom1.pri_genomic.gtf.gz",
    chicken = "imports/genomes_annotations/chicken/GCF_016699485.2_bGalGal1.mat.broiler.GRCg7b_genomic.gtf.gz",
    frog = "imports/genomes_annotations/westernFrog/GCF_000004195.4_UCB_Xtro_10.0_genomic.gtf.gz",
    gar = "imports/genomes_annotations/spottedGar/GCF_040954835.1_fLepOcu1.hap2_genomic.gtf.gz",
    amphioxus = "imports/genomes_annotations/amphioxus/GCF_035083965.1_klBraLanc5.hap2_genomic.gtf.gz",
)
print(gtf_files)

rule all:
    input:
        expand("scratch/gene_annot_with_strand/{org}_genes.tsv", org = gtf_files),
        "scratch/aligned_single_maf/mafa_alignment.txt",
        #"scratch/aligned_single_maf/single_mafs_aligned.xlsx",
        "scratch/aligned_maf_groups/large_alignment.txt",
        "scratch/aligned_maf_groups/small_alignment.txt",
        "scratch/aligned_maf_groups/otx_alignment.txt",
        "scratch/aligned_maf_groups/large_aligned.xlsx",
        "scratch/aligned_maf_groups/small_aligned.xlsx",
        "scratch/aligned_maf_groups/otx_aligned.xlsx",


rule create_annotation_table:
    input:
        lambda wc: gtf_files[wc.org],
    output:
        "scratch/gene_annot_with_strand/{org}_genes.tsv",
    shell:
        """
        zcat {input} | awk '$3 == "gene"' | cut -f1,4,5,7,9 | cut -d';' -f1 | sed -e 's/gene_id //' > {output}
        """

rule align_single_maf:
    input:
        "configs/maf_names.csv",
        "configs/single_maf_chromosome_direction.csv",
        "workflow/scripts/explore_jawed_synteny/align_direction_single_maf_across_species.py",
    output:
        "scratch/aligned_single_maf/{maf}_alignment.txt"
    shell:
        """
        python workflow/scripts/explore_jawed_synteny/align_direction_single_maf_across_species.py {wildcards.maf}
        """
rule create_excel_single_mafs:
    input:
        expand("scratch/aligned_single_maf/{maf}_alignment.txt",
                maf = ["mafa", "mafb", "cmaf", "nrl", "maff", "mafg", "mafk"])
    output:
        "scratch/aligned_single_maf/single_mafs_aligned.xlsx"
    shell:
        """
        python workflow/scripts/explore_jawed_synteny/create_excel_file.py
        """

rule align_maf_group_large:
    input:
        "scratch/aligned_single_maf/mafa_alignment.txt",
        "scratch/aligned_single_maf/mafb_alignment.txt",
        "scratch/aligned_single_maf/cmaf_alignment.txt",
        "scratch/aligned_single_maf/nrl_alignment.txt",
        "workflow/scripts/explore_jawed_synteny/align_multiple_mafs.py",
    output:
        "scratch/aligned_maf_groups/large_prefixes.txt",
        "scratch/aligned_maf_groups/large_alignment.txt",
    shell:
        """
        python workflow/scripts/explore_jawed_synteny/align_multiple_mafs.py large
        """

rule align_maf_group_small:
    input:
        "scratch/aligned_single_maf/maff_alignment.txt",
        "scratch/aligned_single_maf/mafg_alignment.txt",
        "scratch/aligned_single_maf/mafk_alignment.txt",
        "workflow/scripts/explore_jawed_synteny/align_multiple_mafs.py",
    output:
        "scratch/aligned_maf_groups/small_prefixes.txt",
        "scratch/aligned_maf_groups/small_alignment.txt",
    shell:
        """
        python workflow/scripts/explore_jawed_synteny/align_multiple_mafs.py small
        """

rule align_otx_group:
    input:
        "scratch/aligned_single_maf/otx1_alignment.txt",
        "scratch/aligned_single_maf/otx2_alignment.txt",
        "scratch/aligned_single_maf/crx_alignment.txt",
        "workflow/scripts/explore_jawed_synteny/align_multiple_mafs.py",
    output:
        "scratch/aligned_maf_groups/otx_prefixes.txt",
        "scratch/aligned_maf_groups/otx_alignment.txt",
    shell:
        """
        python workflow/scripts/explore_jawed_synteny/align_multiple_mafs.py otx
        """

#rule create_excel_group_mafs:
#    input:
#        expand("scratch/aligned_maf_groups/{maf}_alignment.txt",
#                maf = ["large"]) #, "small"])
#    output:
#        "scratch/aligned_maf_groups/maf_group_aligned.xlsx"
#    shell:
#        """
#        python workflow/scripts/explore_jawed_synteny/create_excel_file_group_alignments.py
#        """

rule create_excel_group_mafs:
    output:
        "scratch/aligned_maf_groups/{maf_group}_aligned.xlsx"
    shell:
        """
        python workflow/scripts/explore_jawed_synteny/create_excel_file_group_alignments.py {wildcards.maf_group}
        """

