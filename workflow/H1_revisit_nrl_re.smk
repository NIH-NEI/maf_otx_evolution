# # -------------------------------------------------------------------------
# Overview
# -------------------------------------------------------------------------
# This Snakemake workflow supports a reanalysis of NRL regulatory elements
# originally described in Kautzmann et al. (2011). That study was based on an
# earlier release of the UCSC mouse genome (mm9, 2007) and comparative genomics
# resources derived from multi-species alignments of 30 vertebrates.
#
# In contrast, the current UCSC human genome assembly (hg38, 2013) provides
# comparative genomics resources based on alignments across approximately
# 100 vertebrate species, enabling a more comprehensive evolutionary analysis.
# We therefore reanalyze these regulatory elements using updated genomic
# assemblies and alignment resources.
#
# Deriving hg38 genomic coordinates for clusters A, B, and C requires an
# indirect mapping strategy. Human DNA sequences corresponding to each cluster
# are first identified by manual inspection of UCSC comparative genomics tracks
# that include alignments to earlier human assemblies (e.g., hg18, 2006).
# These human sequences are then mapped back to hg38 to recover precise
# genomic coordinates.
#
# This process is further complicated by the fact that Figure 3 in Kautzmann
# et al. (2011) depicts alignments that do not correspond exactly to the full
# cluster boundaries defined elsewhere in the paper, but instead represent
# approximate subregions of those clusters. These approximate sequences are
# therefore used as anchors for coordinate recovery and cross-species mapping.
#
# The overall analysis consists of both manual curation steps—guided by figures
# in the original publication and interactive use of the UCSC Genome Browser—
# and automated steps implemented in this workflow.
#
# This file documents all manual and automated steps for reproducibility and
# long-term reference, while the Snakemake rules themselves implement only the
# automated components.

# -------------------------------------------------------------------------
# Step 1 [Manual]: Define approximate mouse cluster sequences
# -------------------------------------------------------------------------
# Begin with the mouse sequences shown in Figure 3 of Kautzmann et al. (2011).
# Importantly, the genomic regions aligned in Figure 3 only approximately
# correspond to the cluster boundaries defined in Figure 1B:
#
#   - Cluster A is defined as [-304, +119] relative to the Nrl TSS in mm9
#     (Figure 1B), whereas Figure 3B shows an alignment for [-117, +137].
#
#   - Cluster B is defined as [-938, -657] (Figure 1B), whereas Figure 3A
#     shows an alignment for [-909, -667].
#
# Consequently, the sequences extracted from Figure 3 represent approximate
# subregions of clusters A and B rather than their full extents.
#
# The corresponding mouse (mm9) DNA sequences are stored in:
#   1. imports/kautzmann_clusters/approx_clusterA_mm9_mouse_dna.txt
#   2. imports/kautzmann_clusters/approx_clusterB_mm9_mouse_dna.txt

# -------------------------------------------------------------------------
# Step 2 [Automated]: Recover approximate mm9 coordinates for clusters A and B
# -------------------------------------------------------------------------
# Using the approximate mouse sequences from Step 1, a Python script queries
# UCSC to recover the precise mm9 genomic coordinates corresponding to each
# cluster.
#
# The resulting coordinates are saved in:
#   3. imports/kautzmann_clusters/clusterA_mm9_coordinates.txt or .psl
#   4. imports/kautzmann_clusters/clusterB_mm9_coordinates.txt or .psl

# -------------------------------------------------------------------------
# Step 3 [Manual]: Create precise mm9 coordinates for all three clusters
# -------------------------------------------------------------------------
# See the comments in the following file for calculation of mm9 co-ordinates
# for all three clusters A, B, C
#   5. imports/kautzmann_clusters/mm9_track_clusters_ABC.bed

# -------------------------------------------------------------------------
# Step 4 [Manual]: Identify aligned human sequences via comparative genomics
# -------------------------------------------------------------------------
# Using the Comparative Genomics track in the mm9 UCSC Genome Browser,
# identify the human DNA sequences aligned to each mouse cluster region.
# Specifically, right click on the multiz track and click "show details for multiz30way..."
# It will show the alignments. Keep the human alignments by using "grep human" and then
# removing intitial part  upto "Human  "
#  eg: "cat alignment.txt | grep -i human | sed "s/Human  /\t/" | cut -f2
#
# The extracted human DNA alignments are saved in:
#   6. imports/kautzmann_clusters/clusterA_mm9_human_alignment.txt
#   7. imports/kautzmann_clusters/clusterB_mm9_human_alignment.txt
#   8. imports/kautzmann_clusters/clusterC_mm9_human_alignment.txt

# -------------------------------------------------------------------------
# Step 5 [Automated]: Recover hg38 coordinates for human cluster regions
# -------------------------------------------------------------------------
# Using the human-aligned sequences from Step 4, a Python script queries UCSC
# to determine the corresponding hg38 genomic coordinates for clusters A, B,
# and C.
#
# The resulting coordinates are saved in:
#   9.  imports/kautzmann_clusters/clusterA_hg38_coordinates.txt
#   10. imports/kautzmann_clusters/clusterB_hg38_coordinates.txt
#   11. imports/kautzmann_clusters/clusterC_hg38_coordinates.txt

# -------------------------------------------------------------------------
# Step 6 [Automated]: Generate a consolidated hg38 BED file
# -------------------------------------------------------------------------
# Finally, the hg38 coordinates for clusters A, B, and C are combined into a
# single BED file for visualization in the UCSC Genome Browser.
#
# The consolidated BED file is saved as:
#   12. imports/kautzmann_clusters/clusters_ABC_hg38_coordinates.bed
#


rule all:
    input:
        #"imports/kautzmann_clusters/clusterA_mm9_coordinates_tmp.txt",
        #"imports/kautzmann_clusters/clusterB_mm9_coordinates_tmp.txt",
        #"imports/kautzmann_clusters/clusterA_hg38_coordinates_tmp.txt",
        #"imports/kautzmann_clusters/clusterB_hg38_coordinates_tmp.txt",
        "imports/kautzmann_clusters/clusterC_hg38_coordinates_tmp.txt",

rule search_dna_ucsc_genome_mm9:
    """
    If the script fails, use local blat 
    """
    input:
        dna="imports/kautzmann_clusters/approx_{cluster}_mm9_mouse_dna.txt",
        scrpt="workflow/scripts/reanalyze_nrl_regulatory_elements/search_dna_in_ucsc_genome.py",
    output:
        "imports/kautzmann_clusters/{cluster}_mm9_coordinates_tmp.txt",
    shell:
        "python {input.scrpt} {input.dna} mm9"

######################
# If the above script fails, you can use local blat program
# 1. Create a single sequence fasta from the sequence
# 2. module load blat
# 3. faToNib imports/kautzmann_clusters/approx_clusterA_mm9_mouse_dna.fa imports/kautzmann_clusters/approx_clusterA_mm9_mouse_dna.nib
# 4. blat /fdb/igenomes/Mus_musculus/UCSC/mm9/Sequence/WholeGenomeFasta/genome.fa cluaterA.nib clusterA.psl
######################

rule search_dna_ucsc_genome_hg38:
    """
    If the script fails, use local blat 
    """
    input:
        dna="imports/kautzmann_clusters/{cluster}_mm9_human_alignment.txt",
        scrpt="workflow/scripts/reanalyze_nrl_regulatory_elements/search_dna_in_ucsc_genome.py",
    output:
        "imports/kautzmann_clusters/{cluster}_hg38_coordinates_tmp.txt",
    shell:
        "python {input.scrpt} {input.dna} hg38"

######################
# If the above script fails, you can use local blat program
# 1. module load blat
# 2. faToNib imports/kautzmann_clusters/clusterA_mm9_human_dna.fa imports/kautzmann_clusters/clusterA_mm9_human_dna.nib
# 3. blat /fdb/igenomes/Homo_sapiens/UCSC/hg38/Sequence/WholeGenomeFasta/genome.fa  \
#        imports/kautzmann_clusters/clusterA_mm9_human_dna.nib \
#        imports/kautzmann_clusters/clusterA_hg38_coordinates.psl

######################
