import pandas as pd
from os.path import basename

IMPORTS_PATH = "imports/genomes_annotations"
URL_COMMON_PREFIX = "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF"

rule all:
    input:
        #"scratch/canonical_peptides/treeShrew.faa",
        #"scratch/canonical_peptides/chipmunk.faa",
        #"scratch/canonical_peptides/grassMouse.faa",
        "scratch/canonical_peptides/brownHagfish.faa",
        "scratch/canonical_peptides/europeanRiverLamprey.faa",
        faa = "scratch/canonical_peptides/europeanBrookLamprey.faa",

rule prepare_canonical_proteome_treeShrew:
    """
    Process pep and gtf files downloaded from http://www.treeshrewdb.org/download.html
    """
    input:
        pep = "imports/genomes_annotations/treeShrew/TS_3.0.pep.fa.gz",
        gtf = "imports/genomes_annotations/treeShrew/TS_3.0.genomeannotation.gtf.gz",
        py = "workflow/scripts/compile_proteome/get_canonical_isoform_for_treeShrew.py",
    output:
        faa = "scratch/canonical_peptides/treeShrew.faa",
        info = "scratch/canonical_peptides/treeShrew_info.tsv",
    shell:
        """
        python {input.py} \
            --organism treeShrew \
            --gtf {input.gtf} \
            --protein-faa {input.pep} \
            --out {output.faa} \
            --info {output.info}
        """

rule prepare_canonical_proteome_brownHagfish:
    """
    Process pep and gtf files downloaded from http://www.treeshrewdb.org/download.html
    """
    input:
        pep = "imports/genomes_annotations/brownHagfish/Pata_ah2p.pep.fa",
        gtf = "imports/genomes_annotations/brownHagfish/Pata_MYAQb_rnm_ah2p.gtf",
        py = "workflow/scripts/compile_proteome/get_canonical_isoform_for_brownHagfish.py",
    output:
        faa = "scratch/canonical_peptides/brownHagfish.faa",
        info = "scratch/canonical_peptides/brownHagfish_info.tsv",
    shell:
        """
        python {input.py} \
            --organism brownHagfish \
            --gtf {input.gtf} \
            --protein-faa {input.pep} \
            --out {output.faa} \
            --info {output.info}
        """

rule prepare_canonical_proteome_europeanBrookLamprey:
    """
    Process pep and gtf files downloaded from http://www.treeshrewdb.org/download.html
    """
    input:
        pep = "imports/genomes_annotations/europeanBrookLamprey/kcLamPlan1.2.hap1.proteins.fa",
        gtf = "imports/genomes_annotations/europeanBrookLamprey/kcLamPlan1.2.hap1.gff",
        py = "workflow/scripts/compile_proteome/get_canonical_isoform_for_europeanRiverLamprey.py",
    output:
        faa = "scratch/canonical_peptides/europeanBrookLamprey.faa",
        info = "scratch/canonical_peptides/europeanBrookLamprey_info.tsv",
    shell:
        """
        python {input.py} \
            --organism europeanBrookLamprey \
            --gtf {input.gtf} \
            --protein-faa {input.pep} \
            --out {output.faa} \
            --info {output.info}
        """

rule prepare_canonical_proteome_europeanRiverLamprey:
    """
    Process pep and gtf files downloaded from http://www.treeshrewdb.org/download.html
    """
    input:
        pep = "imports/genomes_annotations/europeanRiverLamprey/kcLamFluv2.2.hap1.proteins.fa",
        gtf = "imports/genomes_annotations/europeanRiverLamprey/kcLamFluv2.2.hap1.gff",
        py = "workflow/scripts/compile_proteome/get_canonical_isoform_for_europeanRiverLamprey.py",
    output:
        faa = "scratch/canonical_peptides/europeanRiverLamprey.faa",
        info = "scratch/canonical_peptides/europeanRiverLamprey_info.tsv",
    shell:
        """
        python {input.py} \
            --organism europeanRiverLamprey \
            --gtf {input.gtf} \
            --protein-faa {input.pep} \
            --out {output.faa} \
            --info {output.info}
        """

#############################################################################
# For chipmunk (Tamias_sibiricus), the originators kept the cds sequences
# and gff3 annotation file at the following figshare url without providing any
# peptide file. So we need to create the peptide fasta file first.
# https://doi.org/10.6084/m9.figshare.20219664
# ############################################################################

rule prepare_protein_sequences_from_chipmunk_cds:
    input:
        cds = "imports/genomes_annotations/chipmunk/Tamias.representive.cds.fa",
    output:
        pep = "imports/genomes_annotations/chipmunk/Tamias.representive.translated_cds.faa",
    run:
        from Bio import SeqIO
        from Bio.SeqRecord import SeqRecord
        import gzip

        def translate_to_protein(dna_record):
            # Translate in frame 1, stop at first stop codon (recommended for CDS)
            protein_seq = dna_record.seq.translate(to_stop=True)

            # Alternative: include * for stop codons (uncomment if preferred)
            # protein_seq = dna_record.seq.translate(table="Standard", cds=False)

            return SeqRecord(
                seq=protein_seq,
                id=dna_record.id,
                description=dna_record.description  # preserves original description
            )

        # Open input and output as gzip
        with open(input.cds, "r") as infile, open(output.pep, "w") as outfile:
            # Parse input FASTA (text mode 'rt' for gzip)
            input_records = SeqIO.parse(infile, "fasta")

            # Translate each record
            protein_records = (translate_to_protein(rec) for rec in input_records)

            # Write directly to compressed output
            SeqIO.write(protein_records, outfile, "fasta")

        print(f"Conversion complete: {input.cds} -> {output.pep}")

rule fix_chipmunk_gff:
    input:
        gff = "imports/genomes_annotations/chipmunk/Tamias.representive.gff",
    output:
        gff = "imports/genomes_annotations/chipmunk/Tamias.representive_fixed.gff",
    shell:
        """
        module load agat
        agat_convert_sp_gxf2gxf.pl --gff {input.gff} -o {output.gff}
        """

rule prepare_canonical_proteome_chipmunk:
    input:
        pep = "imports/genomes_annotations/chipmunk/Tamias.representive.translated_cds.faa",
        gff = "imports/genomes_annotations/chipmunk/Tamias.representive_fixed.gff",
        py = "workflow/scripts/compile_proteome/get_canonical_isoform_for_chipmunk.py",
    output:
        faa = "scratch/canonical_peptides/chipmunk.faa",
        info = "scratch/canonical_peptides/chipmunk.tsv",
    shell:
        """
        python {input.py} \
            --organism chipmunk \
            --gff {input.gff} \
            --protein-faa {input.pep} \
            --out {output.faa} \
            --info {output.info}
        """

#############################################################################
# For grassMouse (Rhabdomys_pumilio), the originators kept the genome assembly
# and gff3 annotation file at the following figshare url without providing any
# peptide file. So we need to create the peptide fasta file first.
# https://doi.org/10.26188/20321655
# ############################################################################

rule prepare_protein_sequences_from_grassMouse_genome:
    input:
        genome = "imports/genomes_annotations/grassMouse/Rhabdomys_pumilio_Princeton_asm1.0_preNCBI.fasta",
        gff = "imports/genomes_annotations/grassMouse/Rhabdomys_pumilio.mouse_gene_name_final_fixed_trimmed.gff",
    output:
        pep = "imports/genomes_annotations/grassMouse/Rhabdomys_pumilio.mouse_gene_name_final.pep.fa",
    shell:
        """
        module load cufflinks
        gffread -y {output.pep} -g {input.genome} {input.gff}
        """

rule prepare_canonical_proteome_grassMouse:
    """
    Process pep and gtf files downloaded from http://www.treeshrewdb.org/download.html
    """
    input:
        pep = "imports/genomes_annotations/grassMouse/Rhabdomys_pumilio.mouse_gene_name_final.pep.fa",
        gff = "imports/genomes_annotations/grassMouse/Rhabdomys_pumilio.mouse_gene_name_final_fixed_trimmed.gff",
        py = "workflow/scripts/compile_proteome/get_canonical_isoform_for_grassMouse.py",
    output:
        faa = "scratch/canonical_peptides/grassMouse.faa",
        info = "scratch/canonical_peptides/grassMouse_info.tsv",
    shell:
        """
        python {input.py} \
            --organism grassMouse \
            --gff {input.gff} \
            --protein-faa {input.pep} \
            --out {output.faa} \
            --info {output.info}
        """
