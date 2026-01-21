import pandas as pd
from os.path import basename

TRINITY_RES_TOP = "/data/VisionEvo/TrinityOutput/Assemblies"
TRANSDECODER_RES_TOP = "/data/VisionEvo/TransDecoderOutput/OrganismLevel"
TRINOTATE_RES_TOP = "/data/VisionEvo/TrinotateOutput/Reports"
TRINOTATE_DATA_DIR = "/data/pals2/tools/Trinotate_data_dir"

rnaseq = (
    pd.read_csv("configs/species173_table.tsv", sep = "\t")
#    .merge((
#        pd.read_csv("configs/species173_annotations.tsv", sep = "\t")
#        [["OrganismID", "AnnotationSource", "AnnotationName", "AnnotationRelease"]]
#    ), how = "left")
    .query("bulkrna107 == 'Yes'")
)
print(rnaseq)

rule all:
    input:
        expand(TRINOTATE_RES_TOP + "/{org}/{org}_trianotate_report.tsv", org  = rnaseq.OrganismID.unique()),

rule generate_trinty_gene2trans_map:
    input:
        TRINITY_RES_TOP + "/{org}_denovo.fasta",
    output:
        TRINITY_RES_TOP + "/{org}_denovo.fasta.gene_to_trans_map",
    shell:
        """
        module load trinity
        $TRINITY_HOME/util/support_scripts/get_Trinity_gene_to_trans_map.pl {input} > {output}
        """
rule run_trinotate:
    input:
        trans = TRINITY_RES_TOP + "/{org}_denovo.fasta",
        g2t = TRINITY_RES_TOP + "/{org}_denovo.fasta.gene_to_trans_map",
        pep = TRANSDECODER_RES_TOP + "/{org}/{org}_denovo.fasta.transdecoder.pep",
    output:
        db = TRINOTATE_RES_TOP + "/{org}/{org}_trianotate_db.sqlite",
        rep = TRINOTATE_RES_TOP + "/{org}/{org}_trianotate_report.tsv",
    threads: 4
    shell:
        """
        module load trinotate
        cd $(dirname {output.db})
        Trinotate --db {output.db} --create --trinotate_data_dir {TRINOTATE_DATA_DIR}
        Trinotate --db {output.db} --init \
            --gene_trans_map {input.g2t} \
            --transcript_fasta {input.trans} \
            --transdecoder_pep {input.pep}
        Trinotate --db {output.db} --run ALL \
            --trinotate_data_dir {TRINOTATE_DATA_DIR} \
            --transcript_fasta {input.trans} \
            --transdecoder_pep {input.pep} \
            --use_diamond
        Trinotate --db {output.db} --report \
            --incl_pep --incl_trans > {output.rep}
        """
