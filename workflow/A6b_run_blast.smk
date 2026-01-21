from io import StringIO
import pandas as pd
import glob

import pandas as pd
import glob

BLASTDB_DIR = "/data/pals2/tools/blastdb"

rnaseq = (
    pd.read_csv("configs/species173_table.tsv", sep = "\t")
    .query("bulkrna107 == 'Yes'")
)
print(rnaseq)

short2id = rnaseq.set_index("OrganismShortName")["OrganismID"].to_dict()
print(short2id)


localrules: all, consolidate_all_blast_results

rule all:
    input:
        expand("scratch/blast_query_proteins/{short}/best_db_{short}__query_maf_otx.tsv",
                short=rnaseq.OrganismShortName.unique()),


rule run_blast:
    input:
        db = lambda wc: f"scratch/blastdb_per_organism/{short2id[wc.short]}/{short2id[wc.short]}",
        query = "scratch/annotated_proteins/annotated_maf_otx.faa",
    output:
        "scratch/blast_query_proteins/{short}/db_{short}__query_maf_otx.asn",
    threads: 16
    group: "runblast"
    resources:
        runtime = 2*24*60,
        mem_mb = 64*1024,
        disk_mb = 256*1024,
    params:
        dbsize = 44834480755, # from "blastdbcmd -info -db scratch/blastdb/sample/blastdb_mode_sample"
    shell:
        """
        module load blast
        export BLASTDB="{BLASTDB_DIR}"
        db="$(readlink -f {input.db})"
        blastp -db $db -query {input.query} \
            -out {output} -num_threads {threads} \
            -outfmt 11 -dbsize {params.dbsize} \
            -max_target_seqs 1000000
        """

rule format_blast_results:
    input:
        "scratch/blast_query_proteins/{short}/db_{short}__query_maf_otx.asn",
    output:
        "scratch/blast_query_proteins/{short}/db_{short}__query_maf_otx.tsv",
    threads: 16
    group: "formatblast"
    resources:
        runtime = 2*24*60,
        mem_mb = 16*1024,
        disk_mb = 32*1024,
    shell:
        """
        module load blast
        export BLASTDB="{BLASTDB_DIR}"
        blast_formatter -archive {input} -outfmt "6 qseqid qgi qacc qaccver qlen sseqid sallseqid sgi sallgi sacc saccver sallacc slen qstart qend sstart send qseq sseq evalue bitscore score length pident nident mismatch positive gapopen gaps ppos frames qframe sframe btop staxid ssciname scomname sblastname sskingdom staxids sscinames scomnames sblastnames sskingdoms stitle salltitles sstrand qcovs qcovhsp qcovus" > {output}
        """

rule get_best_blast_match:
    input:
        blast_res = "scratch/blast_query_proteins/{short}/db_{short}__query_maf_otx.tsv",
    output:
        "scratch/blast_query_proteins/{short}/best_db_{short}__query_maf_otx.tsv",
    threads: 1
    group: "bestblast"
    resources:
        runtime = 2*24*60,
        mem_mb = 16*1024,
        disk_mb = 32*1024,
    script:
        "scripts/run_blast/filter_best_blast_match.py"
