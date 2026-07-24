rule all:
    input:
        "scratch/grand_list/grand_list_domains.tsv",

rule mark_otx_domains:
    input:
        ips = "scratch/grand_list_interpro/grand_list_interproscan.tsv",
    output:
        "scratch/grand_list/grand_list_domains.tsv",
    script:
        "scripts/analyze_coevolution/mark_maf_otx_domains.py"

