import pandas as pd
from os.path import basename

IMPORTS_PATH = "imports/genomes_annotations"
URL_COMMON_PREFIX = "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF"

suffixes = dict(
    translated_cds = "_translated_cds.faa.gz",
    untranslated_cds = "_cds_from_genomic.fna.gz",
    genomic_gff = "_genomic.gff.gz",
)

urls = (
    pd.read_csv("configs/ncbi_urls.tsv", sep = "\t")
    .assign(acc = lambda tdf: tdf.GCFString.str.split("/").str[-1])
    .assign(url_prefix = lambda tdf: URL_COMMON_PREFIX + "/" + tdf.GCFString)
    .assign(path_prefix = lambda tdf: IMPORTS_PATH + "/" + tdf.ShortName)
    .merge(pd.DataFrame(suffixes.items(), columns = ["file_type", "suffix"]), how = "cross")
    .assign(remote = lambda tdf: tdf.url_prefix + tdf.suffix)
    .assign(local = lambda tdf: tdf.GCFString.apply(lambda x: basename(x)))
    .assign(local = lambda tdf: tdf.path_prefix + "/" + tdf.local + tdf.suffix)
)
print(urls)

path2url = urls.set_index("local")["remote"].to_dict()
org2acc = urls.set_index("ShortName")["acc"].to_dict()


rule all:
    input:
        urls.local,
        expand("scratch/canonical_peptides/{org}.faa", org  = urls.ShortName.unique()),


rule download_ncbi_data:
    output:
        IMPORTS_PATH + "/{org}/{path_prefix}.gz"
    params:
        url = lambda wildcards, output: path2url[output[0]]
    shell:
        """
        wget -O {output} {params.url}
        """

rule map_gene_to_canonical_protein:
    input:
        pep = lambda wc: (
            urls.query("file_type == 'translated_cds'")
            .query("ShortName == @wc.org")
            .local
        ),
        gff = lambda wc: (
            urls.query("file_type == 'genomic_gff'")
            .query("ShortName == @wc.org")
            .local
        ),
        py = "workflow/scripts/prepare_ncbi_proteins/get_canonical_isoform.py",
    params:
        acc = lambda wildcards, output: org2acc[wildcards.org],
    output:
        "scratch/canonical_peptides/{org}.faa",
    shell:
        """
        python {input.py} \
            --organism {wildcards.org} \
            --accession {params.acc} \
            --gff {input.gff} \
            --protein-faa {input.pep} \
            --out {output}
        """

