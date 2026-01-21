import pandas as pd
from os.path import basename

IMPORTS_PATH = "imports/genomes_annotations"
URL_COMMON_PREFIX = "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF"
SPECIES_FILE = "configs/species_table.tsv"

suffixes = dict(
    translated_cds = "_translated_cds.faa.gz",
    untranslated_cds = "_cds_from_genomic.fna.gz",
    genomic_gff = "_genomic.gff.gz",
    genomic_gtf = "_genomic.gtf.gz",
    genomic_fna = "_genomic.fna.gz",
)

urls = (
    pd.read_csv(SPECIES_FILE, sep = "\t")
    .query("AnnotationSource == 'GCF'")
    .assign(acc=lambda tdf: tdf.apply(
        lambda row: (
            f"{row['AnnotationRelease']}_{row['AnnotationName']}"
        ),
        axis = 1,
    ))
    .assign(GCFString=lambda tdf: tdf.apply(
        lambda row: (
            f"{str(row['AnnotationRelease'])[4:7]}/"
            f"{str(row['AnnotationRelease'])[7:10]}/"
            f"{str(row['AnnotationRelease'])[10:13]}/"
            f"{row['acc']}/"
            f"{row['acc']}"
        ),
        axis=1,
      ))
    .assign(url_prefix = lambda tdf: URL_COMMON_PREFIX + "/" + tdf.GCFString)
    .assign(path_prefix = lambda tdf: IMPORTS_PATH + "/" + tdf.OrganismShortName)
    .merge(pd.DataFrame(suffixes.items(), columns = ["file_type", "suffix"]), how = "cross")
    .assign(remote = lambda tdf: tdf.url_prefix + tdf.suffix)
    .assign(local = lambda tdf: tdf.GCFString.apply(lambda x: basename(x)))
    .assign(local = lambda tdf: tdf.path_prefix + "/" + tdf.local + tdf.suffix)
)
print(urls)

path2url = urls.set_index("local")["remote"].to_dict()
org2acc = urls.set_index("OrganismShortName")["acc"].to_dict()


rule all:
    input:
        urls.local,
        expand("scratch/canonical_peptides/{org}.faa", org  = urls.OrganismShortName.unique()),
        expand("scratch/transcripts_annotation/{org}/{org}_cds_pep_gene_info.tsv", org = urls.OrganismShortName.unique()),


rule download_ncbi_data:
    """
    Use this rule separately to download all three files: cds, translated_cds and gff for an organism
    """
    output:
        IMPORTS_PATH + "/{org}/{path_prefix}.gz"
    params:
        url = lambda wildcards, output: path2url[output[0]]
    shell:
        """
        wget -O '{output}' '{params.url}'
        """

def parse_ncbi_fasta_headers(fasta_path, fasta_type="cds", index_name="cds_header"):
    """
    Parse a FASTA file with NCBI-style CDS/PEP headers and extract metadata.
    Returns a pandas DataFrame indexed by transcript ID.
    """
    import pandas as pd
    import gzip
    import re

    records = []

    # Regular expression patterns for key-value pairs in the header
    kv_pattern = re.compile(r'\[([^=]+)=([^]]+)\]')

    with gzip.open(fasta_path, "rt") as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                # Process the header line (without the '>')
                header = line[1:].strip()

                # Split into main ID and the rest
                parts = header.split(' ', 1)
                record_id = parts[0].rstrip()  # e.g., ENSANIT00000013014.1
                if len(parts) < 2:
                    print(f"Warning: No metadata for {record_id}")
                    continue
                annotation = parts[1]

                # Parse all [key=value] fields
                info_dict = {
                    'gene': '',
                    'protein_id': '',
                    'location': '',
                    'db_xref_GeneID': '',
                    'db_xref_CCDS': '',
                    'db_xref_Ensembl': ''
                }

                matches = kv_pattern.findall(annotation)
                for key, value in matches:
                    key = key.strip()
                    value = value.strip()
                    if key == 'gene':
                        info_dict['gene'] = value
                    elif key == 'protein_id':
                        info_dict['protein_id'] = value
                    elif key == 'location':
                        info_dict['location'] = value
                    elif key == 'db_xref':
                        # db_xref can have multiple comma-separated values
                        for ref in [x.strip() for x in value.split(',')]:
                            if ref.startswith('GeneID:'):
                                info_dict['db_xref_GeneID'] = ref[7:]
                            elif ref.startswith('CCDS:'):
                                info_dict['db_xref_CCDS'] = ref[5:]
                            elif ref.startswith('Ensembl:'):
                                info_dict['db_xref_Ensembl'] = ref[8:]

                records.append((record_id, info_dict))

    # Create DataFrame
    if not records:
        print("No records parsed.", file=sys.stderr)
        return pd.DataFrame()

    df = pd.DataFrame.from_records(
        [info for _, info in records],
        index=[tid for tid, _ in records]
    )
    df.index.name = index_name

    # Ensure columns are in desired order
    column_order = [
        'gene',
        'protein_id',
        'location',
        'db_xref_GeneID',
        'db_xref_CCDS',
        'db_xref_Ensembl'
    ]
    df = df[column_order]

    print(df)

    return df

rule extract_cds_info_from_fasta:
    input:
        cds = lambda wc: (
            urls.query("file_type == 'untranslated_cds'")
            .query("OrganismShortName == @wc.org")
            .local
        ),
    output:
        info = "scratch/transcripts_annotation/{org}/{org}_cds_info.tsv",
    run:
        df = parse_ncbi_fasta_headers(input.cds[0])
        df.to_csv(output.info, sep = "\t")


rule extract_pep_info_from_fasta:
    input:
        pep = lambda wc: (
            urls.query("file_type == 'translated_cds'")
            .query("OrganismShortName == @wc.org")
            .local
        ),
    output:
        info = "scratch/transcripts_annotation/{org}/{org}_pep_info.tsv",
    run:
        df = parse_ncbi_fasta_headers(input.pep[0], "pep", "pep_header")
        df.to_csv(output.info, sep = "\t")

rule merge_cds_pep_info:
    input:
        cds_info = "scratch/transcripts_annotation/{org}/{org}_cds_info.tsv",
        pep_info = "scratch/transcripts_annotation/{org}/{org}_pep_info.tsv",
    output:
        "scratch/transcripts_annotation/{org}/{org}_cds_pep_gene_info.tsv",
    run:
        import pandas as pd
        cds_info = (
            pd.read_csv(input.cds_info, sep = "\t")
            [["cds_header", "gene", "protein_id"]]
            .rename(columns = {
                "protein_id": "transcript"
            })
            .assign(common_header = lambda tdf: tdf.cds_header.str.replace("_cds_", "_haha_"))
        )
        print(cds_info)
        pep_info = (
            pd.read_csv(input.pep_info, sep = "\t")
            [["pep_header", "protein_id"]]
            .rename(columns = {
                "protein_id": "protein"
            })
            .assign(common_header = lambda tdf: tdf.pep_header.str.replace("_prot_", "_haha_"))
        )
        print(pep_info)

        assert(cds_info.shape[0] == pep_info.shape[0])

        df = (
            cds_info.merge(pep_info, how = "inner", on="common_header")
            [[
                "cds_header",
                "gene",
                "transcript",
                "protein",
                "pep_header",
            ]]
        )
        print(df)

        assert(df.shape[0] == pep_info.shape[0])

        df.to_csv(output[0], sep ="\t", index = False)


rule map_gene_to_canonical_protein:
    input:
        pep = lambda wc: (
            urls.query("file_type == 'translated_cds'")
            .query("OrganismShortName == @wc.org")
            .local
        ),
        gff = lambda wc: (
            urls.query("file_type == 'genomic_gff'")
            .query("OrganismShortName == @wc.org")
            .local
        ),
        py = "workflow/scripts/prepare_ncbi_proteins/get_canonical_isoform_from_ncbi.py",
    params:
        acc = lambda wildcards, output: org2acc[wildcards.org],
    output:
        faa = "scratch/canonical_peptides/{org}.faa",
        info = "scratch/canonical_peptides/{org}_info.tsv",
    shell:
        """
        python {input.py} \
            --organism {wildcards.org} \
            --accession {params.acc} \
            --gff {input.gff} \
            --protein-faa {input.pep} \
            --out {output.faa} \
            --info {output.info}
        """

