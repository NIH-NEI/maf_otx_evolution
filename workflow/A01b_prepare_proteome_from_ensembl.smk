import pandas as pd
from os.path import basename
import gzip

IMPORTS_PATH = "imports/genomes_annotations"
URL_COMMON_PREFIX = "https://ftp.ensembl.org/pub/release-115"
SPECIES_FILE = "configs/species_table.tsv"

suffixes = {
    "pep" : ".pep.all.fa.gz",
    "cds" : ".cds.all.fa.gz",
    "gff3" : ".115.gff3.gz",
}

def get_url_part(file_type, org_id, acc, suffix):
    if file_type == "gff3":
        return (
            URL_COMMON_PREFIX +
            "/gff3/" +
            org_id.lower() + "/" +
            acc + suffix
        )
    return (
            URL_COMMON_PREFIX +
            "/fasta/" +
            org_id.lower() + "/" +
            file_type + "/" +
            acc + suffix
    )

urls = (
    pd.read_csv(SPECIES_FILE, sep = "\t")
    .query("AnnotationSource == 'Ensembl'")
    .assign(acc=lambda tdf: tdf.apply(
        lambda row: (
            f"{row['OrganismID']}.{row['AnnotationName']}"
        ),
        axis = 1,
    ))
    .assign(GCFString=lambda tdf: tdf.apply(
        lambda row: (
            f"{str(row['OrganismID']).lower()}/"
            f"{str(row['OrganismID']).lower()}."
            f"{str(row['OrganismID'])}"
        ),
        axis=1,
      ))
    .assign(path_prefix = lambda tdf: IMPORTS_PATH + "/" + tdf.OrganismShortName)
    .merge(pd.DataFrame(suffixes.items(), columns = ["file_type", "suffix"]), how = "cross")
    .assign(remote = lambda tdf: tdf.apply(
        lambda row: get_url_part(
            row["file_type"],
            row["OrganismID"],
            row["acc"],
            row["suffix"],
        ),
        axis = 1
    ))
    .assign(local = lambda tdf: tdf.remote.apply(lambda x: basename(x)))
    .assign(local = lambda tdf: tdf.path_prefix + "/" + tdf.local)
)
print(urls)
print(urls.remote.to_list())
print(urls.local.to_list())

path2url = urls.set_index("local")["remote"].to_dict()
org2acc = urls.set_index("OrganismShortName")["acc"].to_dict()


rule all:
    input:
        urls.local,
        expand("scratch/canonical_peptides/{org}.faa", org  = urls.OrganismShortName.unique()),
        expand("scratch/transcripts_annotation/{org}/{org}_cds_pep_gene_info.tsv", org = urls.OrganismShortName.unique()),


rule download_ensembl_data:
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

def parse_ensembl_fasta_headers(fasta_path, fasta_type="cds", index_name="cds_header"):
    """
    Parse a FASTA file with Ensembl-style CDS/PEP headers and extract metadata.
    Returns a pandas DataFrame indexed by transcript ID.
    """
    import pandas as pd
    import re

    records = []

    # Regular expression patterns for key-value pairs in the header
    kv_pattern = re.compile(r'(\w+):([^:\s]+(?:[^\s]*[^\s:])?)')

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
                    print(f"Warning: No metadata for {record_id}", file=sys.stderr)
                    continue
                metadata_str = parts[1]

                # Extract known fixed fields first
                metadata_parts = metadata_str.split(' ', 4 if fasta_type == "cds" else 5)  # Split carefully to preserve description
                info_dict = {
                    'primary_assembly': '',
                    'gene': '',
                    'transcript': '',
                    'gene_biotype': '',
                    'transcript_biotype': '',
                    'description': ''
                }

                # First few fields are positional: cds primary_assembly:... gene:... etc.
                idx = 0
                if metadata_parts[idx] == fasta_type:
                    idx += 1

                if idx < len(metadata_parts) and metadata_parts[idx].startswith('primary_assembly:'):
                    info_dict['primary_assembly'] = metadata_parts[idx][17:]  # remove prefix
                    idx += 1

                # Now parse remaining key:value pairs
                remaining = ' '.join(metadata_parts[idx:])
                matches = kv_pattern.findall(remaining)

                description = ''
                for key, value in matches:
                    if key == 'gene':
                        info_dict['gene'] = value
                    elif key == 'gene_biotype':
                        info_dict['gene_biotype'] = value
                    elif key == 'transcript':
                        info_dict['transcript'] = value
                    elif key == 'transcript_biotype':
                        info_dict['transcript_biotype'] = value
                    elif key == 'description':
                        description = value
                    # Ignore others if needed

                # If there's leftover text after parsing known fields, it might be description
                # This handles cases where description has spaces and comes last
                if 'description:' in remaining:
                    desc_part = remaining.split('description:', 1)[1].strip()
                    info_dict['description'] = desc_part
                else:
                    info_dict['description'] = description

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
    column_order = ['primary_assembly', 'gene', 'gene_biotype', 'transcript', 'transcript_biotype', 'description']
    df = df[column_order]

    print(df)

    return df

rule extract_cds_info_from_fasta:
    input:
        cds = lambda wc: (
            urls.query("file_type == 'cds'")
            .query("OrganismShortName == @wc.org")
            .local
        ),
    output:
        info = "scratch/transcripts_annotation/{org}/{org}_cds_info.tsv",
    run:
        df = parse_ensembl_fasta_headers(input.cds[0])
        df.to_csv(output.info, sep = "\t")


rule extract_pep_info_from_fasta:
    input:
        pep = lambda wc: (
            urls.query("file_type == 'pep'")
            .query("OrganismShortName == @wc.org")
            .local
        ),
    output:
        info = "scratch/transcripts_annotation/{org}/{org}_pep_info.tsv",
    run:
        df = parse_ensembl_fasta_headers(input.pep[0], "pep", "pep_header")
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
            [["cds_header", "gene"]]
            .assign(transcript = lambda tdf: tdf.cds_header)
            # remove version from gene
            .assign(gene = lambda tdf: tdf.gene.apply(lambda s: s.rsplit(".", 1)[0] if s.startswith("ENS") else s))
        )
        pep_info = (
            pd.read_csv(input.pep_info, sep = "\t")
            [["transcript", "pep_header"]]
            .assign(protein = lambda tdf: tdf.pep_header)
        )
        assert(cds_info.shape[0] == pep_info.shape[0])

        df = (
            cds_info.merge(pep_info, how = "inner")
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
            urls.query("file_type == 'pep'")
            .query("OrganismShortName == @wc.org")
            .local
        ),
        gff = lambda wc: (
            urls.query("file_type == 'gff3'")
            .query("OrganismShortName == @wc.org")
            .local
        ),
        py = "workflow/scripts/prepare_ncbi_proteins/get_canonical_isoform_from_ensembl.py",
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

