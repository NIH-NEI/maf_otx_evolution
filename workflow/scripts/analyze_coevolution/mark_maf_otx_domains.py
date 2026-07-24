import pandas as pd
import numpy as np
from collections import defaultdict

# Define domain keywords mapping to regions
DOMAIN_KEYWORDS = {
    "disordered": ["mobidb"],
    "activation": ["PF08383", "IPR013592"],  # Maf transcription factor, N-terminal
    "basic_bzip": ["IPR004827", "IPR046347"],  # Basic leucine zipper and superfamily
    "ehr_basic_bzip": ["IPR004826"],  # Basic leucine zipper domain, Maf-type
    "ehr_basic": ["IPR008917", "IPR001356", "DNA-binding"],  # Skn-1-like and DNA-binding domains
    "maf_tf": ["IPR024874"],  # Transcription factor Maf family
    "homeobox": ["IPR009057", "IPR001356"],
    "otx1_tf": ["PF03529", "IPR013851"],
}
DOMAIN_KEYS = list(DOMAIN_KEYWORDS.keys())

def parse_interproscan_tsv(tsv_file, output_file):
    """
    Parse InterProScan TSV file and generate domain annotations.
    
    Args:
        tsv_file (str): Path to InterProScan TSV file
        output_file (str): Path to output TSV file
    """
    # Read InterProScan TSV file
    df = pd.read_csv(
        tsv_file,
        sep="\t",
        header=None,
        usecols=[0, 3, 4, 5, 6, 7, 11],
        names=["ProteinID", "analysis", "signature", "description", "startpos", "endpos", "details"]
    )

    def classify(interpro_id, description):
        """Classify domains based on InterPro ID or description keywords."""
        results = []
        for region, keywords in DOMAIN_KEYWORDS.items():
            if any(kw.lower() in str(description).lower() or kw in str(interpro_id) for kw in keywords):
                results.append(region)
        return results or "unknown"

    def merge_intervals(group):
        """Merge overlapping intervals in a group."""
        merged = []
        for _, row in group.iterrows():
            if not merged or row['startpos'] > merged[-1][1]:
                merged.append([row['startpos'], row['endpos']])
            else:
                merged[-1][1] = max(merged[-1][1], row['endpos'])
        return pd.DataFrame(merged, columns=['startpos', 'endpos'])

    def get_left(row):
        """Calculate left flanking region between ehr_basic_bzip and basic_bzip."""
        big, small = row['ehr_basic_bzip'], row['basic_bzip']
        if pd.isna(big) or pd.isna(small):
            return np.nan
        try:
            big = [int(x) for x in big.split(":")]
            small = [int(x) for x in small.split(":")]
            #if big[1] < small[1]:
            #    print("left", big, small)
            #    exit(-1)
            #    return np.nan
            return f"{big[0]}:{small[0]-1}" if big[0] < small[0] else np.nan
        except (ValueError, IndexError):
            return np.nan

    def get_right(row, col1, col2):
        """Calculate right flanking region between two columns."""
        big, small = row[col1], row[col2]
        if pd.isna(big) or pd.isna(small):
            return np.nan
        try:
            big = [int(x) for x in big.split(":")]
            small = [int(x) for x in small.split(":")]
            #if big[0] > small[0]:
            #    print("left", big, small)
            #    exit(-1)
            #    return np.nan
            return f"{small[1]+1}:{big[1]}" if big[1] > small[1] else np.nan
        except (ValueError, IndexError):
            return np.nan

    # Process DataFrame
    df = (
        df
        #.query("ProteinID in ['sheep__NRL__XP_042108879.1', 'human__NRL__NP_001341697.1']")
        .assign(description=lambda tdf: tdf.description + " " + tdf.details.fillna(""))
        .assign(domain=lambda tdf: tdf.apply(lambda row: classify(row.signature, row.description), axis=1))
        .explode("domain")
        .query("domain != 'unknown'")
    )
    print(df)
    print(df[["ProteinID", "startpos", "endpos", "domain"]])
    df = (
        df
        .query("domain in @DOMAIN_KEYS and domain != 'disordered'")
        [["ProteinID", "domain", "startpos", "endpos"]]
        .sort_values(["ProteinID", "startpos"])
        .groupby(['ProteinID', 'domain'])
        .apply(merge_intervals, include_groups=False)
        .reset_index()
        .sort_values(['ProteinID', 'domain', "startpos"])
        .assign(pos=lambda tdf: tdf["startpos"].astype(str) + ":" + tdf["endpos"].astype(str))
        .pivot_table(index="ProteinID", columns="domain", values="pos", aggfunc = "first")
        .fillna("")
        .reset_index()
        .assign(ehr=lambda tdf: tdf[["ehr_basic_bzip", "basic_bzip"]].apply(get_left, axis=1))
        .assign(basic=lambda tdf: tdf[["ehr_basic", "ehr"]].apply(get_right, axis=1, args=("ehr_basic", "ehr")))
        .assign(bzip=lambda tdf: tdf[["ehr_basic_bzip", "ehr_basic"]].apply(get_right, axis=1, args=("ehr_basic_bzip", "ehr_basic")))
    )

    # Save results to output file
    df.to_csv(output_file, sep="\t", index=False)
    return df

# Execute parsing
if __name__ == "__main__":
    interproscan_results_file = snakemake.input["ips"]
    maf_domains_file = snakemake.output[0]
    domain_annotations = parse_interproscan_tsv(interproscan_results_file, maf_domains_file)
