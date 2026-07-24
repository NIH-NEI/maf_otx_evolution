import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from utils import *
import sys

# File paths
xlsx_file_path = f"scratch/aligned_maf_groups/maf_group_aligned.xlsx"

# Main part of the script
if len(sys.argv) < 2:
    sys.exit("Error: specify maf-group")

maf_group = sys.argv[1]
if (maf_group == "large"):
    mafs = "mafa mafb cmaf nrl".split()
else:
    mafs = "maff mafg mafk".split()
xlsx_file_path = f"scratch/aligned_maf_groups/{maf_group}_aligned.xlsx"

# Generate a list of colors using matplotlib
def generate_colors(num_colors):
    cmap = plt.get_cmap('tab20')  # You can choose any colormap
    return [cmap(i) for i in range(num_colors)]

# Convert RGBA colors to Excel-compatible hex format
def rgba_to_hex(rgba):
    return '#{:02x}{:02x}{:02x}'.format(int(rgba[0]*255), int(rgba[1]*255), int(rgba[2]*255))

# Function to apply conditional formatting based on prefixes
def apply_conditional_formatting(writer, sheet_name, df, prefixes, colors):
    workbook = writer.book
    worksheet = writer.sheets[sheet_name]

    # Map each prefix to a specific format
    prefix_format_map = {prefix: workbook.add_format({'bg_color': rgba_to_hex(color)}) for prefix, color in zip(prefixes, colors)}
    yellow_format = workbook.add_format({'bg_color': 'yellow'})

    # Apply conditional formatting for each prefix
    for prefix, format in prefix_format_map.items():
        for col_num, col_name in enumerate(df.columns):
            for row_num, value in enumerate(df[col_name]):
                if isinstance(value, str) and value.lower().startswith(prefix.lower()):
                    worksheet.write(row_num + 1 + 3, col_num + 1, value, format)  # row_num + 1 to account for the header

    tmp_df = df.reset_index()
    found_genes = tmp_df[tmp_df.pos ==  0]
#    print(found_genes)
    row_index = found_genes.index[0]
#    print(row_index)

    worksheet.set_row(row_index + 1 + 3, None, yellow_format)
    worksheet.freeze_panes(3, 1)

def get_combined_df(mafs, prefixes):
    df = pd.concat([
        pd.read_hdf(f"scratch/aligned_single_maf/{maf}_alignment.h5", key = maf)
        for maf in mafs
    ], axis = 1, keys = mafs)
    print(df)
    rows_keep = [False for x in range(df.shape[0])]
    print(rows_keep)
    for col_num, col_name in enumerate(df.columns):
        for row_num, value in enumerate(df[col_name]):
            for prefix in mafs + prefixes:
                if isinstance(value, str) and value.upper().startswith(prefix.upper()):
                    rows_keep[row_num] = True
    print(rows_keep)
    df = df.loc[rows_keep,:]
    print(df)
    return df
 
in_file = f"scratch/aligned_maf_groups/{maf_group}_prefixes.txt"
with open(in_file) as f:
    text = f.readlines()
prefixes_left = text[0].rstrip().split(",")
#prefixes_right = text[1].rstrip().split(",")
prefixes = prefixes_left #+ prefixes_right

print(prefixes)

in_file = f"scratch/aligned_maf_groups/{maf_group}_alignment.txt"
genes = {}
pos = {}
with open(in_file) as f:
    text = f.readlines()
    for line in text:
        maf = line.rstrip().split("|")[0].split(",")[1]
        if maf in genes.keys():
            genes[maf] += line.rstrip().split("|")[1].split(",")
        else:
            genes[maf] = line.rstrip().split("|")[1].split(",")
        pos[maf] = [x.split(":")[1] for x in genes[maf]]
        genes[maf] = [x.split(":")[0] for x in genes[maf]]
print(genes)
print(pos)

df = (
    pd.concat([
        pd.DataFrame(dict(maf = maf, gene = candidates, proximity = pos[maf]))
        for maf, candidates in genes.items()
    ])
#    .explode("gene")
#    .drop_duplicates()
    .assign(prefix = lambda tdf: tdf.gene.apply(lambda x: remove_trailing_numbers(x)))
    .groupby(["prefix", "maf"]).agg(dict(gene = "|".join, proximity = min))
    .unstack("maf", fill_value = "") #["gene"]
)

print(df)

df["score"] = (df["gene"] != "").sum(axis=1)
df["pos"] = np.nanmin(df["proximity"].applymap(lambda x: int(x) if x != '' else 10000).values, axis = 1)

print(df)

df = df.sort_values(["score", "pos"], ascending = [False, True])
print(df)

prefixes = df.index.to_list()[0:20]
#prefixes = [
#    "CPNE",
#    "RNF",
#    "TOX",
#    "SLC12A",
#    "NDRG",
#    "JPH",
#    "ADCY",
#    "CBLN",
#    "CHD",
#]
#prefixes = [
#    "KCTD",
#    "SYNGR",
#    "SSTR",
#    "RAC",
#    "CDK",
#    "GGA",
#    "SOX",
#    "SOCS",
#    "CARD",
#    "CYTH",
#    "RBFOX",
#    "ANKRD",
#    "SHISA",
#    "MYH",
#]

# Write DataFrames to an Excel file with each DataFrame in a separate sheet
with pd.ExcelWriter(xlsx_file_path, engine='xlsxwriter') as writer:
    print(prefixes)
    df.to_excel(writer, sheet_name = "summary", index=True)
    detailed =  get_combined_df(mafs, prefixes)
    print(detailed)
    detailed = detailed.loc[:,detailed.columns.get_level_values(2).isin({"gene"})]
    detailed.to_excel(writer, sheet_name = "detailed", index=True)

    # Generate colors
    colors = generate_colors(len(prefixes))

    apply_conditional_formatting(writer, "detailed", detailed, prefixes, colors)

print(f"DataFrames have been written to {xlsx_file_path}")
