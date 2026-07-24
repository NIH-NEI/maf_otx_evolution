import pandas as pd
import matplotlib.pyplot as plt

# File paths
xlsx_file_path = 'scratch/aligned_single_maf/single_mafs_aligned.xlsx'

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
                    worksheet.write(row_num + 1 + 2, col_num + 1, value, format)  # row_num + 1 to account for the header

    tmp_df = df.reset_index()
    found_genes = tmp_df[tmp_df.pos ==  0]
#    print(found_genes)
    row_index = found_genes.index[0]
#    print(row_index)

    worksheet.set_row(row_index + 1 + 2, None, yellow_format)
    worksheet.freeze_panes(2, 1)

# Write DataFrames to an Excel file with each DataFrame in a separate sheet
with pd.ExcelWriter(xlsx_file_path, engine='xlsxwriter') as writer:
    for maf in "mafa mafb cmaf nrl maff mafg mafk".split():
        h5_file = f"scratch/aligned_single_maf/{maf}_alignment.h5"
        df = pd.read_hdf(h5_file, key = maf)
        df.to_excel(writer, sheet_name = maf, index=True)
        in_file = f"scratch/aligned_single_maf/{maf}_alignment.txt"
        with open(in_file) as f:
            text = f.readlines()
        #prefixes_left = text[0].rstrip().split(",")
        #prefixes_right = text[1].rstrip().split(",")
        #prefixes = prefixes_left[0:10] + prefixes_right[0:10]
        prefixes = text[0].rstrip().split(",")[0:20]

        # Generate colors
        colors = generate_colors(len(prefixes))

        apply_conditional_formatting(writer, maf, df, prefixes, colors)

print(f"DataFrames have been written to {xlsx_file_path}")
