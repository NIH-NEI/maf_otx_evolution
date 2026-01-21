import pandas as pd
import re

iqtree_file = "scratch/protein_tree_aln/species50/comparison/MAFL__full__clustalo__fullaln__fast.iqtree"
trees_file = "scratch/protein_tree_aln/species50/comparison/MAFL__full__clustalo__fullaln__fast.trees"

iqtree_file = snakemake.input["result"]
trees_file = iqtree_file.replace(".iqtree", ".trees")
zinfo_file = snakemake.input["zinfo"]
out_file = snakemake.output["table"]

print(iqtree_file)

tree_strings = {}
with open(trees_file, "r") as f:
    lines = f.readlines()

for s in lines:
    # Extract the content within brackets
    m = re.search(r"\[([^\]]+)\](.+)", s)
    bracket_content = m.group(1)
    newick_string = m.group(2)

    # Extract the tree number (after "tree" and before "lh=")
    tree_number = re.search(r"tree\s+(\d+)", bracket_content).group(1)

    tree_strings[tree_number] = newick_string

with open(iqtree_file, "r") as f:
    lines = f.readlines()
    lines = [l.strip() for l in lines]

print(lines)

START_STR = "------------------------------------------------------------------"
END_STR = "deltaL  : logL difference from the maximal logl in the set."

header = ""
data_rows = []
for i in range(1, len(lines)):
    if lines[i].startswith(START_STR):
        header = lines[i-1].split()
        for j in range(i+1, len(lines)):
            if lines[j] == END_STR:
                break
            if lines[j]:
                data_rows.append(lines[j].split())
        break

updated_header = []
for col in header:
    col = col.replace("-", "_")
    updated_header += (
        [col] if col in ["Tree", "logL", "deltaL"] else [col, col + "_sig"]
    )

# Create DataFrame
df = pd.DataFrame(data_rows, columns=updated_header)
print(df)

# Convert numeric columns
for col in df.columns[1:]:
    if col.endswith("_sig"):
        next
    elif col == "Tree":
        df[col] = df[col].astype(str)
    else:
        df[col] = pd.to_numeric(df[col])

df = (
    df.assign(tree_string = lambda tdf: tdf.Tree.apply(lambda x: tree_strings[x]))
)

df2 = pd.read_csv(zinfo_file, sep="\t")
df2["Tree"] = df2["Tree"].astype(str)
df["Tree"] = df["Tree"].astype(str)


df = (
    df.merge(df2, how = "left")
    .set_index(["Tree", "Constraint"])
    .drop(columns = ['tree_string', 'NewickFile', 'TreeString'])
)
print(df)

df.to_csv(out_file, sep = "\t")
