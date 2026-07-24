import pandas as pd
from utils import *

in_file = f"scratch/aligned_maf_groups/large_prefixes.txt"
with open(in_file) as f:
    text = f.readlines()
prefixes_left = text[0].rstrip().split(",")
prefixes_right = text[1].rstrip().split(",")
prefixes = prefixes_left + prefixes_right

print(prefixes)

in_file = f"scratch/aligned_maf_groups/large_alignment.txt"
genes = {}
with open(in_file) as f:
    text = f.readlines()
    for line in text:
        maf = line.rstrip().split(":")[0].split(",")[1]
        if maf in genes.keys():
            genes[maf] += line.rstrip().split(":")[1].split(",")
        else:
            genes[maf] = line.rstrip().split(":")[1].split(",")
print(genes)

df = (
    pd.concat([
        pd.DataFrame(dict(maf = maf, gene = candidates))
        for maf, candidates in genes.items()
    ])
    .explode("gene")
    .assign(prefix = lambda tdf: tdf.gene.apply(lambda x: remove_trailing_numbers(x)))
    .groupby(["prefix", "maf"]).agg(dict(gene = "|".join))
    .unstack("maf", fill_value = "")["gene"]
)

print(df)

df["score"] = (df != "").sum(axis=1)

print(df)

df = df.sort_values("score", ascending = False)

print(df)

