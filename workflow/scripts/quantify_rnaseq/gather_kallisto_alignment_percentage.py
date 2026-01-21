import glob
import pandas as pd
import re

files = [x for x in snakemake.input]
samples = [x.split("/")[2] for x in files]

print(files)
print(samples)

pat = r'\[quant\] processed ([0-9,]+) reads, ([0-9,]+) reads pseudoaligned'

def get_counts(fpath):
    with open(fpath) as f:
        for l in f:
            if res := re.search(pat, l):
              return pd.Series({
                  "total": int(res.group(1).replace(',', '')),
                  "aligned": int(res.group(2).replace(',', '')),
              })

df = pd.DataFrame({
    "sample": samples,
    "file": files
}).set_index("sample")

df = (
    pd.concat([df.file.apply(get_counts)])
    .assign(pct_align = lambda tdf: tdf.aligned*100/tdf.total)
    .sort_values("pct_align")
)
print(df)

df.to_csv(snakemake.output[0], sep="\t")
