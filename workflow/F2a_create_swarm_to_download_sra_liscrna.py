import pandas as pd
import os

swarm_file = "workflow/F2a_download_sra_liscrna.swarm"

samples = (
    pd.read_csv("configs/LiSraRunTable.csv")
    .rename(columns = {"Sample Name": "GSM"})
    [["GSM", "Run", "AvgSpotLen", "Bases"]]
    .sort_values(["Bases"])
)
print(samples)

content = """#SWARM --threads-per-process 28
#SWARM --gb-per-process 240
#SWARM --time 4:00:00
#SWARM --sbatch '--mail-type=FAIL --gres=lscratch:400'
#SWARM --module sratoolkit,pigz
#SWARM --logdir /data/pals2/maf_zenodo/liscrna_swarmlogs
"""

with open(swarm_file, "w") as f:
    f.write(content)

    for idx, row in samples.iterrows():
        gsm = row["GSM"]
        srr = row["Run"]
        folder = f"/data/VisionDevo/LiRetinaScRNA/ByGSM/{gsm}/{srr}"
        os.makedirs(folder, exist_ok=True)
        cmds = [
            f"fasterq-dump {srr} --split-files --threads 28 --temp /lscratch/$SLURM_JOBID --outdir {folder} --progress",
            f"pigz -p 32 {folder}/*.fastq",
        ]
        cmd = "; ".join(cmds)
        f.write(cmd + "\n")

