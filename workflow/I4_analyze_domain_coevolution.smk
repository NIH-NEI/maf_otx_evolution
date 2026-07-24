import sys
print(sys.executable)

import pandas as pd
from Bio import SeqIO

gene_group_subgroups_file = "configs/goi_group_subgroups.csv"

gene_group_subgroups = (
    pd.read_csv(gene_group_subgroups_file)
    .groupby("GeneGroup")
    .agg(list)
    ["GeneSubgroup"]
)
print(gene_group_subgroups)

tfs = ["OTX1", "OTX2", "CRX"]
targets = ["MAFA", "MAFB", "CMAF", "NRL"]
tf_domains = ["homeobox", "otx1_tf"]
target_domains = ["activation", "basic_bzip", "ehr_basic", "ehr_basic_bzip", "maf_tf", "ehr", "basic", "bzip"]


mammals =  (
    pd.read_csv("configs/species_table.tsv", sep = "\t")
    .query("OrganismColor == 'Mammals'")
    .OrganismShortName.tolist()
)
print(mammals)

#tfs = ["OTX1"]
#tf_domains = ["homeobox"]
#targets = ["MAFA"]
#target_domains = ["activation"]

rule all:
    input:
        expand(
            "scratch/coevolution_analysis/coevol_coeff__{tf}__{tfdomain}__{target}__{tgtdomain}.tsv",
            tf = tfs,
            target = targets,
            tfdomain = tf_domains,
            tgtdomain = target_domains,
        ),
        expand(
            "scratch/coevolution_analysis/coevol_summary__OTX__{tfdomain}__MAFL__{tgtdomain}__{spgrp}.tsv",
            tfdomain = tf_domains,
            tgtdomain = target_domains,
            spgrp = ["mammal", "nonmammal", "all"],
        ),
        "scratch/coevolution_analysis/tiled/OTX_otx1_tf_MAFL_tiled.pdf",
        "scratch/coevolution_analysis/tiled/OTX_homeobox_MAFL_tiled.pdf",

rule get_common_vertebrates_for_coevolution_analysis:
    input:
        tf_info = "scratch/final_list/{tf}.tsv",
        target_info = "scratch/final_list/{target}.tsv",
    output:
        "scratch/coevolution_analysis/common_species__{tf}__{target}.tsv"
    run:
        import pandas as pd
        tf_info = (
            pd.read_csv(input.tf_info, sep = "\t")
            [["ProteinID", "OrganismShortName", "GeneGroup", "OrganismColor", "ProteinSeq"]]
            .query("GeneGroup != 'Outgroup'")
        )
        target_info = (
            pd.read_csv(input.target_info, sep = "\t")
            [["ProteinID", "OrganismShortName", "GeneGroup", "OrganismColor", "ProteinSeq"]]
            .query("GeneGroup != 'Outgroup'")
        )

        common_species = list(set(tf_info.OrganismShortName) & set(target_info.OrganismShortName))
        print(common_species)
        print(len(common_species))

        tf_info = (
            tf_info.query("OrganismShortName in @common_species")
        )
        print(tf_info)
        print(len(tf_info))

        target_info = (
            target_info.query("OrganismShortName in @common_species")
        )
        print(target_info)
        print(len(target_info))


        df = pd.concat([tf_info, target_info])

        df[["ProteinID", "OrganismShortName", "GeneGroup", "ProteinSeq"]].to_csv(output[0], sep = "\t", index = False)

rule augment_domains:
    input:
        grand = "scratch/grand_list/grand_list.tsv",
        domains = "scratch/grand_list/grand_list_domains.tsv",
        poi = "scratch/coevolution_analysis/common_species__{tf}__{target}.tsv",
    output:
        df = "scratch/coevolution_analysis/coevol_coeff__{tf}__{tfdomain}__{target}__{tgtdomain}.tsv",
        scatter = "scratch/coevolution_analysis/coevol_coeff_scatter__{tf}__{tfdomain}__{target}__{tgtdomain}.pdf",
        heatmap = "scratch/coevolution_analysis/coevol_coeff_heatmap__{tf}__{tfdomain}__{target}__{tgtdomain}.pdf",
    run:
        import pandas as pd
        df = (
            pd.read_csv(input.poi, sep = "\t")
            .merge(pd.read_csv(input.domains, sep = "\t"), how = "left")
        )
        print(df)

        tf_df = (
            df.query("GeneGroup == @wildcards.tf")
            .assign(startpos = lambda tdf: tdf[wildcards.tfdomain].str.split(":").str[0])
            .assign(endpos = lambda tdf: tdf[wildcards.tfdomain].str.split(":").str[1])
        )
        target_df = (
            df.query("GeneGroup == @wildcards.target")
            .assign(startpos = lambda tdf: tdf[wildcards.tgtdomain].str.split(":").str[0])
            .assign(endpos = lambda tdf: tdf[wildcards.tgtdomain].str.split(":").str[1])
        )

        df = (
            pd.concat([tf_df, target_df])
            [["ProteinID", "ProteinSeq", "GeneGroup", "startpos", "endpos"]]
            .dropna()
            .assign(newid = lambda tdf: tdf.apply(lambda x: f"{x.ProteinID}__{x.startpos}__{x.endpos}", axis=1))
        )
        print(df)

        df["startpos"] = df["startpos"].astype(int)
        df["endpos"] = df["endpos"].astype(int)

        orig_df = (
            df
            .assign(subsequence = lambda tdf: tdf.apply(
                lambda r: r["ProteinSeq"][r["startpos"]-1 : r["endpos"]],
                axis = 1
            ))
            .set_index("newid", drop = False)
        )
        print(orig_df)
        print(orig_df.columns)

        from Bio import pairwise2
        import itertools

        otx_name = wildcards.tf
        mafl_name = wildcards.target

        # Function to compute pairwise sequence similarity (global alignment score)
        def compute_similarity(seq1, seq2):
            print(seq1)
            print(seq2)
            print(type(seq1), getattr(seq1, "shape", None))
            print(type(seq2), getattr(seq2, "shape", None))
            alignments = pairwise2.align.globalxx(seq1, seq2, one_alignment_only=True)
            print(alignments)
            alignment = alignments[0]
            print(alignment)
            score = alignment.score
            print(score)
            alignment_length = len(alignment.seqA.replace("-", ""))  # or len(alignment.seqA) to include gaps
            print(alignment_length)
            normalized_score = score / alignment_length if alignment_length > 0 else 0
            print(normalized_score)
            return normalized_score

        df = orig_df.query("GeneGroup == @wildcards.tf")
        print(df)

        pairs = list(itertools.combinations(df.index, 2))
        data = []
        for i, j in pairs:
            id1 = df.at[i, 'newid']
            id2 = df.at[j, 'newid']
            seq1 = df.at[i, 'subsequence']
            seq2 = df.at[j, 'subsequence']
            print(id1, ":", seq1)
            print(id2, ":", seq2)
            score = compute_similarity(seq1, seq2)
            data.append({'id1': id1, 'id2': id2, 'similarity': score})

        similarity_df = (
            pd.DataFrame(data)
            .assign(org1 = lambda tdf: tdf.id1.str.split("__").str[0])
            .assign(org2 = lambda tdf: tdf.id2.str.split("__").str[0])
            .assign(org_pair = lambda tdf: tdf.apply(
                lambda x: f"{x.org1}__{x.org2}" if x.org1 < x.org2 else f"{x.org2}__{x.org1}",
                axis = 1
            ))
            .query("org1 != org2")
            .groupby(["org_pair"])
            .agg({"similarity": "max"})
            .reset_index()
        )
        print(similarity_df)

        otx_similarity = similarity_df.rename(columns = {"similarity": "otx_similarity"})
        print(otx_similarity)

        df = orig_df.query("GeneGroup == @mafl_name")
        print(df)

        pairs = list(itertools.combinations(df.index, 2))
        data = []
        for i, j in pairs:
            id1 = df.at[i, 'newid']
            id2 = df.at[j, 'newid']
            seq1 = df.at[i, 'subsequence']
            seq2 = df.at[j, 'subsequence']
            print(id1, ":", seq1)
            print(id2, ":", seq2)
            score = compute_similarity(seq1, seq2)
            data.append({'id1': id1, 'id2': id2, 'similarity': score})

        similarity_df = (
            pd.DataFrame(data)
            .assign(org1 = lambda tdf: tdf.id1.str.split("__").str[0])
            .assign(org2 = lambda tdf: tdf.id2.str.split("__").str[0])
            .assign(org_pair = lambda tdf: tdf.apply(
                lambda x: f"{x.org1}__{x.org2}" if x.org1 < x.org2 else f"{x.org2}__{x.org1}",
                axis = 1
            ))
            .query("org1 != org2")
            .groupby(["org_pair"])
            .agg({"similarity": "max"})
            .reset_index()
        )
        print(similarity_df)

        mafl_similarity = similarity_df.rename(columns = {"similarity": "mafl_similarity"})
        print(mafl_similarity)

        
        similarity = (
            otx_similarity.merge(mafl_similarity, how = "inner")
        )
        print(similarity)
        similarity.to_csv(output.df, sep = "\t", index = False)

        final_coevol_score = similarity["otx_similarity"].corr(similarity["mafl_similarity"])
        print(final_coevol_score)

        import pandas as pd
        import seaborn as sns
        import matplotlib
        matplotlib.use("Agg")  # non-interactive backend
        import matplotlib.pyplot as plt

        from scipy.stats import pearsonr
        df = similarity

# Example: Load your DataFrame (replace with your actual DataFrame)
# df = pd.read_csv("your_data.csv")  # or however you're getting the data

# Assuming your DataFrame is already available as `df`

# Calculate Pearson correlation
        r, p_value = pearsonr(df['otx_similarity'], df['mafl_similarity'])

# Create scatter plot with regression line
        plt.figure(figsize=(8, 6))
        sns.regplot(data=df, x='otx_similarity', y='mafl_similarity', ci=None, scatter_kws={'s': 10, 'alpha': 0.7})

# Annotate correlation
        plt.title(f'Pearson r = {r:.3f}, p = {p_value:.3e}', fontsize=14)
        plt.xlabel(f'{otx_name} Similarity')
        plt.ylabel(f'{mafl_name} Similarity')
        plt.grid(True)
        plt.tight_layout()
        #plt.show()
        plt.savefig(output.scatter)


        df = similarity

# Step 1: Split the org_pair into org1 and org2
        df[['org1', 'org2']] = df['org_pair'].str.split('__', expand=True)

        df_reversed = df.copy()
        df_reversed['org1'], df_reversed['org2'] = df['org2'], df['org1']

# Step 3: Concatenate original and reversed to fill both triangles
        df = pd.concat([df, df_reversed], ignore_index=True)
        print(df)


        all_orgs = list(set(list(df.org1) + list(df.org2)))
        diag_df = pd.DataFrame({
            'org1': all_orgs,
            'org2': all_orgs,
            'otx_similarity': 1.0,
            'mafl_similarity': 1.0
        })
        df = pd.concat([df, diag_df], ignore_index=True)

        org_order = (
            pd.DataFrame(dict(ShortName = all_orgs))
            .merge(pd.read_csv("configs/species50_table.tsv", sep = "\t")[["ShortName", "Order50"]])
            .sort_values("Order50")
            ["ShortName"].tolist()
        )
        print(org_order)

# Step 2: Pivot the data to get similarity matrices
        otx_matrix = df.pivot(index='org1', columns='org2', values='otx_similarity').fillna(0)
        mafa_matrix = df.pivot(index='org1', columns='org2', values='mafl_similarity').fillna(0)

# Reorder rows and columns
        otx_matrix = otx_matrix.reindex(index=org_order, columns=org_order)
        mafa_matrix = mafa_matrix.reindex(index=org_order, columns=org_order)

# Step 3: Create and save the heatmaps
        fig, axs = plt.subplots(1, 2, figsize=(24, 10))

        otx_domain = wildcards.tfdomain
        maf_domain = wildcards.tgtdomain

        sns.heatmap(otx_matrix, annot=False, cmap="YlGnBu", ax=axs[0], cbar_kws={'label': f'{otx_name} {otx_domain} Domain Similarity'})
        axs[0].set_title(f"{otx_name} {otx_domain} Domain Similarity Heatmap")
        axs[0].set_xlabel("org2")
        axs[0].set_ylabel("org1")

        sns.heatmap(mafa_matrix, annot=False, cmap="YlOrBr", ax=axs[1], cbar_kws={'label': f'{mafl_name} {maf_domain} Domain Similarity'})
        axs[1].set_title(f"{mafl_name} {maf_domain} Domain Similarity Heatmap")
        axs[1].set_xlabel("org2")
        axs[1].set_ylabel("org1")

        plt.tight_layout()
        #plt.savefig("similarity_heatmaps.png", dpi=300)  # Save to PNG file
        plt.savefig(output.heatmap)
        #plt.show()


def get_coevol_inputs(wc):
    mydict = {
        f"{tf}_{target}": f"scratch/coevolution_analysis/coevol_coeff__{tf}__{{tfdomain}}__{target}__{{tgtdomain}}.tsv"
        for tf in tfs
        for target in targets
    }
    print(mydict)
    return mydict

rule create_consolidated_correlation_table:
    input:
        unpack(get_coevol_inputs),
    output:
        df = "scratch/coevolution_analysis/coevol_summary__OTX__{tfdomain}__MAFL__{tgtdomain}__{spgrp}.tsv",
        heatmap = "scratch/coevolution_analysis/coevol_summary__OTX__{tfdomain}__MAFL__{tgtdomain}__{spgrp}.pdf",
    run:
        import pandas as pd
        import seaborn as sns
        import matplotlib
        matplotlib.use('Agg')  # Non-GUI backend
        import matplotlib.pyplot as plt
        import numpy as np
        from scipy.stats import pearsonr
        values = []
        for cfg, infile in input.items():
            print(cfg, infile)
            df = pd.read_csv(infile, sep = "\t")
            print(df)
            df[['org1', 'org2']] = df['org_pair'].str.split('__', expand=True)
            if wildcards.spgrp == "mammal":
                df = df.query("org1 in @mammals").query("org2 in @mammals") 
            elif wildcards.spgrp == "nonmammal":
                df = df.query("org1 not in @mammals").query("org2 not in @mammals") 
            print(df)
            r, p_value = pearsonr(df['otx_similarity'], df['mafl_similarity'])
            values.append({"config": cfg, "pearsonr": r, "pvalue": p_value})
        df = (
            pd.DataFrame(values)
            .fillna({"pearsonr":0, "pvalue": 1})
        )
        print(df)
        df.to_csv(output.df, sep = "\t", index = False)

# Split config into two columns
        df[['row', 'col']] = df['config'].str.split('_', expand=True)

# Pivot for heatmap
        pval_matrix = df.pivot(index='row', columns='col', values='pvalue')
        corr_matrix = df.pivot(index='row', columns='col', values='pearsonr')

# Log10-transform p-values for better color contrast
        log_pval_matrix = -np.log10(pval_matrix)

        global_min = 0
        global_max = 7

# Plot heatmap
        plt.figure(figsize=(6, 4.5))
        sns.heatmap(log_pval_matrix, annot=corr_matrix.round(2), fmt='', cmap='viridis',
                    cbar_kws={'label': '-log10(p-value)'}, linewidths=0.5, linecolor='gray', vmin=global_min,
                    vmax=global_max,
                    annot_kws={"size": 20})

        otx_domain = wildcards.tfdomain
        maf_domain = wildcards.tgtdomain
        organisms = wildcards.spgrp

        plt.title(f"PCC {organisms} OTX {otx_domain} : MAF {maf_domain}")
        plt.xlabel("MAF Variant")
        plt.ylabel("OTX Variant")
        plt.tight_layout()
        #plt.show()
        plt.savefig(output.heatmap)

def get_tile_inputs(wc):
    # Define the MAFL combinations in order
    MAFL_DOMAINS = [
        "activation",
        "ehr_basic_bzip",
        "ehr",
        "basic",
        "bzip",
        "ehr_basic",
        "basic_bzip",
        "maf_tf",
    ]
    # Define the species variants
    SPECIES = ["all", "mammal", "nonmammal"]

    def get_pdf_path(otx_type, mafl_combo, species):
        """Generate the PDF file path."""
        return (
            f"scratch/coevolution_analysis/"
            f"coevol_summary__OTX__{otx_type}__MAFL__{mafl_combo}__{species}.pdf"
        )

    infiles = []
    for row_idx, mafl_dom in enumerate(MAFL_DOMAINS):
        for col_idx, species in enumerate(SPECIES):
            infiles.append(get_pdf_path(wc.otxdom, mafl_dom, species))

    return infiles


rule tile_pdfs:
    input:
        get_tile_inputs,
    output:
        tex = "scratch/coevolution_analysis/tiled/OTX_{otxdom}_MAFL_tiled.tex",
        pdf = "scratch/coevolution_analysis/tiled/OTX_{otxdom}_MAFL_tiled.pdf",
    params:
        BASE_DIR = "scratch/coevolution_analysis",
    script:
        "scripts/analyze_coevolution/tile_pdfs_using_latex.py"
