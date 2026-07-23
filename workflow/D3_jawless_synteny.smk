import sys
print(sys.executable)
import pandas as pd
import gffutils
from pathlib import Path

annotation = dict(
    amphioxus = "imports/genomes_annotations/amphioxus/GCF_035083965.1_klBraLanc5.hap2_genomic.gtf.gz",
    seaLamprey = "imports/genomes_annotations/seaLamprey/GCF_048934315.1_UKy_Petmar_22M1.pri1.0_genomic.gtf.gz",
    #fareasternBrookLamprey =
    #europeanBrookLamprey =
    #europeanRiverLamprey =
    brownHagfish = "imports/genomes_annotations/brownHagfish/Pata_MYAQb_rnm_ah2p.gtf",
    atlanticHagfish = "imports/genomes_annotations/atlanticHagfish/GCF_040869285.1_UKY_Mglu_1.0_genomic.gtf.gz",
    #inshoreHagfish =
)
print(annotation)

nc2chr = {
    "NC_089723.1": "chr2",
    "NC_089728.1": "chr7",
    "NC_089729.1": "chr8",
    "NC_090422.1": "chr1",
    "NC_090424.1": "chr3",
    "NC_089739.1": "chr18",
    "NC_090423.1": "chr2",
    "NC_090426.1": "chr5",
    "NC_090427.1": "chr6",
    "NC_090429.1": "chr8",
    "NC_090430.1": "chr9",
    "NC_090431.1": "chr10",
    "NC_090432.1": "chr11",
    "NC_090434.1": "chr13",
    "NC_133669.1": "chr3",
    "NC_133678.1": "chr12",
    "NC_133681.1": "chr15",
    "NC_133686.1": "chr20",
    "NC_133689.1": "chr23",
    "NC_133691.1": "chr25",
    "NC_133692.1": "chr26",
    "NC_133703.1": "chr37",
    "NC_133709.1": "chr43",
    "NC_133711.1": "chr45",
    "NC_133720.1": "chr54",
    "NC_133724.1": "chr58",
    "NC_133729.1": "chr63",
    "NC_133735.1": "chr69",
    "NC_133738.1": "chr72",
    "NW_027145746.1": "scaf117",
    "NW_027146371.1": "scaf1738",
    "NW_027423194.1": "scaf222",
    "NW_027423649.1": "scaf633",
    "NC_089725.1": "chr4",
    "NC_089726.1": "chr5",
    "NC_089738.1": "chr17",
    "NC_090425.1": "chr4",
    "NC_090433.1": "chr12",
    "NC_133672.1": "chr6",
    "NC_133683.1": "chr17",
    "NC_133695.1": "chr29",
    "NC_133697.1": "chr31",
    "NC_133699.1": "chr33",
    "NC_133708.1": "chr42",
    "NC_133717.1": "chr51",
    "NC_133719.1": "chr53",
}

rule all:
    input:
        expand("scratch/extracted_annotation/{org}.tsv", org = annotation),
        "scratch/cyclostomes_synteny/largemaf_neighbors_loci.tsv",
        "scratch/cyclostomes_synteny/smallmaf_neighbors_loci.tsv",
        "scratch/cyclostomes_synteny/otx_neighbors_loci.tsv",

rule extract_annotation:
    input:
        lambda wc: annotation[wc.org]
    output:
        "scratch/extracted_annotation/{org}.tsv",
    shell:
        """
        module load agat
        agat_convert_sp_gff2tsv.pl -f {input} -o {output}
        """


rule get_chromosome_locations:
    input:
        oglist = "configs/{grp}_neighbors_in_diverse_orthofinder.txt",
        ogfile = "scratch/orthofinder_diverse/results/Results_Feb11/Orthogroups/Orthogroups.tsv",
    output:
        "scratch/cyclostomes_synteny/{grp}_neighbors_loci.tsv"
    run:
        import pandas as pd

        def get_loci(grp):
            #print(grp)
            organism = grp.name
            #print(organism)
            symbols = grp.GeneSymbol.tolist()
            #print(symbols)
            annot_file = f"scratch/extracted_annotation/{organism}.tsv"
            #print(annot_file)
            symbol_col = dict(
                amphioxus = "gene_id",
                seaLamprey = "gene_id",
                atlanticHagfish = "gene_id",
                brownHagfish = "gene_name",
            )[organism]
            df = (
                pd.read_csv(annot_file, sep = "\t", usecols = [
                    "primary_tag", symbol_col, "seq_id", "start", "end", "strand"
                ])
                .rename(columns = {
                    symbol_col: "GeneSymbol",
                    "seq_id": "Chromosome",
                })
                .query("primary_tag == 'gene'")
                .drop(columns = ["primary_tag"])
                .query("GeneSymbol in @symbols")
                .assign(Chromosome = lambda tdf: tdf.Chromosome.apply(lambda x: nc2chr[x] if x in nc2chr else x))
            )
            #print(df)
            grp = grp.merge(df, how = "left")
            #print(grp)
            return(grp)


        df = (
            pd.read_csv(input.oglist, sep = "\t")
            .groupby("Orthogroup").agg("_".join)
            .reset_index()
            .merge(pd.read_csv(input.ogfile, sep = "\t"), how = "left")
            .melt(id_vars = ["Orthogroup", "Gene"], var_name = "OrganismShortName", value_name = "ProteinID")
            .query("OrganismShortName in @annotation")
            .dropna()
            .assign(ProteinID = lambda pdf: pdf.ProteinID.apply(lambda x: x.split(", ")))
            .explode("ProteinID")
            .assign(GeneSymbol = lambda tdf: tdf.ProteinID.str.split("__").str[1])
            .groupby("OrganismShortName")
            .apply(get_loci)
            .reset_index(drop = True)
            .sort_values(["OrganismShortName", "Chromosome", "start", "end"])
        )
        print(df)
        df.to_csv(output[0], sep = "\t", index = False)


