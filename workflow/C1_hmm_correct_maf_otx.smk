rule all:
    input:
        expand("scratch/hmm_corrected/hmm_corrected_{genegrp}.newick", genegrp = [
 #           "CMAF", "NRL", "MAFA", "MAFB",
 #           "OTX1", "OTX2", "CRX",
 #           "LargeMAF",
            "SmallMAF",
 #           "AllMAF",
 #           "AllOTX",
        ]),

rule hmm_correct_protein:
    input:
        grand_faa = "scratch/grand_list/grand_list.faa",
        train_pred = "configs/train_hmm/hmm_train_pred_{genegrp}.tsv"
    output:
        pep = "scratch/hmm_corrected/hmm_corrected_{genegrp}.faa",
        aln = "scratch/hmm_corrected/hmm_corrected_{genegrp}_aln.fasta",
        tree = "scratch/hmm_corrected/hmm_corrected_{genegrp}.newick",
        pred_list = "scratch/hmm_corrected/hmm_corrected_{genegrp}/pred_list.txt",
        train_list = "scratch/hmm_corrected/hmm_corrected_{genegrp}/train_list.txt",
        pred_faa = "scratch/hmm_corrected/hmm_corrected_{genegrp}/pred.faa",
        train_faa = "scratch/hmm_corrected/hmm_corrected_{genegrp}/train.faa",
        train_aln = "scratch/hmm_corrected/hmm_corrected_{genegrp}/train_aln.fasta",
        train_hmm = "scratch/hmm_corrected/hmm_corrected_{genegrp}/train.hmm",
        pred_sto = "scratch/hmm_corrected/hmm_corrected_{genegrp}/pred.sto",
        masked_faa = "scratch/hmm_corrected/hmm_corrected_{genegrp}/masked.faa",
        masked_sto = "scratch/hmm_corrected/hmm_corrected_{genegrp}/masked.sto",
    shell:
        """
        module load seqkit clustalo FastTree hmmer
        cat {input.train_pred} | grep "pred" | cut -f1 > {output.pred_list}
        cat {input.train_pred} | grep "train" | cut -f1 > {output.train_list}
        seqkit grep -n -f {output.pred_list} {input.grand_faa} -o {output.pred_faa}
        seqkit grep -n -f {output.train_list} {input.grand_faa} -o {output.train_faa}
        clustalo -i {output.train_faa} -o {output.train_aln} --threads 32 --outfmt=fasta
        hmmbuild {output.train_hmm} {output.train_aln} > hmmbuild.log 2>&1
        hmmalign --trim {output.train_hmm} {output.pred_faa} > {output.pred_sto}
        esl-alimask -g --gapthresh 0.5 -o {output.masked_sto} {output.pred_sto}
        esl-reformat fasta {output.masked_sto} > {output.masked_faa}
        cat {output.train_faa} {output.masked_faa} > {output.pep}
        clustalo -i {output.pep} -o {output.aln} --threads=32 --outfmt=fasta
        FastTreeMP < {output.aln} > {output.tree}
        """

