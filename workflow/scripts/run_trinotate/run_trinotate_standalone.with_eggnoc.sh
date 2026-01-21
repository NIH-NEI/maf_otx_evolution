#!/usr/bin/env bash
TRINITY_RES_TOP="/data/VisionEvo/TrinityOutput/Assemblies"
TRANSDECODER_RES_TOP="/data/VisionEvo/TransDecoderOutput/OrganismLevel"
TRINOTATE_RES_TOP="/data/VisionEvo/TrinotateOutput/Reports"
export TRINOTATE_DATA_DIR="/data/pals2/tools/Trinotate_data_dir"

if [ "$#" -ne 1 ]; then
    echo "Usage: $0 <orgid>"
    exit 1
fi

org=$1
trans="${TRINITY_RES_TOP}/${org}_denovo.fasta"
g2t="${TRINITY_RES_TOP}/${org}_denovo.fasta.gene_to_trans_map"
pep="${TRANSDECODER_RES_TOP}/${org}/${org}_denovo.fasta.transdecoder.pep"
db="${TRINOTATE_RES_TOP}/${org}/${org}_trinotate_db.sqlite"
rep="${TRINOTATE_RES_TOP}/${org}/${org}_trinotate_report.tsv"
rundir="${TRINOTATE_RES_TOP}/${org}"

echo "trans=$trans"
echo "g2t=$g2t"
echo "db=$db"
echo "rep=$rep"
echo "rundir=$rundir"

set -euo pipefail

module load trinotate
mkdir -p $rundir
cd $rundir
Trinotate --db ${db} --create --trinotate_data_dir ${TRINOTATE_DATA_DIR}
Trinotate --db ${db} --init \
    --gene_trans_map ${g2t} \
    --transcript_fasta ${trans} \
    --transdecoder_pep ${pep}
Trinotate --db ${db} --run ALL \
    --trinotate_data_dir ${TRINOTATE_DATA_DIR} \
    --transcript_fasta ${trans} \
    --transdecoder_pep ${pep} \
    --use_diamond
Trinotate --db ${db} --report \
    --incl_pep --incl_trans > ${rep}
