#!/usr/bin/env bash

# Usage:
#   ./run_cellranger_auto.sh /home/cellrangetest/fastqs \
#                            /home/cellrangetest/refdata-gex-GRCh38-2024-A \
#                            16

set -euo pipefail  ##exit when encounter an error in sccript

FASTQ_ROOT="/home/cellrangetest/fastqs"
REFERENCE="/home/cellrangetest/refdata-gex-GRCh38-2024-A"
CORES= 16

LOGDIR="cellranger_logs"
mkdir -p "${LOGDIR}"


for SAMPLE_DIR in "${FASTQ_ROOT}"/*; do
  SAMPLE=$(basename "${SAMPLE_DIR}")

  echo " Processing sample: ${SAMPLE}"

 
  cellranger count --id="run_count_test_${SAMPLE}"  --fastqs="${SAMPLE_DIR}" --transcriptome="${REFERENCE}" \
    --create-bam=true \    --localcores="${CORES}" \

  echo "Processed ${SAMPLE}"
  echo
done

