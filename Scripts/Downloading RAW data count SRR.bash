

#!/usr/bin/env bash

 

# Usage: ./download_gse_fastq

# Arguments:

#   $1 = GSE accession

#   $2 = number of threads (default: 4)

 

set -euo pipefail

 

GSE=

THREADS=${2:-4}

 

mkdir -p "${GSE}/fastq"

cd "${GSE}"

 

echo "Fetching SRR accessions for ${GSE}..."

 

SRR_LIST=$(esearch -db gds -query "${GSE}" \

  | elink -target sra \

  | efetch -format runinfo \

  | cut -d',' -f1 \

  | grep SRR)

 

echo "Found $(echo "$SRR_LIST" | wc -l) SRR runs"

 

echo " Downloading FASTQs..."

 

for SRR in $SRR_LIST; do

  echo "Processing ${SRR}"

 

  fasterq-dump "${SRR}" \

    --split-files \

    --threads "${THREADS}" \

    --outdir fastq

 

done

 

echo " .gz formatting FASTQs..."

gzip fastq/*.fastq

 

echo " coverting to Illumina-style filenames..."

 

cd fastq

for f in *_1.fastq.gz; do

  base=$(basename "$f" _1.fastq.gz)

  mv "$f" "${base}_R1_001.fastq.gz"

done

 

for f in *_2.fastq.gz; do

  base=$(basename "$f" _2.fastq.gz)

  mv "$f" "${base}_R2_001.fastq.gz"

done

 

echo "Done!"

