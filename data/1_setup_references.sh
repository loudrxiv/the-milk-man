#!/usr/bin/env bash

#====================================================#
# * Purpose: To download needed refs and other data
# annotations for ATAC-seq matching
#
# * Requires: "the-milk-man" conda env 
#====================================================#

set -ueo pipefail

#===== Setup
printf "(1) Checking if conda is activate...\n"
if [[ $CONDA_DEFAULT_ENV == "the-milk-man" ]]; then
	printf "\t--> The conda env you need is already activate, good...\n\n"
elif [[ $CONDA_DEFAULT_ENV != "the-milk-man" ]]; then
	printf "\t--> You need to activate 'the-milk-man'.\n\n"
else
	printf "\t--> Something is wrong...\n\n"
	exit
fi

printf "(2) Ensuring we have the necessary reference genome...\n"
if [ -d /net/talisker/home/benos/mae117/.local/share/genomes/ARS-UCD2.0 ]; then
	printf "\t--> We have the genome assembey for Bos taurus (domestic cattle).\n\n"
else
	printf "\t--> We are going to use genomepy to get the reference for Bos taurus!\n\n"
	genomepy install --annotation --threads 8 ARS-UCD2.0
fi

#===== Arguments
REF="/net/talisker/home/benos/mae117/.local/share/genomes/ARS-UCD2.0/ARS-UCD2.0.fa"

#===== Pipeline

# normally we would need to generate an index from the ref with an aligner (done via genomepy)
printf "(3) Now we would normally generate an index with an aligner, but genomepy does that for us!\n\n"

# we would also need to generate a fai (but again, instead of using samtools, we used genomepy)
printf "(4) What would follow is the generation of an indexer for the reference, again, handled by genomepy.\n\n"

# now i need to trim the reads
printf "(5) We move on to preprocessing the reads...we use fastp here.\n"
mkdir -p controls/preprocessed
for i in $(cat controls.txt); do PR1="controls/${i}_1.fastq.gz"; PR2="controls/${i}_2.fastq.gz"; fastp -i $PR1 -I $PR2 -o "controls/preprocessed/${i}_1.fastq.gz" -O "controls/preprocessed/${i}_2.fastq.gz" --detect_adapter_for_pe -w 8; done

mkdir -p experiments/peak/preprocessed
for i in $(cat experiments-peak.txt); do PR1="experiments/peak/${i}_1.fastq.gz"; PR2="experiments/peak/${i}_2.fastq.gz"; fastp -i $PR1 -I $PR2 -o "experiments/peak/preprocessed/${i}_1.fastq.gz" -O "experiments/peak/preprocessed/${i}_2.fastq.gz" --detect_adapter_for_pe -w 8; done

mkdir -p experiments/late/preprocessed
for i in $(cat experiments-late.txt); do PR1="experiments/late/${i}_1.fastq.gz"; PR2="experiments/late/${i}_2.fastq.gz"; fastp -i $PR1 -I $PR2 -o "experiments/late/preprocessed/${i}_fastq.gz" -O "experiments/late/preprocessed/${i}_2.fastq.gz" --detect_adapter_for_pe -w 8; done

# now i align the preprocessed reads to the reference
#mkdir -p controls/bams
#for i in $(cat controls.txt); do bwa-mem2 mem -t 16 $REF "controls/preprocessed/${i}_1.fastq.gz" "controls/preprocessed/${i}_2.fastq.gz" | samtools view -b - | samtools sort -@8 -o "controls/bams/${i}.sorted.bam"; samtools index "controls/bam/${i}.sorted.bam"; done

#mkdir -p experiments/peak/bams
#for i in $(cat controls.txt); do bwa-mem2 mem -t 16 $REF "experiments/peak/preprocessed/${i}_1.fastq.gz" "experiments/peak/preprocessed/${i}_2.fastq.gz" | samtools view -b - | samtools sort -@8 -o "experiments/peak/bams/${i}.sorted.bam"; samtools index "experiments/peak/bam/${i}.sorted.bam"; done

#mkdir -p experiments/late/bams
#for i in $(cat controls.txt); do bwa-mem2 mem -t 16 $REF "experiments/late/preprocessed/${i}_1.fastq.gz" "experiments/late/preprocessed/${i}_2.fastq.gz" | samtools view -b - | samtools sort -@8 -o "experiments/late/bams/${i}.sorted.bam"; samtools index "experiments/late/bam/${i}.sorted.bam"; done

#===== Teardown

exit 0
