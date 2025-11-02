#!/usr/bin/env bash

#====================================================#
# * Purpose: To download needed refs and other data
# annotations for ATAC-seq matching
#
# * Requires: "tmm" conda env 
#====================================================#

set -ueo pipefail

#===== Setup (1,2)

printf "(1) Checking if conda is activate...\n"
if [[ $CONDA_DEFAULT_ENV == "tmm" ]]; then
	printf "\t--> The conda env you need is already activate, good...\n\n"
elif [[ $CONDA_DEFAULT_ENV != "tmm" ]]; then
	printf "\t--> You need to activate 'tmm'.\n\n"
else
	printf "\t--> Something is wrong...\n\n"
	exit
fi

printf "(2) Ensuring we have the necessary reference genome...\n"
if [ -d /net/talisker/home/benos/mae117/.local/share/genomes/bosTau9 ]; then
	printf "\t--> We have the genome assembey for Bos taurus (domestic cattle).\n\n"
else
	printf "\t--> We are going to use genomepy to get the reference for Bos taurus!\n\n"
	genomepy install --annotation --threads 8 bosTau9
fi

#===== Arguments
REF="/net/talisker/home/benos/mae117/.local/share/genomes/bosTau9/index/bwa/bosTau9.fa"

#===== Pipeline (3,4,5,6)

#----- Generate aligner indicies

# normally we would need to generate an index from the ref with an aligner (done via genomepy)
printf "(3) Now we would normally generate an index with an aligner, but genomepy does that for us!\n\n"

#----- Index reads to the reference

# we would also need to generate a fai (but again, instead of using samtools, we used genomepy)
printf "(4) What would follow is the generation of an indexer for the reference, again, handled by genomepy.\n\n"

#----- Preprocess the reads

printf "(5) We move on to preprocessing the reads...we use fastp here.\n\n"

# controls
mkdir -p controls/preprocessed
for i in $(cat controls.txt); do PR1="controls/${i}_1.fastq.gz"; PR2="controls/${i}_2.fastq.gz"; fastp -i $PR1 -I $PR2 -o "controls/preprocessed/${i}_1.fastq.gz" -O "controls/preprocessed/${i}_2.fastq.gz" --detect_adapter_for_pe -w 8; done

# peaks
mkdir -p experiments/peak/preprocessed
for i in $(cat experiments-peak.txt); do PR1="experiments/peak/${i}_1.fastq.gz"; PR2="experiments/peak/${i}_2.fastq.gz"; fastp -i $PR1 -I $PR2 -o "experiments/peak/preprocessed/${i}_1.fastq.gz" -O "experiments/peak/preprocessed/${i}_2.fastq.gz" --detect_adapter_for_pe -w 8; done

# lates
mkdir -p experiments/late/preprocessed
for i in $(cat experiments-late.txt); do PR1="experiments/late/${i}_1.fastq.gz"; PR2="experiments/late/${i}_2.fastq.gz"; fastp -i $PR1 -I $PR2 -o "experiments/late/preprocessed/${i}_1.fastq.gz" -O "experiments/late/preprocessed/${i}_2.fastq.gz" --detect_adapter_for_pe -w 8; done

#----- Now I align the preprocessed reads to the reference

printf "(6) We align the processed reads to the reference using bowtie2.\n"

# controls
mkdir -p controls/bostau9/bams
for i in $(cat controls.txt); do 
  bwa-mem2 mem -t 16 $REF "controls/preprocessed/${i}_1.fastq.gz" "controls/preprocessed/${i}_2.fastq.gz" \
  | samtools view -b - \
  | samtools sort -@8 -o "controls/bostau9/bams/${i}.sorted.bam"
  printf "\t-->We are done creating the sorted bam file for: ${i}\n"
  samtools index "controls/bostau9/bams/${i}.sorted.bam"
done

# peaks
mkdir -p experiments/peak/bostau9/bams
for i in $(cat experiments-peak.txt); do 
  bwa-mem2 mem -t 16 $REF "experiments/peak/preprocessed/${i}_1.fastq.gz" "experiments/peak/preprocessed/${i}_2.fastq.gz" \
  | samtools view -b - \
  | samtools sort -@8 -o "experiments/peak/bostau9/bams/${i}.sorted.bam"
  printf "\t-->We are done creating the sorted bam file for: ${i}\n"
  samtools index "experiments/peak/bostau9/bams/${i}.sorted.bam"
done

# lates
mkdir -p experiments/late/bostau9/bams
for i in $(cat experiments-late.txt); do 
  bwa-mem2 mem -t 16 $REF "experiments/late/preprocessed/${i}_1.fastq.gz" "experiments/late/preprocessed/${i}_2.fastq.gz" \
  | samtools view -b - \
  | samtools sort -@8 -o "experiments/late/bostau9/bams/${i}.sorted.bam"
  printf "\t-->We are done creating the sorted bam file for: ${i}\n"
  samtools index "experiments/late/bostau9/bams/${i}.sorted.bam"
done

#----- Deduplicate (ATAC bulk often dedup true)

printf "(7) We now focus on dedup-ing the aligned reads to the reference!\n"

# controls
for i in $(cat controls.txt); do
  picard -Xmx4g AddOrReplaceReadGroups I=controls/bostau9/bams/${i}.sorted.bam O=controls/bostau9/bams/${i}.rg.bam RGID="${i}" RGLB=lib1 RGPL=ILLUMINA RGPU=${i}.1 RGSM=${i} 
  samtools index controls/bostau9/bams/${i}.rg.bam
done
for i in $(cat controls.txt); do
  picard -Xmx8g MarkDuplicates I=controls/bostau9/bams/${i}.rg.bam O=controls/bostau9/bams/${i}.dedup.bam M=controls/bostau9/bams/${i}.dup.txt REMOVE_DUPLICATES=true OPTICAL_DUPLICATE_PIXEL_DISTANCE=2500
  samtools index controls/bostau9/bams/${i}.dedup.bam
done
printf "\t-->Controls are dedupped.\n"

# peaks
for i in $(cat experiments-peak.txt); do
  picard -Xmx4g AddOrReplaceReadGroups I=experiments/peak/bostau9/bams/${i}.sorted.bam O=experiments/peak/bostau9/bams/${i}.rg.bam RGID="${i}" RGLB=lib1 RGPL=ILLUMINA RGPU=${i}.1 RGSM=${i} 
  samtools index experiments/peak/bostau9/bams/${i}.rg.bam
done
for i in $(cat experiments-peak.txt); do
  picard -Xmx8g MarkDuplicates I=experiments/peak/bostau9/bams/${i}.rg.bam O=experiments/peak/bostau9/bams/${i}.dedup.bam M=experiments/peak/bostau9/bams/${i}.dup.txt REMOVE_DUPLICATES=true OPTICAL_DUPLICATE_PIXEL_DISTANCE=2500
  samtools index experiments/peak/bostau9/bams/${i}.dedup.bam
done
printf "\t-->Peaks are dedupped.\n"

# lates
for i in $(cat experiments-late.txt); do
  picard -Xmx4g AddOrReplaceReadGroups I=experiments/late/bostau9/bams/${i}.sorted.bam O=experiments/late/bostau9/bams/${i}.rg.bam RGID="${i}" RGLB=lib1 RGPL=ILLUMINA RGPU=${i}.1 RGSM=${i} 
  samtools index experiments/late/bostau9/bams/${i}.rg.bam
done
for i in $(cat experiments-late.txt); do
  picard -Xmx8g MarkDuplicates I=experiments/late/bostau9/bams/${i}.rg.bam O=experiments/late/bostau9/bams/${i}.dedup.bam M=experiments/late/bostau9/bams/${i}.dup.txt REMOVE_DUPLICATES=true OPTICAL_DUPLICATE_PIXEL_DISTANCE=2500
  samtools index experiments/late/bostau9/bams/${i}.dedup.bam
done
printf "\t-->Lates are dedupped.\n\n"

#----- Filter deduplicated reads and perform additional processing.

printf "(8) Filter deduplicate reads, and additional processing.\n"

# controls
for i in $(cat controls.txt); do
  samtools idxstats "controls/bostau9/bams/${i}.dedup.bam" | cut -f1 | grep -v 'chrM' | grep -v '_' | grep -v '*' > "controls/bostau9/bams/${i}.keep_chrs.txt"
  samtools view -b -f 0x2 -q 30 "controls/bostau9/bams/${i}.dedup.bam" $(<"controls/bostau9/bams/${i}.keep_chrs.txt") > "controls/bostau9/bams/${i}.flt.bam"
  samtools index "controls/bostau9/bams/${i}.flt.bam"
done

# peaks
for i in $(cat experiments-peak.txt); do
  samtools idxstats "experiments/peak/bostau9/bams/${i}.dedup.bam" | cut -f1 | grep -v 'chrM' | grep -v '_' | grep -v '*' > "experiments/peak/bostau9/bams/${i}.keep_chrs.txt"
  samtools view -b -f 0x2 -q 30 "experiments/peak/bostau9/bams/${i}.dedup.bam" $(<"experiments/peak/bostau9/bams/${i}.keep_chrs.txt") > "experiments/peak/bostau9/bams/${i}.flt.bam"
  samtools index "experiments/peak/bostau9/bams/${i}.flt.bam"
done

# lates
for i in $(cat experiments-late.txt); do
  samtools idxstats "experiments/late/bostau9/bams/${i}.dedup.bam" | cut -f1 | grep -v 'chrM' | grep -v '_' | grep -v '*' > "experiments/late/bostau9/bams/${i}.keep_chrs.txt"
  samtools view -b -f 0x2 -q 30 "experiments/late/bostau9/bams/${i}.dedup.bam" $(<"experiments/late/bostau9/bams/${i}.keep_chrs.txt") > "experiments/late/bostau9/bams/${i}.flt.bam"
  samtools index "experiments/late/bostau9/bams/${i}.flt.bam"
done

#----- Make fragments files for scPrinter

printf "(9) We move to the final step: creating the fragment files for scPrinter.\n"

# controls
mkdir -p controls/bostau9/fragments
for i in $(cat controls.txt); do
  inbam="controls/bostau9/bams/${i}.flt.bam"
  nsbam="controls/bostau9/bams/${i}.flt.namesort.bam"
  outgz="controls/bostau9/fragments/${i}.fragments.tsv.gz"

  # Name-sort (required for -bedpe)
  samtools sort -n -@8 -o "$nsbam" "$inbam"

  # BAM -> BEDPE -> fragments (chr start end sample)
  bedtools bamtobed -bedpe -i "$nsbam" 2> "controls/bostau9/fragments/${i}.bamtobed.warn.log" \
  | awk -v OFS='\t' -v BARC="$i" '($1==$4 && $3-$2 < 1000){
        s = ($2<$5 ? $2 : $5);
        e = ($3>$6 ? $3 : $6);
        print $1, s, e, BARC
    }' \
  | sort -k1,1 -k2,2n \
  | bgzip > "$outgz"

  tabix -p bed "$outgz"
done

# peaks
mkdir -p experiments/peak/bostau9/fragments
for i in $(cat experiments-peak.txt); do
  inbam="experiments/peak/bostau9/bams/${i}.flt.bam"
  nsbam="experiments/peak/bostau9/bams/${i}.flt.namesort.bam"
  outgz="experiments/peak/bostau9/fragments/${i}.fragments.tsv.gz"

  # Name-sort (required for -bedpe)
  samtools sort -n -@8 -o "$nsbam" "$inbam"

  # BAM -> BEDPE -> fragments (chr start end sample)
  bedtools bamtobed -bedpe -i "$nsbam" 2> "experiments/peak/bostau9/fragments/${i}.bamtobed.warn.log" \
  | awk -v OFS='\t' -v BARC="$i" '($1==$4 && $3-$2 < 1000){
        s = ($2<$5 ? $2 : $5);
        e = ($3>$6 ? $3 : $6);
        print $1, s, e, BARC
      }' \
  | sort -k1,1 -k2,2n -k3,3n \
  | bgzip > "$outgz"

  tabix -p bed "$outgz"
done

# lates
mkdir -p experiments/late/bostau9/bams
for i in $(cat experiments-late.txt); do
  inbam="experiments/late/bostau9/bams/${i}.flt.bam"
  nsbam="experiments/late/bostau9/bams/${i}.flt.namesort.bam"
  outgz="experiments/late/bostau9/fragments/${i}.fragments.tsv.gz"

  # Name-sort (required for -bedpe)
  samtools sort -n -@8 -o "$nsbam" "$inbam"

  # BAM -> BEDPE -> fragments (chr start end sample)
  bedtools bamtobed -bedpe -i "$nsbam" 2> "experiments/late/bostau9/fragments/${i}.bamtobed.warn.log" \
  | awk -v OFS='\t' -v BARC="$i" '($1==$4 && $3-$2 < 1000){
        s = ($2<$5 ? $2 : $5);
        e = ($3>$6 ? $3 : $6);
        print $1, s, e, BARC
      }' \
  | sort -k1,1 -k2,2n -k3,3n \
  | bgzip > "$outgz"

  tabix -p bed "$outgz"
done

#===== Teardown

exit 0