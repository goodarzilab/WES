INbam=$1
sample=${INbam/.bam/}

# If mean read length is less than 70bp: DO SOMETING ELSE!!
# bwa mem -t 8 -T 0 -R <read_group> <reference> <fastq_1.fq.gz> <fastq_2.fq.gz> | samtools view -Shb -o <output.bam> -

1. correct bam header and write it to a temp bam files 
samtools view -@ 10 -Shb -o ${sample}.temp1.bam $INbam

2. sort bam file 
gatk SortSam CREATE_INDEX=true INPUT=${sample}.temp1.bam OUTPUT=${sample}.sort.bam SORT_ORDER=coordinate VALIDATION_STRINGENCY=STRICT

3. merge bam file 
gatk MergeSamFiles ASSUME_SORTED=false CREATE_INDEX=true MERGE_SEQUENCE_DICTIONARIES=false SORT_ORDER=coordinate USE_THREADING=true VALIDATION_STRINGENCY=STRICT INPUT=${sample}.sort.bam OUTPUT=${sample}.merge.bam

4. MARK DUPLICATES 
gatk MarkDuplicates CREATE_INDEX=true INPUT=${sample}.merge.bam OUTPUT=${sample}.dups.bam METRICS_FILE=${sample}.dups.mtx VALIDATION_STRINGENCY=STRICT

5. Add headers to bam file
gatk AddOrReplaceReadGroups I=${sample}.dups.bam O=${sample}.out.bam RGID=4 RGLB=lib1 RGPL=illumina RGPU=unit1 RGSM=20


# # remove temp files
# rm -v ${sample}.temp1.bam ${sample}.sort.bam ${sample}.merge.bam
# ${sample}.dups.mtx
