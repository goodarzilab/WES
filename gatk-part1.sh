#link to gatk pipeline: https://gencore.bio.nyu.edu/variant-calling-pipeline-gatk4
bam_file=$1
sampleid=$2
genome_fa=$3

mkdir ${sampleid}
# this directory for each sample should have fastq files and bam files after xenofiltering
# all outputs are made in the same folder 

#input is bam files, output is bam file without duplicate reads
gatk MarkDuplicates CREATE_INDEX=true \
  INPUT=${bam_file} \
  OUTPUT=${sampleid}/${sampleid}.dups.bam \
  METRICS_FILE=${sampleid}/${sampleid}.dups.mtx \
  VALIDATION_STRINGENCY=STRICT \
  MAX_FILE_HANDLES=1000

#input is bam file from step above, output is bam file with correct read tags
gatk AddOrReplaceReadGroups \
  I=${sampleid}/${sampleid}.dups.bam \
  O=${sampleid}/${sampleid}.ctags.bam \
  RGID=4 RGLB=lib1 RGPL=illumina RGPU=unit1 RGSM=20

#input is reference and updated bam file, output is alignment file
gatk CollectAlignmentSummaryMetrics  \
  R=$genome_fa \
  I=${sampleid}/${sampleid}.ctags.bam \
  O=${sampleid}/${sampleid}.alignment_metrics.txt

#input is updated bam file, output is metrics file and histogram plots
gatk CollectInsertSizeMetrics \
  I=${sampleid}/${sampleid}.ctags.bam  \
  O=${sampleid}/${sampleid}.insert_metrics.txt \
  H=${sampleid}/${sampleid}.insert_size_histogram.pdf

samtools depth -a ${sampleid}/${sampleid}.ctags.bam  > ${sampleid}/${sampleid}.depth_out.txt
