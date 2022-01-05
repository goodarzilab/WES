#link to gatk pipeline: https://gencore.bio.nyu.edu/variant-calling-pipeline-gatk4
bam_file=$1
sampleid=$2

mkdir ${sampleid}sample
# this directory for each sample should have fastq files and bam files after xenofiltering
# all outputs are made in the same folder 

#input is bam files, output is bam file without duplicate reads
gatk MarkDuplicates CREATE_INDEX=true INPUT=${bam_file} OUTPUT=${sampleid}/${sampleid}.bam METRICS_FILE=${sampleid}/${sampleid}.dups.mtx VALIDATION_STRINGENCY=STRICT

#input is bam file from step above, output is bam file with correct read tags
gatk AddOrReplaceReadGroups I=${sampleid}/${sampleid}.bam O=${sampleid}/${sampleid}out.bam RGID=4 RGLB=lib1 RGPL=illumina RGPU=unit1 RGSM=20

#input is reference and updated bam file, output is alignment file
gatk CollectAlignmentSummaryMetrics  R=hg38.fa I=${sampleid}/${sampleid}out.bam O=${sampleid}/${sampleid}alignment_metrics.txt

#input is updated bam file, output is metrics file and histogram plots
gatk CollectInsertSizeMetrics I=${sampleid}/${sampleid}out.bam  O=${sampleid}/${sampleid}insert_metrics.txt H=${sampleid}/${sampleid}insert_size_histogram.pdf

samtools depth -a ${sampleid}/${sampleid}out.bam  > ${sampleid}/${sampleid}depth_out.txt

#input is reference and updated bam, creates a raw file of variants
gatk HaplotypeCaller -R hg38.fa -I ${sampleid}/${sampleid}out.bam -O ${sampleid}/${sampleid}raw_variants.vcf

#imput is reference and raw variants, output is raw files of snp
gatk SelectVariants -R hg38.fa -V ${sampleid}/${sampleid}raw_variants.vcf  -select-type SNP -O ${sampleid}/${sampleid}raw_snps.vcf

#input is reference and raw variants file, output is raw file for indels
gatk SelectVariants -R hg38.fa -V ${sampleid}/${sampleid}raw_variants.vcf -select-type INDEL -O ${sampleid}/${sampleid}raw_indels.vcf

#input is reference and raw snps, output is the same file filtered with parameters
gatk VariantFiltration -R hg38.fa -V ${sampleid}/${sampleid}raw_snps.vcf -O ${sampleid}/${sampleid}filtered_snps.vcf -filter-name "QD_filter" -filter "QD < 2.0" -filter-name "FS_filter" -filter "FS > 60.0" -filter-name "MQ_filter" -filter "MQ < 40.0" -filter-name "SOR_filter" -filter "SOR > 4.0" -filter-name "MQRankSum_filter" -filter "MQRankSum < -12.5" -filter-name "ReadPosRankSum_filter" -filter "ReadPosRankSum < -8.0"

#input is reference and raw indels file, output is filtered file with parameters
gatk VariantFiltration -R hg38.fa -V ${sampleid}/${sampleid}raw_indels.vcf -O ${sampleid}/${sampleid}filtered_indels.vcf -filter-name "QD_filter" -filter "QD < 2.0" -filter-name "FS_filter" -filter "FS > 200.0" -filter-name "SOR_filter" -filter "SOR > 10.0"

#input is filtered snps file, output is bsqr file for snps
gatk SelectVariants --exclude-filtered -V ${sampleid}/${sampleid}filtered_snps.vcf -O ${sampleid}/${sampleid}bqsr_snps.vcf

#input is filtered indels file, output is bsqr indels 
gatk SelectVariants --exclude-filtered -V ${sampleid}/${sampleid}filtered_indels.vcf -O ${sampleid}/${sampleid}bqsr_indels.vcf

#input is reference and updated bam file and all bsqr files, output is calibrated data table
gatk BaseRecalibrator -R hg38.fa -I ${sampleid}/${sampleid}out.bam --known-sites ${sampleid}/${sampleid}bqsr_snps.vcf --known-sites ${sampleid}/${sampleid}bqsr_indels.vcf -O ${sampleid}/${sampleid}recal_data.table

#input is reference and updated bam and first data table, output is bam file of recalibrated reads
gatk ApplyBQSR -R hg38.fa -I ${sampleid}/${sampleid}out.bam -bqsr ${sampleid}/${sampleid}recal_data.table -O ${sampleid}/${sampleid}recal_reads.bam 

#input is reference, all bsqr files, and calibrated bam file, output is a post data table
gatk BaseRecalibrator -R hg38.fa -I ${sampleid}/${sampleid}recal_reads.bam --known-sites ${sampleid}/${sampleid}bqsr_snps.vcf --known-sites ${sampleid}/${sampleid}bqsr_indels.vcf -O ${sampleid}/${sampleid}post_recal_data.table

#input is both data tables, output is recalibration plots 
gatk AnalyzeCovariates -before ${sampleid}/${sampleid}recal_data.table -after ${sampleid}/${sampleid}post_recal_data.table -plots ${sampleid}/${sampleid}recalibration_plots.pdf

#input is reference and recal bam file, output is recal file of raw variants
gatk HaplotypeCaller -R hg38.fa -I ${sampleid}/${sampleid}recal_reads.bam -O ${sampleid}/${sampleid}raw_variants_recal.vcf

#input is reference and recal variants file, output is recal snps file
gatk SelectVariants -R hg38.fa -V ${sampleid}/${sampleid}raw_variants_recal.vcf -select-type SNP -O ${sampleid}/${sampleid}raw_snps_recal.vcf

#inout is reference and raw variants file, output is recal indels file
gatk SelectVariants -R hg38.fa -V ${sampleid}/${sampleid}raw_variants.vcf -select-type INDEL -O ${sampleid}/${sampleid}raw_indels_recal.vcf

#input is reference and recal snps file, output is filtered snps file with parameters
gatk VariantFiltration -R hg38.fa -V ${sampleid}/${sampleid}raw_snps_recal.vcf -O ${sampleid}/${sampleid}filtered_snps_final.vcf -filter-name "QD_filter" -filter "QD < 2.0" -filter-name "FS_filter" -filter "FS > 60.0" -filter-name "MQ_filter" -filter "MQ < 40.0" -filter-name "SOR_filter" -filter "SOR > 4.0" -filter-name "MQRankSum_filter" -filter "MQRankSum < -12.5" -filter-name "ReadPosRankSum_filter" -filter "ReadPosRankSum < -8.0"

#input is reference and recal indels file, output is filtered indels file with parameters
gatk VariantFiltration -R hg38.fa -V ${sampleid}/${sampleid}raw_indels_recal.vcf -O ${sampleid}/${sampleid}filtered_indels_final.vcf -filter-name "QD_filter" -filter "QD < 2.0" -filter-name "FS_filter" -filter "FS > 200.0" -filter-name "SOR_filter" -filter "SOR > 10.0"

# download snpEff and create path in directory before running this command
# select the correct reference genome database to match human
java -jar ./snpEff/snpEff.jar -v GRCh38.99 ./dna/filtered_snps_final.vcf > ${sampleid}/${sampleid}filtered_snps_final.ann.vcf

parse_metrics.sh ${sampleid} > ${sampleid}/${sampleid}_report.csv









