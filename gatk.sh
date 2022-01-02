mkdir {sampleid}sample
# this directory for each sample should have fastq files and bam files after xenofiltering
# all outputs are made in the same folder 

gatk MarkDuplicates CREATE_INDEX=true INPUT={sampleid}.processed_Filtered.bam OUTPUT={sampleid}.bam METRICS_FILE={sampleid}.dups.mtx VALIDATION_STRINGENCY=STRICT

gatk AddOrReplaceReadGroups I={sampleid}.bam O={sampleid}out.bam RGID=4 RGLB=lib1 RGPL=illumina RGPU=unit1 RGSM=20

gatk CollectAlignmentSummaryMetrics  R=hg38.fa I={sampleid}out.bam O={sampleid}alignment_metrics.txt

gatk CollectInsertSizeMetrics I={sampleid}out.bam  O={sampleid}insert_metrics.txt H={sampleid}insert_size_histogram.pdf

samtools depth -a {sampleid}out.bam  > {sampleid}depth_out.txt

gatk CreateSequenceDictionary R=hg38.fa O=hg38.dict
# creates necessary reference files for next steps

gatk HaplotypeCaller -R hg38.fa -I {sampleid}out.bam -O {sampleid}raw_variants.vcf

gatk SelectVariants -R hg38.fa -V {sampleid}raw_variants.vcf  -select-type SNP -O {sampleid}raw_snps.vcf

gatk SelectVariants -R hg38.fa -V {sampleid}raw_variants.vcf -select-type INDEL -O {sampleid}raw_indels.vcf

gatk VariantFiltration -R hg38.fa -V {sampleid}raw_snps.vcf -O {sampleid}filtered_snps.vcf -filter-name "QD_filter" -filter "QD < 2.0" -filter-name "FS_filter" -filter "FS > 60.0" -filter-name "MQ_filter" -filter "MQ < 40.0" -filter-name "SOR_filter" -filter "SOR > 4.0" -filter-name "MQRankSum_filter" -filter "MQRankSum < -12.5" -filter-name "ReadPosRankSum_filter" -filter "ReadPosRankSum < -8.0"

gatk VariantFiltration -R hg38.fa -V {sampleid}raw_indels.vcf -O {sampleid}filtered_indels.vcf -filter-name "QD_filter" -filter "QD < 2.0" -filter-name "FS_filter" -filter "FS > 200.0" -filter-name "SOR_filter" -filter "SOR > 10.0"

gatk SelectVariants --exclude-filtered -V {sampleid}filtered_snps.vcf -O {sampleid}bqsr_snps.vcf

gatk SelectVariants --exclude-filtered -V {sampleid}filtered_indels.vcf -O {sampleid}bqsr_indels.vcf

gatk BaseRecalibrator -R hg38.fa -I {sampleid}out.bam --known-sites {sampleid}bqsr_snps.vcf --known-sites {sampleid}bqsr_indels.vcf -O {sampleid}recal_data.table

gatk ApplyBQSR -R hg38.fa -I {sampleid}out.bam -bqsr {sampleid}recal_data.table -O {sampleid}recal_reads.bam 

gatk BaseRecalibrator -R hg38.fa -I {sampleid}recal_reads.bam --known-sites {sampleid}bqsr_snps.vcf --known-sites {sampleid}bqsr_indels.vcf -O {sampleid}post_recal_data.table

gatk AnalyzeCovariates -before {sampleid}recal_data.table -after {sampleid}post_recal_data.table -plots {sampleid}recalibration_plots.pdf

gatk HaplotypeCaller -R hg38.fa -I {sampleid}recal_reads.bam -O {sampleid}raw_variants_recal.vcf

gatk SelectVariants -R hg38.fa -V {sampleid}raw_variants_recal.vcf -select-type SNP -O {sampleid}raw_snps_recal.vcf

gatk SelectVariants -R hg38.fa -V {sampleid}raw_variants.vcf -select-type INDEL -O {sampleid}raw_indels_recal.vcf

gatk VariantFiltration -R hg38.fa -V {sampleid}raw_snps_recal.vcf -O {sampleid}filtered_snps_final.vcf -filter-name "QD_filter" -filter "QD < 2.0" -filter-name "FS_filter" -filter "FS > 60.0" -filter-name "MQ_filter" -filter "MQ < 40.0" -filter-name "SOR_filter" -filter "SOR > 4.0" -filter-name "MQRankSum_filter" -filter "MQRankSum < -12.5" -filter-name "ReadPosRankSum_filter" -filter "ReadPosRankSum < -8.0"

gatk VariantFiltration -R hg38.fa -V {sampleid}raw_indels_recal.vcf -O {sampleid}filtered_indels_final.vcf -filter-name "QD_filter" -filter "QD < 2.0" -filter-name "FS_filter" -filter "FS > 200.0" -filter-name "SOR_filter" -filter "SOR > 10.0"

# download snpEff and create path in directory before running this command
# select the correct reference genome database to match human
java -jar ./snpEff/snpEff.jar -v GRCh38.99 ./dna/filtered_snps_final.vcf > {sampleid}filtered_snps_final.ann.vcf

parse_metrics.sh {sampleid} > {sampleid}_report.csv









