#link to gatk pipeline: https://gencore.bio.nyu.edu/variant-calling-pipeline-gatk4
sampleid=$1
genome_fa=$2

#input is reference and updated bam, creates a raw file of variants
gatk HaplotypeCaller \
  -R $genome_fa \
  -I ${sampleid}/${sampleid}.ctags.bam \
  -O ${sampleid}/${sampleid}.raw_variants.vcf

#imput is reference and raw variants, output is raw files of snp
gatk SelectVariants \
  -R $genome_fa \
  -select-type SNP \
  -V ${sampleid}/${sampleid}.raw_variants.vcf  \
  -O ${sampleid}/${sampleid}.raw_snps.vcf

#input is reference and raw variants file, output is raw file for indels
gatk SelectVariants \
  -R $genome_fa \
  -select-type INDEL \
  -V ${sampleid}/${sampleid}.raw_variants.vcf \
  -O ${sampleid}/${sampleid}.raw_indels.vcf

#input is reference and raw snps, output is the same file filtered with parameters
gatk VariantFiltration \
  -R $genome_fa \
  -filter-name "QD_filter" -filter "QD < 2.0" \
  -filter-name "FS_filter" -filter "FS > 60.0" \
  -filter-name "MQ_filter" -filter "MQ < 40.0" \
  -filter-name "SOR_filter" -filter "SOR > 4.0" \
  -filter-name "MQRankSum_filter" -filter "MQRankSum < -12.5" \
  -filter-name "ReadPosRankSum_filter" -filter "ReadPosRankSum < -8.0" \
  -V ${sampleid}/${sampleid}.raw_snps.vcf \
  -O ${sampleid}/${sampleid}.filtered_snps.vcf

#input is reference and raw indels file, output is filtered file with parameters
gatk VariantFiltration \
  -R $genome_fa \
  -filter-name "QD_filter" -filter "QD < 2.0" \
  -filter-name "FS_filter" -filter "FS > 200.0" \
  -filter-name "SOR_filter" -filter "SOR > 10.0" \
  -V ${sampleid}/${sampleid}.raw_indels.vcf \
  -O ${sampleid}/${sampleid}.filtered_indels.vcf 


### ApplyBQSR
#input is filtered snps file, output is bsqr file for snps
gatk SelectVariants \
  --exclude-filtered \
  -V ${sampleid}/${sampleid}.filtered_snps.vcf \
  -O ${sampleid}/${sampleid}.bqsr_snps.vcf

#input is filtered indels file, output is bsqr indels 
gatk SelectVariants \
  --exclude-filtered \
  -V ${sampleid}/${sampleid}.filtered_indels.vcf \
  -O ${sampleid}/${sampleid}.bqsr_indels.vcf

#input is reference and updated bam file and all bsqr files, output is calibrated data table
gatk BaseRecalibrator \
  -R $genome_fa \
  -I ${sampleid}/${sampleid}.ctags.bam \
  --known-sites ${sampleid}/${sampleid}.bqsr_snps.vcf \
  --known-sites ${sampleid}/${sampleid}.bqsr_indels.vcf \
  -O ${sampleid}/${sampleid}.recal_data.table

#input is reference and updated bam and first data table, output is bam file of recalibrated reads
gatk ApplyBQSR \
  -R $genome_fa \
  -I ${sampleid}/${sampleid}.ctags.bam \
  -bqsr ${sampleid}/${sampleid}.recal_data.table \
  -O ${sampleid}/${sampleid}.recal_reads.bam 

#input is reference, all bsqr files, and calibrated bam file, output is a post data table
gatk BaseRecalibrator \
  -R $genome_fa \
  -I ${sampleid}/${sampleid}.recal_reads.bam \
  --known-sites ${sampleid}/${sampleid}.bqsr_snps.vcf \
  --known-sites ${sampleid}/${sampleid}.bqsr_indels.vcf \
  -O ${sampleid}/${sampleid}.post_recal_data.table

#input is both data tables, output is recalibration plots 
gatk AnalyzeCovariates \
  -before ${sampleid}/${sampleid}.recal_data.table \
  -after ${sampleid}/${sampleid}.post_recal_data.table \
  -plots ${sampleid}/${sampleid}.recalibration_plots.pdf

#input is reference and recal bam file, output is recal file of raw variants
gatk HaplotypeCaller \
  -R $genome_fa \
  -I ${sampleid}/${sampleid}.recal_reads.bam \
  -O ${sampleid}/${sampleid}.raw_variants_recal.vcf

#input is reference and recal variants file, output is recal snps file
gatk SelectVariants \
  -R $genome_fa \
  -select-type SNP \
  -V ${sampleid}/${sampleid}.raw_variants_recal.vcf \
  -O ${sampleid}/${sampleid}.raw_snps_recal.vcf

#inout is reference and raw variants file, output is recal indels file
gatk SelectVariants \
  -R $genome_fa \
  -select-type INDEL \
  -V ${sampleid}/${sampleid}.raw_variants.vcf \
  -O ${sampleid}/${sampleid}.raw_indels_recal.vcf

#input is reference and recal snps file, output is filtered snps file with parameters
gatk VariantFiltration \
  -R $genome_fa \
  -filter-name "QD_filter" -filter "QD < 2.0" \
  -filter-name "FS_filter" -filter "FS > 60.0" \
  -filter-name "MQ_filter" -filter "MQ < 40.0" \
  -filter-name "SOR_filter" -filter "SOR > 4.0" \
  -filter-name "MQRankSum_filter" -filter "MQRankSum < -12.5" \
  -filter-name "ReadPosRankSum_filter" -filter "ReadPosRankSum < -8.0" \
  -V ${sampleid}/${sampleid}.raw_snps_recal.vcf \
  -O ${sampleid}/${sampleid}.filtered_snps_final.vcf

#input is reference and recal indels file, output is filtered indels file with parameters
gatk VariantFiltration \
  -R $genome_fa \
  -V ${sampleid}/${sampleid}.raw_indels_recal.vcf \
  -O ${sampleid}/${sampleid}.filtered_indels_final.vcf \
  -filter-name "QD_filter" -filter "QD < 2.0" \
  -filter-name "FS_filter" -filter "FS > 200.0" \
  -filter-name "SOR_filter" -filter "SOR > 10.0"
