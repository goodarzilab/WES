### Tumor-Only Variant Call Command-Line Parameters 
## 1. Generate OXOG metrics:
gatk CollectSequencingArtifactMetrics -I 31.processed_Filtered.dups.bam -O 31.output.txt -R hg38.fa

## 2. Generate pileup summaries on tumor sample:
gatk GetPileupSummaries \
-I Tumor_Sample_Alignment.bam \
-O <job_identifier>.targeted_sequencing.table \
-V af-only-gnomad-common-biallelic.grch38.main.vcf.gz \ # Germline reference from gnomad
-L intervals.bed \ ## Only chr1-22 + XYM
-R GRCh38.d1.vd1.fa

## 3. Calculate contamination on tumor sample
gatk CalculateContamination \
-I <job_identifier>.targeted_sequencing.table \ # From step 2
-O <job_identifier>.targeted_sequencing.contamination.table

## 4. Find tumor sample name from BAM
gatk GetSampleName \
-I Tumor_Sample_Alignment.bam \
-O <job_identifier>.targeted_sequencing.sample_name

## 5. Run MuTect2 using only tumor sample on chromosome level (25 commands with different intervals)
gatk Mutect2 \
-R GRCh38.d1.vd1.fa \
-L chr4:1-190214555 \ # Specify chromosome
-I Tumor_Sample_Alignment.bam \
-O 3.mt2.vcf \
-tumor <tumor_sample_name> \ # From step 4
--af-of-alleles-not-in-resource 2.5e-06 \
--germline-resource af-only-gnomad.hg38.vcf.gz \ # Germline reference from gnomad
-pon gatk4_mutect2_4136_pon.vcf.gz # New panel of normal created by 4136 TCGA curated normal samples, using GATK4

# After this step, all chromosome level VCFs are merged into one.

## 6. Sort VCF with Picard
gatk SortVcf \
SEQUENCE_DICTIONARY=GRCh38.d1.vd1.dict \
OUTPUT=<job_identifier>.targeted_sequencing.mutect2.tumor_only.sorted.vcf.gz \
I=merged_multi_gatk4_mutect2_tumor_only_calling.vcf \ # From step 5
CREATE_INDEX=true

## 7. Filter variant calls from MuTect
gatk FilterMutectCalls \
-O <job_identifier>.targeted_sequencing.mutect2.tumor_only.contFiltered.vcf.gz \
-V <job_identifier>.targeted_sequencing.mutect2.tumor_only.sorted.vcf.gz \ # From step 6
--contamination-table <job_identifier>.targeted_sequencing.contamination.table \ # From step 3
-L intervals.bed

## 8. Filter variants by orientation bias
gatk FilterByOrientationBias \
-O <job_identifier>.targeted_sequencing.tumor_only.gatk4_mutect2.raw_somatic_mutation.vcf.gz \ # final output
-P <job_identifier>.pre_adapter_detail_metrics.txt \ # From step 1
-V <job_identifier>.targeted_sequencing.mutect2.tumor_only.contFiltered.vcf.gz \ # From step 7
-L intervals.bed \
-R GRCh38.d1.vd1.fa \
-AM G/T \
-AM C/T
