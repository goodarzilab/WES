export bam_file=__
read_counts=${bam_file/.bam/.read_counts.tsv}


gatk PreprocessIntervals \
        -R /rumi/shams/genomes/hg38/hg38.fa \
        --padding 0 \
        -imr OVERLAPPING_ONLY \
        -O hg38.preprocessed.interval_list
        
        
gatk AnnotateIntervals \
        -L hg38.preprocessed.interval_list \
        -R /rumi/shams/genomes/hg38/hg38.fa \
        -imr OVERLAPPING_ONLY \
        -O hg38.annotated.tsv

gatk CollectReadCounts \
        -L hg38.preprocessed.interval_list \
        -R /rumi/shams/genomes/hg38/hg38.fa \
        -imr OVERLAPPING_ONLY \
        -I $bam_file \
        --format TSV \
        -O $read_counts
        
        
gatk FilterIntervals \
        -L hg38.preprocessed.interval_list \
        --annotated-intervals hg38.annotated.tsv \
        -I $read_counts
        -imr OVERLAPPING_ONLY \
        -O hg38.sample.gc.filtered.interval_list
        
 gatk DetermineGermlineContigPloidy \
        -L hg38.sample.gc.filtered.interval_list \
        --interval-merging-rule OVERLAPPING_ONLY \
        -I $read_counts
        --contig-ploidy-priors hg38_contig_ploidy_priors.tsv \
        --output . \
        --output-prefix ploidy \
        --verbosity DEBUG
        
gatk DetermineGermlineContigPloidy \
        --model cohort-23wgs-20190213-contig-ploidy-model \
        -I $read_counts \
        -O . \
        --output-prefix ploidy-case \
        --verbosity DEBUG
        
gatk GermlineCNVCaller \
        --run-mode COHORT \
        -L scatter-sm/twelve_1of2.interval_list \
        -I $read_counts 
        --contig-ploidy-calls ploidy-calls \
        --annotated-intervals twelveregions.annotated.tsv \
        --interval-merging-rule OVERLAPPING_ONLY \
        --output cohort24-twelve \
        --output-prefix cohort24-twelve_1of2 \
        --verbosity DEBUG
