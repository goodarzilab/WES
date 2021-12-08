PDIR=$1
bamDIR=$2

cd $PDIR

for INbam in ${bamDIR}/*.bam; do
    fq_base=`basename $fq_file`;
    
    sample_id=${INbam/.bam/}
    
    mkdir $sample_id

    echo '--------------' $sample_id '--------------'
    # If mean read length is less than 70bp: DO SOMETING ELSE!!
    # bwa mem -t 8 -T 0 -R <read_group> <reference> <fastq_1.fq.gz> <fastq_2.fq.gz> | samtools view -Shb -o <output.bam> -
    # bam file from xenofiltering output should be the initial bam for this script

    # 1. correct bam header and write it to a temp bam files 
    samtools view -@ 10 -Shb -o ${sample_id}/{sample_id}.temp1.bam $INbam

    # 2. sort bam file 
    gatk SortSam CREATE_INDEX=true INPUT=${sample_id}/{sample_id}.temp1.bam OUTPUT=${sample_id}/{sample_id}.sort.bam SORT_ORDER=coordinate VALIDATION_STRINGENCY=STRICT

    # 3. merge bam file 
    gatk MergeSamFiles ASSUME_SORTED=false CREATE_INDEX=true MERGE_SEQUENCE_DICTIONARIES=false SORT_ORDER=coordinate USE_THREADING=true VALIDATION_STRINGENCY=STRICT INPUT=${sample_id}/{sample_id}.sort.bam OUTPUT=${sample_id}/{sample_id}.merge.bam

    # 4. MARK DUPLICATES 
    gatk MarkDuplicates CREATE_INDEX=true INPUT=${sample_id}/{sample_id}.merge.bam OUTPUT=${sample_id}/{sample_id}.dups.bam METRICS_FILE=${sample_id}/{sample_id}.dups.mtx VALIDATION_STRINGENCY=STRICT

    # 5. Add headers to bam file
    gatk AddOrReplaceReadGroups I=$${sample_id}/{sample_id}.dups.bam O=${sample_id}/{sample_id}.out.bam RGID=4 RGLB=lib1 RGPL=illumina RGPU=unit1 RGSM=20
    
done

multiqc ${sample_id}/ -n mutiqc-${sample_id}

# remove temp files
rm -v */*.temp1.bam
rm -v */*.sort.bam
rm -v */*.merge.bam
rm -v */*.dups.bam
