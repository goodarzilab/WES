genomeDir=$1
genomeRef=$2 # hg38 or mm10
fastqDir=$3
bamDir=$4
JOBS=$5

mkdir $bamDir

for fq1 in $fastqDir/*_1.fastq.gz; do
	base=`basename $fq1`
	fq2=${fq1/1.fastq.gz/_2.fastq.gz}; # the R1 pattern might be different!
	out=${base/_1.fastq.gz/};
	echo $out;
	bwa mem -t $JOBS ${genomeDir}/${genomeRef}.fa $fq1 $fq2 | samtools view -bS - &> bam/${out}.${genomeRef}.bam;
	#cd $bamDir
	#samtools sort -@ $JOBS -m 4G -o ${out}.${genomeRef}.srt.bam ${out}.${genomeRef}.bam
	#samtools index ${out}.${genomeRef}.srt.bam;
	#cd ../
done
