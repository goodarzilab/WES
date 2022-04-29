hg38='/sadra/goodarzilab/abe/People/Seema/hg38.fa'
# alignment
mkdir bam
for fq1 in fastq/*_1.fastq.gz; do
	base=`basename $fq1`
	fq2=${fq1/_1.fast/_2.fast};
	out=${base/_1.fastq.gz/};
	echo $out;
	bwa mem -t 16 $hg38 $fq1 $fq2 | samtools view -bS - &> bam/${out}.bam;
done
