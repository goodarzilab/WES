hg38=/path/to/reference/genome
# Sets the path to human genome for reference

# bwa aligner is used to index reference genome
bwa index hg38

# script to take paired end fastq files and convert to sam files
 
for f1 in *1.fastq.gz; do
f2=${f1/1.fastq.gz/2.fastq.gz}
out=${f1/_1.fastq.gz/.human.sam}
echo "bwa mem -t 10 $hg38 $f1 $f2 > $out"
bwa mem -t 10 $hg38 $f1 $f2 > $out
done

# same process is done for mouse genome
mm10=/path/to/mouse/genome

bwa index mm10

# aligment to mouse genome is done for same fastq files
for f1 in *1.fastq.gz; do
f2=${f1/1.fastq.gz/2.fastq.gz}
out=${f1/_1.fastq.gz/.mouse.sam}
echo "bwa mem -t 10 $mm10 $f1 $f2 > $out"
bwa mem -t 10 $mm10 $f1 $f2 > $out
done

# sam files are convert to bam files to reduce file size using samtools
for f in *.sam; do
out=${f/.sam/.bam}
samtools view -@ 10 -b $f > $out
done

# bam files are processed and indexed before xenofiltering
for f in *.bam; do
out=${f/.bam/.srt.bam};
echo "samtools sort -o $out $f";
samtools sort -@ 10 -m 4G -o $out $f
samtools index $out;
done
