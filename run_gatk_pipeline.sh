genome="/rumi/shams/genomes/hg38/hg38.fa"

mkdir gatk-results
cd gatk-results

for bam in ../bam/*; do
sample=${bam/.bam/}
sample=`basename $sample`
echo $sample
echo `date`
echo "Start!"
bash $WES/gatk-part1.sh $bam $sample $genome &> ${sample}.log;
wait;
echo `date`
echo "Part 1, done!"

bash $WES/gatk-part2.sh $sample $genome &>> ${sample}.log;
wait;
echo `date`
echo "Part 2, done!"

bash $WES/gatk-part3.sh $sample &>> ${sample}.log;
wait;
echo `date`
echo "Part 3, done!"
echo "------------------------------------"

done
