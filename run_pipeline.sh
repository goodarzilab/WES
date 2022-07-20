genomeDir=$1
genomeRef=$2 # hg38 or mm10
bamDir=$3
wildcard=$4

genome=${genomeDir}/${genomeRef}.fa
genome=`realpath $genome`

mkdir gatk-results
cd gatk-results

for bam in ../${bamDir}/${wildcard}; do
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
