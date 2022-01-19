sampleID=$1

echo $sampleID;

java -jar \
	${WES}tools/snpEff/snpEff.jar \
	-v GRCh38.99 \
	${sampleID}/${sampleID}.filtered_snps_final.vcf &> ${sampleID}/${sampleID}.filtered_snps_final.ann.vcf;
wait;

mv -v snpEff_genes.txt ${sampleID}/${sampleID}.snpEff_genes.txt
mv -v snpEff_summary.html ${sampleID}/${sampleID}.snpEff_summary.html
