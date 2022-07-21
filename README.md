# WES
A Whole Exome Sequencing (WES) analysis workflow

# Quick start
## [Set up tools](tools/README.md)
Follow above link to install and set up all related tools. 
## Alignment

```bash 
nohup bash $WES/alignment.sh ../genome mm10 fastq-2 align_mouse 5 &> align_mouse.out
```
```bash 
nohup bash $WES/alignment.sh ../genome mm10 fastq-2 align_mouse 5 &> align_mouse.out
```

## Xenofilter [optional]
```bash
nohup Rscript $WES/xenofilter.R "align_human" "*.hg38" "align_mouse" "*.mm10" 30 &> xenofilter.out & 
```
```bash
mv Filtered_bams/XenofilteR.log xenofilter.log

mkdir align_xenofilter
for f in Filtered_bams/*; 
	do b=`basename $f`; o=${b/_Filtered.bam/xenofilter.bam}; 
	mv -v $f align_xenofilter/$o; 
done

rm -rv Filtered_bams/
```

## Run GATK pipeline
```bash
nohup bash $WES/run_pipeline.sh ../genome hg38 bam "*" > run_pipeline.out
```

___ 
This pipeline developed by Sonakshi Bhalla (bhallasonakshi@gmail.com) and Abolfazl Arab (abarbiology@gmail.com) at Goodarzi Lab, UCSF – 2021-2022. 
