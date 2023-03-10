#!/bin/bash

#change directory to the folder and run as: sh S5_microhaplotyper_pipeline.sh 

for f in *.fastq; do bwa mem VISAGE_ET_microhaplotyper_refgenome.fa $f >  ${f/.fastq/_microhaplotyper.sam}; done

for f in *.sam; do  samtools view -h -L VISAGE_ET_microhaplotyper_targets.bed $f -q 30 |  awk '/^@/ || length($10) >= 80' | samtools sort -o ${f/_microhaplotyper.sam/_filtered.sam}; done

for f in *_filtered.sam; do  samtools view -u $f | samtools sort -o ${f/_filtered.sam/_microhaplotyper.bam}; done

ls -1 | grep -E "_microhaplotyper.bam" >bamlist.txt

samtools mpileup --redo-BAQ --min-BQ 30 --per-sample-mF -l VISAGE_ET_microhaplotyper_hotspots.bed --output-tags DP,AD -f VISAGE_ET_microhaplotyper_refgenome.fa --BCF -b bamlist.txt | bcftools call --multiallelic-caller --skip-variants indels -Ov -Ob| bcftools norm -Ov --check-ref w -f VISAGE_ET_microhaplotyper_refgenome.fa > CALLS.vcf

ls -1 | grep -E "_filtered.sam" | while read line; do printf $line'\t'$line'\t'NA'\n'; done >label.txt


find . -name \*_microhaplotyper.bam -delete
find . -name \*_microhaplotyper.sam -delete
find . -name \*bamlist.txt -delete
