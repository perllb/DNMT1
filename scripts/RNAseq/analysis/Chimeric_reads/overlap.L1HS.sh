## Get reads that overlap TE and NCBI gene 

#!/bin/bash

module load GCC/5.4.0-2.26  OpenMPI/1.10.3
module load BEDTools/2.26.0

file=$1
echo $file
sample=$(echo $file | cut -f2 -d"_" )
echo $sample


## Get reads overlapping L1HS

bedtools intersect -bed -split -wo -abam $file -b /projects/fs1/medpvb/no_backup/genomicData/hg38/RepeatMasker/L1/L1HS.hg38.uniqID.bed > ${sample}_L1HS.bed

## cut out info
awk ' { print $1,$2,$3,$16"-"$4,$5,$6,$7,$8,$9,$10,$11,$12 } ' OFS="\t" ${sample}_L1HS.bed  > ${sample}_L1HS.cut.bed

## Get reads overlapping NCBI
bedtools intersect -split -wo -a ${sample}_L1HS.cut.bed -b /projects/fs1/medpvb/no_backup/genomicData/hg38/NCBI_iGenomes/genes.gtf > ${sample}_L1HS_NCBIgenes.bed

## remove duplicates (unique in terms of reads)
sort -u -k1,4 ${sample}_L1HS_NCBIgenes.bed > tmp.txt
mv tmp.txt ${sample}_L1HS_NCBIgenes.bed

## Get only pairs where read is sense to NCBI
# if read 1, then it is the opposite strand
# if read 2, then it is correct in file
cut -f4 ${sample}_L1HS_NCBIgenes.bed | sed 's/\//\t/g' | cut -f2 > ${sample}_L1HS_NCBIgenes.READ.txt
paste ${sample}_L1HS_NCBIgenes.bed ${sample}_L1HS_NCBIgenes.READ.txt > ${sample}_L1HS_NCBIgenes.read.txt
rm ${sample}_L1HS_NCBIgenes.READ.txt

# if (read 1 AND opposite strand) OR (read 2 AND same strand)
awk -F $'\t' ' { if($23==1 && $6!=$19 || $23==2 && $6==$19) print $0 } ' OFS="\t"  ${sample}_L1HS_NCBIgenes.read.txt > ${sample}_L1HS_NCBIgenes.sameStrand.txt



## Get number of intersecting reads of each gene
cut -f21 ${sample}_L1HS_NCBIgenes.bed | sed 's/;/\t/g' | cut -f1 | sort | uniq -c
cut -f21 ${sample}_L1HS_NCBIgenes.sameStrand.txt | sed 's/;/\t/g' | cut -f1 | sort | uniq -c

