cd "/Users/perludvik/Dropbox (MN)/Per/PhD/Projects/DNAmeth/hNES DNMT1 KO/2018.08.30_DNMT1L1_figures_scripts/WGBS/data"

## make bedGraph from coverage
cut -f1,2,3,4 WT_bismark_bt2_pe.deduplicated.bismark.cov  > WT_bismark_bt2_pe.deduplicated.bismark.bedGraph
cut -f1,2,3,4 DNMT1_KO_bismark_bt2_pe.deduplicated.bismark.cov  > DNMT1_KO_bismark_bt2_pe.deduplicated.bismark.bedGraph

## hg38 chrom sizes > bed
awk ' { print $1,"0",$2 } ' OFS="\t" hg38.chrom.sizes.primary.txt  > hg38.chrom.sizes.primary.bed

## Chop hg38 into 1000bp bins..
bedops --chop 1000 hg38.chrom.sizes.primary.bed > hg38.chrom.sizes.primary.1000bin.bed

## make bed format
awk ' { print $1,$2-1,$3,".",$4 } ' OFS="\t" WT_bismark_bt2_pe.deduplicated.bismark.bedGraph > WT_bismark_bt2_pe.deduplicated.bed
awk ' { print $1,$2-1,$3,".",$4 } ' OFS="\t" DNMT1_KO_bismark_bt2_pe.deduplicated.bismark.bedGraph > DNMT1_KO_bismark_bt2_pe.deduplicated.bed

## sort bed
sortBed -i WT_bismark_bt2_pe.deduplicated.bed > WT_bismark_bt2_pe.deduplicated.sorted.bed

sortBed -i DNMT1_KO_bismark_bt2_pe.deduplicated.bed > DNMT1_KO_bismark_bt2_pe.deduplicated.sorted.bed

## get mean of 1000bins
bedmap --mean --delim "\t" --echo hg38.chrom.sizes.primary.1000bin.bed WT_bismark_bt2_pe.deduplicated.sorted.bed > WT_bismark_1000bin.mean_bed
bedmap --mean --delim "\t" --echo hg38.chrom.sizes.primary.1000bin.bed DNMT1_KO_bismark_bt2_pe.deduplicated.sorted.bed > DNMT1KO_bismark_1000bin.mean_bed

cut -f1  WT_bismark_1000bin.mean_bed | grep -v NAN >  WT_bismark_1000bin.mean_only.bed
cut -f1  DNMT1KO_bismark_1000bin.mean_bed | grep -v NAN >  DNMT1KO_bismark_1000bin.mean_only.bed


#### 3x coverage
mkdir 3x

## Get only 3x coverage bases

awk ' { if ( $6+$5 > 2 ) { print $0 } } ' DNMT1_KO_bismark_bt2_pe.deduplicated.bismark.cov > 3x/DNMT1_KO_bismark_bt2_pe.deduplicated.bismark.3x.cov
awk ' { if ( $6+$5 > 2 ) { print $0 } } ' WT_bismark_bt2_pe.deduplicated.bismark.cov > 3x/WT_bismark_bt2_pe.deduplicated.bismark.3x.cov

# make bedGraph
cut -f1,2,3,4 3x/DNMT1_KO_bismark_bt2_pe.deduplicated.bismark.3x.cov > 3x/DNMT1_KO_bismark_bt2_pe.deduplicated.bismark.3x.bedGraph
cut -f1,2,3,4 3x/WT_bismark_bt2_pe.deduplicated.bismark.3x.cov > 3x/WT_bismark_bt2_pe.deduplicated.bismark.3x.bedGraph

# make bed
awk ' { print $1,$2-1,$3,".",$4 } ' OFS="\t" 3x/DNMT1_KO_bismark_bt2_pe.deduplicated.bismark.3x.bedGraph > 3x/DNMT1_KO_bismark_bt2_pe.deduplicated.bismark.3x.bed
awk ' { print $1,$2-1,$3,".",$4 } ' OFS="\t" 3x/WT_bismark_bt2_pe.deduplicated.bismark.3x.bedGraph > 3x/WT_bismark_bt2_pe.deduplicated.bismark.3x.bed

# sort bed
sortBed -i 3x/DNMT1_KO_bismark_bt2_pe.deduplicated.bismark.3x.bed > 3x/tmp.txt
mv 3x/tmp.txt 3x/DNMT1_KO_bismark_bt2_pe.deduplicated.bismark.3x.bed

sortBed -i 3x/WT_bismark_bt2_pe.deduplicated.bismark.3x.bed > 3x/tmp.txt
mv 3x/tmp.txt 3x/WT_bismark_bt2_pe.deduplicated.bismark.3x.bed

## get mean of 1000bins
bedmap --mean --delim "\t" --echo hg38.chrom.sizes.primary.1000bin.bed 3x/DNMT1_KO_bismark_bt2_pe.deduplicated.bismark.3x.bed > 3x/DNMT1ko_bismark_1000bin.mean.3x.bed
bedmap --mean --delim "\t" --echo hg38.chrom.sizes.primary.1000bin.bed 3x/WT_bismark_bt2_pe.deduplicated.bismark.3x.bed > 3x/WT_bismark_1000bin.mean.3x.bed

cut -f1  3x/DNMT1ko_bismark_1000bin.mean.3x.bed | grep -v NAN >  3x/DNMT1ko_bismark_1000bin.mean_only.3x.bed
cut -f1  3x/WT_bismark_1000bin.mean.3x.bed | grep -v NAN >  3x/WT_bismark_1000bin.mean_only.3x.bed
