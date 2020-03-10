################################################################################################
#### This script compute the mean CpG ratios (mCpG/CpG) over given features.####################
################################################################################################

### L1 promoters ###
cd "/Users/perludvik/Dropbox (MN)/Per/PhD/Projects/DNAmeth/hNES DNMT1 KO/2018.08.30_DNMT1L1_figures_scripts/WGBS"

sort-bed features/L1/L1HS.FL6000.promoter.bed > tmp.txt
mv tmp.txt features/L1/L1HS.FL6000.promoter.bed

## get mean of each feature
bedmap --mean --count --delim "\t" --echo features/L1/L1HS.FL6000.promoter.bed data/3x/DNMT1_KO_bismark_bt2_pe.deduplicated.bismark.3x.bed > analysis/compFeatures/3x/L1/DNMT1ko_L1HS.FL6000mean.promoter.3x.bed
cut -f1 analysis/compFeatures/3x/L1/DNMT1ko_L1HS.FL6000mean.promoter.3x.bed | grep -v NAN > analysis/compFeatures/3x/L1/DNMT1ko_L1HS.FL6000.promoter_mean.3x_rmNAN.bed
cut -f2 analysis/compFeatures/3x/L1/DNMT1ko_L1HS.FL6000mean.promoter.3x.bed > analysis/compFeatures/3x/L1/DNMT1ko_L1HS.FL6000.promoter.3x_count.bed

bedmap --mean --count --delim "\t" --echo features/L1/L1HS.FL6000.promoter.bed data/3x/WT_bismark_bt2_pe.deduplicated.bismark.3x.bed > analysis/compFeatures/3x/L1/WT_L1HS.FL6000.promoter.3x.bed
cut -f1  analysis/compFeatures/3x/L1/WT_L1HS.FL6000.promoter.3x.bed | grep -v NAN > analysis/compFeatures/3x/L1/WT_L1HS.FL6000.promoter_mean.3x_rmNAN.bed
cut -f2  analysis/compFeatures/3x/L1/WT_L1HS.FL6000.promoter.3x.bed  > analysis/compFeatures/3x/L1/WT_L1HS.FL6000.promoter.3x_count.bed

### L1 FL6000 ###
for l1 in L1HS L1PA2 L1PA3 L1PA4
do
  echo $l1

  sort-bed features/L1/$l1.FL6000.bed > tmp.txt
  mv tmp.txt features/L1/$l1.FL6000.bed

## KO
  bedmap --mean --count --delim "\t" --echo features/L1/$l1.FL6000.bed data/3x/DNMT1_KO_bismark_bt2_pe.deduplicated.bismark.3x.bed > analysis/compFeatures/3x/outputTE/beds/DNMT1ko_$l1.FL6000mean.3x.bed
  cut -f1 analysis/compFeatures/3x/outputTE/beds/DNMT1ko_$l1.FL6000mean.3x.bed | grep -v NAN > analysis/compFeatures/3x/outputTE/beds/DNMT1ko_$l1.FL6000mean.3x_rmNAN.bed
  cut -f2 analysis/compFeatures/3x/outputTE/beds/DNMT1ko_$l1.FL6000mean.3x.bed  > analysis/compFeatures/3x/outputTE/beds/DNMT1ko_$l1.FL6000mean.3x_count.bed

## WT
  bedmap --mean --count --delim "\t" --echo features/L1/$l1.FL6000.bed data/3x/WT_bismark_bt2_pe.deduplicated.bismark.3x.bed > analysis/compFeatures/3x/outputTE/beds/WT_$l1.FL6000mean.3x.bed
  cut -f1 analysis/compFeatures/3x/outputTE/beds/WT_$l1.FL6000mean.3x.bed | grep -v NAN > analysis/compFeatures/3x/outputTE/beds/WT_$l1.FL6000mean.3x_rmNAN.bed
  cut -f2 analysis/compFeatures/3x/outputTE/beds/WT_$l1.FL6000mean.3x.bed  > analysis/compFeatures/3x/outputTE/beds/WT_$l1.FL6000mean.3x_count.bed
done

### TE classes ###
for class in LINE SINE LTR SVA
do
  echo $class

  sort-bed features/classes/hg38.fa.out.$class.bed > tmp.txt
  mv tmp.txt features/classes/hg38.fa.out.$class.bed

  bedmap --mean --count --delim "\t" --echo features/classes/hg38.fa.out.$class.bed data/3x/DNMT1_KO_bismark_bt2_pe.deduplicated.bismark.3x.bed > analysis/compFeatures/3x/outputTE/beds/DNMT1ko_$class.mean.3x.bed
  cut -f1 analysis/compFeatures/3x/outputTE/beds/DNMT1ko_$class.mean.3x.bed | grep -v NAN > analysis/compFeatures/3x/outputTE/beds/DNMT1ko_$class.mean.3x_rmNAN.bed
  cut -f2 analysis/compFeatures/3x/outputTE/beds/DNMT1ko_$class.mean.3x.bed > analysis/compFeatures/3x/outputTE/beds/DNMT1ko_$class.mean.3x_count.bed

  bedmap --mean --count --delim "\t" --echo features/classes/hg38.fa.out.$class.bed data/3x/WT_bismark_bt2_pe.deduplicated.bismark.3x.bed > analysis/compFeatures/3x/outputTE/beds/WT_$class.mean.3x.bed
  cut -f1 analysis/compFeatures/3x/outputTE/beds/WT_$class.mean.3x.bed | grep -v NAN > analysis/compFeatures/3x/outputTE/beds/WT_$class.mean.3x_rmNAN.bed
  cut -f2 analysis/compFeatures/3x/outputTE/beds/WT_$class.mean.3x.bed > analysis/compFeatures/3x/outputTE/beds/WT_$class.mean.3x_count.bed
done

### FLI L1 ###
class=FLI-L1

sort-bed features/L1/l1base2.bed > tmp.txt
mv tmp.txt features/L1/l1base2.bed

bedmap --mean --count --delim "\t" --echo features/L1/l1base2.bed data/3x/DNMT1_KO_bismark_bt2_pe.deduplicated.bismark.3x.bed > analysis/compFeatures/3x/outputTE/beds/DNMT1ko_$class.mean.3x.bed
cut -f1 analysis/compFeatures/3x/outputTE/beds/DNMT1ko_$class.mean.3x.bed | grep -v NAN > analysis/compFeatures/3x/outputTE/beds/DNMT1ko_$class.mean.3x_rmNAN.bed
cut -f2 analysis/compFeatures/3x/outputTE/beds/DNMT1ko_$class.mean.3x.bed > analysis/compFeatures/3x/outputTE/beds/DNMT1ko_$class.3x_count.bed

bedmap --mean --count --delim "\t" --echo features/L1/l1base2.bed data/3x/WT_bismark_bt2_pe.deduplicated.bismark.3x.bed > analysis/compFeatures/3x/outputTE/beds/WT_$class.mean.3x.bed
cut -f1 analysis/compFeatures/3x/outputTE/beds/WT_$class.mean.3x.bed | grep -v NAN > analysis/compFeatures/3x/outputTE/beds/WT_$class.mean.3x_rmNAN.bed
cut -f2 analysis/compFeatures/3x/outputTE/beds/WT_$class.mean.3x.bed > analysis/compFeatures/3x/outputTE/beds/WT_$class.3x_count.bed
