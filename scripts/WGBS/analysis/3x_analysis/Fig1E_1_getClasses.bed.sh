## Get LINE SINE SVA LTR classes bed

cd "/Users/perludvik/Documents/bioinformatics/genomicData/hg38/Repeats"

awk ' $1 ~ /^chr/' hg38.fa.out.FamClass.bed | awk ' { if ($5 == "-") { print $1,$3,$2,$4,$5,$6 } else { print $0 } } ' OFS="\t" > tmp.txt
mv tmp.txt hg38.fa.out.FamClass.bed

grep LINE hg38.fa.out.FamClass.bed > hg38.fa.out.LINE.bed
grep SINE hg38.fa.out.FamClass.bed > hg38.fa.out.SINE.bed
grep LTR hg38.fa.out.FamClass.bed > hg38.fa.out.LTR.bed
grep SVA hg38.fa.out.FamClass.bed > hg38.fa.out.SVA.bed

cd "/Users/perludvik/Dropbox (MN)/Per/PhD/Projects/DNAmeth/hNES DNMT1 KO/2018.08.30_DNMT1L1_figures_scripts/WGBS/features/classes"

dir="/Users/perludvik/Documents/bioinformatics/genomicData/hg38/Repeats"
cp $dir/hg38.fa.out.SVA.bed .
cp $dir/hg38.fa.out.SINE.bed .
cp $dir/hg38.fa.out.LINE.bed .
cp $dir/hg38.fa.out.LTR.bed .
