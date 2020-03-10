cd "/Users/perludvik/Dropbox (MN)/Per/PhD/Projects/DNAmeth/hNES DNMT1 KO/2018.08.30_DNMT1L1_figures_scripts/NearbyGenes/deeptools_neargenes"

d=50000

mkdir $d
cd $d
mkdir output
mkdir output/"WGBS.rnaseq.k27ac.pol2"

chipdir="../../../ChIP-integration/chip.means"
featdir="../../../NearbyGenes/output/"$d

### ALL TEs
computeMatrix scale-regions -S $chipdir/WT_bismark_bt2_pe.deduplicated.bismark.bw \
  $chipdir/DNMT1_KO_bismark_bt2_pe.deduplicated.bismark.bw \
  $chipdir/DNMT1KO_hg38.unique_WT_REV.FWD.sum.bw \
  $chipdir/DNMT1KO_hg38.unique_KO_REV.FWD.sum.bw \
  $chipdir/chip_k27ac_WT_3.4.mean.bw \
  $chipdir/chip_k27ac_KO_1.2.mean.bw  \
  $chipdir/chip_pol_WT_7.8.mean.bw \
  $chipdir/chip_pol_KO_5.6.mean.bw \
  -R $featdir/L1HS.close.$d.bed \
  $featdir/L1PA2.close.$d.bed \
  $featdir/L1PA3.close.$d.bed \
  $featdir/HERV.close.$d.bed \
  $featdir/SVA.close.$d.bed \
  -a 10000 -b 10000 --regionBodyLength 10000 -bs 500 -p 4  \
  -out output/WGBS.rnaseq.k27ac.pol2/l1hs.pa2.pa3.herv.sva_nearGenes_WGBS_RNAseq_k27ac_polII_scaleRegion.mat


plotHeatmap -m output/WGBS.rnaseq.k27ac.pol2/l1hs.pa2.pa3.herv.sva_nearGenes_WGBS_RNAseq_k27ac_polII_scaleRegion.mat -o output/WGBS.rnaseq.k27ac.pol2/l1hs.pa2.pa3.herv.sva_nearGenes_WGBS_RNAseq_k27ac_polII_scaleRegion.pdf --missingDataColor 1  --sortUsingSamples 6   --zMin 0 --yMin 0 --yMax 100 100 200 200 12 12 17 1.4  --zMax 100 100 3 3 10 10 12 1   --colorMap Greys Greys Greens Greens Reds Reds Blues Blues --regionsLabel "L1HS" "L1PA2" "L1PA3" "HERV" "SVA" --samplesLabel "WGBS WT" "WGBS KO" "RNAseq_WT" "RNAseq_KO" "H3K27ac WT" "H3K27ac KO" "polII WT" "polII KO"  --startLabel TE --endLabel gene

plotProfile -m output/WGBS.rnaseq.k27ac.pol2/l1hs.pa2.pa3.herv.sva_nearGenes_WGBS_RNAseq_k27ac_polII_scaleRegion.mat -o output/WGBS.rnaseq.k27ac.pol2/l1hs.pa2.pa3.herv.sva_nearGenes_WGBS_RNAseq_k27ac_polII_scaleRegion.profile.pdf --yMin 0 --yMax 100 100 200 200 12 12 17 1.5  --colors black green red blue orange --regionsLabel "L1HS" "L1PA2" "L1PA3" "HERV" "SVA"  --samplesLabel "WGBS WT" "WGBS KO" "RNAseq_WT" "RNAseq_KO" "H3K27ac WT" "H3K27ac KO" "polII WT" "polII KO" --startLabel TE --endLabel gene

## L1s
computeMatrix scale-regions -S $chipdir/WT_bismark_bt2_pe.deduplicated.bismark.bw \
  $chipdir/DNMT1_KO_bismark_bt2_pe.deduplicated.bismark.bw \
  $chipdir/DNMT1KO_hg38.unique_WT_REV.FWD.sum.bw \
  $chipdir/DNMT1KO_hg38.unique_KO_REV.FWD.sum.bw \
  $chipdir/chip_k27ac_WT_3.4.mean.bw \
  $chipdir/chip_k27ac_KO_1.2.mean.bw  \
  $chipdir/chip_pol_WT_7.8.mean.bw \
  $chipdir/chip_pol_KO_5.6.mean.bw \
  -R $featdir/L1HS.close.$d.bed \
  $featdir/l1pa2.close.$d.bed \
  $featdir/l1pa2.close.$d.bed \
  -a 10000 -b 10000 --regionBodyLength 10000 -bs 500 -p 4  \
  -out output/WGBS.rnaseq.k27ac.pol2/l1hs.pa2.pa3_nearGenes_WGBS_RNAseq_k27ac_polII_scaleRegion.mat


plotHeatmap -m output/WGBS.rnaseq.k27ac.pol2/l1hs.pa2.pa3_nearGenes_WGBS_RNAseq_k27ac_polII_scaleRegion.mat -o output/WGBS.rnaseq.k27ac.pol2/l1hs.pa2.pa3_nearGenes_WGBS_RNAseq_k27ac_polII_scaleRegion.pdf --missingDataColor 1  --sortUsingSamples 6   --zMin 0 --yMin 0 --yMax 100 100 200 200 12 12 17 1.4  --zMax 100 100 3 3 10 10 12 1   --colorMap Greys Greys Greens Greens Reds Reds Blues Blues --regionsLabel "L1HS" "L1PA2" "L1PA3" --samplesLabel "WGBS WT" "WGBS KO" "RNAseq_WT" "RNAseq_KO" "H3K27ac WT" "H3K27ac KO" "polII WT" "polII KO"  --startLabel TE --endLabel gene

plotProfile -m output/WGBS.rnaseq.k27ac.pol2/l1hs.pa2.pa3_nearGenes_WGBS_RNAseq_k27ac_polII_scaleRegion.mat -o output/WGBS.rnaseq.k27ac.pol2/l1hs.pa2.pa3_nearGenes_WGBS_RNAseq_k27ac_polII_scaleRegion.profile.pdf --yMin 0 --yMax 100 100 200 200 12 12 17 1.5  --colors black green red blue orange --regionsLabel "L1HS" "L1PA2" "L1PA3"  --samplesLabel "WGBS WT" "WGBS KO" "RNAseq_WT" "RNAseq_KO" "H3K27ac WT" "H3K27ac KO" "polII WT" "polII KO" --startLabel TE --endLabel gene


## HERV SVA
computeMatrix scale-regions -S $chipdir/WT_bismark_bt2_pe.deduplicated.bismark.bw \
  $chipdir/DNMT1_KO_bismark_bt2_pe.deduplicated.bismark.bw \
  $chipdir/DNMT1KO_hg38.unique_WT_REV.FWD.sum.bw \
  $chipdir/DNMT1KO_hg38.unique_KO_REV.FWD.sum.bw \
  $chipdir/chip_k27ac_WT_3.4.mean.bw \
  $chipdir/chip_k27ac_KO_1.2.mean.bw  \
  $chipdir/chip_pol_WT_7.8.mean.bw \
  $chipdir/chip_pol_KO_5.6.mean.bw \
  -R $featdir/HERV.close.$d.bed \
  $featdir/SVA.close.$d.bed \
  -a 10000 -b 10000 --regionBodyLength 10000 -bs 500 -p 4  \
  -out output/WGBS.rnaseq.k27ac.pol2/herv.sva_nearGenes_WGBS_RNAseq_k27ac_polII_scaleRegion.mat


plotHeatmap -m output/WGBS.rnaseq.k27ac.pol2/herv.sva_nearGenes_WGBS_RNAseq_k27ac_polII_scaleRegion.mat -o output/WGBS.rnaseq.k27ac.pol2/herv.sva_nearGenes_WGBS_RNAseq_k27ac_polII_scaleRegion.pdf --missingDataColor 1  --sortUsingSamples 6   --zMin 0 --yMin 0 --yMax 100 100 200 200 12 12 17 1.4  --zMax 100 100 3 3 10 10 12 1   --colorMap Greys Greys Greens Greens Reds Reds Blues Blues --regionsLabel "HERV" "SVA" --samplesLabel "WGBS WT" "WGBS KO" "RNAseq_WT" "RNAseq_KO" "H3K27ac WT" "H3K27ac KO" "polII WT" "polII KO"  --startLabel TE --endLabel gene

plotProfile -m output/WGBS.rnaseq.k27ac.pol2/herv.sva_nearGenes_WGBS_RNAseq_k27ac_polII_scaleRegion.mat -o output/WGBS.rnaseq.k27ac.pol2/herv.sva_nearGenes_WGBS_RNAseq_k27ac_polII_scaleRegion.profile.pdf --yMin 0 --yMax 100 100 200 200 12 12 17 1.5  --colors black green red blue orange --regionsLabel "HERV" "SVA"  --samplesLabel "WGBS WT" "WGBS KO" "RNAseq_WT" "RNAseq_KO" "H3K27ac WT" "H3K27ac KO" "polII WT" "polII KO" --startLabel TE --endLabel gene
