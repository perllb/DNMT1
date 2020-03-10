#### Plot L1 FL6000, SVA (1000), PJ HERV (score300) ####

### WT: K27me3, k9me3, trim28 ###

cd "/Users/perludvik/Dropbox (MN)/Per/PhD/Projects/DNAmeth/hNES DNMT1 KO/2018.08.30_DNMT1L1_figures_scripts/ChIP-integration/"

##  L1 HS,PA2-PA7 ##

computeMatrix scale-regions -S chip.means/trim28.log2ratioIN.bw \
  chip.means/h3k9me3.log2ratioIN.bw \
  chip.means/chip_k27m3_WT_11.12.mean.bw \
  -R L1/features/l1.beds/L1HS.FL6000.bed \
  L1/features/l1.beds/L1PA2.FL6000.bed \
  L1/features/l1.beds/L1PA3.FL6000.bed \
  L1/features/l1.beds/L1PA4.FL6000.bed \
  L1/features/l1.beds/L1PA5.FL6000.bed \
  L1/features/l1.beds/L1PA6.FL6000.bed \
  L1/features/l1.beds/L1PA7.FL6000.bed \
  -a 5000 -b 5000 --regionBodyLength 3000 -bs 200 -p 8 \
  -out finalGroups/output/k27me3.k9me3.trim28/l1hs.pa7_wt_k27me3.k9me3.trim28_scaleRegion.mat

plotHeatmap -m finalGroups/output/k27me3.k9me3.trim28/l1hs.pa7_wt_k27me3.k9me3.trim28_scaleRegion.mat -o finalGroups/output/k27me3.k9me3.trim28/l1hs.pa7_wt_k27me3.k9me3.trim28_scaleRegion.pdf --missingDataColor 1  --sortUsingSamples 1   --zMin 0 --yMin -0.3 --yMax 1.2 2 15 --zMax 1 1.9 15   --colorMap Blues Greens Greys --regionsLabel "L1HS" "L1PA2" "L1PA3" "L1PA4" "L1PA5" "L1PA6" "L1PA7"  --samplesLabel "TRIM28 WT" "H3K9me3 WT" "H3K27me3 WT"

plotProfile -m finalGroups/output/k27me3.k9me3.trim28/l1hs.pa7_wt_k27me3.k9me3.trim28_scaleRegion.mat -o finalGroups/output/k27me3.k9me3.trim28/l1hs.pa7_wt_k27me3.k9me3.trim28_scaleRegion.profile.pdf --yMin -0.3 --yMax 1.2 2 15  --colors purple red blue green yellow orange grey black --regionsLabel "L1HS" "L1PA2" "L1PA3" "L1PA4" "L1PA5" "L1PA6" "L1PA7"  --samplesLabel "TRIM28 WT" "H3K9me3 WT" "H3K27me3 WT"


## PJ HERV - 300
# WGBS + RNA seq + K27ac + PolII with HERV
computeMatrix scale-regions -S chip.means/trim28.log2ratioIN.bw \
  chip.means/h3k9me3.log2ratioIN.bw \
  chip.means/chip_k27m3_WT_11.12.mean.bw \
  -R HERV/features/hg38.HERV.score300.bed \
  -a 5000 -b 5000 --regionBodyLength 5000 -bs 200 -p 8 \
  -out finalGroups/output/k27me3.k9me3.trim28/HERV300_wt_k27me3.k9me3.trim28_scaleRegion.mat

plotHeatmap -m finalGroups/output/k27me3.k9me3.trim28/HERV300_wt_k27me3.k9me3.trim28_scaleRegion.mat -o finalGroups/output/k27me3.k9me3.trim28/HERV300_wt_k27me3.k9me3.trim28_scaleRegion.pdf --missingDataColor 1  --sortUsingSamples 1   --zMin 0 --yMin -0.3 --yMax 1.2 2 15 --zMax 1 1.9 15   --colorMap Blues Greens Greys  --regionsLabel "HERV"  --samplesLabel "TRIM28 WT" "H3K9me3 WT" "H3K27me3 WT"

plotProfile -m finalGroups/output/k27me3.k9me3.trim28/HERV300_wt_k27me3.k9me3.trim28_scaleRegion.mat -o finalGroups/output/k27me3.k9me3.trim28/HERV300_wt_k27me3.k9me3.trim28_scaleRegion.profile.pdf --yMin -0.3 --yMax 1.2 2 15  --colors black red blue --regionsLabel "HERV"  --samplesLabel "TRIM28 WT" "H3K9me3 WT" "H3K27me3 WT"

## SVA 1000bp

computeMatrix scale-regions -S chip.means/trim28.log2ratioIN.bw \
  chip.means/h3k9me3.log2ratioIN.bw \
  chip.means/chip_k27m3_WT_11.12.mean.bw \
  -R SVA/features/SVA_F.hg38_FL1000.bed \
  SVA/features/SVA_E.hg38_FL1000.bed \
  SVA/features/SVA_D.hg38_FL1000.bed \
  SVA/features/SVA_C.hg38_FL1000.bed \
  SVA/features/SVA_B.hg38_FL1000.bed \
  SVA/features/SVA_A.hg38_FL1000.bed \
  -a 500 -b 500 --regionBodyLength 500 -bs 10 -p 8 \
  -out finalGroups/output/k27me3.k9me3.trim28/SVA.FL1000_wt_k27me3.k9me3.trim28_scaleRegion.mat

plotHeatmap -m finalGroups/output/k27me3.k9me3.trim28/SVA.FL1000_wt_k27me3.k9me3.trim28_scaleRegion.mat -o finalGroups/output/k27me3.k9me3.trim28/SVA.FL1000_wt_k27me3.k9me3.trim28_scaleRegion.pdf --missingDataColor 1  --sortUsingSamples 1   --zMin 0 --yMin -0.3 --yMax 2 4 15 --zMax 2 4 15   --colorMap Blues Greens Greys  --regionsLabel "SVA_F" "SVA_E" "SVA_D" "SVA_C" "SVA_B" "SVA_A" --samplesLabel "TRIM28 WT" "H3K9me3 WT" "H3K27me3 WT"

plotProfile -m finalGroups/output/k27me3.k9me3.trim28/SVA.FL1000_wt_k27me3.k9me3.trim28_scaleRegion.mat -o finalGroups/output/k27me3.k9me3.trim28/SVA.FL1000_wt_k27me3.k9me3.trim28_scaleRegion.profile.pdf --yMin -0.3 --yMax 1.2 2 15  --colors purple red blue green yellow orange grey black --regionsLabel "SVA_F" "SVA_E" "SVA_D" "SVA_C" "SVA_B" "SVA_A" --samplesLabel "TRIM28 WT" "H3K9me3 WT" "H3K27me3 WT"




##  FLI -LA

computeMatrix scale-regions -S chip.means/trim28.log2ratioIN.bw \
  chip.means/h3k9me3.log2ratioIN.bw \
  chip.means/chip_k27m3_WT_11.12.mean.bw \
  -R L1/features/l1.beds/l1base2.bed \
  -a 5000 -b 5000 --regionBodyLength 3000 -bs 200 -p 8 \
  -out finalGroups/output/k27me3.k9me3.trim28/l1base2_wt_k27me3.k9me3.trim28_scaleRegion.mat

plotHeatmap -m finalGroups/output/k27me3.k9me3.trim28/l1base2_wt_k27me3.k9me3.trim28_scaleRegion.mat -o finalGroups/output/k27me3.k9me3.trim28/l1base2_wt_k27me3.k9me3.trim28_scaleRegion.pdf --missingDataColor 1  --sortUsingSamples 1   --zMin 0 --yMin -0.3 --yMax 1.2 2 15 --zMax 1 1.9 15   --colorMap Blues Greens Greys --regionsLabel "FLI-L1"  --samplesLabel "TRIM28 WT" "H3K9me3 WT" "H3K27me3 WT"


## LTR5hs

computeMatrix scale-regions -S chip.means/trim28.log2ratioIN.bw \
  chip.means/h3k9me3.log2ratioIN.bw \
  chip.means/chip_k27m3_WT_11.12.mean.bw \
  -R LTR5_Hs/LTR5_Hs.hg38.uniqID.FL900.bed \
  -a 1000 -b 1000 --regionBodyLength 500 -bs 50 -p 4 \
  -out LTR5_Hs/l1base2_wt_k27me3.k9me3.trim28_scaleRegion.mat

plotHeatmap -m LTR5_Hs/l1base2_wt_k27me3.k9me3.trim28_scaleRegion.mat -o LTR5_Hs/l1base2_wt_k27me3.k9me3.trim28_scaleRegion.pdf --missingDataColor 1  --sortUsingSamples 1   --zMin 0 --yMin -0.3 --yMax 1.2 2 15 --zMax 1 1.9 15   --colorMap Blues Greens Greys --regionsLabel "LTR5_Hs"  --samplesLabel "TRIM28 WT" "H3K9me3 WT" "H3K27me3 WT"
