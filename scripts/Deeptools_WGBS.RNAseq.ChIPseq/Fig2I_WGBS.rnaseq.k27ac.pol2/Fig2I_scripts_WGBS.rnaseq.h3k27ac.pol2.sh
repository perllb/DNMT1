#### Plot L1 FL6000, SVA (1000), PJ HERV (score300) ####


### WGBS - RNAseq - K27ac - PolII ###

cd "/Users/perludvik/Dropbox (MN)/Per/PhD/Projects/DNAmeth/hNES DNMT1 KO/2018.08.30_DNMT1L1_figures_scripts/ChIP-integration/"

##  L1 HS,PA2-PA7 ##

# get FL L1s
grep L1HS ~/Documents/bioinformatics/genomicData/hg38/Retro.hg38.uniqID.bed | awk ' { if($3-$2>6000 || $2-$3>6000) { print $0; } } ' > L1HS.FL6000.bed
grep L1PA2 ~/Documents/bioinformatics/genomicData/hg38/Retro.hg38.uniqID.bed | awk ' { if($3-$2>6000 || $2-$3>6000) { print $0; } } ' > L1PA2.FL6000.bed
grep L1PA3 ~/Documents/bioinformatics/genomicData/hg38/Retro.hg38.uniqID.bed | awk ' { if($3-$2>6000 || $2-$3>6000) { print $0; } } ' > L1PA3.FL6000.bed
grep L1PA4 ~/Documents/bioinformatics/genomicData/hg38/Retro.hg38.uniqID.bed | awk ' { if($3-$2>6000 || $2-$3>6000) { print $0; } } ' > L1PA4.FL6000.bed
grep L1PA5 ~/Documents/bioinformatics/genomicData/hg38/Retro.hg38.uniqID.bed | awk ' { if($3-$2>6000 || $2-$3>6000) { print $0; } } ' > L1PA5.FL6000.bed
grep L1PA6 ~/Documents/bioinformatics/genomicData/hg38/Retro.hg38.uniqID.bed | awk ' { if($3-$2>6000 || $2-$3>6000) { print $0; } } ' > L1PA6.FL6000.bed
grep L1PA7 ~/Documents/bioinformatics/genomicData/hg38/Retro.hg38.uniqID.bed | awk ' { if($3-$2>6000 || $2-$3>6000) { print $0; } } ' > L1PA7.FL6000.bed


# WGBS + RNA seq + K27ac + PolII with L1HS-PA3

### MissingDataZero
computeMatrix scale-regions -S chip.means/WT_bismark_bt2_pe.deduplicated.bismark.bw \
  chip.means/DNMT1_KO_bismark_bt2_pe.deduplicated.bismark.bw \
  chip.means/DNMT1KO_hg38.unique_WT_REV.FWD.sum.bw \
  chip.means/DNMT1KO_hg38.unique_KO_REV.FWD.sum.bw \
  chip.means/chip_k27ac_WT_3.4.mean.bw \
  chip.means/chip_k27ac_KO_1.2.mean.bw  \
  chip.means/chip_pol_WT_7.8.mean.bw \
  chip.means/chip_pol_KO_5.6.mean.bw \
  -R L1/features/l1.beds/L1HS.FL6000.bed \
  L1/features/l1.beds/L1PA2.FL6000.bed \
  L1/features/l1.beds/L1PA3.FL6000.bed \
  L1/features/l1.beds/L1PA4.FL6000.bed \
  L1/features/l1.beds/L1PA5.FL6000.bed \
  L1/features/l1.beds/L1PA6.FL6000.bed \
  L1/features/l1.beds/L1PA7.FL6000.bed \
  -a 5000 -b 5000 --regionBodyLength 3000 -bs 200 -p 8 --missingDataAsZero \
  -out finalGroups/output/WGBS.rnaseq.k27ac.pol2/MissingDataZero/l1hs.pa7_WGBS_RNAseq_k27ac_polII_scaleRegion.MissingDataZero.mat

plotProfile -m finalGroups/output/WGBS.rnaseq.k27ac.pol2/MissingDataZero/l1hs.pa7_WGBS_RNAseq_k27ac_polII_scaleRegion.MissingDataZero.mat -o finalGroups/output/WGBS.rnaseq.k27ac.pol2/MissingDataZerol1hs.pa7_WGBS_RNAseq_k27ac_polII_scaleRegion.MissingDataZero.profile.pdf --yMin 0 --yMax 10 10 10 10 10 10 12 1.4  --colors purple red blue green yellow orange grey black --regionsLabel "L1HS" "L1PA2" "L1PA3" "L1PA4" "L1PA5" "L1PA6" "L1PA7"  --samplesLabel "WGBS WT" "WGBS KO" "RNAseq_WT" "RNAseq_KO" "H3K27ac WT" "H3K27ac KO" "polII WT" "polII KO"

### MissingDataZero - TSS - TES only
computeMatrix scale-regions -S chip.means/WT_bismark_bt2_pe.deduplicated.bismark.bw \
  chip.means/DNMT1_KO_bismark_bt2_pe.deduplicated.bismark.bw \
  chip.means/DNMT1KO_hg38.unique_WT_REV.FWD.sum.bw \
  chip.means/DNMT1KO_hg38.unique_KO_REV.FWD.sum.bw \
  chip.means/chip_k27ac_WT_3.4.mean.bw \
  chip.means/chip_k27ac_KO_1.2.mean.bw  \
  chip.means/chip_pol_WT_7.8.mean.bw \
  chip.means/chip_pol_KO_5.6.mean.bw \
  -R L1/features/l1.beds/L1HS.FL6000.bed \
  L1/features/l1.beds/L1PA2.FL6000.bed \
  L1/features/l1.beds/L1PA3.FL6000.bed \
  L1/features/l1.beds/L1PA4.FL6000.bed \
  L1/features/l1.beds/L1PA5.FL6000.bed \
  L1/features/l1.beds/L1PA6.FL6000.bed \
  L1/features/l1.beds/L1PA7.FL6000.bed \
  -a 200 -b 200 --regionBodyLength 3000 -bs 200 -p 8 --missingDataAsZero \
  -out finalGroups/output/WGBS.rnaseq.k27ac.pol2/MissingDataZero/l1hs.pa7_WGBS_RNAseq_k27ac_polII_scaleRegion.MissingDataZero_geneBody.mat

plotProfile -m finalGroups/output/WGBS.rnaseq.k27ac.pol2/MissingDataZero/l1hs.pa7_WGBS_RNAseq_k27ac_polII_scaleRegion.MissingDataZero_geneBody.mat -o finalGroups/output/WGBS.rnaseq.k27ac.pol2/MissingDataZero/l1hs.pa7_WGBS_RNAseq_k27ac_polII_scaleRegion.MissingDataZero_geneBody.profile.pdf --yMin 0 --yMax 10 10 3 3 10 10 12 1.4  --colors purple red blue green yellow orange grey black --regionsLabel "L1HS" "L1PA2" "L1PA3" "L1PA4" "L1PA5" "L1PA6" "L1PA7"  --samplesLabel "WGBS WT" "WGBS KO" "RNAseq_WT" "RNAseq_KO" "H3K27ac WT" "H3K27ac KO" "polII WT" "polII KO"

### Normal - main!
computeMatrix scale-regions -S chip.means/WT_bismark_bt2_pe.deduplicated.bismark.bw \
  chip.means/DNMT1_KO_bismark_bt2_pe.deduplicated.bismark.bw \
  chip.means/DNMT1KO_hg38.unique_WT_REV.FWD.sum.bw \
  chip.means/DNMT1KO_hg38.unique_KO_REV.FWD.sum.bw \
  chip.means/chip_k27ac_WT_3.4.mean.bw \
  chip.means/chip_k27ac_KO_1.2.mean.bw  \
  chip.means/chip_pol_WT_7.8.mean.bw \
  chip.means/chip_pol_KO_5.6.mean.bw \
  -R L1/features/l1.beds/L1HS.FL6000.bed \
  L1/features/l1.beds/L1PA2.FL6000.bed \
  L1/features/l1.beds/L1PA3.FL6000.bed \
  L1/features/l1.beds/L1PA4.FL6000.bed \
  L1/features/l1.beds/L1PA5.FL6000.bed \
  L1/features/l1.beds/L1PA6.FL6000.bed \
  L1/features/l1.beds/L1PA7.FL6000.bed \
  -a 5000 -b 5000 --regionBodyLength 3000 -bs 200 -p 8  \
  -out finalGroups/output/WGBS.rnaseq.k27ac.pol2/l1hs.pa7_WGBS_RNAseq_k27ac_polII_scaleRegion.mat

plotProfile -m finalGroups/output/WGBS.rnaseq.k27ac.pol2/l1hs.pa7_WGBS_RNAseq_k27ac_polII_scaleRegion.mat -o finalGroups/output/WGBS.rnaseq.k27ac.pol2/l1hs.pa7_WGBS_RNAseq_k27ac_polII_scaleRegion.profile.pdf --yMin 0 --yMax 100 100 200 200 10 10 12 1.4  --colors purple red blue green yellow orange grey black --regionsLabel "L1HS" "L1PA2" "L1PA3" "L1PA4" "L1PA5" "L1PA6" "L1PA7"  --samplesLabel "WGBS WT" "WGBS KO" "RNAseq_WT" "RNAseq_KO" "H3K27ac WT" "H3K27ac KO" "polII WT" "polII KO"

plotHeatmap -m finalGroups/output/WGBS.rnaseq.k27ac.pol2/l1hs.pa7_WGBS_RNAseq_k27ac_polII_scaleRegion.mat -o finalGroups/output/WGBS.rnaseq.k27ac.pol2/l1hs.pa7_WGBS_RNAseq_k27ac_polII_scaleRegion.pdf --missingDataColor 1  --sortUsingSamples 6   --zMin 0 --yMin 0 --yMax 100 100 200 200 12 12 12 1.4  --zMax 100 100 3 3 10 10 12 1   --colorMap Greys Greys Greens Greens Reds Reds Blues Blues --regionsLabel "L1HS" "L1PA2" "L1PA3" "L1PA4" "L1PA5" "L1PA6" "L1PA7"  --samplesLabel "WGBS WT" "WGBS KO" "RNAseq_WT" "RNAseq_KO" "H3K27ac WT" "H3K27ac KO" "polII WT" "polII KO"


## PJ HERV - 300
# WGBS + RNA seq + K27ac + PolII with HERV

## missingDataZero
computeMatrix scale-regions -S chip.means/WT_bismark_bt2_pe.deduplicated.bismark.bw  \
                                chip.means/DNMT1_KO_bismark_bt2_pe.deduplicated.bismark.bw \
                                 chip.means/DNMT1KO_hg38.unique_WT_REV.FWD.sum.bw \
                                 chip.means/DNMT1KO_hg38.unique_KO_REV.FWD.sum.bw \
                                 chip.means/chip_k27ac_WT_3.4.mean.bw \
                                 chip.means/chip_k27ac_KO_1.2.mean.bw  \
                                 chip.means/chip_pol_WT_7.8.mean.bw \
                                 chip.means/chip_pol_KO_5.6.mean.bw \
                                 -R HERV/features/hg38.HERV.score300.bed \
                                  -a 5000 -b 5000 --regionBodyLength 5000 -bs 200 -p 8 --missingDataAsZero \
                                  -out finalGroups/output/WGBS.rnaseq.k27ac.pol2/MissingDataZero/HERV300_WGBS_RNAseq_k27ac_polII_scaleRegion.missingDataZero.mat

plotProfile -m finalGroups/output/WGBS.rnaseq.k27ac.pol2/MissingDataZero/HERV300_WGBS_RNAseq_k27ac_polII_scaleRegion.missingDataZero.mat -o finalGroups/output/WGBS.rnaseq.k27ac.pol2/MissingDataZero/HERV300_WGBS_RNAseq_k27ac_polII_scaleRegion.missingDataZero.profile.pdf --yMin 0 --yMax 10 10 40 40 10 10 12 1.4  --colors purple red blue green yellow orange grey black --regionsLabel "HERV_300s"  --samplesLabel "WGBS WT" "WGBS KO" "RNAseq_WT" "RNAseq_KO" "H3K27ac WT" "H3K27ac KO" "polII WT" "polII KO"

computeMatrix scale-regions -S chip.means/WT_bismark_bt2_pe.deduplicated.bismark.bw  \
                                chip.means/DNMT1_KO_bismark_bt2_pe.deduplicated.bismark.bw \
                                chip.means/DNMT1KO_hg38.unique_WT_REV.FWD.sum.bw \
                                chip.means/DNMT1KO_hg38.unique_KO_REV.FWD.sum.bw \
                                chip.means/chip_k27ac_WT_3.4.mean.bw \
                                chip.means/chip_k27ac_KO_1.2.mean.bw  \
                                chip.means/chip_pol_WT_7.8.mean.bw \
                                chip.means/chip_pol_KO_5.6.mean.bw \
                                -R HERV/features/hg38.HERV.score300.bed \
                                -a 5000 -b 5000 --regionBodyLength 5000 -bs 200 -p 8 \
                                 -out finalGroups/output/WGBS.rnaseq.k27ac.pol2/HERV300_WGBS_RNAseq_k27ac_polII_scaleRegion.mat

plotProfile -m finalGroups/output/WGBS.rnaseq.k27ac.pol2/HERV300_WGBS_RNAseq_k27ac_polII_scaleRegion.mat -o finalGroups/output/WGBS.rnaseq.k27ac.pol2/HERV300_WGBS_RNAseq_k27ac_polII_scaleRegion.profile.pdf --yMin 0 --yMax 100 100 200 200 10 10 12 1.4  --colors purple red blue green yellow orange grey black --regionsLabel "HERV_300s"  --samplesLabel "WGBS WT" "WGBS KO" "RNAseq_WT" "RNAseq_KO" "H3K27ac WT" "H3K27ac KO" "polII WT" "polII KO"

plotHeatmap -m finalGroups/output/WGBS.rnaseq.k27ac.pol2/HERV300_WGBS_RNAseq_k27ac_polII_scaleRegion.mat -o finalGroups/output/WGBS.rnaseq.k27ac.pol2/HERV300_WGBS_RNAseq_k27ac_polII_scaleRegion.pdf --missingDataColor 1  --sortUsingSamples 6   --zMin 0 --yMin 0 --yMax 100 100 200 200 12 12 12 1.4  --zMax 100 100  3 3  20 20 20 10   --colorMap Greys Greys Greens Greens Reds Reds Blues Blues  --regionsLabel "HERV"  --samplesLabel "WGBS WT" "WGBS KO" "RNAseq_WT" "RNAseq_KO" "H3K27ac WT" "H3K27ac KO" "polII WT" "polII KO"


## SVA 1000bp
cd "/Users/perludvik/Dropbox (MN)/Per/PhD/Projects/DNAmeth/hNES DNMT1 KO/2018.08.30_DNMT1L1_figures_scripts/ChIP-integration/SVA/features"

awk ' { if($3-$2>1000 || $2-$3>1000) { print $0 } } ' SVA_A.hg38.bed > SVA_A.hg38_FL1000.bed
awk ' { if($3-$2>1000 || $2-$3>1000) { print $0 } } ' SVA_B.hg38.bed > SVA_B.hg38_FL1000.bed
awk ' { if($3-$2>1000 || $2-$3>1000) { print $0 } } ' SVA_C.hg38.bed > SVA_C.hg38_FL1000.bed
awk ' { if($3-$2>1000 || $2-$3>1000) { print $0 } } ' SVA_D.hg38.bed > SVA_D.hg38_FL1000.bed
awk ' { if($3-$2>1000 || $2-$3>1000) { print $0 } } ' SVA_E.hg38.bed > SVA_E.hg38_FL1000.bed
awk ' { if($3-$2>1000 || $2-$3>1000) { print $0 } } ' SVA_F.hg38.bed > SVA_F.hg38_FL1000.bed

cd .. ;  cd ..

# WGBS + RNA seq + K27ac + PolII with SVA
## MissingDataZero
computeMatrix scale-regions -S chip.means/WT_bismark_bt2_pe.deduplicated.bismark.bw  \
                                chip.means/DNMT1_KO_bismark_bt2_pe.deduplicated.bismark.bw \
                                chip.means/DNMT1KO_hg38.unique_WT_REV.FWD.sum.bw \
                                chip.means/DNMT1KO_hg38.unique_KO_REV.FWD.sum.bw \
                                chip.means/chip_k27ac_WT_3.4.mean.bw \
                                chip.means/chip_k27ac_KO_1.2.mean.bw  \
                                chip.means/chip_pol_WT_7.8.mean.bw \
                                chip.means/chip_pol_KO_5.6.mean.bw \
                                -R SVA/features/SVA_F.hg38_FL1000.bed \
                                SVA/features/SVA_E.hg38_FL1000.bed \
                                SVA/features/SVA_D.hg38_FL1000.bed \
                                SVA/features/SVA_C.hg38_FL1000.bed \
                                SVA/features/SVA_B.hg38_FL1000.bed \
                                SVA/features/SVA_A.hg38_FL1000.bed \
                                -a 500 -b 500 --regionBodyLength 500 -bs 20 -p 8 --missingDataAsZero \
                                -out finalGroups/output/WGBS.rnaseq.k27ac.pol2/MissingDataZero/SVA.FL1000_WGBS_RNAseq_k27ac_polII_scaleRegion.missingDataZero.mat

plotProfile -m finalGroups/output/WGBS.rnaseq.k27ac.pol2/MissingDataZero/SVA.FL1000_WGBS_RNAseq_k27ac_polII_scaleRegion.missingDataZero.mat -o finalGroups/output/WGBS.rnaseq.k27ac.pol2/MissingDataZero/SVA.FL1000_WGBS_RNAseq_k27ac_polII_scaleRegion.MissingDataZero.profile.pdf --yMin 0 --yMax 15 15 40 40 10 10 12 1.4  --colors purple red blue green yellow orange grey black --regionsLabel "SVA_F" "SVA_E" "SVA_D" "SVA_C" "SVA_B" "SVA_A"  --samplesLabel "WGBS WT" "WGBS KO" "RNAseq_WT" "RNAseq_KO" "H3K27ac WT" "H3K27ac KO" "polII WT" "polII KO"

# Missing Data NA
computeMatrix scale-regions -S chip.means/WT_bismark_bt2_pe.deduplicated.bismark.bw  \
                                chip.means/DNMT1_KO_bismark_bt2_pe.deduplicated.bismark.bw \
                                chip.means/DNMT1KO_hg38.unique_WT_REV.FWD.sum.bw \
                                chip.means/DNMT1KO_hg38.unique_KO_REV.FWD.sum.bw \
                                chip.means/chip_k27ac_WT_3.4.mean.bw \
                                chip.means/chip_k27ac_KO_1.2.mean.bw  \
                                chip.means/chip_pol_WT_7.8.mean.bw \
                                chip.means/chip_pol_KO_5.6.mean.bw \
                                -R SVA/features/SVA_F.hg38_FL1000.bed \
                                SVA/features/SVA_E.hg38_FL1000.bed \
                                SVA/features/SVA_D.hg38_FL1000.bed \
                                SVA/features/SVA_C.hg38_FL1000.bed \
                                SVA/features/SVA_B.hg38_FL1000.bed \
                                SVA/features/SVA_A.hg38_FL1000.bed \
                                -a 500 -b 500 --regionBodyLength 500 -bs 10 -p 8 \
                                -out finalGroups/output/WGBS.rnaseq.k27ac.pol2/SVA.FL1000_WGBS_RNAseq_k27ac_polII_scaleRegion.mat

plotHeatmap -m finalGroups/output/WGBS.rnaseq.k27ac.pol2/SVA.FL1000_WGBS_RNAseq_k27ac_polII_scaleRegion.mat -o finalGroups/output/WGBS.rnaseq.k27ac.pol2/SVA.FL1000_WGBS_RNAseq_k27ac_polII_scaleRegion.pdf --missingDataColor 1  --sortUsingSamples 6   --zMin 0 --yMin 0 --yMax 100 100 200 200 12 12 12 1.4  --zMax 100 100 20 20 20 20 20 12   --colorMap Greys Greys Greens Greens Reds Reds Blues Blues --regionsLabel "SVA_F" "SVA_E" "SVA_D" "SVA_C" "SVA_B" "SVA_A"  --samplesLabel "WGBS WT" "WGBS KO" "RNAseq_WT" "RNAseq_KO" "H3K27ac WT" "H3K27ac KO" "polII WT" "polII KO"

plotProfile -m finalGroups/output/WGBS.rnaseq.k27ac.pol2/SVA.FL1000_WGBS_RNAseq_k27ac_polII_scaleRegion.mat -o finalGroups/output/WGBS.rnaseq.k27ac.pol2/SVA.FL1000_WGBS_RNAseq_k27ac_polII_scaleRegion.profile.pdf --yMin 0 --yMax 100 100 200 200 10 10 12 1.4  --colors purple red blue green yellow orange grey black --regionsLabel "SVA_F" "SVA_E" "SVA_D" "SVA_C" "SVA_B" "SVA_A"  --samplesLabel "WGBS WT" "WGBS KO" "RNAseq_WT" "RNAseq_KO" "H3K27ac WT" "H3K27ac KO" "polII WT" "polII KO"


#### FLI-L1

computeMatrix scale-regions -S chip.means/WT_bismark_bt2_pe.deduplicated.bismark.bw  \
  chip.means/DNMT1_KO_bismark_bt2_pe.deduplicated.bismark.bw \
  chip.means/DNMT1KO_hg38.unique_WT_REV.FWD.sum.bw \
  chip.means/DNMT1KO_hg38.unique_KO_REV.FWD.sum.bw \
  chip.means/chip_k27ac_WT_3.4.mean.bw \
  chip.means/chip_k27ac_KO_1.2.mean.bw \
  chip.means/chip_pol_WT_7.8.mean.bw  \
  chip.means/chip_pol_KO_5.6.mean.bw \
  -R L1/features/l1.beds/l1base2.bed \
   -a 5000 -b 5000 --regionBodyLength 3000 -bs 200 -p 4 \
   -out finalGroups/output/WGBS.rnaseq.k27ac.pol2/l1base2_wgbs.rnaseq.k27ac.pol2_scaleRegions.mat

plotHeatmap -m finalGroups/output/WGBS.rnaseq.k27ac.pol2/l1base2_wgbs.rnaseq.k27ac.pol2_scaleRegions.mat  -o finalGroups/output/WGBS.rnaseq.k27ac.pol2/l1base2_wgbs.rnaseq.k27ac.pol2_scaleRegions.pdf --missingDataColor 1  --sortUsingSamples 6  --zMin 0 --yMin 0 --yMax 100 100 200 200 12 5 12 0.5 --zMax 100 100 3 3 10 10 12 1   --colorMap Greys Greys Greens Greens Reds Reds Blues Blues --regionsLabel "FLI-L1"  --samplesLabel "WGBS WT" "WGBS KO" "RNAseq_WT" "RNAseq_KO" "H3K27ac WT" "H3K27ac KO" "polII WT" "polII KO"


### LTR5_Hs


computeMatrix scale-regions -S chip.means/WT_bismark_bt2_pe.deduplicated.bismark.bw  \
  chip.means/DNMT1_KO_bismark_bt2_pe.deduplicated.bismark.bw \
  chip.means/DNMT1KO_hg38.unique_WT_REV.FWD.sum.bw \
  chip.means/DNMT1KO_hg38.unique_KO_REV.FWD.sum.bw \
  chip.means/chip_k27ac_WT_3.4.mean.bw \
  chip.means/chip_k27ac_KO_1.2.mean.bw \
  chip.means/chip_pol_WT_7.8.mean.bw  \
  chip.means/chip_pol_KO_5.6.mean.bw \
  -R LTR5_Hs/LTR5_Hs.hg38.uniqID.FL900.bed \
   -a 1000 -b 1000 --regionBodyLength 500 -bs 50 -p 4 \
   -out LTR5_Hs/LTR5hs_wgbs.rnaseq.k27ac.pol2_scaleRegions.mat

plotHeatmap -m LTR5_Hs/LTR5hs_wgbs.rnaseq.k27ac.pol2_scaleRegions.mat  -o LTR5_Hs/LTR5hs_wgbs.rnaseq.k27ac.pol2_scaleRegions.pdf --missingDataColor 1  --sortUsingSamples 6  --zMin 0 --yMin 0 --yMax 100 100 200 200 12 5 12 0.5 --zMax 100 100 3 3 10 10 12 1   --colorMap Greys Greys Greens Greens Reds Reds Blues Blues --regionsLabel "LTR5_Hs"  --samplesLabel "WGBS WT" "WGBS KO" "RNAseq_WT" "RNAseq_KO" "H3K27ac WT" "H3K27ac KO" "polII WT" "polII KO"
