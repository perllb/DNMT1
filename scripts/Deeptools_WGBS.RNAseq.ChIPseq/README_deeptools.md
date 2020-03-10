## Visualize WGBS, RNA-seq and ChIP-seq data with Deeptools software

* ### Heatmaps were generated in 3 pipelines:
##### 1. ChIP-seq data for H3K27me3, H3K9me3 and TRIM28 (all only WT) at Transposons
##### 2. WGBS (WT/KO), RNA-seq (WT/KO), H3K27ac (WT/KO) and pol2 (WT/KO) at Transposons
##### 3. WGBS (WT/KO), RNA-seq (WT/KO), H3K27ac (WT/KO) and pol2 (WT/KO) at L1-gene fusion pairs


* #### Fig 2I:
    ##### 1. Get coordinates of TEs and create heatmap:
    > Fig2I_WGBS.rnaseq.k27ac.pol2/Fig2I_scripts_WGBS.rnaseq.h3k27ac.pol2.sh

* #### Fig 2J:
    ##### 1. Get coordinates of TEs and create heatmap:
    > Fig2J_k27me3.k9me3.trim28/Fig2J_scripts_wt_k27me3.k9me3.trim28.sh

* #### Fig 3C:
    ##### 1. Get genes close to upregulated L1s
    > Fig3C_NearbyGenes/Fig3C_get.near.genes.R
   
    ##### 2. Create heatmap Fig3C:
    > Fig3C_NearbyGenes/Fig3C_scripts_deeptools.neargenes.sh

