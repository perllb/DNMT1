## Analyse WGBS data 

* ### Generation of figures with WGBS data

    #### Fig 1C:
    1. Compute global levels (1000bp bins): 
        > Process_WGBS.sh
    
    2. Fig1C violin plots: 
        > Fig1C_vioplot.binSize_1000.R
    
    #### Fig 1E and 2E:
    
    1. Get TE Classes coordinates from RepeatMasker hg38 bed-file
        > Fig1E_1_getClasses.bed.sh

    2. Get L1-promoter regions
        > Fig1E_2_L1.promoter.sh
        
    3. Compute mean mCpG/CpG ratio over features 
        > Fig1E_3_compute.features.3x.sh
        
    4. Fig1E violin plots:
        > Fig1E_4_TE_analysis.3x.WGBS.R
        
    #### Heatmaps (deeptools) are described under Deeptools_WGBS.RNAseq.ChIPseq folder
