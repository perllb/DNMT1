## RNA-seq pipelines

* #### STAR alignment scripts (unique and multimapping)
    > Mapping.sh

* #### Figures:
* ##### Fig 2A,B,C: upregulation of TEs:
    > analysis/TE_activation/Fig2A.B.C_scatter all TE.R
* ##### Fig 2F: FLI-L1 sense activation:
    > analysis/TE_activation/Fig2F_FLI-L1 sense anti activation.R
* ##### Fig 3A.C: Genes expression changes, genes close to L1 up and SVA / HERV
    > analysis/mRNA/NearbyGenes/Fig3A.C_near.genes.R
* ##### Fig 3B: maPlot coding genes
    > analysis/mRNA/Fig3B_mRNA scatter.R
* ##### Fig 3E: Heatmap TE-gene fusion reads
    1. Get reads that overlap NCBI gene AND L1
    > analysis/mRNA/Chimeric_reads/*
    
    2. Plot heatmap with reads  
    > analysis/mRNA/Chimeric_reads/Fig3E_FusionReads.heatmap.R

* ##### Fig 3F: Heatmap fusion genes hNES diff
    > analysis/mRNA/Fig3F_heatmap_hNESdiff_fusionGenes.R
