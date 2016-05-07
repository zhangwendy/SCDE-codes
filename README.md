# SCDE-codes

The R scripts implements a set of statistical methods from SCDE for analyzing single-cell RNA-seq data. 

1. STAR.sh   The job submitting bash file for STAR alignment with STAR parameters 
2. mergy_STAR_stat.py   This python script will mergy read count from STAR for each gene each sample to make a big matrix for SCDE

3. scde_pagoda.r   The R script implemented in the scde resolves multiple, potentially overlapping aspects of transcriptional heterogeneity by identifying known pathways that show significant excess of coordinated variability among the measured cells
4. scde_de.r  The R script implements routines for fitting individual error models and call differently expressed genes for single-cell RNA-seq measurements
5. scde.sh   The job submitting bash file for SCDE
