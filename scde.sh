#!/bin/bash
#$ -cwd
#$ -j y
#$ -l h_vmem=2G
#$ -pe shm 20
#$ -l h_rt=12:00:00
#$ -A baas_wu
#$ -M yuezhang@stanford.edu
date
hostname
module load r/3.2.2
module load gcc/5.2.0

# go terms overdispersion analysis
Rscript  scde_pagoda.r
# call DE genes
Rscript  scde_de.r
