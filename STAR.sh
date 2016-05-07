#!/bin/bash
#$ -cwd
#$ -j y
#$ -l h_vmem=40G
#$ -l h_rt=120:00:00
#$ -A baas_wu
#$ -M yuezhang@stanford.edu
date
hostname
module load STAR/2.5.1b
STAR --genomeDir $STAR_HG19_GENOME --genomeLoad LoadAndExit
STAR --genomeDir $STAR_HG19_GENOME --genomeLoad LoadAndKeep  --readFilesCommand zcat  --quantMode GeneCounts TranscriptomeSAM  --runThreadN 8 --outReadsUnmapped Fastx --outFilterType BySJout --outSAMattributes NH HI AS nM NM MD XS --outSAMstrandField intronMotif  --outSAMtype BAM Unsorted  --alignSoftClipAtReferenceEnds No  --outFilterIntronMotifs RemoveNoncanonicalUnannotated --readFilesIn  /srv/gsfs0/projects/wu/Ioannis/IK-1-TAAGGCGA-ATAGAGAG_S1_R1_001.fastq.gz  /srv/gsfs0/projects/wu/Ioannis/IK-1-TAAGGCGA-ATAGAGAG_S1_R2_001.fastq.gz  --outFileNamePrefix IK-1_
STAR --genomeDir $STAR_HG19_GENOME --genomeLoad LoadAndKeep  --readFilesCommand zcat  --quantMode GeneCounts TranscriptomeSAM  --runThreadN 8 --outReadsUnmapped Fastx --outFilterType BySJout --outSAMattributes NH HI AS nM NM MD XS --outSAMstrandField intronMotif  --outSAMtype BAM Unsorted  --alignSoftClipAtReferenceEnds No  --outFilterIntronMotifs RemoveNoncanonicalUnannotated --readFilesIn  /srv/gsfs0/projects/wu/Ioannis/IK-10-CGAGGCTG-ATAGAGAG_S10_R1_001.fastq.gz  /srv/gsfs0/projects/wu/Ioannis/IK-10-CGAGGCTG-ATAGAGAG_S10_R2_001.fastq.gz  --outFileNamePrefix IK-10_
