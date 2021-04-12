#!/bin/bash -l
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=00:20:00
#SBATCH --mem=5000M
#SBATCH --mail-type=all
#SBATCH --mail-user=lweis@princeton.edu
#SBATCH -o H3K27acbedtoolsout.%j
#SBATCH -e H3K27acbedtoolserr.%j

#Set Path for necessary programs
export PATH=$PATH:/tigress/AUTISM/Programs/bedtools-master/bin

file1=16CRShyperDMRs_ordered.bed
file2=52CRShypoDMRs_ordered.bed

# CTCF intersection
bedtools intersect -wa -a $file1 -b ../CTCF-ENCFF715PXI.bed > bedintersect_output/ctcf_$file1
bedtools intersect -wa -a $file2 -b ../CTCF-ENCFF715PXI.bed > bedintersect_output/ctcf_$file2

# H3K27ac intersection
bedtools intersect -wa -a $file1 -b ../H3K27ac-ENCFF792TCD.bed > bedintersect_output/H3K27ac_$file1
bedtools intersect -wa -a $file2 -b ../H3K27ac-ENCFF792TCD.bed > bedintersect_output/H3K27ac_$file2

# H3K4me1 intersection
bedtools intersect -wa -a $file1 -b ../H3K4me1-ENCFF801KBC.bed > bedintersect_output/H3K4me1_$file1
bedtools intersect -wa -a $file2 -b ../H3K4me1-ENCFF801KBC.bed > bedintersect_output/H3K4me1_$file2

# H3K4me3 intersection
bedtools intersect -wa -a $file1 -b ../H3K4me3-ENCFF953HUV.bed > bedintersect_output/H3K4me3_$file1
bedtools intersect -wa -a $file2 -b ../H3K4me3-ENCFF953HUV.bed > bedintersect_output/H3K4me3_$file2