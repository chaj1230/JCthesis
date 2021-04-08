#!/bin/bash -l
#SBATCH --nodes=1
#SBATCH --ntasks=5
#SBATCH --time=20:00:00
#SBATCH --mem=50G
#SBATCH --mail-type=all
#SBATCH --mail-user=jaeyoonc@princeton.edu
#SBATCH -o bismarkextractrr16out.%j
#SBATCH -e bismarkextractrr16err.%j

#Set Path for necessary programs
export PATH=$PATH:/tigress/AUTISM/Programs/FastQC:/tigress/AUTISM/Programs/pythonenv/bin/:/tigress/AUTISM/Programs/TrimGalore-0.6.6:/tigress/AUTISM/Programs/bowtie2-2.2.8:/tigress/AUTISM/Programs/Bismark-0.22.3:/tigress/AUTISM/Programs/samtools-1.2

#Trim using trim galore
module load anaconda3/5.3.1
source activate /tigress/AUTISM/Programs/pythonenv

cd /tigress/AUTISM/ssm/fastq/combinedfastq
trim_galore --2colour 20 --illumina --keep --rrbs --paired -o ../../trimmed G695_rr16_R1_001.fastq.gz G695_rr16_R2_001.fastq.gz


# BISMARK MAPPER
cd ../../trimmed/
mkdir ../aligned/rr16combined
bismark --genome /tigress/AUTISM/refGenomes/mouse/mouse_GRCm38 --parallel 5 -1 G695_rr16_R1_001_val_1.fq.gz -2 G695_rr16_R2_001_val_2.fq.gz -N 0 --temp_dir tmp --rg_tag --rg_id rr16 --rg_sample rr16 --unmapped --ambiguous -o ../aligned/rr16combined --non_bs_mm --nucleotide_coverage
#
#
#
## EXTRACTOR
cd ../aligned/rr16combined
bismark_methylation_extractor --gzip --bedGraph --buffer_size 2G --paired-end --comprehensive -merge_non_CpG --remove_spaces --ignore 3 --ignore_r2 2 --cytosine_report --genome_folder /tigress/AUTISM/refGenomes/mouse/mouse_GRCm38 --multicore 1 G695_rr16_R1_001_val_1_bismark_bt2_pe.bam
#
bismark2report
#
samtools sort -O bam -o G695_rr16_R1_001_val_1_bismark_bt2_pe.sorted.bam -T tmp G695_rr16_R1_001_val_1_bismark_bt2_pe.bam
samtools index G695_rr16_R1_001_val_1_bismark_bt2_pe.sorted.bam
#
##Add samtools flagstat step
samtools flagstat G695_rr16_R1_001_val_1_bismark_bt2_pe.sorted.bam > G695_rr16_R1_001_val_1_bismark_bt2_pe.sorted.metrics.txt
#
##Add check insert size step
cd /tigress/AUTISM/Programs/picard-tools-2.5.0
java -jar picard.jar CollectInsertSizeMetrics I=/tigress/AUTISM/ssm/aligned/rr16combined/G695_rr16_R1_001_val_1_bismark_bt2_pe.sorted.bam O=/tigress/AUTISM/ssm/aligned/rr16combined/G695_rr16_R1_001_insert_size_metrics.txt H=/tigress/AUTISM/ssm/aligned/rr16combined/G695_rr16_R1_001_insert_size_histogram.pdf M=0.5
