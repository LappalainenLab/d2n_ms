#!/bin/bash
#SBATCH --job-name=Bedtools_K562
#SBATCH --cpus-per-task=2
#SBATCH --mem=16g
#SBATCH --output=../../logs/stdout/bedtools
#SBATCH --error=../../logs/stderr/bedtools

mkdir -p ../../logs/stdout/bedtools ../../logs/stderr/bedtools ../../processed_data/01-library_design/bedtools/

bedtools sort -i ../../processed_data/01-library_design/bedtools/D2N_GDT3_hg19.bed | bedtools intersect -v -a - -b /gpfs/commons/groups/lappalainen_lab/jdomingo/data/genomes/K562/ENCODE/ENCFF538YDL.vcf.gz >  ../../processed_data/01-library_design/bedtools/D2N_GDT3_hg19_intK562.bed

