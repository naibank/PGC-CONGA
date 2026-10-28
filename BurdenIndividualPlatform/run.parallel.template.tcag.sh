#!/bin/bash
#SBATCH -N 1 -c 8
#SBATCH --job-name MGRB
#SBATCH --mem=36G
#SBATCH --time=1:00:00
#SBATCH --partition=tcag 
#SBATCH --reservation=tcag250gb
#SBATCH --constraint=AlmaLinux9
module load R/4.5.2;module load openssl/1.0

Rscript ../burdenASD.R
