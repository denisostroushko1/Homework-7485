#!/bin/bash
#SBATCH -A wksh0001
#SBATCH --time=12:00:00
#SBATCH --ntasks=8
#SBATCH --mem=64g
#SBATCH --tmp=64g
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=wu001125@umn.edu

cd /home/wksh0001/wu001125/7485/project
module load R/4.3.0-openblas
Rscript final12.R




