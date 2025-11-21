#!/bin/bash -l
#SBATCH --nodes 1
#SBATCH --ntasks 8
#SBATCH --cpus-per-task 1
#SBATCH --mem=1G
#SBATCH --time=0
#SBATCH --mail-user=your_mailj@whatever.com
#SBATCH --error=script.err

# $1 basename of braker's gtf file

get_longest_isoform.py -g $1.gtf -o $1_longest.gtf
