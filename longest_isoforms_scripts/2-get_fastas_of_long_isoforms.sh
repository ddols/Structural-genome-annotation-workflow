#!/bin/bash -l
#SBATCH --nodes 1
#SBATCH --ntasks 8
#SBATCH --cpus-per-task 1
#SBATCH --mem=1G
#SBATCH --time=0
#SBATCH --mail-user=your_mail@whatever.com
#SBATCH --error=script.err

# $1 masked genome from /results/1_MaskRepeat/RepeatMasker/
# $2 gtf file obtained from running 1-get_longest_isoforms.sh

module load augustus

/aplic/augustus/3.4.0/scripts/getAnnoFastaFromJoingenes.py -g $1 \
 -o longest_isoforms -f $2
