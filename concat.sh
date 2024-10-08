#!/bin/bash
#SBATCH --account=priority-davidnidever
#SBATCH --job-name=concat
#SBATCH --output=/home/x25h971/canfind_dr2/outfiles/concat.out
#SBATCH --error=/home/x25h971/canfind_dr2/outfiles/concat.err
#SBATCH --partition=priority
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=6000
#SBATCH --time=00:60:00
module load Anaconda3
source activate $HOME/condaenv/
python /home/x25h971/canfind_dr2/files/tracklet_concat.py 2 2

