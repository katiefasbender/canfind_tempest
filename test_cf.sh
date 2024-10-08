#!/bin/bash
#SBATCH --account=priority-davidnidever
#SBATCH --job-name=cf_dr2_1645749692.9475627_6058
#SBATCH --output=/home/x25h971/canfind_dr2/outfiles/cf_dr2_1645749692.9475627_6058.out
#SBATCH --error=/home/x25h971/canfind_dr2/outfiles/cf_dr2_1645749692.9475627_6058.err
#SBATCH --partition=priority
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2000
#SBATCH --time=00:60:00
module load Anaconda3/2020.07
source activate $HOME/condaenv/
python /home/x25h971/canfind_dr2/files/canfind_v2.py 144310 1

