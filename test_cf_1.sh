#!/bin/bash
#SBATCH --account=priority-davidnidever
#SBATCH --job-name=cf_dr2_1645749692.9475627_6059
#SBATCH --output=/home/x25h971/canfind_dr2/outfiles/cf_dr2_1645749692.9475627_6059.out
#SBATCH --error=/home/x25h971/canfind_dr2/outfiles/cf_dr2_1645749692.9475627_6059.err
#SBATCH --partition=priority
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --cpus-per-task=1
#SBATCH --mem=2000
#SBATCH --time=00:60:00
module load Anaconda3/2020.07
source activate $HOME/condaenv/
srun --ntasks=1 python /home/x25h971/canfind_dr2/files/canfind_v2.py 165290 1 &
srun --ntasks=1 python /home/x25h971/canfind_dr2/files/canfind_v2.py 182433 1 &
wait
