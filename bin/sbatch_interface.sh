#!/bin/env bash

#SBATCH --cpus-per-task=1
#SBATCH --job-name=GZ
#SBATCH -t 0:20:00
#SBATCH -o /dev/null
#SBATCH -p atlas_all
#SBATCH --mem=1400

# parameters:
#------------
# 1) start_pos_x
# 2) start_pos_y
# 3) start_pos_z
# 4) Z
# 5) energy
# 6) save filepath

cd ../python
python -c "import Solve; Solve.trajectory([$1, $2, $3], $4, $5, savefile=\"$6\")"
