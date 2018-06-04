#!/bin/env bash

#SBATCH --cpus-per-task=1
#SBATCH --job-name=GZ
#SBATCH -t 8:00:00
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
while read job; do
    param=( $job )
    python -c "import Solve; Solve.trajectory([${param[0]}, ${param[1]}, ${param[2]}], ${param[3]}, ${param[4]}, savefile=\"${param[5]}\")"
done