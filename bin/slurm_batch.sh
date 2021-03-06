#!/bin/env bash

#SBATCH --cpus-per-task=1
#SBATCH --job-name=GZ
#SBATCH -t 8:00:00
#SBATCH -o /beegfs/DATA/atlas/ealbin/GZ/jobs/dump/slurm-%j.out
#SBATCH -p atlas_all
#SBATCH --mem=5000

# intended use:
# $ sbatch slurm_batch.sh job_XXXXXX.txt

cd ../python
while read job; do
    echo $job
    param=( $job )
    # parameters:
    #------------
    # 0) start_pos_x
    # 1) start_pos_y
    # 2) start_pos_z
    # 3) Z
    # 4) energy
    # 5) save filepath
    python -c "import Solve; Solve.trajectory([${param[0]}, ${param[1]}, ${param[2]}], ${param[3]}, ${param[4]}, savefile=\"${param[5]}\")"
done < "$1"

rm $1
