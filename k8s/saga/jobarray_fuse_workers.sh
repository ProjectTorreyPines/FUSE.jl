#!/bin/bash
#SBATCH --job-name=fuse_workers
#SBATCH --nodes=1 --ntasks=1
#SBATCH --time=4:00:00
#SBATCH --output=fuse_worker_%a.out
#SBATCH --error=fuse_worker_%a.err
#SBATCH --array=1-120

echo "FUSE WORKER $SLURM_ARRAY_TASK_ID"
echo "BROKER $FUSE_BROKER"

julia -e 'import FUSE; FUSE.worker_start()'

echo "DONE"