#!/bin/bash

#SBATCH --partition=bgmp        ### Partition (like a queue in PBS)
#SBATCH --job-name=dedup_%j      ### Job Name
#SBATCH --output=dedup_%j.out         ### File in which to store job output
#SBATCH --error=dedup_%j.err          ### File in which to store job error messages
#SBATCH --time=0-00:30:00       ### Wall clock time limit in Days-HH:MM:SS
#SBATCH --nodes=1               ### Number of nodes needed for the job
#SBATCH --ntasks-per-node=1     ### Number of tasks to be launched per Node
#SBATCH --cpus-per-task=1     ### Number of cpus per task
#SBATCH --account=bgmp      ### Account used for job submission

/usr/bin/time -v ./Ferina_deduper.py -f /projects/bgmp/shared/deduper/C1_SE_uniqAlign.sam -o deduplicated.sam  -u STL96.txt