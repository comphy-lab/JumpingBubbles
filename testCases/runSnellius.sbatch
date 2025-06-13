#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks=192
#SBATCH --time=100:00:00
#SBATCH --partition=genoa
#SBATCH --mail-type=ALL
#SBATCH --mail-user=v.sanjay@utwente.nl

srun --mpi=pmi2 JumpingBubbles
