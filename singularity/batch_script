#!/bin/bash
#SBATCH -n 4        # Number of cores
#SBATCH -t 0-00:45  # Runtime in D-HH:MM
#SBATCH --job-name=singularity
#SBATCH -p regular1


#load singularity and MPI modules
module load singularity/3.4.1
module load intel/2021.2
module load openmpi3


echo "Starting singularity on host $HOSTNAME"

#this binds mpi, so we do not need to 
#srun -n 4 singularity exec ithacafv.sif /bin/bash Of.sh

#without MPI
singularity exec ithicafv.sif /bin/bash Of.sh

#Other valid 

#mpi run not recomeneded
#mpirun -n 4 singularity exec ithacafv.sif /bin/bash Of.sh

#mpiexec -n 4 singularity exec ithacafv.sif /bin/bash Of.sh

#singularity exec ithicafv.sif /bin/bash Of.sh

echo "Completed singularity on host $HOSTNAME"


