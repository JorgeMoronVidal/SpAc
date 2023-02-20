#!/bin/sh
#PJM -L  "node=8"                          # Assign node 1 node
#PJM -L  "elapse=01:00:00"                 # Elapsed time limit 1 hour
#PJM -g  hp220345                          # group name
#PJM -x  PJM_LLIO_GFSCACHE=/vol0002        # volume names that job uses
#PJM -s                                    # Statistical information output
#PJM --mpi "max-proc-per-node=12"           # Upper limit of number of MPI process created at 1 node

export PARALLEL=4                          # Specify the number of thread
export OMP_NUM_THREADS=${PARALLEL}         # Specify the number of thread

echo "start"
mpiexec ./main_parallel
echo "end"
~                         