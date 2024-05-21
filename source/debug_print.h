#pragma once

#ifdef DEBUG_PRINT_ON
#define DEBUG_PRINTF(...)   {printf(__VA_ARGS__);fflush(stdout);}
#define DEBUG_BARRIER(mpi_comm)   MPI_Barrier(mpi_comm);
#else
//nothing to do
#define DEBUG_PRINTF(...)
#define DEBUG_BARRIER(mpi_comm)
#endif
