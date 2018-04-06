/*
 *  Copyright (C) 2012, 2013, 2014, 2015, 2016, 2017, 2018 Michael Hofmann, Chemnitz University of Technology
 *  
 *  This file is part of the ZMPI All-to-all Specific Library.
 *  
 *  The ZMPI-ATASP is free software: you can redistribute it and/or
 *  modify it under the terms of the GNU Lesser Public License as published
 *  by the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *  
 *  The ZMPI-ATASP is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Lesser Public License for more details.
 *  
 *  You should have received a copy of the GNU Lesser Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */


#include <stdio.h>
#include <stdlib.h>

#include <mpi.h>

#include "zmpi_atasp.h"


#define TPROC_TYPE   0
#define TPROC_EXDEF  0

#define ATASP_TYPE   0


/* distribution function variant #1: one or none target process per data element */
int tproc_int(void *b,
#if MPI_VERSION >= 3
  MPI_Count
#else
  ZMPI_Count
#endif
  x, void *data)
{
  /* return target process rank or MPI_PROC_NULL */
  return (((int *) b)[x] < 0)?MPI_PROC_NULL:((int *) b)[x];
}


/* distribution function variant #2: one or none target process per data element and an individual copy for each process */
int tproc_mod_int(void *b,
#if MPI_VERSION >= 3
  MPI_Count
#else
  ZMPI_Count
#endif
  x, void *data, void *mod)
{
  /* create an individual copy of the data element */
  if (mod) ((int *) mod)[0] = ((int *) b)[x];

  /* return target process rank or MPI_PROC_NULL */
  return (((int *) b)[x] < 0)?MPI_PROC_NULL:((int *) b)[x];
}


/* distribution function variant #3: several target processes per data element */
void tprocs_int(void *b,
#if MPI_VERSION >= 3
  MPI_Count
#else
  ZMPI_Count
#endif
  x, void *data, int *nprocs, int *procs)
{
  /* no target process */
  if (((int *) b)[x] < 0)
  {
    /* set number of target processes */
    *nprocs = 0;
    return;
  }

  /* set number of target processes */
  *nprocs = 1;

  /* set target process */
  procs[0] = ((int *) b)[x];
}


/* distribution function variant #4: several target processes per data element and an individual copy for each process */
void tprocs_mod_int(void *b,
#if MPI_VERSION >= 3
  MPI_Count
#else
  ZMPI_Count
#endif
  x, void *data, int *nprocs, int *procs, void *mods)
{
  /* no target process */
  if (((int *) b)[x] < 0)
  {
    /* set number of target processes */
    *nprocs = 0;
    return;
  }

  /* set number of target processes */
  *nprocs = 1;

  /* set target process */
  procs[0] = ((int *) b)[x];

  /* create an individual copy of the data element */
  if (mods) ((int *) mods)[0] = ((int *) b)[x];
}


#if TPROC_EXDEF
ZMPI_TPROC_EXDEF_DEFINE_TPROC(exdef_tproc_int, tproc_int);
ZMPI_TPROC_EXDEF_DEFINE_TPROC_MOD(exdef_tproc_mod_int, tproc_mod_int);
ZMPI_TPROC_EXDEF_DEFINE_TPROCS(exdef_tprocs_int, tprocs_int);
ZMPI_TPROC_EXDEF_DEFINE_TPROCS_MOD(exdef_tprocs_mod_int, tprocs_mod_int);
#endif


int main(int argc, char **argv)
{
  const int nlocal = 10000;
  const int nlocal_max = nlocal * 1.1;

  int i, received, total_received;
  int *sbuf, *rbuf;

#if MPI_VERSION >= 3
  MPI_Status status;
#else
  ZMPI_Status status;
#endif

  MPI_Comm comm;
  int comm_size, comm_rank;
  
  ZMPI_Tproc tproc;


  MPI_Init(&argc, &argv);

  comm = MPI_COMM_WORLD;

  MPI_Comm_size(comm, &comm_size);
  MPI_Comm_rank(comm, &comm_rank);

  sbuf = malloc(nlocal * sizeof(int));
  rbuf = malloc(nlocal_max * sizeof(int));

  srandom((comm_rank + 1) * 2501);

  for (i = 0; i < nlocal; ++i) sbuf[i] = (random() % (comm_size + 1)) - 1;

  switch (TPROC_TYPE)
  {
    case 0:
#if TPROC_EXDEF
      ZMPI_Tproc_create_tproc(&tproc, tproc_int, NULL, exdef_tproc_int);
#else
      ZMPI_Tproc_create_tproc(&tproc, tproc_int, NULL, ZMPI_TPROC_EXDEF_NULL);
#endif
      break;
    case 1:
#if TPROC_EXDEF
      ZMPI_Tproc_create_tproc_mod(&tproc, tproc_mod_int, NULL, exdef_tproc_mod_int);
#else
      ZMPI_Tproc_create_tproc_mod(&tproc, tproc_mod_int, NULL, ZMPI_TPROC_EXDEF_NULL);
#endif
      break;
    case 2:
#if TPROC_EXDEF
      ZMPI_Tproc_create_tprocs(&tproc, 1, tprocs_int, NULL, exdef_tprocs_int);
#else
      ZMPI_Tproc_create_tprocs(&tproc, 1, tprocs_int, NULL, ZMPI_TPROC_EXDEF_NULL);
#endif
      break;
    case 3:
#if TPROC_EXDEF
      ZMPI_Tproc_create_tprocs_mod(&tproc, 1, tprocs_mod_int, NULL, exdef_tprocs_mod_int);
#else
      ZMPI_Tproc_create_tprocs_mod(&tproc, 1, tprocs_mod_int, NULL, ZMPI_TPROC_EXDEF_NULL);
#endif
      break;
  }

  switch (ATASP_TYPE)
  {
    case 0:
      ZMPI_Alltoall_specific_type = ZMPI_ALLTOALL_SPECIFIC_TYPE_ALLTOALLV;
      break;
    case 1:
      ZMPI_Alltoall_specific_type = ZMPI_ALLTOALL_SPECIFIC_TYPE_ALLTOALLW;
      break;
    case 2:
      ZMPI_Alltoall_specific_type = ZMPI_ALLTOALL_SPECIFIC_TYPE_PUT;
      break;
    case 3:
      ZMPI_Alltoall_specific_type = ZMPI_ALLTOALL_SPECIFIC_TYPE_SENDRECV;
      break;
    default:
      ZMPI_Alltoall_specific_type = ZMPI_ALLTOALL_SPECIFIC_TYPE_DEFAULT;
      break;
  }

  ZMPI_Alltoall_specific(sbuf, nlocal, MPI_INT, rbuf, nlocal_max, MPI_INT, tproc, NULL, comm, &status);

#if MPI_VERSION >= 3
  MPI_Get_elements(&status, MPI_INT, &received);
#else
  ZMPI_Get_elements(&status, MPI_INT, &received);
#endif

  ZMPI_Tproc_free(&tproc);

  printf("%d: received = %d\n", comm_rank, received);

  MPI_Reduce(&received, &total_received, 1, MPI_INT, MPI_SUM, 0, comm);
  if (comm_rank == 0) printf("%d: total received = %d\n", comm_rank, total_received);

  free(sbuf);
  free(rbuf);

  MPI_Finalize();

  return 0;
}
