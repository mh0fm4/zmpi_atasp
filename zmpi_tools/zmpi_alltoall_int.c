/*
 *  Copyright (C) 2012, 2013, 2014, 2015 Michael Hofmann, Chemnitz University of Technology
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

#include "z_pack.h"

#include "zmpi_tools.h"


#define DEFAULT_INT  0


int ZMPI_Alltoall_int_alltoall(int *sendbuf, int *recvbuf, MPI_Comm comm) /* zmpi_func ZMPI_Alltoall_int_alltoall */
{
  return MPI_Alltoall(sendbuf, 1, MPI_INT, recvbuf, 1, MPI_INT, comm);
}


int ZMPI_Alltoall_int_2step(int *sendbuf, int *recvbuf, MPI_Comm comm) /* zmpi_func ZMPI_Alltoall_int_2step */
{
  return ZMPI_Alltoall_2step_int(sendbuf, 1, MPI_INT, recvbuf, 1, MPI_INT, comm);
}


static int _ZMPI_Alltoall_int_proclists_put(int alloc_mem, int nphases, int *sendbuf, int nsprocs, int *sprocs, int *recvbuf, int nrprocs, int *rprocs, MPI_Comm comm)
{
  int i, p, size, rank, *rcounts_put;

  MPI_Win win;


  MPI_Comm_size(comm, &size);
  MPI_Comm_rank(comm, &rank);

  if (alloc_mem) MPI_Alloc_mem(size * sizeof(int), MPI_INFO_NULL, &rcounts_put);
  else rcounts_put = recvbuf;

  if (nrprocs >= 0)
    for (i = 0; i < nrprocs; ++i) rcounts_put[rprocs[i]] = DEFAULT_INT;
  else
    for (i = 0; i < size; ++i) rcounts_put[i] = DEFAULT_INT;

  MPI_Win_create(rcounts_put, size * sizeof(int), sizeof(int), MPI_INFO_NULL, comm, &win);
  MPI_Win_fence(MPI_MODE_NOSTORE|MPI_MODE_NOPRECEDE, win);

  for (p = 0; p < nphases; ++p)
  {
/*    printf("%d: phase = %d of %d\n", rank, p, nphases);*/
  
    if (rank % nphases == p)
    {
      if (nsprocs >= 0)
      {
        for (i = 0; i < nsprocs; ++i)
          if (sendbuf[sprocs[i]] != DEFAULT_INT) MPI_Put(&sendbuf[sprocs[i]], 1, MPI_INT, sprocs[i], rank, 1, MPI_INT, win);

      } else
      {
        for (i = 0; i < size; ++i)
          if (sendbuf[i] != DEFAULT_INT) MPI_Put(&sendbuf[i], 1, MPI_INT, i, rank, 1, MPI_INT, win);
      }
    }

    if (p < nphases - 1) MPI_Win_fence(0, win);
  }

  MPI_Win_fence(MPI_MODE_NOPUT|MPI_MODE_NOSUCCEED, win);
  MPI_Win_free(&win);

  if (alloc_mem)
  {
    if (nrprocs >= 0)
      for (i = 0; i < nrprocs; ++i) recvbuf[rprocs[i]] = rcounts_put[rprocs[i]];
    else
      for (i = 0; i < size; ++i) recvbuf[i] = rcounts_put[i];

    MPI_Free_mem(rcounts_put);    
  }

  return MPI_SUCCESS;
}


int ZMPI_Alltoall_int_put(int *sendbuf, int *recvbuf, MPI_Comm comm) /* zmpi_func ZMPI_Alltoall_int_put */
{
  return _ZMPI_Alltoall_int_proclists_put(0, 1, sendbuf, -1, NULL, recvbuf, -1, NULL, comm);
}


int ZMPI_Alltoall_int_put_alloc(int *sendbuf, int *recvbuf, MPI_Comm comm) /* zmpi_func ZMPI_Alltoall_int_put_alloc */
{
  return _ZMPI_Alltoall_int_proclists_put(1, 1, sendbuf, -1, NULL, recvbuf, -1, NULL, comm);
}


int ZMPI_Alltoall_int_put_2phases(int *sendbuf, int *recvbuf, MPI_Comm comm) /* zmpi_func ZMPI_Alltoall_int_put_2phases */
{
  return _ZMPI_Alltoall_int_proclists_put(0, 2, sendbuf, -1, NULL, recvbuf, -1, NULL, comm);
}


int ZMPI_Alltoall_int_put_2phases_alloc(int *sendbuf, int *recvbuf, MPI_Comm comm) /* zmpi_func ZMPI_Alltoall_int_put_2phases_alloc */
{
  return _ZMPI_Alltoall_int_proclists_put(1, 2, sendbuf, -1, NULL, recvbuf, -1, NULL, comm);
}


int ZMPI_Alltoall_int_put_3phases(int *sendbuf, int *recvbuf, MPI_Comm comm) /* zmpi_func ZMPI_Alltoall_int_put_3phases */
{
  return _ZMPI_Alltoall_int_proclists_put(0, 3, sendbuf, -1, NULL, recvbuf, -1, NULL, comm);
}


int ZMPI_Alltoall_int_put_3phases_alloc(int *sendbuf, int *recvbuf, MPI_Comm comm) /* zmpi_func ZMPI_Alltoall_int_put_3phases_alloc */
{
  return _ZMPI_Alltoall_int_proclists_put(1, 3, sendbuf, -1, NULL, recvbuf, -1, NULL, comm);
}


int ZMPI_Alltoall_int_proclists_isendirecv(int *sendbuf, int nsprocs, int *sprocs, int *recvbuf, int nrprocs, int *rprocs, MPI_Comm comm) /* zmpi_func ZMPI_Alltoall_int_proclists_isendirecv */
{
  return ZMPI_Alltoall_proclists_isendirecv(sendbuf, 1, MPI_INT, nsprocs, sprocs, recvbuf, 1, MPI_INT, nrprocs, rprocs, comm);
}


int ZMPI_Alltoall_int_proclists_alltoallv(int *sendbuf, int nsprocs, int *sprocs, int *recvbuf, int nrprocs, int *rprocs, MPI_Comm comm) /* zmpi_func ZMPI_Alltoall_int_proclists_alltoallv */
{
  int i, size;

  int *scounts2, *sdispls2, *rcounts2, *rdispls2;


  MPI_Comm_size(comm, &size);

  scounts2 = z_alloc(4 * size, sizeof(int));
  sdispls2 = scounts2 + 1 * size;
  rcounts2 = scounts2 + 2 * size;
  rdispls2 = scounts2 + 3 * size;

  for (i = 0; i < size; ++i)
  {
    scounts2[i] = rcounts2[i] = DEFAULT_INT;
    sdispls2[i] = rdispls2[i] = i;
    recvbuf[i] = 0;
  }

  for (i = 0; i < nsprocs; ++i) scounts2[sprocs[i]] = 1;
  for (i = 0; i < nrprocs; ++i) rcounts2[rprocs[i]] = 1;

  MPI_Alltoallv(sendbuf, scounts2, sdispls2, MPI_INT, recvbuf, rcounts2, rdispls2, MPI_INT, comm);

  z_free(scounts2);

  return MPI_SUCCESS;
}


int ZMPI_Alltoall_int_proclists_put(int *sendbuf, int nsprocs, int *sprocs, int *recvbuf, int nrprocs, int *rprocs, MPI_Comm comm) /* zmpi_func ZMPI_Alltoall_int_proclists_put */
{
  return _ZMPI_Alltoall_int_proclists_put(0, 1, sendbuf, nsprocs, sprocs, recvbuf, nrprocs, rprocs, comm);
}


int ZMPI_Alltoall_int_proclists_put_alloc(int *sendbuf, int nsprocs, int *sprocs, int *recvbuf, int nrprocs, int *rprocs, MPI_Comm comm) /* zmpi_func ZMPI_Alltoall_int_proclists_put_alloc */
{
  return _ZMPI_Alltoall_int_proclists_put(1, 1, sendbuf, nsprocs, sprocs, recvbuf, nrprocs, rprocs, comm);
}


int ZMPI_Alltoall_int_proclists_put_2phases(int *sendbuf, int nsprocs, int *sprocs, int *recvbuf, int nrprocs, int *rprocs, MPI_Comm comm) /* zmpi_func ZMPI_Alltoall_int_proclists_put_2phases */
{
  return _ZMPI_Alltoall_int_proclists_put(0, 2, sendbuf, nsprocs, sprocs, recvbuf, nrprocs, rprocs, comm);
}


int ZMPI_Alltoall_int_proclists_put_2phases_alloc(int *sendbuf, int nsprocs, int *sprocs, int *recvbuf, int nrprocs, int *rprocs, MPI_Comm comm) /* zmpi_func ZMPI_Alltoall_int_proclists_put_2phases_alloc */
{
  return _ZMPI_Alltoall_int_proclists_put(1, 2, sendbuf, nsprocs, sprocs, recvbuf, nrprocs, rprocs, comm);
}
