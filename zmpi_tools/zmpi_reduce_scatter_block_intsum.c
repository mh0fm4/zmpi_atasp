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


#include <string.h>
#include <mpi.h>

#include "z_pack.h"

#include "zmpi_tools.h"


#define DEFAULT_INT  0


static int _ZMPI_Reduce_scatter_block_intsum_accumulate(const int *sendbuf, int nsendprocs, int *sendprocs, int *recvbuf, int recvcount, int nrecvprocs, int *recvprocs, MPI_Comm comm)
{
  int i, j, size, rank;

  MPI_Win win;


  MPI_Comm_size(comm, &size);
  MPI_Comm_rank(comm, &rank);

  for (i = 0; i < recvcount; ++i) recvbuf[i] = DEFAULT_INT;

  MPI_Win_create(recvbuf, recvcount * sizeof(int), sizeof(int), MPI_INFO_NULL, comm, &win);
  MPI_Win_fence(MPI_MODE_NOSTORE|MPI_MODE_NOPRECEDE, win);

  if (nsendprocs >= 0)
  {
    for (j = 0; j < nsendprocs; ++j)
    {
      for (i = 0; i < recvcount; ++i) if (sendbuf[sendprocs[j] * recvcount + i] != DEFAULT_INT) break;

      if (i < recvcount) MPI_Accumulate((void *) &sendbuf[sendprocs[j] * recvcount], recvcount, MPI_INT, sendprocs[j], 0, recvcount, MPI_INT, MPI_SUM, win);
    }

  } else
  {
    for (j = 0; j < size; ++j)
    {
      for (i = 0; i < recvcount; ++i) if (sendbuf[j * recvcount + i] != DEFAULT_INT) break;

      if (i < recvcount) MPI_Accumulate((void *) &sendbuf[j * recvcount], recvcount, MPI_INT, j, 0, recvcount, MPI_INT, MPI_SUM, win);
    }
  }

  MPI_Win_fence(MPI_MODE_NOPUT|MPI_MODE_NOSUCCEED, win);
  MPI_Win_free(&win);

  return MPI_SUCCESS;
}


int ZMPI_Reduce_scatter_block_intsum_accumulate(const int *sendbuf, int *recvbuf, int recvcount, MPI_Comm comm) /* zmpi_func ZMPI_Reduce_scatter_block_intsum_accumulate */
{
  return _ZMPI_Reduce_scatter_block_intsum_accumulate(sendbuf, -1, NULL, recvbuf, recvcount, -1, NULL, comm);
}


int ZMPI_Reduce_scatter_block_intsum_proclists_isendirecv(const int *sendbuf, int nsendprocs, int *sendprocs, int *recvbuf, int recvcount, int nrecvprocs, int *recvprocs, MPI_Comm comm) /* zmpi_func ZMPI_Reduce_scatter_block_intsum_proclists_isendirecv */
{
  int i, j;

  int *recvbuf_full;
  MPI_Request *reqs;


  recvbuf_full = z_alloc(nrecvprocs * recvcount, sizeof(int));

  reqs = z_alloc(nsendprocs + nrecvprocs, sizeof(MPI_Request));

  for (j = 0; j < nrecvprocs; ++j) MPI_Irecv(&recvbuf_full[j * recvcount], recvcount, MPI_INT, recvprocs[j], 0, comm, &reqs[j]);

  for (j = 0; j < nsendprocs; ++j) MPI_Isend((void *) &sendbuf[sendprocs[j] * recvcount], recvcount, MPI_INT, sendprocs[j], 0, comm, &reqs[nrecvprocs + j]);

  MPI_Waitall(nsendprocs + nrecvprocs, reqs, MPI_STATUSES_IGNORE);

  for (i = 0; i < recvcount; ++i) recvbuf[i] = DEFAULT_INT;

  for (j = 0; j < nrecvprocs; ++j)
    for (i = 0; i < recvcount; ++i) recvbuf[i] += recvbuf_full[j * recvcount + i];

  z_free(reqs);

  z_free(recvbuf_full);

  return MPI_SUCCESS;
}


int ZMPI_Reduce_scatter_block_intsum_proclists_alltoallv(const int *sendbuf, int nsendprocs, int *sendprocs, int *recvbuf, int recvcount, int nrecvprocs, int *recvprocs, MPI_Comm comm) /* zmpi_func ZMPI_Reduce_scatter_block_intsum_proclists_alltoallv */
{
  int i, j, size, rank;

  int *recvbuf_full;
  int *scounts, *sdispls, *rcounts, *rdispls;


  MPI_Comm_size(comm, &size);
  MPI_Comm_rank(comm, &rank);

  recvbuf_full = z_alloc(nrecvprocs * recvcount, sizeof(int));

  scounts = z_alloc(4 * size, sizeof(int));
  sdispls = scounts + 1 * size;
  rcounts = scounts + 2 * size;
  rdispls = scounts + 3 * size;

  memset(scounts, 0, 4 * size * sizeof(int));

  for (j = 0; j < nrecvprocs; ++j)
  {
    rcounts[recvprocs[j]] = recvcount;
    rdispls[recvprocs[j]] = j * recvcount;
  }

  for (j = 0; j < nsendprocs; ++j)
  {
    scounts[sendprocs[j]] = recvcount;
    sdispls[sendprocs[j]] = sendprocs[j] * recvcount;
  }

  MPI_Alltoallv((void *) sendbuf, scounts, sdispls, MPI_INT, recvbuf_full, rcounts, rdispls, MPI_INT, comm);

  for (i = 0; i < recvcount; ++i) recvbuf[i] = DEFAULT_INT;

  for (j = 0; j < nrecvprocs; ++j)
    for (i = 0; i < recvcount; ++i) recvbuf[i] += recvbuf_full[j * recvcount + i];

  z_free(scounts);

  z_free(recvbuf_full);

  return MPI_SUCCESS;
}


int ZMPI_Reduce_scatter_block_intsum_proclists_accumulate(const int *sendbuf, int nsendprocs, int *sendprocs, int *recvbuf, int recvcount, int nrecvprocs, int *recvprocs, MPI_Comm comm) /* zmpi_func ZMPI_Reduce_scatter_block_intsum_proclists_accumulate */
{
  return _ZMPI_Reduce_scatter_block_intsum_accumulate(sendbuf, nsendprocs, sendprocs, recvbuf, recvcount, nrecvprocs, recvprocs, comm);
}
