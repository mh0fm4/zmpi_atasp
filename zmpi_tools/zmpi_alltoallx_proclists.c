/*
 *  Copyright (C) 2012, 2013, 2014, 2015, 2016 Michael Hofmann, Chemnitz University of Technology
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


int ZMPI_Alltoall_proclists_isendirecv(void *sendbuf, int sendcount, MPI_Datatype sendtype, int nsendprocs, int *sendprocs, void *recvbuf, int recvcount, MPI_Datatype recvtype, int nrecvprocs, int *recvprocs, MPI_Comm comm) /* zmpi_func ZMPI_Alltoall_proclists_isendirecv */
{
  int i, j;

  const int tag = 0;

  int nreqs;
  MPI_Request *reqs;
  MPI_Status *stats;

  MPI_Aint sendtype_lb, sendtype_extent, recvtype_lb, recvtype_extent;


  reqs = z_alloc(nrecvprocs + nsendprocs, sizeof(MPI_Request) + sizeof(MPI_Status));
  stats = (MPI_Status *) (reqs + nrecvprocs + nsendprocs);

  MPI_Type_get_extent(sendtype, &sendtype_lb, &sendtype_extent);
  MPI_Type_get_extent(recvtype, &recvtype_lb, &recvtype_extent);

  nreqs = 0;

  for (i = 0; i < nrecvprocs; ++i)
  {
    j = recvprocs[i];
    MPI_Irecv(((char *) recvbuf) + (j * recvcount * recvtype_extent), recvcount, recvtype, j, tag, comm, &reqs[nreqs]);
    ++nreqs;
  }

  for (i = 0; i < nsendprocs; ++i)
  {
    j = sendprocs[i];
    MPI_Isend(((char *) sendbuf) + (j * sendcount * sendtype_extent), sendcount, sendtype, j, tag, comm, &reqs[nreqs]);
    ++nreqs;
  }

  MPI_Waitall(nreqs, reqs, stats);

  z_free(reqs);

  return MPI_SUCCESS;
}


int ZMPI_Alltoall_proclists(void *sendbuf, int sendcount, MPI_Datatype sendtype, int nsendprocs, int *sendprocs, void *recvbuf, int recvcount, MPI_Datatype recvtype, int nrecvprocs, int *recvprocs, MPI_Comm comm) /* zmpi_func ZMPI_Alltoall_proclists */
{
  return ZMPI_Alltoall_proclists_isendirecv(sendbuf, sendcount, sendtype, nsendprocs, sendprocs, recvbuf, recvcount, recvtype, nrecvprocs, recvprocs, comm);
}


int ZMPI_Alltoallv_proclists_isendirecv(void *sendbuf, int *sendcounts, int *senddispls, MPI_Datatype sendtype, int nsendprocs, int *sendprocs, void *recvbuf, int *recvcounts, int *recvdispls, MPI_Datatype recvtype, int nrecvprocs, int *recvprocs, MPI_Comm comm) /* zmpi_func ZMPI_Alltoallv_proclists_isendirecv */
{
  int i, j;

  const int tag = 0;

  int nreqs;
  MPI_Request *reqs;
  MPI_Status *stats;

  MPI_Aint sendtype_lb, sendtype_extent, recvtype_lb, recvtype_extent;


  reqs = z_alloc(nrecvprocs + nsendprocs, sizeof(MPI_Request) + sizeof(MPI_Status));
  stats = (MPI_Status *) (reqs + nrecvprocs + nsendprocs);

  MPI_Type_get_extent(sendtype, &sendtype_lb, &sendtype_extent);
  MPI_Type_get_extent(recvtype, &recvtype_lb, &recvtype_extent);

  nreqs = 0;

  for (i = 0; i < nrecvprocs; ++i)
  {
    j = recvprocs[i];
    if (recvcounts[j] > 0)
    {
      MPI_Irecv(((char *) recvbuf) + (recvdispls[j] * recvtype_extent), recvcounts[j], recvtype, j, tag, comm, &reqs[nreqs]);
      ++nreqs;
    }
  }

  for (i = 0; i < nsendprocs; ++i)
  {
    j = sendprocs[i];
    if (sendcounts[j] > 0)
    {
      MPI_Isend(((char *) sendbuf) + (senddispls[j] * sendtype_extent), sendcounts[j], sendtype, j, tag, comm, &reqs[nreqs]);
      ++nreqs;
    }
  }

  MPI_Waitall(nreqs, reqs, stats);

  z_free(reqs);

  return MPI_SUCCESS;
}


int ZMPI_Alltoallv_proclists(void *sendbuf, int *sendcounts, int *senddispls, MPI_Datatype sendtype, int nsendprocs, int *sendprocs, void *recvbuf, int *recvcounts, int *recvdispls, MPI_Datatype recvtype, int nrecvprocs, int *recvprocs, MPI_Comm comm) /* zmpi_func ZMPI_Alltoallv_proclists */
{
  return ZMPI_Alltoallv_proclists_isendirecv(sendbuf, sendcounts, senddispls, sendtype, nsendprocs, sendprocs, recvbuf, recvcounts, recvdispls, recvtype, nrecvprocs, recvprocs, comm);
}


int ZMPI_Alltoallw_proclists_isendirecv(void *sendbuf, int *sendcounts, int *senddispls, MPI_Datatype *sendtypes, int nsendprocs, int *sendprocs, void *recvbuf, int *recvcounts, int *recvdispls, MPI_Datatype *recvtypes, int nrecvprocs, int *recvprocs, MPI_Comm comm) /* zmpi_func ZMPI_Alltoallw_proclists_isendirecv */
{
  int i, j;

  const int tag = 0;

  int nreqs;
  MPI_Request *reqs;
  MPI_Status *stats;

  MPI_Aint sendtype_lb, sendtype_extent, recvtype_lb, recvtype_extent;


  reqs = z_alloc(nrecvprocs + nsendprocs, sizeof(MPI_Request) + sizeof(MPI_Status));
  stats = (MPI_Status *) (reqs + nrecvprocs + nsendprocs);

  nreqs = 0;

  for (i = 0; i < nrecvprocs; ++i)
  {
    j = recvprocs[i];
    if (recvcounts[j] > 0)
    {
      MPI_Type_get_extent(recvtypes[j], &recvtype_lb, &recvtype_extent);
      MPI_Irecv(((char *) recvbuf) + recvdispls[j], recvcounts[j], recvtypes[j], j, tag, comm, &reqs[nreqs]);
      ++nreqs;
    }
  }

  for (i = 0; i < nsendprocs; ++i)
  {
    j = sendprocs[i];
    if (sendcounts[j] > 0)
    {
      MPI_Type_get_extent(sendtypes[j], &sendtype_lb, &sendtype_extent);
      MPI_Isend(((char *) sendbuf) + senddispls[j], sendcounts[j], sendtypes[j], j, tag, comm, &reqs[nreqs]);
      ++nreqs;
    }
  }

  MPI_Waitall(nreqs, reqs, stats);

  z_free(reqs);

  return MPI_SUCCESS;
}


int ZMPI_Alltoallw_proclists(void *sendbuf, int *sendcounts, int *senddispls, MPI_Datatype *sendtypes, int nsendprocs, int *sendprocs, void *recvbuf, int *recvcounts, int *recvdispls, MPI_Datatype *recvtypes, int nrecvprocs, int *recvprocs, MPI_Comm comm) /* zmpi_func ZMPI_Alltoallw_proclists */
{
  return ZMPI_Alltoallw_proclists_isendirecv(sendbuf, sendcounts, senddispls, sendtypes, nsendprocs, sendprocs, recvbuf, recvcounts, recvdispls, recvtypes, nrecvprocs, recvprocs, comm);
}
