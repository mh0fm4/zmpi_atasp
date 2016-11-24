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


int ZMPI_Alltoallx_isendirecv_max_procs = 0;


int ZMPI_Alltoall_isendirecv(void* sendbuf, int sendcount, MPI_Datatype sendtype, void* recvbuf, int recvcount, MPI_Datatype recvtype, MPI_Comm comm)
{
  int i, j, k;
  int comm_size, comm_rank;
  const int tag = 0;

  int nreqs;
  MPI_Request *reqs;
  MPI_Status *stats;

  MPI_Aint sendtype_lb, sendtype_extent, recvtype_lb, recvtype_extent;

  int sub_comm_size, sub_comm_rank, sub_comm_base, ext_comm_size;


  MPI_Comm_size(comm, &comm_size);
  MPI_Comm_rank(comm, &comm_rank);

  if (ZMPI_Alltoallx_isendirecv_max_procs < 0) sub_comm_size = comm_size / -ZMPI_Alltoallx_isendirecv_max_procs;
  else if (ZMPI_Alltoallx_isendirecv_max_procs > 0) sub_comm_size = ZMPI_Alltoallx_isendirecv_max_procs;
  else sub_comm_size = comm_size;
  
  sub_comm_size = z_minmax(0, sub_comm_size, comm_size);

  sub_comm_rank = comm_rank % sub_comm_size;
  sub_comm_base = comm_rank - sub_comm_rank;

  ext_comm_size = ((double) comm_size / (double) sub_comm_size) + 0.5;
  ext_comm_size *= sub_comm_size;

  reqs = z_alloc(2 * sub_comm_size, sizeof(MPI_Request));
  stats = z_alloc(2 * sub_comm_size, sizeof(MPI_Status));

  MPI_Type_get_extent(sendtype, &sendtype_lb, &sendtype_extent);
  MPI_Type_get_extent(recvtype, &recvtype_lb, &recvtype_extent);

  k = 0;
  while (k < ext_comm_size)
  {
    nreqs = 0;

    for (i = 0; i < sub_comm_size; ++i)
    {
      j = (((sub_comm_rank + i) % sub_comm_size) + sub_comm_base - k + ext_comm_size) % ext_comm_size;
      if (j < comm_size)
      {
        MPI_Irecv(((char *) recvbuf) + (j * recvcount * recvtype_extent), recvcount, recvtype, j, tag, comm, &reqs[nreqs]);
        ++nreqs;
      }
    }

    for (i = 0; i < sub_comm_size; ++i)
    {
      j = (((sub_comm_rank + i) % sub_comm_size) + sub_comm_base + k + ext_comm_size) % ext_comm_size;
      if (j < comm_size)
      {
        MPI_Isend(((char *) sendbuf) + (j * sendcount * sendtype_extent), sendcount, sendtype, j, tag, comm, &reqs[nreqs]);
        ++nreqs;
      }
    }

    MPI_Waitall(nreqs, reqs, stats);

    k += sub_comm_size;
  }

  z_free(reqs);
  z_free(stats);

  return MPI_SUCCESS;
}


int ZMPI_Alltoallv_isendirecv(void* sendbuf, int *sendcounts, int *sdispls, MPI_Datatype sendtype, void* recvbuf, int *recvcounts, int *rdispls, MPI_Datatype recvtype, MPI_Comm comm)
{
  int i, j, k;
  int comm_size, comm_rank;
  const int tag = 0;

  int nreqs;
  MPI_Request *reqs;
  MPI_Status *stats;

  MPI_Aint sendtype_lb, sendtype_extent, recvtype_lb, recvtype_extent;

  int sub_comm_size, sub_comm_rank, sub_comm_base, ext_comm_size;


  MPI_Comm_size(comm, &comm_size);
  MPI_Comm_rank(comm, &comm_rank);

  MPI_Type_get_extent(sendtype, &sendtype_lb, &sendtype_extent);
  MPI_Type_get_extent(recvtype, &recvtype_lb, &recvtype_extent);

  if (ZMPI_Alltoallx_isendirecv_max_procs < 0) sub_comm_size = comm_size / -ZMPI_Alltoallx_isendirecv_max_procs;
  else if (ZMPI_Alltoallx_isendirecv_max_procs > 0) sub_comm_size = ZMPI_Alltoallx_isendirecv_max_procs;
  else sub_comm_size = comm_size;
  
  sub_comm_size = z_minmax(0, sub_comm_size, comm_size);

  sub_comm_rank = comm_rank % sub_comm_size;
  sub_comm_base = comm_rank - sub_comm_rank;

  ext_comm_size = ((double) comm_size / (double) sub_comm_size) + 0.5;
  ext_comm_size *= sub_comm_size;

  reqs = z_alloc(2 * sub_comm_size, sizeof(MPI_Request));
  stats = z_alloc(2 * sub_comm_size, sizeof(MPI_Status));

  k = 0;
  while (k < ext_comm_size)
  {
    nreqs = 0;

    for (i = 0; i < sub_comm_size; ++i)
    {
      j = (((sub_comm_rank + i) % sub_comm_size) + sub_comm_base - k + ext_comm_size) % ext_comm_size;
      if (j < comm_size && recvcounts[j] > 0)
      {
        MPI_Irecv(((char *) recvbuf) + (rdispls[j] * recvtype_extent), recvcounts[j], recvtype, j, tag, comm, &reqs[nreqs]);
        ++nreqs;
      }
    }

    for (i = 0; i < sub_comm_size; ++i)
    {
      j = (((sub_comm_rank + i) % sub_comm_size) + sub_comm_base + k + ext_comm_size) % ext_comm_size;
      if (j < comm_size && sendcounts[j] > 0)
      {
        MPI_Isend(((char *) sendbuf) + (sdispls[j] * sendtype_extent), sendcounts[j], sendtype, j, tag, comm, &reqs[nreqs]);
        ++nreqs;
      }
    }

    MPI_Waitall(nreqs, reqs, stats);

    k += sub_comm_size;
  }

  z_free(reqs);
  z_free(stats);

  return MPI_SUCCESS;
}


int ZMPI_Alltoallw_isendirecv(void* sendbuf, int sendcounts[], int sdispls[], MPI_Datatype sendtypes[], void *recvbuf, int recvcounts[], int rdispls[], MPI_Datatype recvtypes[], MPI_Comm comm)
{
  int i, j, k;
  int comm_size, comm_rank;
  const int tag = 0;

  int nreqs;
  MPI_Request *reqs;
  MPI_Status *stats;

  int sub_comm_size, sub_comm_rank, sub_comm_base, ext_comm_size;


  MPI_Comm_size(comm, &comm_size);
  MPI_Comm_rank(comm, &comm_rank);

  if (ZMPI_Alltoallx_isendirecv_max_procs < 0) sub_comm_size = comm_size / -ZMPI_Alltoallx_isendirecv_max_procs;
  else if (ZMPI_Alltoallx_isendirecv_max_procs > 0) sub_comm_size = ZMPI_Alltoallx_isendirecv_max_procs;
  else sub_comm_size = comm_size;
  
  sub_comm_size = z_minmax(0, sub_comm_size, comm_size);

  sub_comm_rank = comm_rank % sub_comm_size;
  sub_comm_base = comm_rank - sub_comm_rank;

  ext_comm_size = ((double) comm_size / (double) sub_comm_size) + 0.5;
  ext_comm_size *= sub_comm_size;

  reqs = z_alloc(2 * sub_comm_size, sizeof(MPI_Request));
  stats = z_alloc(2 * sub_comm_size, sizeof(MPI_Status));

  k = 0;
  while (k < ext_comm_size)
  {
    nreqs = 0;

    for (i = 0; i < sub_comm_size; ++i)
    {
      j = (((sub_comm_rank + i) % sub_comm_size) + sub_comm_base - k + ext_comm_size) % ext_comm_size;
      if (j < comm_size && recvcounts[j] > 0)
      {
        MPI_Irecv(((char *) recvbuf) + rdispls[j], recvcounts[j], recvtypes[j], j, tag, comm, &reqs[nreqs]);
        ++nreqs;
      }
    }

    for (i = 0; i < sub_comm_size; ++i)
    {
      j = (((sub_comm_rank + i) % sub_comm_size) + sub_comm_base + k + ext_comm_size) % ext_comm_size;
      if (j < comm_size && sendcounts[j] > 0)
      {
        MPI_Isend(((char *) sendbuf) + sdispls[j], sendcounts[j], sendtypes[j], j, tag, comm, &reqs[nreqs]);
        ++nreqs;
      }
    }

    MPI_Waitall(nreqs, reqs, stats);

    k += sub_comm_size;
  }

  z_free(reqs);
  z_free(stats);

  return MPI_SUCCESS;
}
