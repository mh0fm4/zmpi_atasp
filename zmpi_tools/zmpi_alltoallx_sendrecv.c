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


int ZMPI_Alltoall_sendrecv(void* sendbuf, int sendcount, MPI_Datatype sendtype, void* recvbuf, int recvcount, MPI_Datatype recvtype, MPI_Comm comm)
{
  int s, r;
  int comm_size, comm_rank;
  MPI_Status stat;
  const int tag = 0;

  MPI_Aint sendtype_lb, sendtype_extent, recvtype_lb, recvtype_extent;


  MPI_Type_get_extent(sendtype, &sendtype_lb, &sendtype_extent);
  MPI_Type_get_extent(recvtype, &recvtype_lb, &recvtype_extent);

  MPI_Comm_size(comm, &comm_size);
  MPI_Comm_rank(comm, &comm_rank);

  for (s = 0; s < comm_size; ++s)
  for (r = 0; r < comm_size; ++r)
  MPI_Sendrecv(((char *) sendbuf) + (sendcount * r * sendtype_extent), sendcount, sendtype, (comm_rank == s)?r:MPI_PROC_NULL, tag,
               ((char *) recvbuf) + (recvcount * s * recvtype_extent), recvcount, recvtype, (comm_rank == r)?s:MPI_PROC_NULL, tag, comm, &stat);

  return MPI_SUCCESS;
}


int ZMPI_Alltoallv_sendrecv(void* sendbuf, int *sendcounts, int *sdispls, MPI_Datatype sendtype, void* recvbuf, int *recvcounts, int *rdispls, MPI_Datatype recvtype, MPI_Comm comm)
{
  int s, r;
  int comm_size, comm_rank;
  MPI_Status stat;
  const int tag = 0;

  MPI_Aint sendtype_lb, sendtype_extent, recvtype_lb, recvtype_extent;


  MPI_Type_get_extent(sendtype, &sendtype_lb, &sendtype_extent);
  MPI_Type_get_extent(recvtype, &recvtype_lb, &recvtype_extent);

  MPI_Comm_size(comm, &comm_size);
  MPI_Comm_rank(comm, &comm_rank);

  for (s = 0; s < comm_size; ++s)
  for (r = 0; r < comm_size; ++r)
  MPI_Sendrecv(((char *) sendbuf) + (sdispls[r] * sendtype_extent), sendcounts[r], sendtype, (comm_rank == s)?r:MPI_PROC_NULL, tag,
               ((char *) recvbuf) + (rdispls[s] * recvtype_extent), recvcounts[s], recvtype, (comm_rank == r)?s:MPI_PROC_NULL, tag, comm, &stat);

  return MPI_SUCCESS;
}


int ZMPI_Alltoallw_sendrecv(void* sendbuf, int sendcounts[], int sdispls[], MPI_Datatype sendtypes[], void *recvbuf, int recvcounts[], int rdispls[], MPI_Datatype recvtypes[], MPI_Comm comm)
{
  int s, r;
  int comm_size, comm_rank;
  MPI_Status stat;
  const int tag = 0;


  MPI_Comm_size(comm, &comm_size);
  MPI_Comm_rank(comm, &comm_rank);

  for (s = 0; s < comm_size; ++s)
  for (r = 0; r < comm_size; ++r)
  MPI_Sendrecv(((char *) sendbuf) + sdispls[r], sendcounts[r], sendtypes[r], (comm_rank == s)?r:MPI_PROC_NULL, tag,
               ((char *) recvbuf) + rdispls[s], recvcounts[s], recvtypes[s], (comm_rank == r)?s:MPI_PROC_NULL, tag, comm, &stat);

  return MPI_SUCCESS;
}
