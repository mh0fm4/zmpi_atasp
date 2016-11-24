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


#include <mpi.h>

#include "z_pack.h"

#include "zmpi_tools.h"


int ZMPI_Reduce_scatter_block(const void *sendbuf, void *recvbuf, int recvcount, MPI_Datatype datatype, MPI_Op op, MPI_Comm comm) /* zmpi_func ZMPI_Reduce_scatter_block */
{
#if MPI_VERSION >= 2 && MPI_SUBVERSION >= 2

  return MPI_Reduce_scatter_block((void *) sendbuf, recvbuf, recvcount, datatype,op, comm);

#else

  int comm_size, *recvcounts, i, exit_code;


  MPI_Comm_size(comm, &comm_size);

  recvcounts = z_alloc(comm_size, sizeof(int));

  for (i = 0; i < comm_size; ++i) recvcounts[i] = recvcount;

  exit_code = MPI_Reduce_scatter((void *) sendbuf, recvbuf, recvcounts, datatype, op, comm);

  z_free(recvcounts);

  return exit_code;
#endif
}
