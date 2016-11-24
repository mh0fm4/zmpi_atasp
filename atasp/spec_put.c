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

#include "spec_core.h"
#include "spec_common.h"


/*#define SPEC_PRINT*/


#if MPI_VERSION >= 3
# define HAVE_MPI_RPUT  1
#else
# define HAVE_MPI_RPUT  0
#endif


static spint_t _spec_put_db(spint_t nphases, spec_elem_t *sb, spec_elem_t *rb, spec_elem_t *xb, spec_tproc_t tproc, spec_tproc_data_t tproc_data, int size, int rank, MPI_Comm comm)
{
  spint_t exit_code = SPEC_EXIT_SUCCESS;

  spint_t x, i, p;

  struct { spec_elem_index_t i; spec_proc_t p; } spec0pd;
  struct { spec_elem_index_t i; spec_int_t j, n; } spec2pd;

  spec_proc_t *procs = NULL;
  spec_elem_t *mods = NULL;

  int *scounts, *pdispls, rcount;

  spint_t rtotal;

  MPI_Aint lb, extent;
  MPI_Win *wins;

#if HAVE_MPI_RPUT
  MPI_Request *reqs = NULL;
#endif

#ifdef Z_PACK_TIMING
  double t[6] = { 0, 0, 0, 0, 0, 0 };
#endif


  Z_TIMING_SYNC(comm); Z_TIMING_START(t[0]);

  /* check supported tproc functions */
  if (!spec_check_tproc_support(tproc, 1, HAVE_MPI_RPUT, 1, HAVE_MPI_RPUT, rank, "spec_put_db"))
  {
    spec_elem_set_n(rb, 0);
    exit_code = SPEC_EXIT_FAILED;
    goto exit;
  }

  scounts = z_alloc(2 * size, sizeof(int));
  pdispls = scounts + 1 * size;

  /* setup tproc buffers */
  spec_tproc_setup(tproc, sb, &procs, &mods);

#if HAVE_MPI_RPUT
  if (mods) reqs = z_alloc(tproc->max_tprocs * spec_elem_x, sizeof(MPI_Request));
#endif

#ifdef SPEC_PRINT
  printf("input\n");
  spec_print(tproc, tproc_data, sb);
#endif

  /* reset tproc if necessary */
  if (tproc->reset) tproc->reset(tproc_data);

  /* make local send counts */
  Z_TIMING_SYNC(comm); Z_TIMING_START(t[1]);

  spec_make_counts(tproc, tproc_data, sb, 0, size, scounts, procs);

  Z_TIMING_SYNC(comm); Z_TIMING_STOP(t[1]);

  /* redistribute and reduce local send counts to put displacements */
  Z_TIMING_SYNC(comm); Z_TIMING_START(t[2]);

  spec_prefix_counts(scounts, pdispls, size,
#ifdef SPEC_PROCLISTS
    tproc->nsend_procs, tproc->send_procs, tproc->nrecv_procs, tproc->recv_procs,
#endif
    size, rank, comm);

  Z_TIMING_SYNC(comm); Z_TIMING_STOP(t[2]);

  MPI_Scatter(pdispls, 1, MPI_INT, &rcount, 1, MPI_INT, size - 1, comm);

  for (i = 0; i < size; ++i) pdispls[i] -= scounts[i];

  rtotal = rcount;

  /* check size of receive buffer */
  if (!spec_check_buffer_size(rb, rtotal, 1, comm, rank, "spec_put_db", "receive"))
  {
    spec_elem_set_n(rb, rtotal);
    exit_code = SPEC_EXIT_FAILED;
    goto free_and_exit;
  }

  /* redistribute with puts */
  wins = z_alloc(spec_elem_x, sizeof(MPI_Win));

  /* setup RMA windows */
  Z_TIMING_SYNC(comm); Z_TIMING_START(t[3]);

  for (x = 0; x < spec_elem_x; ++x)
  {
    MPI_Type_get_extent(spec_elem_x_type_mpi(rb, x), &lb, &extent);
    MPI_Win_create(spec_elem_x_get_buf(rb, x), spec_elem_x_nextent(rb, x, spec_elem_get_nmax(rb)), extent, MPI_INFO_NULL, comm, &wins[x]);
    MPI_Win_fence(MPI_MODE_NOSTORE|MPI_MODE_NOPRECEDE, wins[x]);
  }

  Z_TIMING_SYNC(comm); Z_TIMING_STOP(t[3]);

  /* reset tproc if necessary */
  if (tproc->reset) tproc->reset(tproc_data);

  Z_TIMING_SYNC(comm); Z_TIMING_START(t[4]);

  for (p = 0; p < nphases; ++p)
  {
    if (rank % nphases == p)
    {
      if (tproc->tproc)
      {
        for (spec0pd.i = 0; spec0pd.i < spec_elem_get_n(sb); ++spec0pd.i)
        {
          spec0pd.p = tproc->tproc(spec_elem_get_buf(sb), spec0pd.i, tproc_data);

          if (spec0pd.p == SPEC_PROC_NONE) continue;

          for (x = 0; x < spec_elem_x; ++x)
          {
            MPI_Put(spec_elem_x_at(sb, x, spec0pd.i), 1, spec_elem_x_type_mpi(sb, x), spec0pd.p, pdispls[spec0pd.p], 1, spec_elem_x_type_mpi(rb, x), wins[x]);
          }

          ++pdispls[spec0pd.p];
        }

      } else if (tproc->tproc_mod)
      {
#if HAVE_MPI_RPUT
        for (spec0pd.i = 0; spec0pd.i < spec_elem_get_n(sb); ++spec0pd.i)
        {
          spec0pd.p = tproc->tproc_mod(spec_elem_get_buf(sb), spec0pd.i, tproc_data, mods);

          if (spec0pd.p == SPEC_PROC_NONE) continue;

          for (x = 0; x < spec_elem_x; ++x)
          {
            MPI_Rput(spec_elem_x_at(mods, x, 0), 1, spec_elem_x_type_mpi(mods, x), spec0pd.p, pdispls[spec0pd.p], 1, spec_elem_x_type_mpi(rb, x), wins[x], &reqs[x]);
          }
          MPI_Waitall(spec_elem_x, reqs, MPI_STATUSES_IGNORE);

          ++pdispls[spec0pd.p];
        }
#endif
      } else if (tproc->tprocs)
      {
        for (spec2pd.i = 0; spec2pd.i < spec_elem_get_n(sb); ++spec2pd.i)
        {
          tproc->tprocs(spec_elem_get_buf(sb), spec2pd.i, tproc_data, &spec2pd.n, procs);

          for (spec2pd.j = 0; spec2pd.j < spec2pd.n; ++spec2pd.j)
          {
            for (x = 0; x < spec_elem_x; ++x)
            {
              MPI_Put(spec_elem_x_at(sb, x, spec2pd.i), 1, spec_elem_x_type_mpi(sb, x), procs[spec2pd.j], pdispls[procs[spec2pd.j]], 1, spec_elem_x_type_mpi(rb, x), wins[x]);
            }

            ++pdispls[procs[spec2pd.j]];
          }
        }

      } else if (tproc->tprocs_mod)
      {
#if HAVE_MPI_RPUT
        for (spec2pd.i = 0; spec2pd.i < spec_elem_get_n(sb); ++spec2pd.i)
        {
          tproc->tprocs(spec_elem_get_buf(sb), spec2pd.i, tproc_data, &spec2pd.n, procs);

          for (spec2pd.j = 0; spec2pd.j < spec2pd.n; ++spec2pd.j)
          {
            for (x = 0; x < spec_elem_x; ++x)
            {
              MPI_Rput(spec_elem_x_at(mods, x, spec2pd.j), 1, spec_elem_x_type_mpi(mods, x), procs[spec2pd.j], pdispls[procs[spec2pd.j]], 1, spec_elem_x_type_mpi(rb, x), wins[x], &reqs[spec2pd.j * spec_elem_x + x]);
            }

            ++pdispls[procs[spec2pd.j]];
          }
          MPI_Waitall(spec2pd.n * spec_elem_x, reqs, MPI_STATUSES_IGNORE);
        }
#endif
      }
    }

    if (p < nphases - 1)
    {
      for (x = 0; x < spec_elem_x; ++x) MPI_Win_fence(0, wins[x]);

      /* reset tproc if necessary */
      if (tproc->reset) tproc->reset(tproc_data);
    }
  }

  Z_TIMING_SYNC(comm); Z_TIMING_STOP(t[4]);

  Z_TIMING_SYNC(comm); Z_TIMING_START(t[5]);

  /* release RMA windows */
  for (x = 0; x < spec_elem_x; ++x)
  {
    MPI_Win_fence(MPI_MODE_NOPUT|MPI_MODE_NOSUCCEED, wins[x]);
    MPI_Win_free(&wins[x]);
  }

  Z_TIMING_SYNC(comm); Z_TIMING_STOP(t[5]);

  spec_elem_set_n(rb, rtotal);

  z_free(wins);

free_and_exit:

#if HAVE_MPI_RPUT
  if (reqs) z_free(reqs);
#endif

  /* free tproc buffers */
  spec_tproc_release(&procs, &mods);

  z_free(scounts);

exit:
  Z_TIMING_SYNC(comm); Z_TIMING_STOP(t[0]);

  Z_TIMING_PRINT(0, __func__, sizeof(t) / sizeof(t[0]), t, rank);

#if defined(Z_PACK_TIMING) && defined(SPEC_TIMING)
  if (spec_timing)
    for (i = 0; i < sizeof(t) / sizeof(t[0]); ++i) spec_timing[i] = t[i];
#endif

  return exit_code;
}


spint_t spec_put_db(spec_elem_t *sb, spec_elem_t *rb, spec_elem_t *xb, spec_tproc_t tproc, spec_tproc_data_t tproc_data, int size, int rank, MPI_Comm comm) /* sp_func spec_put_db */
{
  return _spec_put_db(1, sb, rb, xb, tproc, tproc_data, size, rank, comm);
}


spint_t spec_put_2phases_db(spec_elem_t *sb, spec_elem_t *rb, spec_elem_t *xb, spec_tproc_t tproc, spec_tproc_data_t tproc_data, int size, int rank, MPI_Comm comm) /* sp_func spec_put_2phases_db */
{
  return _spec_put_db(2, sb, rb, xb, tproc, tproc_data, size, rank, comm);
}
