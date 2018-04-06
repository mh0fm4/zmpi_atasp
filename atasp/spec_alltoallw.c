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

#include "spec_core.h"
#include "spec_common.h"


/*#define SPEC_PRINT*/


spint_t spec_alltoallw_db(spec_elem_t *sb, spec_elem_t *rb, spec_elem_t *xb, spec_tproc_t tproc, spec_tproc_data_t tproc_data, int size, int rank, MPI_Comm comm) /* sp_func spec_alltoallw_db */
{
  spint_t exit_code = SPEC_EXIT_SUCCESS;

  spint_t i;

  spec_proc_t *procs = NULL;
  spec_elem_t *mods = NULL;

  int *scounts, *sdispls, *rcounts, *rdispls;

  spint_t stotal, rtotal;

  int *icounts, *idispls, *indices;
  MPI_Datatype *stypes, *rtypes;

  spec_elem_t _xb;

  SPEC_DECLARE_TPROC_INDICES_DB
  SPEC_DECLARE_TPROC_MOD_INDICES_DB
  SPEC_DECLARE_TPROCS_INDICES_DB
  SPEC_DECLARE_TPROCS_MOD_INDICES_DB

  DECLARE_SPEC_ELEM_ALLTOALLW_INDEXED_DB

#ifdef Z_PACK_TIMING
  double t[5] = { 0, 0, 0, 0, 0 };
#endif


  Z_TIMING_SYNC(comm); Z_TIMING_START(t[0]);

  /* check supported tproc functions */
  if (!spec_check_tproc_support(tproc, 1, 1, 1, 1, rank, "spec_alltoallw_db"))
  {
    spec_elem_set_n(rb, 0);
    exit_code = SPEC_EXIT_FAILED;
    goto exit;
  }

  scounts = z_alloc(4 * size, sizeof(int));
  sdispls = scounts + 1 * size;
  rcounts = scounts + 2 * size;
  rdispls = scounts + 3 * size;

  /* setup tproc buffers */
  spec_tproc_setup(tproc, sb, &procs, &mods);

#ifdef SPEC_PRINT
  printf("input\n");
  spec_print(tproc, tproc_data, sb);
#endif

  /* reset tproc if necessary */
  if (tproc->reset) tproc->reset(tproc_data);

  /* make local counts */
  Z_TIMING_SYNC(comm); Z_TIMING_START(t[1]);

  spec_make_counts(tproc, tproc_data, sb, 0, size, scounts, procs);

  Z_TIMING_SYNC(comm); Z_TIMING_STOP(t[1]);

#ifdef SPEC_PRINT
  printf("after count\n");
  spec_print(tproc, tproc_data, sb);
#endif

  /* redistribute local counts */
  Z_TIMING_SYNC(comm); Z_TIMING_START(t[2]);

  spec_redistribute_counts(tproc, scounts, rcounts,
#ifdef SPEC_PROCLISTS
    tproc->nsend_procs, tproc->send_procs, tproc->nrecv_procs, tproc->recv_procs,
#endif
    size, rank, comm);

  Z_TIMING_SYNC(comm); Z_TIMING_STOP(t[2]);

  sdispls[0] = rdispls[0] = 0;
  for (i = 1; i < size; ++i)
  {
    sdispls[i] = sdispls[i - 1] + scounts[i - 1];
    rdispls[i] = rdispls[i - 1] + rcounts[i - 1];
  }

  stotal = sdispls[size - 1] + scounts[size - 1];
  rtotal = rdispls[size - 1] + rcounts[size - 1];

  /* check size of receive buffer */
  if (!spec_check_buffer_size(rb, rtotal, 1, comm, rank, "spec_alltoallw_db", "receive"))
  {
    spec_elem_set_n(rb, rtotal);
    exit_code = SPEC_EXIT_FAILED;
    goto free_and_exit;
  }

  /* reset tproc if necessary */
  if (tproc->reset) tproc->reset(tproc_data);

  /* index creation */
  Z_TIMING_SYNC(comm); Z_TIMING_START(t[3]);

  icounts = z_alloc(2 * size + stotal, sizeof(int));
  idispls = icounts + 1 * size;
  indices = icounts + 2 * size;

  for (i = 0; i < size; ++i)
  {
    icounts[i] = scounts[i];
    idispls[i] = sdispls[i];
  }

  if (tproc->tproc_mod || tproc->tprocs_mod)
  {
    if (xb == NULL || spec_elem_get_nmax(xb) < stotal)
    {
      spec_elem_copy_type(sb, &_xb);
      spec_elem_alloc_tmp(&_xb, stotal);

      xb = &_xb;
    }

    spec_elem_copy_type(sb, xb);
  }

  if (tproc->tproc)
  {
    if (tproc->tproc_ext.indices_db) tproc->tproc_ext.indices_db(tproc_data, sb, indices, sdispls);
    else SPEC_DO_TPROC_INDICES_DB(tproc->tproc, tproc_data, sb, indices, sdispls);

  } else if (tproc->tproc_mod)
  {
    if (tproc->tproc_mod_ext.indices_db) tproc->tproc_mod_ext.indices_db(tproc_data, sb, indices, sdispls, mods, xb);
    else SPEC_DO_TPROC_MOD_INDICES_DB(tproc->tproc_mod, tproc_data, sb, indices, sdispls, mods, xb);

    sb = xb;

  } else if (tproc->tprocs)
  {
    if (tproc->tprocs_ext.indices_db) tproc->tprocs_ext.indices_db(tproc_data, sb, indices, sdispls, procs);
    else SPEC_DO_TPROCS_INDICES_DB(tproc->tprocs, tproc_data, sb, indices, sdispls, procs);

  } else if (tproc->tprocs_mod)
  {
    if (tproc->tprocs_mod_ext.indices_db) tproc->tprocs_mod_ext.indices_db(tproc_data, sb, indices, sdispls, procs, mods, xb);
    else SPEC_DO_TPROCS_MOD_INDICES_DB(tproc->tprocs_mod, tproc_data, sb, indices, sdispls, procs, mods, xb);

    sb = xb;
  }

  stypes = z_alloc(2 * size, sizeof(MPI_Datatype));
  rtypes = stypes + size;

  for (i = 0; i < size; ++i)
  {
    scounts[i] = (scounts[i] > 0)?1:0;
    sdispls[i] = 0;
  }

  Z_TIMING_SYNC(comm); Z_TIMING_STOP(t[3]);

  /* redistribute with alltoallw */
  Z_TIMING_SYNC(comm); Z_TIMING_START(t[4]);

#ifdef SPEC_PROCLISTS
# ifdef spec_elem_alltoallw_indexed_proclists_db
  if (tproc->nsend_procs >= 0 || tproc->nrecv_procs >= 0)
  {
    spec_elem_alltoallw_indexed_proclists_db(sb, scounts, sdispls, stypes, tproc->nsend_procs, tproc->send_procs, rb, rcounts, rdispls, rtypes, tproc->nrecv_procs, tproc->recv_procs, indices, icounts, idispls, size, rank, comm);
  }
  else
# endif
#endif
  {
    spec_elem_alltoallw_indexed_db(sb, scounts, sdispls, stypes, rb, rcounts, rdispls, rtypes, indices, icounts, idispls, size, rank, comm);
  }

  Z_TIMING_SYNC(comm); Z_TIMING_STOP(t[4]);

  spec_elem_set_n(rb, rtotal);

  z_free(stypes);

  if (tproc->tproc_mod || tproc->tprocs_mod)
  {
    if (xb == &_xb) spec_elem_free_tmp(&_xb);
  }

  z_free(icounts);

free_and_exit:

  /* free tproc buffers */
  spec_tproc_release(&procs, &mods);

  z_free(scounts);

exit:
  Z_TIMING_SYNC(comm); Z_TIMING_STOP(t[0]);

  Z_TIMING_PRINT(0, __func__, sizeof(t) / sizeof(t[0]), t, rank);

#if defined(Z_PACK_TIMING) && defined(SPEC_TIMING)
  nspec_timings = sizeof(t) / sizeof(t[0]);
  if (spec_timing)
    for (i = 0; i < nspec_timings; ++i) spec_timing[i] = t[i];
#endif

  return exit_code;
}
