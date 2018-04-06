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


#ifndef SPEC_ALLTOALLV_TRACE_IF
# define SPEC_ALLTOALLV_TRACE_IF  (rank == -1)
#endif


/*#define SPEC_PRINT*/


spint_t spec_alltoallv_db(spec_elem_t *sb, spec_elem_t *rb, spec_elem_t *xb, spec_tproc_t tproc, spec_tproc_data_t tproc_data, int size, int rank, MPI_Comm comm) /* sp_func spec_alltoallv_db */
{
  spint_t exit_code = SPEC_EXIT_SUCCESS;

  spint_t i;

  spec_proc_t *procs = NULL;
  spec_elem_t *mods = NULL;

  int *scounts, *sdispls, *rcounts, *rdispls;

  spec_elem_t _xb;

  spint_t stotal, rtotal;

  SPEC_DECLARE_TPROC_REARRANGE_DB
  SPEC_DECLARE_TPROC_MOD_REARRANGE_DB
  SPEC_DECLARE_TPROCS_REARRANGE_DB
  SPEC_DECLARE_TPROCS_MOD_REARRANGE_DB

#ifdef Z_PACK_TIMING
  double t[5] = { 0, 0, 0, 0, 0 };
#endif


  Z_TRACE_IF(SPEC_ALLTOALLV_TRACE_IF, "spec_alltoallv_db");

  Z_TIMING_SYNC(comm); Z_TIMING_START(t[0]);

  /* check supported tproc functions */
  if (!spec_check_tproc_support(tproc, 1, 1, 1, 1, rank, "spec_alltoallv_db"))
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
  if (!spec_check_buffer_size(rb, rtotal, 1, comm, rank, "spec_alltoallv_db", "receive"))
  {
    spec_elem_set_n(rb, rtotal);
    exit_code = SPEC_EXIT_FAILED;
    goto free_and_exit;
  }

  if (xb == NULL || spec_elem_get_nmax(xb) < stotal)
  {
    spec_elem_copy_type(sb, &_xb);
    spec_elem_alloc_tmp(&_xb, stotal);

    xb = &_xb;
  }

  spec_elem_copy_type(sb, xb);

  /* reset tproc if necessary */
  if (tproc->reset) tproc->reset(tproc_data);

  /* local rearrange */
  Z_TIMING_SYNC(comm); Z_TIMING_START(t[3]);

  if (tproc->tproc)
  {
    if (tproc->tproc_ext.rearrange_db) tproc->tproc_ext.rearrange_db(tproc_data, sb, xb, sdispls);
    else SPEC_DO_TPROC_REARRANGE_DB(tproc->tproc, tproc_data, sb, xb, sdispls);

  } else if (tproc->tproc_mod)
  {
    if (tproc->tproc_mod_ext.rearrange_db) tproc->tproc_mod_ext.rearrange_db(tproc_data, sb, xb, sdispls, mods);
    else SPEC_DO_TPROC_MOD_REARRANGE_DB(tproc->tproc_mod, tproc_data, sb, xb, sdispls, mods);

  } else if (tproc->tprocs)
  {
    if (tproc->tprocs_ext.rearrange_db) tproc->tprocs_ext.rearrange_db(tproc_data, sb, xb, sdispls, procs);
    else SPEC_DO_TPROCS_REARRANGE_DB(tproc->tprocs, tproc_data, sb, xb, sdispls, procs);

  } else if (tproc->tprocs_mod)
  {
    if (tproc->tprocs_mod_ext.rearrange_db) tproc->tprocs_mod_ext.rearrange_db(tproc_data, sb, xb, sdispls, procs, mods);
    else SPEC_DO_TPROCS_MOD_REARRANGE_DB(tproc->tprocs_mod, tproc_data, sb, xb, sdispls, procs, mods);
  }

  Z_TIMING_SYNC(comm); Z_TIMING_STOP(t[3]);

/*#ifdef SPEC_PRINT
  printf("after rearrange\n");
  spec_print(tproc, tproc_data, xb);
#endif*/

  /* redistribute with alltoallv */
  Z_TIMING_SYNC(comm); Z_TIMING_START(t[4]);

  sdispls[0] = 0;
  for (i = 1; i < size; ++i) sdispls[i] = sdispls[i - 1] + scounts[i - 1];

#ifdef SPEC_PROCLISTS
# ifdef spec_elem_alltoallv_proclists_db
  if (tproc->nsend_procs >= 0 || tproc->nrecv_procs >= 0)
  {
    spec_elem_alltoallv_proclists_db(xb, scounts, sdispls, tproc->nsend_procs, tproc->send_procs, rb, rcounts, rdispls, tproc->nrecv_procs, tproc->recv_procs, size, rank, comm);

  } else
# endif
#endif
  {
    spec_elem_alltoallv_db(xb, scounts, sdispls, rb, rcounts, rdispls, size, rank, comm);
  }

  Z_TIMING_SYNC(comm); Z_TIMING_STOP(t[4]);

  spec_elem_set_n(rb, rtotal);

  if (xb == &_xb) spec_elem_free_tmp(&_xb);

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


spint_t spec_alltoallv_ip(spec_elem_t *b, spec_elem_t *xb, spec_tproc_t tproc, spec_tproc_data_t tproc_data, int size, int rank, MPI_Comm comm) /* sp_func spec_alltoallv_ip */
{
  spint_t exit_code = SPEC_EXIT_SUCCESS;

  spint_t i;

  spec_proc_t *procs = NULL;
  spec_elem_t *mods = NULL;

  int *scounts, *sdispls, *rcounts, *rdispls;

  spec_elem_t _xb, *tb;

  spint_t stotal, rtotal;

  spint_t local_alltoallv_inplace, global_alltoallv_inplace;

  SPEC_DECLARE_TPROC_REARRANGE_IP
  SPEC_DECLARE_TPROC_REARRANGE_DB
  SPEC_DECLARE_TPROC_MOD_REARRANGE_IP
  SPEC_DECLARE_TPROC_MOD_REARRANGE_DB
  SPEC_DECLARE_TPROCS_REARRANGE_DB
  SPEC_DECLARE_TPROCS_REARRANGE_IP
  SPEC_DECLARE_TPROCS_MOD_REARRANGE_DB
  SPEC_DECLARE_TPROCS_MOD_REARRANGE_IP

#ifdef Z_PACK_TIMING
  double t[5] = { 0, 0, 0, 0, 0 };
#endif


  Z_TRACE_IF(SPEC_ALLTOALLV_TRACE_IF, "spec_alltoallv_ip");

  Z_TIMING_SYNC(comm); Z_TIMING_START(t[0]);

  /* check supported tproc functions */
  if (!spec_check_tproc_support(tproc, 1, 1, 1, 1, rank, "spec_alltoallv_ip"))
  {
    spec_elem_set_n(b, 0);
    exit_code = SPEC_EXIT_FAILED;
    goto exit;
  }

  scounts = z_alloc(4 * size, sizeof(int));
  sdispls = scounts + 1 * size;
  rcounts = scounts + 2 * size;
  rdispls = scounts + 3 * size;

  /* setup tproc buffers */
  spec_tproc_setup(tproc, b, &procs, &mods);

#ifdef SPEC_PRINT
  printf("input\n");
  spec_print(tproc, tproc_data, b);
#endif

  /* reset tproc if necessary */
  if (tproc->reset) tproc->reset(tproc_data);

  /* make local counts */
  Z_TIMING_SYNC(comm); Z_TIMING_START(t[1]);

  spec_make_counts(tproc, tproc_data, b, 1, size, scounts, procs);

  Z_TIMING_SYNC(comm); Z_TIMING_STOP(t[1]);

#ifdef SPEC_PRINT
  printf("after count\n");
  spec_print(tproc, tproc_data, b);
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

  /* check size of buffer */
  if (!spec_check_buffer_size(b, z_max(stotal, rtotal), 0, comm, rank, "spec_alltoallv_ip", ""))
  {
    spec_elem_set_n(b, z_max(stotal, rtotal));
    exit_code = SPEC_EXIT_FAILED;
    goto free_and_exit;
  }

  if (xb == NULL)
  {
    spec_elem_copy_type(b, &_xb);
    spec_elem_alloc_tmp(&_xb, 1);

    xb = &_xb;
  }

  tb = NULL;
  if (spec_elem_get_nmax(xb) >= stotal) tb = xb;

  spec_elem_copy_type(b, xb);

  /* reset tproc if necessary */
  if (tproc->reset) tproc->reset(tproc_data);

  /* local rearrange */
  Z_TIMING_SYNC(comm); Z_TIMING_START(t[3]);

  if (tproc->tproc)
  {
    if (tb)
    {
      if (tproc->tproc_ext.rearrange_db) tproc->tproc_ext.rearrange_db(tproc_data, b, tb, sdispls);
      else SPEC_DO_TPROC_REARRANGE_DB(tproc->tproc, tproc_data, b, tb, sdispls);

    } else
    {
      if (tproc->tproc_ext.rearrange_ip) tproc->tproc_ext.rearrange_ip(tproc_data, b, xb, sdispls, scounts, size);
      else SPEC_DO_TPROC_REARRANGE_IP(tproc->tproc, tproc_data, b, xb, sdispls, scounts, size);
    }

  } else if (tproc->tproc_mod)
  {
    if (tb)
    {
      if (tproc->tproc_mod_ext.rearrange_db) tproc->tproc_mod_ext.rearrange_db(tproc_data, b, tb, sdispls, mods);
      else SPEC_DO_TPROC_MOD_REARRANGE_DB(tproc->tproc_mod, tproc_data, b, tb, sdispls, mods);

    } else
    {
      if (tproc->tproc_mod_ext.rearrange_ip) tproc->tproc_mod_ext.rearrange_ip(tproc_data, b, xb, sdispls, scounts, size, mods);
      else SPEC_DO_TPROC_MOD_REARRANGE_IP(tproc->tproc_mod, tproc_data, b, xb, sdispls, scounts, size, mods);
    }

  } else if (tproc->tprocs)
  {
    if (tb)
    {
      if (tproc->tprocs_ext.rearrange_db) tproc->tprocs_ext.rearrange_db(tproc_data, b, tb, sdispls, procs);
      else SPEC_DO_TPROCS_REARRANGE_DB(tproc->tprocs, tproc_data, b, tb, sdispls, procs);

    } else
    {
      if (tproc->tprocs_ext.rearrange_ip) tproc->tprocs_ext.rearrange_ip(tproc_data, b, xb, sdispls, scounts, size, procs);
      else SPEC_DO_TPROCS_REARRANGE_IP(tproc->tprocs, tproc_data, b, xb, sdispls, scounts, size, procs);
    }

  } else if (tproc->tprocs_mod)
  {
    if (tb)
    {
      if (tproc->tprocs_mod_ext.rearrange_db) tproc->tprocs_mod_ext.rearrange_db(tproc_data, b, tb, sdispls, procs, mods);
      else SPEC_DO_TPROCS_MOD_REARRANGE_DB(tproc->tprocs_mod, tproc_data, b, tb, sdispls, procs, mods);

    } else
    {
      if (tproc->tprocs_mod_ext.rearrange_ip) tproc->tprocs_mod_ext.rearrange_ip(tproc_data, b, xb, sdispls, scounts, size, procs, mods);
      else SPEC_DO_TPROCS_MOD_REARRANGE_IP(tproc->tprocs_mod, tproc_data, b, xb, sdispls, scounts, size, procs, mods);
    }
  }

  Z_TIMING_SYNC(comm); Z_TIMING_STOP(t[3]);

#ifdef SPEC_PRINT
  printf("after rearrange\n");
  spec_print(tproc, tproc_data, b);
#endif

  /* redistribute with alltoallv */
  Z_TIMING_SYNC(comm); Z_TIMING_START(t[4]);

  sdispls[0] = 0;
  for (i = 1; i < size; ++i) sdispls[i] = sdispls[i - 1] + scounts[i - 1];

  if (tb) local_alltoallv_inplace = 0;
  else
  {
    if (spec_elem_get_nmax(xb) >= rtotal) local_alltoallv_inplace = 0;
    else local_alltoallv_inplace = 1;
  }

  MPI_Allreduce(&local_alltoallv_inplace, &global_alltoallv_inplace, 1, MPI_SPINT, MPI_SUM, comm);

  /* no process wants inplace? */
  if (global_alltoallv_inplace == 0)
  {
    /* double-buffered alltoallv  */
    if (tb)
    {
#ifdef SPEC_PROCLISTS
# ifdef spec_elem_alltoallv_proclists_db
      if (tproc->nsend_procs >= 0 || tproc->nrecv_procs >= 0)
      {
        spec_elem_alltoallv_proclists_db(tb, scounts, sdispls, tproc->nsend_procs, tproc->send_procs, b, rcounts, rdispls, tproc->nrecv_procs, tproc->recv_procs, size, rank, comm);

      } else
# endif
#endif
      {
        spec_elem_alltoallv_db(tb, scounts, sdispls, b, rcounts, rdispls, size, rank, comm);
      }
    } else
    {
#ifdef SPEC_PROCLISTS
# ifdef spec_elem_alltoallv_proclists_db
      if (tproc->nsend_procs >= 0 || tproc->nrecv_procs >= 0)
      {
        spec_elem_alltoallv_proclists_db(b, scounts, sdispls, tproc->nsend_procs, tproc->send_procs, xb, rcounts, rdispls, tproc->nrecv_procs, tproc->recv_procs, size, rank, comm);

      } else
# endif
#endif
      {
        spec_elem_alltoallv_db(b, scounts, sdispls, xb, rcounts, rdispls, size, rank, comm);
      }
      spec_elem_ncopy_at(xb, 0, b, 0, rtotal);
    }

  } else
  {
    if (tb) spec_elem_ncopy_at(tb, 0, b, 0, stotal);
#ifdef SPEC_PROCLISTS
# ifdef spec_elem_alltoallv_proclists_ip
    if (tproc->nsend_procs >= 0 || tproc->nrecv_procs >= 0)
    {
      spec_elem_alltoallv_proclists_ip(b, xb, scounts, sdispls, tproc->nsend_procs, tproc->send_procs, rcounts, rdispls, tproc->nrecv_procs, tproc->recv_procs, size, rank, comm);

    } else
# endif
#endif
    {
#ifdef spec_elem_alltoallv_ip
      spec_elem_alltoallv_ip(b, xb, scounts, sdispls, rcounts, rdispls, size, rank, comm);
#else
#ifdef SPEC_ERROR_FILE
      fprintf(SPEC_ERROR_FILE, "%d: spec_alltoallv_ip: error: required in-place support for alltoallv not available", rank);
#endif
      exit_code = SPEC_EXIT_FAILED;
      goto free_and_exit;
#endif
    }
  }

  Z_TIMING_SYNC(comm); Z_TIMING_STOP(t[4]);

  spec_elem_set_n(b, rtotal);

  if (xb == &_xb) spec_elem_free_tmp(&_xb);

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
