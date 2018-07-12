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
#include <string.h>

#include <mpi.h>

#ifdef HAVE_ZMPI_TOOLS_H
# include "zmpi_tools.h"
#endif

#include "spec_core.h"
#include "spec_common.h"


#ifndef SPEC_COMMON_TRACE_IF
# define SPEC_COMMON_TRACE_IF  (rank == -1)
#endif


#ifdef SPEC_ERROR_FILE
static const char *spec_tproc_name(spint_t id)
{
  switch (id)
  {
    case 0: return "none";
    case 1: return "tproc";
    case 2: return "tproc_mod";
    case 3: return "tprocs";
    case 4: return "tprocs_mod";
  }

  return "unsupported";
}
#endif


static spint_t spec_tproc_supported(spec_tproc_t tproc, spint_t have_tproc, spint_t have_tproc_mod, spint_t have_tprocs, spint_t have_tprocs_mod, spint_t *id)
{
  if (tproc != SPEC_TPROC_NULL)
  {
    if (tproc->tproc)
    {
      *id = 1;
      return have_tproc;

    } else if (tproc->tproc_mod)
    {
      *id = 2;
      return have_tproc_mod;

    } else if (tproc->tprocs)
    {
      *id = 3;
      return have_tprocs;

    } else if (tproc->tprocs_mod)
    {
      *id = 4;
      return have_tprocs_mod;
    }
  }

  *id = 0;
  return 0;
}


spint_t spec_check_tproc_support(spec_tproc_t tproc, spint_t have_tproc, spint_t have_tproc_mod, spint_t have_tprocs, spint_t have_tprocs_mod, int rank, const char *name) /* sp_func spec_check_tproc_support */
{
  spint_t id;

  if (!spec_tproc_supported(tproc, have_tproc, have_tproc_mod, have_tprocs, have_tprocs_mod, &id))
  {
#ifdef SPEC_ERROR_FILE
    fprintf(SPEC_ERROR_FILE, "%d: spec_alltoallv_db: error: target process function type '%s' is not supported\n", rank, spec_tproc_name(id));
#endif
    return 0;
  }

  return 1;
}


spint_t spec_check_buffer_size(spec_elem_t *b, spint_t min_size, spint_t allocatable, MPI_Comm comm, int rank, const char *name, const char *buf_name) /* sp_func spec_check_buffer_size */
{
  spint_t local_buffer_exit, global_buffer_exit;


  local_buffer_exit = (min_size > spec_elem_get_nmax(b))?1:0;

#ifdef spec_elem_alloc_rbuf
  if (local_buffer_exit > 0 && spec_elem_alloc_rbuf(b) && allocatable)
  {
    spec_elem_free_buf(b);
    spec_elem_alloc_buf(b, min_size);
    local_buffer_exit = 0;
  }
#endif

#ifdef SPEC_GLOBAL_EXIT_ON_ERROR
  MPI_Allreduce(&local_buffer_exit, &global_buffer_exit, 1, MPI_SPINT, MPI_SUM, comm);
#else
  global_buffer_exit = local_buffer_exit;
#endif

  if (global_buffer_exit > 0)
  {
#ifdef SPEC_ERROR_FILE
    fprintf(SPEC_ERROR_FILE,
# ifdef SPEC_GLOBAL_EXIT_ON_ERROR
      "%d: %s: error: one %s%sbuffer is too small (local: %" spint_fmt " vs. %" spint_fmt")\n",
# else
      "%d: %s: error: local %s%sbuffer too small (%" spint_fmt " vs. %" spint_fmt")\n",
# endif
      rank, name, buf_name, ((strlen(buf_name) > 0)?" ":""), (spint_t) spec_elem_get_nmax(b), min_size);
#endif
    return 0;
  }

  return 1;
}


void spec_tproc_setup(spec_tproc_t tproc, spec_elem_t *b, spec_proc_t **procs, spec_elem_t **mods) /* sp_func spec_tproc_setup */
{
  spint_t tprocs_max;
  spint_t alloc_mods;


  tprocs_max = 0;
  alloc_mods = 0;

  if (tproc->tproc)
  {
    /* nothing */

  } else if (tproc->tproc_mod)
  {
    alloc_mods = 1;

  } else if (tproc->tprocs)
  {
    tprocs_max = tproc->max_tprocs;

  } else if (tproc->tprocs_mod)
  {
    tprocs_max = tproc->max_tprocs;
    alloc_mods = 1;
  }

  if (procs)
  {
    if (tprocs_max) *procs = z_alloc(tprocs_max, sizeof(spec_proc_t));
    else *procs = NULL;
  }

  if (mods)
  {
    if (alloc_mods)
    {
      *mods = z_alloc(1, sizeof(spec_elem_t));

      spec_elem_copy_type(b, *mods);
      spec_elem_alloc_buf(*mods, 2 * z_max(1, tprocs_max));

    } else *mods = NULL;
  }
}


void spec_tproc_release(spec_proc_t **procs, spec_elem_t **mods) /* sp_func spec_tproc_release */
{
  if (procs && *procs) z_free(*procs);

  if (mods && *mods)
  {
    spec_elem_free_buf(*mods);
    z_free(*mods);
  }
}


spint_t spec_make_counts(spec_tproc_t tproc, spec_tproc_data_t tproc_data, spec_elem_t *b, int ip, int *counts, spec_proc_t *procs, int size, int rank) /* sp_func spec_make_counts */
{
  spint_t i;

  SPEC_DECLARE_TPROC_COUNT_DB
  SPEC_DECLARE_TPROC_MOD_COUNT_DB
  SPEC_DECLARE_TPROCS_COUNT_DB
  SPEC_DECLARE_TPROCS_MOD_COUNT_DB

  SPEC_DECLARE_TPROC_COUNT_IP
  SPEC_DECLARE_TPROC_MOD_COUNT_IP
  SPEC_DECLARE_TPROCS_COUNT_IP
  SPEC_DECLARE_TPROCS_MOD_COUNT_IP

  Z_TRACE_IF(SPEC_COMMON_TRACE_IF, "spec_make_counts");

#ifdef SPEC_COUNTS
  if (tproc->send_counts)
  {
    memcpy(counts, tproc->send_counts, size * sizeof(int));
    return 0;
  }
#endif /* SPEC_COUNTS */

  for (i = 0; i < size; ++i) counts[i] = 0;

  if (!ip)
  {
    if (tproc->tproc)
    {
      if (tproc->tproc_ext.count_db) tproc->tproc_ext.count_db(tproc_data, b, counts);
      else SPEC_DO_TPROC_COUNT_DB(tproc->tproc, tproc_data, b, counts);

    } else if (tproc->tproc_mod)
    {
      if (tproc->tproc_mod_ext.count_db) tproc->tproc_mod_ext.count_db(tproc_data, b, counts);
      else SPEC_DO_TPROC_MOD_COUNT_DB(tproc->tproc_mod, tproc_data, b, counts);

    } else if (tproc->tprocs)
    {
      if (tproc->tprocs_ext.count_db) tproc->tprocs_ext.count_db(tproc_data, b, counts, procs);
      else SPEC_DO_TPROCS_COUNT_DB(tproc->tprocs, tproc_data, b, counts, procs);

    } else if (tproc->tprocs_mod)
    {
      if (tproc->tprocs_mod_ext.count_db) tproc->tprocs_mod_ext.count_db(tproc_data, b, counts, procs);
      else SPEC_DO_TPROCS_MOD_COUNT_DB(tproc->tprocs_mod, tproc_data, b, counts, procs);
    }

  } else {

    if (tproc->tproc)
    {
      if (tproc->tproc_ext.count_ip) tproc->tproc_ext.count_ip(tproc_data, b, counts);
      else SPEC_DO_TPROC_COUNT_IP(tproc->tproc, tproc_data, b, counts);

    } else if (tproc->tproc_mod)
    {
      if (tproc->tproc_mod_ext.count_ip) tproc->tproc_mod_ext.count_ip(tproc_data, b, counts);
      else SPEC_DO_TPROC_MOD_COUNT_IP(tproc->tproc_mod, tproc_data, b, counts);

    } else if (tproc->tprocs)
    {
      if (tproc->tprocs_ext.count_ip) tproc->tprocs_ext.count_ip(tproc_data, b, counts, procs);
      else SPEC_DO_TPROCS_COUNT_IP(tproc->tprocs, tproc_data, b, counts, procs);

    } else if (tproc->tprocs_mod)
    {
      if (tproc->tprocs_mod_ext.count_ip) tproc->tprocs_mod_ext.count_ip(tproc_data, b, counts, procs);
      else SPEC_DO_TPROCS_MOD_COUNT_IP(tproc->tprocs_mod, tproc_data, b, counts, procs);
    }
  }

  return 0;
}


/* sp_var spec_redistribute_counts_type spec_redistribute_counts_proclists_type */
spint_t spec_redistribute_counts_type = SPEC_REDISTRIBUTE_COUNTS_DEFAULT;
spint_t spec_redistribute_counts_proclists_type = SPEC_REDISTRIBUTE_COUNTS_PROCLISTS_DEFAULT;


spint_t spec_redistribute_counts(spec_tproc_t tproc, int *scounts, int *rcounts,
#ifdef SPEC_PROCLISTS
  spint_t nsend_procs, sproc_t *send_procs, spint_t nrecv_procs, sproc_t *recv_procs,
#endif
  int size, int rank, MPI_Comm comm) /* sp_func spec_redistribute_counts */
{
  Z_TRACE_IF(SPEC_COMMON_TRACE_IF, "spec_redistribute_counts");

#ifdef SPEC_COUNTS
  if (tproc->recv_counts)
  {
    memcpy(rcounts, tproc->recv_counts, size * sizeof(int));
    return 0;
  }
#endif /* SPEC_COUNTS */

#ifdef HAVE_ZMPI_ALLTOALL_INT

#ifdef SPEC_PROCLISTS
  if ((nsend_procs >= 0 || nrecv_procs >= 0) && spec_redistribute_counts_proclists_type != SPEC_REDISTRIBUTE_COUNTS_PROCLISTS_IGNORE)
  {
    memset(rcounts, 0, size * sizeof(int));

    switch (spec_redistribute_counts_proclists_type)
    {
      case SPEC_REDISTRIBUTE_COUNTS_PROCLISTS_DEFAULT:
        Z_TRACE_IF(SPEC_COMMON_TRACE_IF, "SPEC_REDISTRIBUTE_COUNTS_PROCLISTS_DEFAULT -> isendirecv");
      case SPEC_REDISTRIBUTE_COUNTS_PROCLISTS_ISENDIRECV:
        Z_TRACE_IF(SPEC_COMMON_TRACE_IF, "SPEC_REDISTRIBUTE_COUNTS_PROCLISTS_ISENDIRECV");
        ZMPI_Alltoall_int_proclists_isendirecv(scounts, nsend_procs, send_procs, rcounts, nrecv_procs, recv_procs, comm);
        break;
      case SPEC_REDISTRIBUTE_COUNTS_PROCLISTS_ALLTOALLV:
        Z_TRACE_IF(SPEC_COMMON_TRACE_IF, "SPEC_REDISTRIBUTE_COUNTS_PROCLISTS_ALLTOALLV");
        ZMPI_Alltoall_int_proclists_alltoallv(scounts, nsend_procs, send_procs, rcounts, nrecv_procs, recv_procs, comm);
        break;
      case SPEC_REDISTRIBUTE_COUNTS_PROCLISTS_PUT:
        Z_TRACE_IF(SPEC_COMMON_TRACE_IF, "SPEC_REDISTRIBUTE_COUNTS_PROCLISTS_PUT");
        ZMPI_Alltoall_int_proclists_put(scounts, nsend_procs, send_procs, rcounts, nrecv_procs, recv_procs, comm);
        break;
      case SPEC_REDISTRIBUTE_COUNTS_PROCLISTS_PUT_ALLOC:
        Z_TRACE_IF(SPEC_COMMON_TRACE_IF, "SPEC_REDISTRIBUTE_COUNTS_PROCLISTS_PUT_ALLOC");
        ZMPI_Alltoall_int_proclists_put_alloc(scounts, nsend_procs, send_procs, rcounts, nrecv_procs, recv_procs, comm);
        break;
      case SPEC_REDISTRIBUTE_COUNTS_PROCLISTS_PUT_2PHASES:
        Z_TRACE_IF(SPEC_COMMON_TRACE_IF, "SPEC_REDISTRIBUTE_COUNTS_PROCLISTS_PUT_2PHASES");
        ZMPI_Alltoall_int_proclists_put_2phases(scounts, nsend_procs, send_procs, rcounts, nrecv_procs, recv_procs, comm);
        break;
      case SPEC_REDISTRIBUTE_COUNTS_PROCLISTS_PUT_2PHASES_ALLOC:
        Z_TRACE_IF(SPEC_COMMON_TRACE_IF, "SPEC_REDISTRIBUTE_COUNTS_PROCLISTS_PUT_2PHASES_ALLOC");
        ZMPI_Alltoall_int_proclists_put_2phases_alloc(scounts, nsend_procs, send_procs, rcounts, nrecv_procs, recv_procs, comm);
        break;
    }

  } else
#endif
  {
    switch (spec_redistribute_counts_type)
    {
      case SPEC_REDISTRIBUTE_COUNTS_DEFAULT:
#ifdef SPEC_REDISTRIBUTE_COUNTS_2STEP_THRESHOLD
        if (size >= SPEC_REDISTRIBUTE_COUNTS_2STEP_THRESHOLD)
        {
          Z_TRACE_IF(SPEC_COMMON_TRACE_IF, "SPEC_REDISTRIBUTE_COUNTS_DEFAULT -> 2step");
          ZMPI_Alltoall_int_2step(scounts, rcounts, comm);

        } else
#endif
        {
          Z_TRACE_IF(SPEC_COMMON_TRACE_IF, "SPEC_REDISTRIBUTE_COUNTS_DEFAULT -> alltoall");
          ZMPI_Alltoall_int_alltoall(scounts, rcounts, comm);
        }
        break;
      case SPEC_REDISTRIBUTE_COUNTS_ALLTOALL:
        Z_TRACE_IF(SPEC_COMMON_TRACE_IF, "SPEC_REDISTRIBUTE_COUNTS_ALLTOALL");
        ZMPI_Alltoall_int_alltoall(scounts, rcounts, comm);
        break;
      case SPEC_REDISTRIBUTE_COUNTS_2STEP:
        Z_TRACE_IF(SPEC_COMMON_TRACE_IF, "SPEC_REDISTRIBUTE_COUNTS_2STEP");
        ZMPI_Alltoall_int_2step(scounts, rcounts, comm);
        break;
      case SPEC_REDISTRIBUTE_COUNTS_PUT:
        Z_TRACE_IF(SPEC_COMMON_TRACE_IF, "SPEC_REDISTRIBUTE_COUNTS_PUT");
        ZMPI_Alltoall_int_put(scounts, rcounts, comm);
        break;
      case SPEC_REDISTRIBUTE_COUNTS_PUT_ALLOC:
        Z_TRACE_IF(SPEC_COMMON_TRACE_IF, "SPEC_REDISTRIBUTE_COUNTS_PUT_ALLOC");
        ZMPI_Alltoall_int_put_alloc(scounts, rcounts, comm);
        break;
      case SPEC_REDISTRIBUTE_COUNTS_PUT_2PHASES:
        Z_TRACE_IF(SPEC_COMMON_TRACE_IF, "SPEC_REDISTRIBUTE_COUNTS_PUT_2PHASES");
        ZMPI_Alltoall_int_put_2phases(scounts, rcounts, comm);
        break;
      case SPEC_REDISTRIBUTE_COUNTS_PUT_2PHASES_ALLOC:
        Z_TRACE_IF(SPEC_COMMON_TRACE_IF, "SPEC_REDISTRIBUTE_COUNTS_PUT_2PHASES_ALLOC");
        ZMPI_Alltoall_int_put_2phases_alloc(scounts, rcounts, comm);
        break;
      case SPEC_REDISTRIBUTE_COUNTS_PUT_3PHASES:
        Z_TRACE_IF(SPEC_COMMON_TRACE_IF, "SPEC_REDISTRIBUTE_COUNTS_PUT_3PHASES");
        ZMPI_Alltoall_int_put_3phases(scounts, rcounts, comm);
        break;
      case SPEC_REDISTRIBUTE_COUNTS_PUT_3PHASES_ALLOC:
        Z_TRACE_IF(SPEC_COMMON_TRACE_IF, "SPEC_REDISTRIBUTE_COUNTS_PUT_3PHASES_ALLOC");
        ZMPI_Alltoall_int_put_3phases_alloc(scounts, rcounts, comm);
        break;
    }
  }

#else /* HAVE_ZMPI_ALLTOALL_INT */

  Z_TRACE_IF(SPEC_COMMON_TRACE_IF, "MPI_Alltoall");
  MPI_Alltoall(scounts, 1, MPI_INT, rcounts, 1, MPI_INT, comm);

#endif /* HAVE_ZMPI_ALLTOALL_INT */

  return 0;
}


/* sp_var spec_reduce_scatter_counts_type spec_reduce_scatter_counts_proclists_type */
spint_t spec_reduce_scatter_counts_type = SPEC_REDUCE_SCATTER_COUNTS_DEFAULT;
spint_t spec_reduce_scatter_counts_proclists_type = SPEC_REDUCE_SCATTER_COUNTS_PROCLISTS_DEFAULT;


spint_t spec_reduce_scatter_counts(int *scounts, int *rcounts, int ncounts,
#ifdef SPEC_PROCLISTS
  spint_t nsend_procs, sproc_t *send_procs, spint_t nrecv_procs, sproc_t *recv_procs,
#endif
  int size, int rank, MPI_Comm comm) /* sp_func spec_reduce_scatter_counts */
{
  Z_TRACE_IF(SPEC_COMMON_TRACE_IF, "spec_reduce_scatter_counts");

#ifdef HAVE_ZMPI_REDUCE_SCATTER_BLOCK_INTSUM

#ifdef SPEC_PROCLISTS
  if ((nsend_procs >= 0 || nrecv_procs >= 0) && spec_reduce_scatter_counts_proclists_type != SPEC_REDUCE_SCATTER_COUNTS_PROCLISTS_IGNORE)
  {
    memset(rcounts, 0, ncounts * sizeof(int));

    switch (spec_reduce_scatter_counts_proclists_type)
    {
      case SPEC_REDUCE_SCATTER_COUNTS_PROCLISTS_DEFAULT:
        Z_TRACE_IF(SPEC_COMMON_TRACE_IF, "SPEC_REDUCE_SCATTER_COUNTS_PROCLISTS_DEFAULT -> isendirecv");
      case SPEC_REDUCE_SCATTER_COUNTS_PROCLISTS_ISENDIRECV:
        Z_TRACE_IF(SPEC_COMMON_TRACE_IF, "SPEC_REDUCE_SCATTER_COUNTS_PROCLISTS_ISENDIRECV");
        ZMPI_Reduce_scatter_block_intsum_proclists_isendirecv(scounts, nsend_procs, send_procs, rcounts, ncounts, nrecv_procs, recv_procs, comm);
        break;
      case SPEC_REDUCE_SCATTER_COUNTS_PROCLISTS_ALLTOALLV:
        Z_TRACE_IF(SPEC_COMMON_TRACE_IF, "SPEC_REDUCE_SCATTER_COUNTS_PROCLISTS_ALLTOALLV");
        ZMPI_Reduce_scatter_block_intsum_proclists_alltoallv(scounts, nsend_procs, send_procs, rcounts, ncounts, nrecv_procs, recv_procs, comm);
        break;
      case SPEC_REDUCE_SCATTER_COUNTS_PROCLISTS_ACCUMULATE:
        Z_TRACE_IF(SPEC_COMMON_TRACE_IF, "SPEC_REDUCE_SCATTER_COUNTS_PROCLISTS_ACCUMULATE");
        ZMPI_Reduce_scatter_block_intsum_proclists_accumulate(scounts, nsend_procs, send_procs, rcounts, ncounts, nrecv_procs, recv_procs, comm);
        break;
    }

  } else
#endif
  {
    switch (spec_reduce_scatter_counts_type)
    {
      case SPEC_REDUCE_SCATTER_COUNTS_DEFAULT:
        Z_TRACE_IF(SPEC_COMMON_TRACE_IF, "SPEC_REDUCE_SCATTER_COUNTS_DEFAULT -> redscat");
      case SPEC_REDUCE_SCATTER_COUNTS_REDSCAT:
        Z_TRACE_IF(SPEC_COMMON_TRACE_IF, "SPEC_REDUCE_SCATTER_COUNTS_REDSCAT");
        ZMPI_Reduce_scatter_block(scounts, rcounts, ncounts, MPI_INT, MPI_SUM, comm);
        break;
      case SPEC_REDUCE_SCATTER_COUNTS_ACCUMULATE:
        Z_TRACE_IF(SPEC_COMMON_TRACE_IF, "SPEC_REDUCE_SCATTER_COUNTS_ACCUMULATE");
        ZMPI_Reduce_scatter_block_intsum_accumulate(scounts, rcounts, ncounts, comm);
        break;
    }
  }

#elif defined(HAVE_ZMPI_REDUCE_SCATTER_BLOCK)

  Z_TRACE_IF(SPEC_COMMON_TRACE_IF, "ZMPI_Reduce_scatter_block");
  ZMPI_Reduce_scatter_block(scounts, rcounts, ncounts, MPI_INT, MPI_SUM, comm);

#else

#ifdef SPEC_ERROR_FILE
  fprintf(SPEC_ERROR_FILE, "%d: spec_reduce_scatter_counts: error: no implementation available\n", rank);
#endif

#endif

  return 0;
}


/* sp_var spec_prefix_counts_type spec_prefix_counts_proclists_type */
spint_t spec_prefix_counts_type = SPEC_PREFIX_COUNTS_DEFAULT;
spint_t spec_prefix_counts_proclists_type = SPEC_PREFIX_COUNTS_PROCLISTS_DEFAULT;


spint_t spec_prefix_counts(int *scounts, int *rcounts, int ncounts,
#ifdef SPEC_PROCLISTS
  spint_t nsend_procs, sproc_t *send_procs, spint_t nrecv_procs, sproc_t *recv_procs,
#endif
  int size, int rank, MPI_Comm comm) /* sp_func spec_prefix_counts */
{
  Z_TRACE_IF(SPEC_COMMON_TRACE_IF, "spec_prefix_counts");

#ifdef SPEC_PROCLISTS
  if ((nsend_procs >= 0 || nrecv_procs >= 0) && spec_prefix_counts_proclists_type != SPEC_PREFIX_COUNTS_PROCLISTS_IGNORE)
  {
    memset(rcounts, 0, ncounts * sizeof(int));

    switch (spec_prefix_counts_proclists_type)
    {
    }

  } else
#endif
  {
    switch (spec_prefix_counts_type)
    {
      case SPEC_PREFIX_COUNTS_SCAN:
        Z_TRACE_IF(SPEC_COMMON_TRACE_IF, "SPEC_PREFIX_COUNTS_SCAN");
        MPI_Scan(scounts, rcounts, ncounts, MPI_INT, MPI_SUM, comm);
        break;
    }
  }

  return 0;
}
