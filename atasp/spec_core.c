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

#include "spec_core.h"
#include "spec_common.h"


#if defined(Z_PACK_TIMING) && defined(SPEC_TIMING)
int nspec_timings = 0; /* sp_var nspec_timings */
double *spec_timing = NULL; /* sp_var spec_timing */
#endif


/*#define DUMMY_TEST*/

#ifdef DUMMY_TEST

static int dummy_spec_tproc(spec_elem_buf_t b, spec_elem_index_t x, spec_tproc_data_t tproc_data)
{
  return SPEC_PROC_NONE;
}

SPEC_DEFINE_TPROC(dummy, dummy_spec_tproc);

static int dummy_spec_tproc_mod(spec_elem_buf_t b, spec_elem_index_t x, spec_tproc_data_t tproc_data, spec_elem_buf_t mod)
{
  return SPEC_PROC_NONE;
}

SPEC_DEFINE_TPROC_MOD(dummy, dummy_spec_tproc_mod);

static int dummy_spec_tprocs(spec_elem_buf_t b, spec_elem_index_t x, spec_tproc_data_t tproc_data, int *procs)
{
  return 0;
}

SPEC_DEFINE_TPROCS(dummy, dummy_spec_tprocs);

static int dummy_spec_tprocs_mod(spec_elem_buf_t b, spec_elem_index_t x, spec_tproc_data_t tproc_data, int *procs, spec_elem_buf_t mod)
{
  return 0;
}

SPEC_DEFINE_TPROCS_MOD(dummy, dummy_spec_tprocs_mod);

#endif


static void spec_tproc_unset_tproc(spec_tproc_t tproc)
{
  tproc->tproc = NULL;
}


static void spec_tproc_unset_ext_tproc(spec_tproc_t tproc)
{
  spec_tproc_ext_t tproc_ext = SPEC_EXT_PARAM_TPROC_NULL;

  tproc->tproc_ext = tproc_ext;
}


static void spec_tproc_unset_tproc_mod(spec_tproc_t tproc)
{
  tproc->tproc_mod = NULL;
}


static void spec_tproc_unset_ext_tproc_mod(spec_tproc_t tproc)
{
  spec_tproc_mod_ext_t tproc_mod_ext = SPEC_EXT_PARAM_TPROC_MOD_NULL;

  tproc->tproc_mod_ext = tproc_mod_ext;
}


static void spec_tproc_unset_tprocs(spec_tproc_t tproc)
{
  tproc->tprocs = NULL;
}


static void spec_tproc_unset_ext_tprocs(spec_tproc_t tproc)
{
  spec_tprocs_ext_t tprocs_ext = SPEC_EXT_PARAM_TPROCS_NULL;

  tproc->tprocs_ext = tprocs_ext;
}


static void spec_tproc_unset_tprocs_mod(spec_tproc_t tproc)
{
  tproc->tprocs_mod = NULL;
}


static void spec_tproc_unset_ext_tprocs_mod(spec_tproc_t tproc)
{
  spec_tprocs_mod_ext_t tprocs_mod_ext = SPEC_EXT_PARAM_TPROCS_MOD_NULL;

  tproc->tprocs_mod_ext = tprocs_mod_ext;
}


spint_t spec_tproc_create(spec_tproc_t *tproc, spec_tproc_f *func, spec_tproc_mod_f *func_mod, spec_tprocs_f *func_s, spec_tprocs_mod_f *func_s_mod, spint_t max_tprocs) /* sp_func spec_tproc_create */
{
  *tproc = z_alloc(1, sizeof(struct _spec_tproc_t));

  (*tproc)->max_tprocs = max_tprocs;

  (*tproc)->tproc = func;
  spec_tproc_unset_ext_tproc(*tproc);

  (*tproc)->tproc_mod = func_mod;
  spec_tproc_unset_ext_tproc_mod(*tproc);

  (*tproc)->tprocs = func_s;
  spec_tproc_unset_ext_tprocs(*tproc);

  (*tproc)->tprocs_mod = func_s_mod;
  spec_tproc_unset_ext_tprocs_mod(*tproc);

  (*tproc)->reset = NULL;

#ifdef SPEC_PROCLISTS
  (*tproc)->nsend_procs = -1;
  (*tproc)->send_procs = NULL;
  (*tproc)->nrecv_procs = -1;
  (*tproc)->recv_procs = NULL;
#endif

#ifdef SPEC_COUNTS
  (*tproc)->send_counts = NULL;
  (*tproc)->recv_counts = NULL;
#endif

  return SPEC_EXIT_SUCCESS;
}


spint_t spec_tproc_destroy(spec_tproc_t *tproc) /* sp_func spec_tproc_destroy */
{
#ifdef SPEC_PROCLISTS
  spec_tproc_set_proclists(*tproc, -1, NULL, -1, NULL, 0, -1, MPI_COMM_NULL);
#endif

#ifdef SPEC_COUNTS
  spec_tproc_set_counts(*tproc, NULL, NULL, 0, -1, MPI_COMM_NULL);
#endif

  z_free(*tproc);

  *tproc = NULL;

  return SPEC_EXIT_SUCCESS;
}


spint_t spec_tproc_duplicate(spec_tproc_t tproc, spec_tproc_t *newtproc) /* sp_func spec_tproc_duplicate */
{
  *newtproc = z_alloc(1, sizeof(struct _spec_tproc_t));

  memcpy(*newtproc, tproc, sizeof(struct _spec_tproc_t));

  return SPEC_EXIT_SUCCESS;
}


spint_t spec_tproc_set_tproc(spec_tproc_t tproc, spec_tproc_f *func) /* sp_func spec_tproc_set_tproc */
{
  spec_tproc_unset_tproc(tproc);
  spec_tproc_unset_tproc_mod(tproc);
  spec_tproc_unset_tprocs(tproc);
  spec_tproc_unset_tprocs_mod(tproc);

  tproc->tproc = func;

  return SPEC_EXIT_SUCCESS;
}


spint_t spec_tproc_set_ext_tproc(spec_tproc_t tproc, const spec_tproc_ext_t *tproc_ext) /* sp_func spec_tproc_set_ext_tproc */
{
  tproc->tproc_ext = *tproc_ext;

  return SPEC_EXIT_SUCCESS;
}


spint_t spec_tproc_set_tproc_mod(spec_tproc_t tproc, spec_tproc_mod_f *func_mod) /* sp_func spec_tproc_set_tproc_mod */
{
  spec_tproc_unset_tproc(tproc);
  spec_tproc_unset_tproc_mod(tproc);
  spec_tproc_unset_tprocs(tproc);
  spec_tproc_unset_tprocs_mod(tproc);

  tproc->tproc_mod = func_mod;

  return SPEC_EXIT_SUCCESS;
}


spint_t spec_tproc_set_ext_tproc_mod(spec_tproc_t tproc, const spec_tproc_mod_ext_t *tproc_mod_ext) /* sp_func spec_tproc_set_ext_tproc_mod */
{
  tproc->tproc_mod_ext = *tproc_mod_ext;

  return SPEC_EXIT_SUCCESS;
}


spint_t spec_tproc_set_tprocs(spec_tproc_t tproc, spec_tprocs_f *func_s, spint_t max_tprocs) /* sp_func spec_tproc_set_tprocs */
{
  spec_tproc_unset_tproc(tproc);
  spec_tproc_unset_tproc_mod(tproc);
  spec_tproc_unset_tprocs(tproc);
  spec_tproc_unset_tprocs_mod(tproc);

  tproc->max_tprocs = max_tprocs;

  tproc->tprocs = func_s;

  return SPEC_EXIT_SUCCESS;
}


spint_t spec_tproc_set_ext_tprocs(spec_tproc_t tproc, const spec_tprocs_ext_t *tprocs_ext) /* sp_func spec_tproc_set_ext_tprocs */
{
  tproc->tprocs_ext = *tprocs_ext;

  return SPEC_EXIT_SUCCESS;
}


spint_t spec_tproc_set_tprocs_mod(spec_tproc_t tproc, spec_tprocs_mod_f *func_s_mod, spint_t max_tprocs) /* sp_func spec_tproc_set_tprocs_mod */
{
  spec_tproc_unset_tproc(tproc);
  spec_tproc_unset_tproc_mod(tproc);
  spec_tproc_unset_tprocs(tproc);
  spec_tproc_unset_tprocs_mod(tproc);

  tproc->max_tprocs = max_tprocs;

  tproc->tprocs_mod = func_s_mod;

  return SPEC_EXIT_SUCCESS;
}


spint_t spec_tproc_set_ext_tprocs_mod(spec_tproc_t tproc, const spec_tprocs_mod_ext_t *tprocs_mod_ext) /* sp_func spec_tproc_set_ext_tprocs_mod */
{
  tproc->tprocs_mod_ext = *tprocs_mod_ext;

  return SPEC_EXIT_SUCCESS;
}


spint_t spec_tproc_set_reset(spec_tproc_t tproc, spec_tproc_reset_f *reset) /* sp_func spec_tproc_set_reset */
{
  tproc->reset = reset;

  return SPEC_EXIT_SUCCESS;
}


#ifdef SPEC_PROCLISTS

void spec_make_recv_proclist(spint_t nsend_procs, sproc_t *send_procs, spint_t *nrecv_procs, sproc_t **recv_procs, int size, int rank, MPI_Comm comm) /* sp_func spec_make_recv_proclist */
{
  spint_t i, j;
  int *s, *r;


  s = z_alloc(2 * size, sizeof(int));
  r = s + 1 * size;

  for (i = 0; i < size; ++i) s[i] = 0;
  for (i = 0; i < nsend_procs; ++i) s[send_procs[i]] = 1;

  MPI_Alltoall(s, 1, MPI_INT, r, 1, MPI_INT, comm);

  for (j = 0, i = 0; i < size; ++i) j += r[i];

  *nrecv_procs = j;
  *recv_procs = z_alloc(*nrecv_procs, sizeof(sproc_t));

  for (j = 0, i = 0; i < size; ++i)
  if (r[i])
  {
    (*recv_procs)[j] = i;
    ++j;
  }

  z_free(s);
}


spint_t spec_tproc_set_proclists(spec_tproc_t tproc, spint_t nsend_procs, sproc_t *send_procs, spint_t nrecv_procs, sproc_t *recv_procs, int size, int rank, MPI_Comm comm) /* sp_func spec_tproc_set_proclists */
{
  spint_t i;

  if (tproc->send_procs) z_free(tproc->send_procs);
  if (tproc->recv_procs) z_free(tproc->recv_procs);

  tproc->nsend_procs = -1;
  tproc->send_procs = NULL;
  tproc->nrecv_procs = -1;
  tproc->recv_procs = NULL;

  if (nsend_procs >= 0)
  {
    tproc->nsend_procs = nsend_procs;
    tproc->send_procs = z_alloc(nsend_procs, sizeof(sproc_t));

    for (i = 0; i < nsend_procs; ++i) tproc->send_procs[i] = send_procs[i];
  }

  if (nrecv_procs >= 0)
  {
    tproc->nrecv_procs = nrecv_procs;
    tproc->recv_procs = z_alloc(nrecv_procs, sizeof(sproc_t));

    for (i = 0; i < nrecv_procs; ++i) tproc->recv_procs[i] = recv_procs[i];

  } else if (nsend_procs >= 0)
  {
    if (comm == MPI_COMM_NULL) return 1;

    spec_make_recv_proclist(nsend_procs, send_procs, &tproc->nrecv_procs, &tproc->recv_procs, size, rank, comm);
  }

/*  printf("%d: send_procs (%" spint_fmt ") = ", rank, tproc->nsend_procs);
  for (i = 0; i < tproc->nsend_procs; ++i) printf("  %" sproc_fmt, tproc->send_procs[i]);
  printf("\n");

  printf("%d: recv_procs (%" spint_fmt ") = ", rank, tproc->nrecv_procs);
  for (i = 0; i < tproc->nrecv_procs; ++i) printf("  %" sproc_fmt, tproc->recv_procs[i]);
  printf("\n");*/

  return SPEC_EXIT_SUCCESS;
}

#endif

#ifdef SPEC_COUNTS

spint_t spec_tproc_set_counts(spec_tproc_t tproc, int *send_counts, int *recv_counts, int size, int rank, MPI_Comm comm)
{
  if (tproc->send_counts) z_free(tproc->send_counts);
  if (tproc->recv_counts) z_free(tproc->recv_counts);

  tproc->send_counts = NULL;
  tproc->recv_counts = NULL;

  if (send_counts)
  {
    tproc->send_counts = z_alloc(size, sizeof(int));
    memcpy(tproc->send_counts, send_counts, size * sizeof(int));
  }

  if (recv_counts)
  {
    tproc->recv_counts = z_alloc(size, sizeof(int));
    memcpy(tproc->recv_counts, recv_counts, size * sizeof(int));
  }

  return SPEC_EXIT_SUCCESS;
}

#endif /* SPEC_COUNTS */


spint_t spec_print(spec_tproc_t tproc, spec_tproc_data_t tproc_data, spec_elem_t *b) /* sp_func spec_print */
{
  spint_t i, j, n;
  spec_proc_t p;

  spec_proc_t *procs = NULL;
  spec_elem_t *mods = NULL;


  /* setup tproc buffers */
  spec_tproc_setup(tproc, b, &procs, &mods);

  /* reset tproc if necessary */
  if (tproc->reset) tproc->reset(tproc_data);


  if (tproc->tproc)
  {
    for (i = 0; i < spec_elem_get_n(b); ++i)
    {
      p = tproc->tproc(spec_elem_get_buf(b), i, tproc_data);
      printf("%" spint_fmt ": %" spec_proc_fmt "\n", i, p);
    }

  } else if (tproc->tproc_mod)
  {
    for (i = 0; i < spec_elem_get_n(b); ++i)
    {
      p = tproc->tproc_mod(spec_elem_get_buf(b), i, tproc_data, spec_elem_get_buf(mods));
      printf("%" spint_fmt ": %" spec_proc_fmt "\n", i, p);
    }

  } else if (tproc->tprocs)
  {
    for (i = 0; i < spec_elem_get_n(b); ++i)
    {
      tproc->tprocs(spec_elem_get_buf(b), i, tproc_data, &n, procs);
      printf("%" spint_fmt ":", i);
      for (j = 0; j < n; ++j) printf(" %" spec_proc_fmt, procs[j]);
      printf("\n");
    }

  } else if (tproc->tprocs_mod)
  {
    for (i = 0; i < spec_elem_get_n(b); ++i)
    {
      tproc->tprocs_mod(spec_elem_get_buf(b), i, tproc_data, &n, procs, spec_elem_get_buf(mods));
      printf("%" spint_fmt ":", i);
      for (j = 0; j < n; ++j) printf(" %" spec_proc_fmt, procs[j]);
      printf("\n");
    }
  }

  /* free tproc buffers */
  spec_tproc_release(&procs, &mods);

  return 0;
}
