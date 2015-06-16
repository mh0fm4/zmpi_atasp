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

#include "spec_core.h"
#include "spec_common.h"


#ifndef SPEC_SENDRECV_TRACE_IF
# define SPEC_SENDRECV_TRACE_IF  (rank == -1)
#endif


/*#define SPEC_PRINT*/


/* sp_var spec_sendrecv_aux spec_sendrecv_aux_size spec_sendrecv_send_requests spec_sendrecv_receive_requests */
void *spec_sendrecv_aux = NULL;
spint_t spec_sendrecv_aux_size = -16*1024;
spint_t spec_sendrecv_send_requests = 10;
spint_t spec_sendrecv_receive_requests = 10;


#define SPEC_SENDRECV_SEND_ALL_AT_ONCE  0


spint_t spec_sendrecv_db(spec_elem_t *sb, spec_elem_t *rb, spec_elem_t *xb, spec_tproc_t tproc, spec_tproc_data_t tproc_data, int size, int rank, MPI_Comm comm) /* sp_func spec_sendrecv_db */
{
  spint_t exit_code = SPEC_EXIT_SUCCESS;

  spint_t i, j, nprocs;
  spec_proc_t p;

  spec_proc_t *procs = NULL;
  spec_elem_t *mods = NULL;

  int *scounts, rcounts[3], ii;

  spint_t scount, sdispl, rdispl, rdispl_old, rdisplpart;
  spint_t stotal, rtotal, sdone, rdone, nfullrecvs, npartrecvs;

  const spint_t sreqs_max = spec_sendrecv_send_requests;
  const spint_t rreqs_max = spec_sendrecv_receive_requests;

/*  printf("%d: requests: %d / %d\n", rank, (int) sreqs_max, (int) rreqs_max);*/

  MPI_Request reqs[(rreqs_max + sreqs_max) * spec_elem_x];
  MPI_Status stats[(rreqs_max + sreqs_max) * spec_elem_x];
  int inds[(rreqs_max + sreqs_max) * spec_elem_x], completed_nreqs;
  spint_t completed_nsends, completed_nrecvs;
  spint_t reqs_i[(rreqs_max + sreqs_max) * 2], reqs_tmp;
  spint_t rreqs_nfree, rreqs_free[rreqs_max], sreqs_nfree, sreqs_free[sreqs_max];
  const spint_t rreqs_base = 0;
  const spint_t sreqs_base = rreqs_max;

  const spint_t aux_size_min = 1;
  const spint_t aux_size_max = -1;
  spint_t aux_size, aux_n, *aux_bases, *aux_displs, *aux_queue, aux_queue_size, aux_queue_first, aux_queue_next, aux_done, aux_tmp, aux_alloc;
  spec_elem_t aux;

  SPEC_DECLARE_TPROC_SENDRECV_DB

#define _MAX_ROUNDS 10

#ifdef MAX_ROUNDS
  spint_t max_rounds = MAX_ROUNDS;
#endif

#ifdef Z_PACK_TIMING
  double tt, t[12] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
#endif


  Z_TIMING_SYNC(comm); Z_TIMING_START(t[0]);

  /* check supported tproc functions */
  if (!spec_check_tproc_support(tproc, 1, 1, 1, 1, rank, "spec_sendrecv_buffer_db"))
  {
    spec_elem_set_n(rb, 0);
    exit_code = SPEC_EXIT_FAILED;
    goto exit;
  }

  scounts = z_alloc(size * 3, sizeof(int));

  /* setup tproc buffers */
  spec_tproc_setup(tproc, sb, &procs, &mods);

#ifdef SPEC_PRINT
  printf("input\n");
  spec_print(&tproc, tproc_data, sb);
#endif

  /* reset tproc if necessary */
  if (tproc->reset) tproc->reset(tproc_data);

  /* make local counts */
  Z_TIMING_SYNC(comm); Z_TIMING_START(t[1]);

  spec_make_counts(tproc, tproc_data, sb, 0, size, scounts, procs);

  Z_TIMING_SYNC(comm); Z_TIMING_STOP(t[1]);

  /* make total, full, and partial receives */

#ifdef SPEC_PROCLISTS
  if (tproc->nsend_procs >= 0 || tproc->nrecv_procs >= 0)
  {
    aux_n = tproc->nsend_procs;

  } else
#endif
  {
    aux_n = size - 1;
  }

  if (spec_sendrecv_aux_size > 0) aux_size = spec_elem_sizefor(sb, spec_sendrecv_aux_size) / z_max(1, aux_n);
  else aux_size = spec_elem_sizefor(sb, -spec_sendrecv_aux_size);

  if (aux_size_min > 0) aux_size = z_max(aux_size, aux_size_min);
  if (aux_size_max > 0) aux_size = z_min(aux_size, aux_size_max);

  Z_TRACE_IF(SPEC_SENDRECV_TRACE_IF, "aux_n: %" spint_fmt ", aux_size: %" spint_fmt, aux_n, aux_size);

  scount = spec_elem_get_n(sb);

  rdisplpart = scounts[rank];

  stotal = 0;
  for (i = size - 1; i >= 0; --i)
  {
    stotal += scounts[i];
    scounts[3 * i + 0] = scounts[i];
    scounts[3 * i + 1] = scounts[i] / aux_size;
    scounts[3 * i + 2] = (scounts[i] % aux_size > 0)?1:0;

    Z_TRACE_IF(SPEC_SENDRECV_TRACE_IF, "scount[%" spint_fmt "]: %d,%d,%d", i, scounts[3 * i + 0], scounts[3 * i + 1], scounts[3 * i + 2]);
  }

  Z_TIMING_SYNC(comm); Z_TIMING_START(t[2]);

  spec_reduce_scatter_counts(scounts, rcounts, 3,
#ifdef SPEC_PROCLISTS
    tproc->nsend_procs, tproc->send_procs, tproc->nrecv_procs, tproc->recv_procs,
#endif
    size, rank, comm);

  Z_TIMING_SYNC(comm); Z_TIMING_STOP(t[2]);

  rtotal = rcounts[0];
  nfullrecvs = rcounts[1] - scounts[3 * rank + 1];
  npartrecvs = rcounts[2] - scounts[3 * rank + 2];

  rdisplpart += nfullrecvs * aux_size;

  /* check size of receive buffer */
  if (!spec_check_buffer_size(rb, rtotal, 1, comm, rank, "spec_sendrecv_buffer_db", "receive"))
  {
    spec_elem_set_n(rb, rtotal);
    exit_code = SPEC_EXIT_FAILED;
    goto free_and_exit;
  }

  /* reset tproc if necessary */
  if (tproc->reset) tproc->reset(tproc_data);

  Z_TRACE_IF(SPEC_SENDRECV_TRACE_IF, "rtotal = %" spint_fmt ", nfullrecvs = %" spint_fmt ", npartrecvs = %" spint_fmt, rtotal, nfullrecvs, npartrecvs);

  /* redistribute */
  Z_TIMING_SYNC(comm); Z_TIMING_START(t[3]);

#define TAG_FULL  0
#define TAG_PART  spec_elem_x

#define AUX_BASE(_p_)         aux_bases[_p_]
#define AUX_DISPL(_p_)        aux_displs[_p_]
#define AUX_DISPL_BEGIN(_p_)  AUX_BASE(_p_)
#define AUX_DISPL_END(_p_)    AUX_BASE((_p_) + 1)
#define AUX_DISPL_SET(_p_, _d_)  aux_displs[_p_] = (_d_)
#define AUX_DISPL_INC(_p_)    ++aux_displs[_p_]
#define AUX_SIZE(_p_)         AUX_DISPL(_p_) - AUX_DISPL_BEGIN(_p_)
#define AUX_ENQUEUE(_p_)      Z_MOP( \
  Z_TRACE_IF(SPEC_SENDRECV_TRACE_IF, "enqueue aux send of %" spint_fmt " to process %" spec_proc_fmt, AUX_SIZE(_p_), (_p_)); \
  aux_queue[aux_queue_next] = (_p_); ++aux_queue_next; aux_queue_next %= aux_queue_size;)
#define AUX_DEQUEUE()         (aux_tmp = aux_queue[aux_queue_first], ++aux_queue_first, aux_queue_first %= aux_queue_size, aux_tmp)
#define AUX_QUEUED()          ((aux_queue_next + aux_queue_size - aux_queue_first) % aux_queue_size)

#define REQS                          reqs
#define REQS_N                        (rreqs_max + sreqs_max) * spec_elem_x
#define REQS_PTR(_r_, _x_)            &reqs[(_r_) * spec_elem_x + (_x_)]

#define REQS_COMPLETE(_r_)            reqs_i[2 * (_r_) + 1]
#define REQS_COMPLETE_RESET(_r_)      reqs_i[2 * (_r_) + 1] = spec_elem_x
#define REQS_COMPLETE_ONE(_r_)        --reqs_i[2 * (_r_) + 1]
#define REQS_COMPLETED(_r_)           (REQS_COMPLETE(_r_) == 0)
#define REQS_COMPLETED_FIRST(_r_)     (REQS_COMPLETE(_r_) == spec_elem_x - 1)

#define REQS_NCUR()                   (sreqs_max - sreqs_nfree + rreqs_max - rreqs_nfree)

#define REQS_SEND_DST_SET(_r_, _p_)   reqs_i[2 * (_r_) + 0] = (_p_)
#define REQS_SEND_DST(_r_)            reqs_i[2 * (_r_) + 0]
#define REQS_SEND_SIZE(_r_)           AUX_SIZE(REQS_SEND_DST(_r_))
#define REQS_SEND_NFREE()             sreqs_nfree
#define REQS_SEND_FREE(_r_)           Z_MOP(sreqs_free[sreqs_nfree] = (_r_); ++sreqs_nfree; reqs_i[2 * (_r_) + 0] = reqs_i[2 * (_r_) + 1] = -2501;)

#define REQS_RECV(_r_)                ((_r_) < rreqs_max)
#define REQS_RECV_SIZE(_r_)           reqs_i[2 * (_r_) + 0]
#define REQS_RECV_SIZE_SET(_r_, _s_)  reqs_i[2 * (_r_) + 0] = (_s_)
#define REQS_RECV_NFREE()             rreqs_nfree
#define REQS_RECV_FREE(_r_)           Z_MOP(rreqs_free[rreqs_nfree] = (_r_); ++rreqs_nfree; reqs_i[2 * (_r_) + 0] = reqs_i[2 * (_r_) + 1] = -2501;)
#define REQS_RECV_FULL_FIRST(_r_)     (REQS_COMPLETE(_r_) == 1)

#define SEND_AUX_FIRST(_p_)  Z_MOP( \
  Z_TRACE_IF(SPEC_SENDRECV_TRACE_IF, "send_aux_first: process: %" spec_proc_fmt ", size: %" spint_fmt ", displ: %" spint_fmt, (_p_), AUX_SIZE(_p_), AUX_BASE(_p_)); \
  --sreqs_nfree; reqs_tmp = sreqs_free[sreqs_nfree]; \
  spec_elem_isend_first(&aux, AUX_BASE(_p_), AUX_SIZE(_p_), _p_, (AUX_SIZE(_p_) == aux_size)?TAG_FULL:TAG_PART, REQS_PTR(reqs_tmp, 0), size, rank, comm); \
  REQS_COMPLETE_RESET(reqs_tmp); \
  REQS_SEND_DST_SET(reqs_tmp, _p_); \
)

#define SEND_AUX_NEXT(_p_, _r_)  Z_MOP( \
  Z_TRACE_IF(SPEC_SENDRECV_TRACE_IF, "send_aux_next: process: %" spec_proc_fmt ", req: %" spint_fmt ", size: %" spint_fmt ", displ: %" spint_fmt, (_p_), (_r_), AUX_SIZE(_p_), AUX_BASE(_p_)); \
  spec_elem_isend_next(&aux, AUX_BASE(_p_), AUX_SIZE(_p_), _p_, (AUX_SIZE(_p_) == aux_size)?TAG_FULL:TAG_PART, REQS_PTR(_r_, 0), size, rank, comm); \
)

#define RECV_FULL_FIRST()  Z_MOP( \
  Z_TRACE_IF(SPEC_SENDRECV_TRACE_IF, "recv_full_first: displ: %" spint_fmt, rdispl); \
  --rreqs_nfree; reqs_tmp = rreqs_free[rreqs_nfree]; \
  spec_elem_irecv_first(rb, rdispl, aux_size, MPI_ANY_SOURCE, TAG_FULL, REQS_PTR(reqs_tmp, 0), size, rank, comm); \
  REQS_COMPLETE_RESET(reqs_tmp); \
  REQS_RECV_SIZE_SET(reqs_tmp, rdispl); \
  rdispl += aux_size; \
)

#define RECV_FULL_NEXT(_p_, _n_, _r_)  Z_MOP( \
  Z_TRACE_IF(SPEC_SENDRECV_TRACE_IF, "recv_full_next: process: %" spec_proc_fmt ", req: %" spint_fmt ", size: %" spint_fmt ", displ: %" spint_fmt, (_p_), (_r_), (_n_), REQS_RECV_SIZE(_r_)); \
  spec_elem_irecv_next(rb, REQS_RECV_SIZE(_r_), _n_, _p_, TAG_FULL, REQS_PTR(_r_, 0), size, rank, comm); \
  REQS_RECV_SIZE_SET(_r_, _n_); \
)

#define RECV_PART_FIRST()  Z_MOP( \
  Z_TRACE_IF(SPEC_SENDRECV_TRACE_IF, "recv_part_first: displ: %" spint_fmt, rdisplpart); \
  --rreqs_nfree; reqs_tmp = rreqs_free[rreqs_nfree]; \
  spec_elem_irecv_first(rb, rdisplpart, aux_size, MPI_ANY_SOURCE, TAG_PART, REQS_PTR(reqs_tmp, 0), size, rank, comm); \
  REQS_COMPLETE_RESET(reqs_tmp); \
)

#define RECV_PART_NEXT(_p_, _n_, _r_)  Z_MOP( \
  Z_TRACE_IF(SPEC_SENDRECV_TRACE_IF, "recv_part_next: process: %" spec_proc_fmt ", req: %" spint_fmt ", size: %" spint_fmt ", displ: %" spint_fmt, (_p_), (_r_), (_n_), rdisplpart); \
  spec_elem_irecv_next(rb, rdisplpart, _n_, _p_, TAG_PART, REQS_PTR(_r_, 0), size, rank, comm); \
  REQS_RECV_SIZE_SET(_r_, _n_); \
  rdisplpart += _n_; \
)

  p = SPEC_PROC_NONE;
  nprocs = -1;

  sdone = rdone = 0;
  sdispl = rdispl = 0;

  rreqs_nfree = rreqs_max;
  for (i = 0; i < rreqs_max; ++i)
  {
    rreqs_free[i] = rreqs_base + i;
    for (j = 0; j < spec_elem_x; ++j) reqs[rreqs_free[i] * spec_elem_x + j] = MPI_REQUEST_NULL;
    reqs_i[2 * rreqs_free[i] + 0] = reqs_i[2 * rreqs_free[i] + 1] = -2501;
  }
  sreqs_nfree = sreqs_max;
  for (i = 0; i < sreqs_max; ++i)
  {
    sreqs_free[i] = sreqs_base + i;
    for (j = 0; j < spec_elem_x; ++j) reqs[sreqs_free[i] * spec_elem_x + j] = MPI_REQUEST_NULL;
    reqs_i[2 * sreqs_free[i] + 0] = reqs_i[2 * sreqs_free[i] + 1] = -2501;
  }

  aux_queue_size = aux_n + 1;  /* aux queue is a circular buffer implemented with only two indices 'next' (to write) and 'first' (to read), thus the queue has to be one slot larger than max. number of required slots to prevent a full queue */
  aux_bases = z_alloc(2 * size + 1 + aux_queue_size, sizeof(spint_t));
  aux_displs = aux_bases + size + 1;

  for (i = 0; i < size; ++i) aux_bases[i] = 0;
#ifdef SPEC_PROCLISTS
  if (tproc->nsend_procs >= 0 || tproc->nrecv_procs >= 0)
  {
    for (i = 0; i < tproc->nsend_procs; ++i) aux_bases[tproc->send_procs[i]] = aux_size;

  } else
#endif
  {
    for (i = 0; i < size; ++i) aux_bases[i] = aux_size;
  }
  aux_bases[rank] = 0;
  j = 0;
  for (i = 0; i < size; ++i)
  {
    aux_displs[i] = j;
    j += aux_bases[i];
    aux_bases[i] = aux_displs[i];
  }
  aux_bases[size] = j;

  aux_queue = aux_displs + size;
  aux_queue_first = aux_queue_next = 0;
  aux_done = 0;

  spec_elem_unset(&aux);

  spec_elem_copy_type(sb, &aux);

  spec_elem_set_nmax(&aux, 0);
#ifdef spec_elem_alloc_tmp_from_block
  if (spec_sendrecv_aux)
  {
    i = (spec_sendrecv_aux_size > 0)?spec_sendrecv_aux_size:(-spec_sendrecv_aux_size * aux_n);
    spec_elem_alloc_tmp_from_block(&aux, spec_sendrecv_aux, i);
  }
#endif

  aux_alloc = 0;
  if (spec_elem_get_nmax(&aux) < aux_size * aux_n)
  {
    aux_alloc = 1;
    spec_elem_alloc_tmp(&aux, aux_n * aux_size);
  }

  Z_TIMING_SYNC(comm); Z_TIMING_STOP(t[3]);

  Z_TIMING_SYNC(comm); Z_TIMING_START(t[4]);

  Z_TIMING_DECL(const int ottmps = 5;);

  while (1)
  {
    Z_TRACE_IF(SPEC_SENDRECV_TRACE_IF, "send: %" spint_fmt " < %" spint_fmt ", recv: %" spint_fmt " < %" spint_fmt "", sdone, stotal, rdone, rtotal);

    if (!(sdone < stotal || rdone < rtotal)) break;

#ifdef MAX_ROUNDS
    if (max_rounds-- <= 0) break;
#endif

    /* if there are no more data elements, then skip processing loop */
    if (sdispl >= scount) goto do_send;

    Z_TRACE_IF(SPEC_SENDRECV_TRACE_IF, "loop over data elements");

    Z_TIMING_START(tt);

    rdispl_old = rdispl;

    if (tproc->tproc)
    {
#if 0
      if (tproc->tproc_ext.sendrecv_db) p = tproc->tproc_ext.sendrecv_db(tproc_data, sb, rb, scount, &sdispl, &rdispl, &aux, aux_displs, aux_size, aux_queue, &aux_queue_next, aux_queue_size, rank, p);
      else SPEC_DO_TPROC_SENDRECV_DB(tproc->tproc, tproc_data, sb, rb, scount, &sdispl, &rdispl, &aux, aux_displs, aux_size, aux_queue, &aux_queue_next, aux_queue_size, rank, p);
#else
      while (sdispl < scount)
      {
        if (p == SPEC_PROC_NONE) p = tproc->tproc(spec_elem_get_buf(sb), sdispl, tproc_data);

        if (p != SPEC_PROC_NONE)
        {
          if (p == rank)
          {
            spec_elem_copy_at(sb, sdispl, rb, rdispl);
            ++rdispl;

          } else
          {
            if (AUX_DISPL(p) >= AUX_DISPL_END(p)) break;

            spec_elem_copy_at(sb, sdispl, &aux, AUX_DISPL(p));

            AUX_DISPL_INC(p);

            if (AUX_DISPL(p) >= AUX_DISPL_END(p)) AUX_ENQUEUE(p);
          }
        }

        p = SPEC_PROC_NONE;
        ++sdispl;
      }
#endif

    } else if (tproc->tproc_mod)
    {
      while (sdispl < scount)
      {
        if (p == SPEC_PROC_NONE) p = tproc->tproc_mod(spec_elem_get_buf(sb), sdispl, tproc_data, spec_elem_get_buf(mods));

        if (p != SPEC_PROC_NONE)
        {
          if (p == rank)
          {
            spec_elem_copy_at(mods, 0, rb, rdispl);
            ++rdispl;

          } else
          {
            if (AUX_DISPL(p) >= AUX_DISPL_END(p)) break;

            spec_elem_copy_at(mods, 0, &aux, AUX_DISPL(p));

            AUX_DISPL_INC(p);

            if (AUX_DISPL(p) >= AUX_DISPL_END(p)) AUX_ENQUEUE(p);
          }
        }

        p = SPEC_PROC_NONE;
        ++sdispl;
      }

    } else if (tproc->tprocs)
    {
      while (sdispl < scount)
      {
        if (nprocs < 0) { tproc->tprocs(spec_elem_get_buf(sb), sdispl, tproc_data, &nprocs, procs); --nprocs; }

        while (nprocs >= 0)
        {
          if (procs[nprocs] == rank)
          {
            spec_elem_copy_at(sb, sdispl, rb, rdispl);
            ++rdispl;

          } else
          {
            if (AUX_DISPL(procs[nprocs]) >= AUX_DISPL_END(procs[nprocs])) break;

            spec_elem_copy_at(sb, sdispl, &aux, AUX_DISPL(procs[nprocs]));

            AUX_DISPL_INC(procs[nprocs]);

            if (AUX_DISPL(procs[nprocs]) >= AUX_DISPL_END(procs[nprocs])) AUX_ENQUEUE(procs[nprocs]);
          }

          --nprocs;
        }

        if (nprocs >= 0) break;

        ++sdispl;
      }

    } else if (tproc->tprocs_mod)
    {
      while (sdispl < scount)
      {
        if (nprocs < 0) { tproc->tprocs_mod(spec_elem_get_buf(sb), sdispl, tproc_data, &nprocs, procs, spec_elem_get_buf(mods)); --nprocs; }

        while (nprocs >= 0)
        {
          if (procs[nprocs] == rank)
          {
            spec_elem_copy_at(mods, nprocs, rb, rdispl);
            ++rdispl;

          } else
          {
            if (AUX_DISPL(procs[nprocs]) >= AUX_DISPL_END(procs[nprocs])) break;

            spec_elem_copy_at(mods, nprocs, &aux, AUX_DISPL(procs[nprocs]));

            AUX_DISPL_INC(procs[nprocs]);

            if (AUX_DISPL(procs[nprocs]) >= AUX_DISPL_END(procs[nprocs])) AUX_ENQUEUE(procs[nprocs]);
          }

          --nprocs;
        }

        if (nprocs >= 0) break;

        ++sdispl;
      }
    }

    sdone += rdispl - rdispl_old;
    rdone += rdispl - rdispl_old;

    Z_TIMING_STOP_ADD(tt, t[ottmps + 0]);

    Z_TRACE_IF(SPEC_SENDRECV_TRACE_IF, "scount: %" spint_fmt ", sdispl: %" spint_fmt ", rdispl: %" spint_fmt ", rdisplpart: %" spint_fmt, scount, sdispl, rdispl, rdisplpart);
    Z_TRACE_IF(SPEC_SENDRECV_TRACE_IF, "stotal: %" spint_fmt ", sdone: %" spint_fmt ", rtotal: %" spint_fmt ", rdone: %" spint_fmt, stotal, sdone, rtotal, rdone);
    Z_TRACE_IF(SPEC_SENDRECV_TRACE_IF, "aux_queue: first: %" spint_fmt ", next: %" spint_fmt, aux_queue_first, aux_queue_next);

    Z_TIMING_START(tt);

    if (!aux_done && sdispl >= scount)
    {
      Z_TRACE_IF(SPEC_SENDRECV_TRACE_IF, "no more local data elements, enqueueing all remaining aux sends");

      for (i = 0; i < size; ++i)
      {
        j = (rank + i) % size;
        if (AUX_SIZE(j) > 0 && AUX_SIZE(j) < aux_size) AUX_ENQUEUE((spec_proc_t) j);
      }

      aux_done = 1;

      Z_TRACE_IF(SPEC_SENDRECV_TRACE_IF, "aux_queue: first: %" spint_fmt ", next: %" spint_fmt, aux_queue_first, aux_queue_next);
    }

    Z_TIMING_STOP_ADD(tt, t[ottmps + 1]);

do_send:
    Z_TRACE_IF(SPEC_SENDRECV_TRACE_IF, "A: sreqs_nfree: %" spint_fmt ", aux_queued: %" spint_fmt "", REQS_SEND_NFREE(), AUX_QUEUED());

    Z_TIMING_START(tt);

    /* initiate queued sends */
    while (REQS_SEND_NFREE() > 0 && AUX_QUEUED() > 0)
    {
      i = AUX_DEQUEUE();
      SEND_AUX_FIRST((spec_proc_t) i);
#if SPEC_SENDRECV_SEND_ALL_AT_ONCE
      SEND_AUX_NEXT((spec_proc_t) i, reqs_tmp);
#endif
    }

    Z_TIMING_STOP_ADD(tt, t[ottmps + 2]);

do_recv:
    Z_TRACE_IF(SPEC_SENDRECV_TRACE_IF, "B: sreqs_nfree: %" spint_fmt ", rreqs_nfree: %" spint_fmt "", REQS_SEND_NFREE(), REQS_RECV_NFREE());

    Z_TIMING_START(tt);

    /* initiate partial recv */
    if (REQS_RECV_NFREE() > 0 && npartrecvs > 0)
    {
      Z_TRACE_IF(SPEC_SENDRECV_TRACE_IF, "initiate part first recv");
      RECV_PART_FIRST();
      --npartrecvs;
      npartrecvs *= -1; /* negative value signals that a part recv is active */
    }

    /* initiate remaining full recvs */
    while (REQS_RECV_NFREE() > 0 && nfullrecvs > 0)
    {
      Z_TRACE_IF(SPEC_SENDRECV_TRACE_IF, "initiate full first recv");
      RECV_FULL_FIRST();
      --nfullrecvs;
    }

    Z_TIMING_STOP_ADD(tt, t[ottmps + 3]);

/*    Z_TRACE_IF(SPEC_SENDRECV_TRACE_IF, "reqs:");
    for (i = 0; i < rreqs_max + sreqs_max; ++i) Z_TRACE_IF(SPEC_SENDRECV_TRACE_IF, "  (%" spint_fmt " / %" spint_fmt ")", reqs_i[2 * i + 0], reqs_i[2 * i + 1]);*/

do_check_requests:
    Z_TRACE_IF(SPEC_SENDRECV_TRACE_IF, "C: sreqs_nfree: %" spint_fmt ", rreqs_nfree: %" spint_fmt "", REQS_SEND_NFREE(), REQS_RECV_NFREE());

    if (REQS_NCUR() == 0)
    {
      Z_TRACE_IF(SPEC_SENDRECV_TRACE_IF, "no requests, continue loop");
      continue;
    }

    Z_TIMING_START(tt);

#if 0
    Z_TRACE_IF(SPEC_SENDRECV_TRACE_IF, "waitany:");
    MPI_Waitany(REQS_N, REQS, &inds[0], &stats[0]);
    completed_nreqs = (inds[0] == MPI_UNDEFINED)?0:1;
#else
    Z_TRACE_IF(SPEC_SENDRECV_TRACE_IF, "waitsome:");
    MPI_Waitsome(REQS_N, REQS, &completed_nreqs, inds, stats);
#endif

    Z_TIMING_STOP_ADD(tt, t[ottmps + 4]);
    Z_TIMING_CMD(++t[ottmps + 5];);

    Z_TRACE_IF(SPEC_SENDRECV_TRACE_IF, "completed requests: %d", completed_nreqs);

    Z_TIMING_START(tt);

    completed_nsends = completed_nrecvs = 0;

    for (j = 0; j < completed_nreqs; ++j)
    {
      Z_TRACE_IF(SPEC_SENDRECV_TRACE_IF, "check request #%" spint_fmt ": index: %d", j, inds[j]);

      if (inds[j] == MPI_UNDEFINED) continue;

      i = inds[j] / spec_elem_x;

      REQS_COMPLETE_ONE(i);

      Z_TRACE_IF(SPEC_SENDRECV_TRACE_IF, "%s-request %" spint_fmt ": (%" spint_fmt " / %" spint_fmt ")", (REQS_RECV(i))?"recv":"send", i, reqs_i[2 * i + 0], reqs_i[2 * i + 1]);

      if (REQS_RECV(i)) /* recv-req done */
      {
        Z_TRACE_IF(SPEC_SENDRECV_TRACE_IF, "recv-request %" spint_fmt ": tag: %d, first: %s, completed: %s", i, stats[j].MPI_TAG, REQS_COMPLETED_FIRST(i)?"yes":"no", REQS_COMPLETED(i)?"yes":"no");

        /* if first component of receive, then receive next */
        if (REQS_COMPLETED_FIRST(i))
        {
          spec_elem_get_recv_count(rb, &stats[j], &ii);
          if (stats[j].MPI_TAG == TAG_FULL) RECV_FULL_NEXT(stats[j].MPI_SOURCE, (spint_t) ii, i);
          else RECV_PART_NEXT(stats[j].MPI_SOURCE, (spint_t) ii, i);
        }

        /* nothing was complete, thus continue with next completed request */
        if (!REQS_COMPLETED(i)) continue;

        /* count completed receives */
        ++completed_nrecvs;

        /* register recv */
        rdone += REQS_RECV_SIZE(i);

        /* free recv-req */
        REQS_RECV_FREE(i);

        /* set part recv inactive */
        if (stats[j].MPI_TAG >= TAG_PART) npartrecvs *= -1;

      } else /* send-req done */
      {
        Z_TRACE_IF(SPEC_SENDRECV_TRACE_IF, "send-request %" spint_fmt " first: %s, completed: %s", i, REQS_COMPLETED_FIRST(i)?"yes":"no", REQS_COMPLETED(i)?"yes":"no");

#if !(SPEC_SENDRECV_SEND_ALL_AT_ONCE)
        /* if first component of send, then send next */
        if (REQS_COMPLETED_FIRST(i)) SEND_AUX_NEXT((spec_proc_t) REQS_SEND_DST(i), i);
#endif

        /* nothing was complete, thus continue with next completed request */
        if (!REQS_COMPLETED(i)) continue;

        /* count completed sends */
        ++completed_nsends;

        /* register send */
        sdone += REQS_SEND_SIZE(i);

        /* reset aux */
        AUX_DISPL_SET(REQS_SEND_DST(i), AUX_DISPL_BEGIN(REQS_SEND_DST(i)));

        /* free send-req */
        REQS_SEND_FREE(i);

        Z_TRACE_IF(SPEC_SENDRECV_TRACE_IF, "sreqs_nfree: %" spint_fmt "", sreqs_nfree);
      }
    }

    Z_TIMING_STOP_ADD(tt, t[ottmps + 6]);

    Z_TRACE_IF(SPEC_SENDRECV_TRACE_IF, "completed sends: %" spint_fmt ", completed receives: %" spint_fmt "", completed_nsends, completed_nrecvs);

    /* if nothing was completed, then wait for next requests */
    if (completed_nsends == 0 && completed_nrecvs == 0) goto do_check_requests;

    /* if only receives were completed, then start new receives */
    if (completed_nsends == 0 && completed_nrecvs > 0) goto do_recv;
  }

  Z_TIMING_SYNC(comm); Z_TIMING_STOP(t[4]);

  if (aux_alloc) spec_elem_free_tmp(&aux);

  z_free(aux_bases);

#if 1
  Z_TIMING_CMD(
    const int nttmps = (sizeof(t) / sizeof(double)) - ottmps;
    double ttmps[nttmps];
    MPI_Allreduce(&t[ottmps], ttmps, nttmps, MPI_DOUBLE, MPI_MAX, comm);
    for (i = 0; i < nttmps; ++i) t[ottmps + i] = ttmps[i];
  );
#endif

  spec_elem_set_n(rb, rtotal);

free_and_exit:

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

#undef TAG_FULL
#undef TAG_PART

#undef AUX_BASE
#undef AUX_DISPL
#undef AUX_DISPL_BEGIN
#undef AUX_DISPL_END
#undef AUX_DISPL_SET
#undef AUX_DISPL_INC
#undef AUX_SIZE
#undef AUX_ENQUEUE
#undef AUX_DEQUEUE
#undef AUX_QUEUED

#undef REQS
#undef REQS_N
#undef REQS_PTR

#undef REQS_COMPLETE
#undef REQS_COMPLETE_RESET
#undef REQS_COMPLETE_ONE
#undef REQS_COMPLETED

#undef REQS_NCUR

#undef REQS_SEND_DST_SET
#undef REQS_SEND_DST
#undef REQS_SEND_SIZE
#undef REQS_SEND_NFREE
#undef REQS_SEND_FREE

#undef REQS_RECV
#undef REQS_RECV_SIZE
#undef REQS_RECV_SIZE_SET
#undef REQS_RECV_NFREE
#undef REQS_RECV_FREE
#undef REQS_RECV_FULL_FIRST

#undef SEND_AUX
#undef RECV_FULL_FIRST
#undef RECV_FULL_NEXT
#undef RECV_PART_FIRST
#undef RECV_PART_NEXT
