/*
  Copyright (c) 2015 Annika Seidel <annika.seidel@studium.uni-hamburg.de>
  Copyright (c) 2015 Center for Bioinformatics, University of Hamburg

  Permission to use, copy, modify, and distribute this software for any
  purpose with or without fee is hereby granted, provided that the above
  copyright notice and this permission notice appear in all copies.

  THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES
  WITH REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF
  MERCHANTABILITY AND FITNESS. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR
  ANY SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
  WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN
  ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF
  OR IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.
*/

#ifndef LINSPACE_MANAGEMENT_H
#define LINSPACE_MANAGEMENT_H

#include "core/types_api.h"

/* The <GtLinspaceManagement> interface. All useful uitlities for different
   algorithms, which work in linear space */

/* The following type is to organize all allocated space */
typedef struct GtLinspaceManagement GtLinspaceManagement;

/* Return a <GtLinspaceManagement> object. */
GtLinspaceManagement* gt_linspace_management_new();

/* Delete the given <spacemenager>. */
void          gt_linspace_management_delete(GtLinspaceManagement *spacemanager);
/* Return bytes of allocated space for linear space algorithms of the given
   <spacemanager>*/
size_t        gt_linspace_management_get_spacepeak(const GtLinspaceManagement
                                                                 *spacemanager);
/* Check if enough space is allocated for the three tabs, whichs are stored in
   the given <spacemanager> for global alignments and resize if necessary. The
   required space depends on the size of one entry, given by <valuesize>,
   <rtabsize>, <crosspointsize> and the count of values given by sequence
   lengths <ulen> and <vlen>. */
void          gt_linspace_management_check(GtLinspaceManagement *spacemanager,
                                           GtUword ulen, GtUword vlen,
                                           size_t valuesize,
                                           size_t rtabsize,
                                           size_t crosspointsize);
/* Check if enough space is allocated for the two tabs, which are stored in
   the given <spacemanager> for local alignments and resize if necessary. The
   required space depends on of the size of one entry, given by <valuesize>,
   <rtabsize> and the count of values given by sequence lengths <ulen>
   and <vlen>. */
void          gt_linspace_management_check_local(GtLinspaceManagement
                                                 *spacemanager,
                                                 GtUword ulen, GtUword vlen,
                                                 size_t valuesize,
                                                 size_t rstabsize);

/* Check if enough space to use square space functions in global case is
   allocated in the given <spacemanager>. To use this function <valuesize>,
   <rsize> and the sequence lengths <ulen> and <vlen> are required. If a
   timesquarefactor is set, the function checks also the relation with this
   factor and resize the space of <spacemanager> if is it necessary.
   */
bool          gt_linspace_management_checksquare(GtLinspaceManagement
                                                 *spacemanager,
                                                 GtUword ulen, GtUword vlen,
                                                 size_t valuesize,
                                                 size_t rsize);

/* Check if enough space to use square space functions in local case is
   allocated in the given <spacemanager>. To use this function <valuesize>,
   <rsize> and the sequence lengths <ulen> and <vlen> are required. If a
   timesquarefactor is set, the function checks also the relation with this
   factor and resize the space of <spacemanager> if is it necessary. */
bool          gt_linspace_management_checksquare_local(GtLinspaceManagement
                                                       *spacemanager,
                                                       GtUword ulen,
                                                       GtUword vlen,
                                                       size_t valuesize,
                                                       size_t rsize);
/* Set sequence length <ulen> for the given <spacemanager>. */
void          gt_linspace_management_set_ulen(GtLinspaceManagement
                                              *spacemanager,
                                              GtUword ulen);
/* Return pointer to valueTab space of the given <spacemanager>. */
void*         gt_linspace_management_get_valueTabspace(
                                                      const GtLinspaceManagement
                                                      *spacemanager);
/* Return pointer to rTab space of the given <spacemanager>. */
void*         gt_linspace_management_get_rTabspace(const GtLinspaceManagement
                                                   *spacemanager);
/* Return pointer to crosspointTab space of the given <spacemanager>. */
void*         gt_linspace_management_get_crosspointTabspace(
                                                      const GtLinspaceManagement
                                                      *spacemanager);
/* Return pointer to Gtmaxcoordvalue space of the given <spacemanager>. */
void*         gt_linspace_management_get_maxspace(const GtLinspaceManagement
                                                  *spacemanager);

/* Change allocated linear space of the of given <spacemanager> in 2dim matrix
   of size (<ulen>+1)*(<vlen>+1) and return pointer to these space. */
GtUword**     gt_linspace_management_change_to_square(GtLinspaceManagement
                                                      *spacemanager,
                                                      GtUword ulen,
                                                      GtUword vlen);
/* Return size of valueTab space of the given <spacemanager>. */
size_t        gt_linspace_management_get_valueTabsize(const GtLinspaceManagement
                                                      *spacemanager);
/* Set <timesquarefactor> for the given <spacemanager>. */
void          gt_linspace_management_set_TSfactor(GtLinspaceManagement
                                                  *spacemanager,
                                                  GtUword timesquarefactor);

#define add_safe(val1, val2, exception) (((val1) != (exception))\
                                           ? (val1) + (val2)\
                                           : (exception))

#define add_safe_max(val1, val2) add_safe(val1,val2,GT_WORD_MAX)

#define add_safe_min(val1, val2) add_safe(val1,val2,GT_WORD_MIN)

#define add_safe_umax(val1, val2) add_safe(val1,val2,GT_UWORD_MAX)

#endif
