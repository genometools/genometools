/*
  Copyright (c) 2015 Annika <annika.seidel@studium.uni-hamburg.de>
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

#ifndef LINSPACEMANAGEMENT_H
#define LINSPACEMANAGEMENT_H
#include "core/types_api.h"

/* file with all useful uitlities for different algorithms,
 * which work in linear space */

/* struct to organize all space allocated */
typedef struct LinspaceManagement LinspaceManagement;

LinspaceManagement* gt_linspaceManagement_new();

void gt_linspaceManagement_delete(LinspaceManagement *spacemanager);

size_t gt_linspaceManagement_get_spacepeak(const LinspaceManagement
                                                                 *spacemanager);

/* checks if enough space is allocated in global case and resize if necessary */
void gt_linspaceManagement_check(LinspaceManagement *spacemanager,
                                 GtUword ulen, GtUword vlen,
                                 size_t valuesize,
                                 size_t rtabsize,
                                 size_t crosspointsize);

/* checks if enough space is allocated in local case and resize if necessary */
void  gt_linspaceManagement_check_local(LinspaceManagement *spacemanager,
                                        GtUword ulen, GtUword vlen,
                                        size_t valuesize,
                                        size_t rstabsize);

/* checks if enough space to use square space functions is allocated
 * in global case*/
bool gt_linspaceManagement_checksquare(LinspaceManagement *spacemanager,
                                       GtUword ulen, GtUword vlen,
                                       size_t valuesize,
                                       size_t rsize);

/* checks if enough space to use square space functions is allocated
 * in local case*/
bool gt_linspaceManagement_checksquare_local(LinspaceManagement *spacemanager,
                                             GtUword ulen, GtUword vlen,
                                             size_t valuesize,
                                             size_t rsize);

void gt_linspaceManagement_set_ulen(LinspaceManagement *spacemanager,
                                    GtUword ulen);

void *gt_linspaceManagement_get_valueTabspace(const LinspaceManagement
                                                                 *spacemanager);

GtUword gt_linspaceManagement_get_valueTabLen(const LinspaceManagement
                                                                 *spacemanager);

void *gt_linspaceManagement_get_rTabspace(const LinspaceManagement
                                                                 *spacemanager);

void *gt_linspaceManagement_get_crosspointTabspace(const LinspaceManagement
                                                                 *spacemanager);

void *gt_linspaceManagement_get_maxspace(const LinspaceManagement
                                                                 *spacemanager);

/* change allocated linear space in 2dim format */
GtUword **gt_linspaceManagement_change_to_square(LinspaceManagement
                                                 *spacemanager,
                                                 GtUword ulen, GtUword vlen);

void gt_linspaceManagement_set_TSfactor(LinspaceManagement *spacemanager,
                                        GtUword timesquarefactor);

inline GtWord add_safe_max(GtWord val1, GtWord val2);
inline GtWord add_safe_min(GtWord val1, GtWord val2);
inline GtUword add_safe_umax(GtUword val1, GtUword val2);
#endif
