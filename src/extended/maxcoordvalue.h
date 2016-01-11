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
#ifndef MAXCOORDVALUE_H
#define MAXCOORDVALUE_H
#include "core/types_api.h"

/* The <GtMaxcoordvalue> interface. */

/* The following type is to manage maximal score value and (i,j) coordinates
   of this value in DP-matrix. */
typedef struct GtMaxcoordvalue GtMaxcoordvalue;
/* Return a <GtMaxcoordvalue> object. */
GtMaxcoordvalue* gt_maxcoordvalue_new(void);
/* Delete the given <max> object. */
void             gt_maxcoordvalue_delete(GtMaxcoordvalue *max);
/* Return the stored maximal value of the given <max>. */
GtWord           gt_maxcoordvalue_get_value(const GtMaxcoordvalue *max);
/* Set start coordinates <starta> and <startb> for the given <max>. */
void             gt_maxcoordvalue_set_start(GtMaxcoordvalue *max,
                                            GtUword starta,
                                            GtUword startb);
/* Return the start coordinates for the given <max>. */
GtUwordPair      gt_maxcoordvalue_get_start(const GtMaxcoordvalue *max);
/* Set end coordinates <end> for the given <max>. */
void             gt_maxcoordvalue_set_end_with_pair(GtMaxcoordvalue *max,
                                                    GtUwordPair end);
/* Return the end coordinates for the given <max>. */
GtUwordPair      gt_maxcoordvalue_get_end(const GtMaxcoordvalue *max);
/* Set a new maximum <value> and update coordinates for the given <max> by
   calling set-functions for start coordinates <start> and end coordinates
   <enda> and <endb>. */
void             gt_maxcoordvalue_coord_update(GtMaxcoordvalue *max,
                                               GtWord value,
                                               GtUwordPair start,
                                               GtUword enda,
                                               GtUword endb);
/* Set a new maximum <value> and update coordinates for the given <max> by
   calling set-functions for end coordinates <enda> and <endb>. */
void             gt_maxcoordvalue_coord_update_without_start(
                                                           GtMaxcoordvalue *max,
                                                           GtWord value,
                                                           GtUword enda,
                                                           GtUword endb);
/* Return the difference between row indices of start and end coordinates of the
   given <max>. */
GtUword          gt_maxcoordvalue_get_row_length(const GtMaxcoordvalue *max);
/* Return the difference between column indices of start and end coordinates of
   the given <max>. */
GtUword          gt_maxcoordvalue_get_col_length(const GtMaxcoordvalue *max);
/* Check coordinate values of the given <max>. Return false if start coordinates
 * are equal to end coordinates, otherwise false.*/
bool             gt_maxcoordvalue_get_length_safe(const GtMaxcoordvalue *max);
/* Reset the given <max>. */
void             gt_maxcoordvalue_reset(GtMaxcoordvalue *max);
#endif
