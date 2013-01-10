/*
  Copyright (c) 2012-2013 Dirk Willrodt <willrodt@zbh.uni-hamburg.de>
  Copyright (c) 2012-2013 Center for Bioinformatics, University of Hamburg

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

#ifndef POPCOUNT_TAB_H
#define POPCOUNT_TAB_H

#include "core/error_api.h"

/* The <GtPopcountTab> class represents a table of values of a fixed bit
   width sorted by their popcount (number of bits set to 1). Values are sorted
   by increasing value within one popcount class (set of values with the
   same popcount). */
typedef struct GtPopcountTab GtPopcountTab;

/* Returns <GtPopcountTab> object with tables for unsigned values of
   <blocksize> bit width. */
GtPopcountTab* gt_popcount_tab_new(unsigned int blocksize);

/* Returns the class (i.e. the popcount) of a given <block>. Value of <block>
   must be smaller than 2^<blocksize> */
unsigned int   gt_popcount_tab_class(unsigned long block,
                                     unsigned int blocksize);

/* Returns the bits needed to store the offsets for a given <class> (i.e.
   popcount) and <blocksize>. <class> <= <blocksize> and might overflow for
   large values, (<blocksize> choose <class>) */
unsigned int   gt_popcount_tab_offset_bits(unsigned int blocksize,
                                           unsigned int class);

/* Return the <i>-th value from <popcount_tab> with <popcount_c> bits set, <i>
   has to be in range 0..(blocksize choose popcount). */
unsigned long  gt_popcount_tab_get(GtPopcountTab *popcount_tab,
                                   unsigned int popcount_c,
                                   unsigned long i);

/* Return rank of 1s in <i>-th block given for <popcount_c> bits set
   up to and including <pos>, which is a relative bit position within that
   block. Note that <pos> <= blocksize of <popcount_tab>. */
unsigned int   gt_popcount_tab_rank_1(GtPopcountTab *popcount_tab,
                                      unsigned int popcount_c,
                                      unsigned long i,
                                      unsigned int pos);

/* Return rank of 0s in <i>-th block given for <popcount_c> bits set
   up to and including <pos>, which is a relative bit position within that
   block.  Note that <pos> <= blocksize of <popcount_tab>. */
unsigned int   gt_popcount_tab_rank_0(GtPopcountTab *popcount_tab,
                                      unsigned int popcount_c,
                                      unsigned long i,
                                      unsigned int pos);

/* Return the offset of a <block> within its class. Note that
   <block> < 2^blocksize of <popcount_tab> */
unsigned long  gt_popcount_tab_get_offset_for_block(GtPopcountTab *popcount_tab,
                                                    unsigned long block);

/* Return size of a <GtPopcountTab> with <blocksize> in bytes. */
size_t         gt_popcount_tab_calculate_size(unsigned int blocksize);

/* Deletes <popcount_tab> and frees all associated memory. */
void           gt_popcount_tab_delete(GtPopcountTab *popcount_tab);

int            gt_popcount_tab_unit_test(GtError *err);

#endif
