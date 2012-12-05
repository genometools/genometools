/*
  Copyright (c) 2012 Dirk Willrodt <willrodt@zbh.uni-hamburg.de>
  Copyright (c) 2012 Center for Bioinformatics, University of Hamburg

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

/* Class <GtPopcountTab> stores a table of unsigned values of fixed bit width
   sorted by popcount value. Values are sorted by increasing value within one
   popcount class. */
typedef struct GtPopcountTab GtPopcountTab;

/* Returns <GtPopcountTab> object with tables for unsigned values of <blocksize>
   bit width */
GtPopcountTab *gt_popcount_tab_new(unsigned char blocksize) ;

/* Deletes <popcount_tab> and frees all associated memory */
void gt_popcount_tab_delete(GtPopcountTab *popcount_tab);

/* Return unsigned long with <popcount> bits, there might be more than one value
   with that <popcount>, <offset> defines which to return, <offset> = 0 beeing
   the smallest value. <offset> has to be in range
   0..(blocksize choose popcount) */
unsigned long gt_popcount_tab_get(GtPopcountTab *popcount_tab,
                                  unsigned char popcount,
                                  unsigned long offset);

int gt_popcount_tab_unit_test(GtError *err);
#endif
