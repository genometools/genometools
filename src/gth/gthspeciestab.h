/*
  Copyright (c) 2003-2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2003-2007 Center for Bioinformatics, University of Hamburg

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

#ifndef GTHSPECIESTAB_H
#define GTHSPECIESTAB_H
#include "core/unused_api.h"

#define NUMOFSPECIES (sizeof (speciestab)/sizeof (speciestab[0]))

/*
  The following defines the table of species. For each species the
  corresponding index in the table is the species number in the GthCallInfo
  structure.
*/

GT_UNUSED static char *speciestab[] =
{
  "human",         /*  0 */
  "mouse",         /*  1 */
  "rat",           /*  2 */
  "chicken",       /*  3 */
  "drosophila",    /*  4 */
  "nematode",      /*  5 */
  "fission_yeast", /*  6 */
  "aspergillus",   /*  7 */
  "arabidopsis",   /*  8 (old) */
  "maize",         /*  9 (old) last species stored in header file */
  "rice",
  "medicago"
};

#endif
