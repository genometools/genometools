/*
  Copyright (c) 2004-2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2004-2007 Center for Bioinformatics, University of Hamburg

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

#ifndef GTHERROR_H
#define GTHERROR_H

/* This file contains all error return values used in GenomeThreader.
   If something unspecific goes wrong in a function, usually -1 is returned.
   For specific events, which should cause specific reactions in the calling
   functions the following codes are used. Error values have to be negative. */

#define GTH_ERROR_CUTOUT_NOT_IN_INTRON            (-10)
#define GTH_ERROR_DP_PARAMETER_ALLOCATION_FAILED  (-11)
#define GTH_ERROR_MATRIX_ALLOCATION_FAILED        (-12)
#define GTH_ERROR_SA_COULD_NOT_BE_DETERMINED      (-13)

#endif
