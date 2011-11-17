/*
  Copyright (c) 2003-2009 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2003-2008 Center for Bioinformatics, University of Hamburg

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

#ifndef ALIGN_COMMON_H
#define ALIGN_COMMON_H

/* This modules contains the common functions which are used in the dna and the
   protein DP. */

#include "gth/bssm_param.h"
#include "gth/dp_param.h"
#include "gth/stat.h"

#define SCORE(T,N,M)       dpm->core.score[T][N][M]
#define GTH_MINUSINFINITY  ((GthFlt) -99999.0)
#define DASH               (SEPARATOR-3)
#define UNSET              (SEPARATOR-4)

/*
  The following updates maxvalue if it is smaller than value.
  The variable <retrace> is defined accordingly.
*/

#define UPDATEMAX(C)\
        if (maxvalue < value)\
        {\
          maxvalue = value;\
          retrace = (GthPath) C;\
        }

typedef unsigned char GthPath;

#endif
