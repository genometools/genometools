/*
  Copyright (c) 2009 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
  Copyright (c) 2009 Center for Bioinformatics, University of Hamburg

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

#ifndef IDXLOCALISW_H
#define IDXLOCALISW_H
#include <inttypes.h>
#include "encodedsequence.h"
#include "procmatch.h"

typedef long Scoretype;

typedef struct
{
  Scoretype matchscore,   /* must be positive */
            mismatchscore,/* must be negative */
            gapstart,     /* must be negative */
            gapextend;    /* must be negative */
} Scorevalues;

#define REPLACEMENTSCORE(SV,A,B) (((A) != (B) || ISSPECIAL(A))\
                                   ? (SV)->mismatchscore\
                                   : (SV)->matchscore)

typedef struct SWdpresource SWdpresource;

SWdpresource *newSWdpresource(Scoretype matchscore,
                              Scoretype mismatchscore,
                              Scoretype gapextend,
                              unsigned long scorethreshold,
                              bool showalignment,
                              Processmatch processmatch,
                              void *processmatchinfo);

void multiapplysmithwaterman(SWdpresource *dpresource,
                             const GtEncodedsequence *encseq,
                             const GtUchar *query,
                             unsigned long querylen);

void freeSWdpresource(SWdpresource *swdpresource);

#endif
