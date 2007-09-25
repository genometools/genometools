/*
  Copyright (c) 2007 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
  Copyright (c) 2007 Center for Bioinformatics, University of Hamburg

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

#include <math.h>
#include "libgtcore/chardef.h"
#include "divmodmul.h"
#include "fmindex.h"

#include "fmi-occ.gen"

unsigned long skfmuniqueforward (const Fmindex *fmindex,
                                 const Uchar *qstart,
                                 const Uchar *qend)
{
  Uchar cc;
  const Uchar *qptr;
  Bwtbound bwtbound;

  qptr = qstart;
  cc = *qptr++;
#undef mydebug
#ifdef mydebug
  printf("# start cc=%u\n",cc);
#endif
  if (ISSPECIAL(cc))
  {
    return 0;
  }
  bwtbound.lbound = fmindex->tfreq[cc];
  bwtbound.ubound = fmindex->tfreq[cc+1];
#ifdef mydebug
  printf("# bounds=" FormatSeqpos "," FormatSeqpos " = " FormatSeqos
          "occurrences\n",
         PRINTSeqposcast(bwtbound.lbound),
         PRINTSeqposcast(bwtbound.ubound),
         PRINTSeqposcast(bwtbound.ubound - bwtbound.lbound));
#endif
  while (qptr < qend && bwtbound.lbound + 1 < bwtbound.ubound)
  {
    cc = *qptr;
#ifdef mydebug
    printf("# cc=%u\n",cc);
#endif
    if (ISSPECIAL (cc))
    {
      return 0;
    }
    bwtbound.lbound = fmindex->tfreq[cc] +
                      fmoccurrence (fmindex, cc, bwtbound.lbound);
    bwtbound.ubound = fmindex->tfreq[cc] +
                      fmoccurrence (fmindex, cc, bwtbound.ubound);
#ifdef mydebug
    printf("# bounds=" FormatSeqpos "," FormatSeqpos " = " FormatSeqos
            "occurrences\n",
           PRINTSeqposcast(bwtbound.lbound),
           PRINTSeqposcast(bwtbound.ubound),
           PRINTSeqposcast(bwtbound.ubound - bwtbound.lbound));
#endif
    qptr++;
  }
  if (bwtbound.lbound + 1 == bwtbound.ubound)
  {
    return (unsigned long) (qptr - qstart);
  }
  return 0;
}
