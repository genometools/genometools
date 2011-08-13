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

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <inttypes.h>
#include "core/error.h"
#include "core/minmax.h"
#include "core/codetype.h"
#include "core/logger.h"
#include "bcktab.h"
#include "sfx-apfxlen.h"
#include "initbasepower.h"

/*
  We need \texttt{prefixlenbits} bits to store the length of
  a matching prefix. So we can store the following maximal value
  in the remaining bits.
*/

#define MAXREMAININGAFTERPREFIXLEN(PFXLENBITS)\
        ((((GtCodetype) 1) << (32-(PFXLENBITS))) - 1)

/*
  We allow to choose the prefixlength \(l\) in such a way that the size of
  table bcktab never exceeds the
  $\texttt{GT_MAXMULTIPLIEROFTOTALLENGTH}\cdot n$,
  where \(n\) is the total length of the input sequence.
*/

#define GT_MAXMULTIPLIEROFTOTALLENGTH   4.0
#define GT_MAXVALUEWITHBITS(BITNUM)     ((1U << (BITNUM)) - 1)

static unsigned int prefixlengthwithmaxspace(unsigned int numofchars,
                                             unsigned long maxbytes,
                                             double factor,
                                             unsigned long maxvalue,
                                             bool withspecialsuffixes)
{
  unsigned int prefixlength;
  uint64_t sizeofrep;

#ifdef WITHINFO
  printf("maxbytes = %lu\n",(unsigned long) (maxbytes * factor));
#endif
  for (prefixlength = 1U; /* Nothing */; prefixlength++)
  {
    sizeofrep = gt_bcktab_sizeoftable(numofchars,prefixlength,
                                      maxvalue,withspecialsuffixes);
#ifdef WITHINFO
    printf("sizeofrep = %lu, after divide %lu\n",(unsigned long) sizeofrep,
                                        (unsigned long) (sizeofrep/factor));
#endif
    if (sizeofrep/factor > (uint64_t) maxbytes)
    {
#ifdef WITHINFO
      printf("prefixlengthwithmaxspace = %u\n",prefixlength-1);
#endif
      return prefixlength-1;
    }
  }
  /*@ignore@*/
  return 1U;
  /*@end@*/
}

unsigned int gt_recommendedprefixlength(unsigned int numofchars,
                                        unsigned long totallength,
                                        double recommendedmultiplier,
                                        bool withspecialsuffixes)
{
  unsigned int prefixlength;

  prefixlength = prefixlengthwithmaxspace(numofchars,totallength,
                                          recommendedmultiplier,
                                          totallength+1,
                                          withspecialsuffixes);
  if (prefixlength == 0)
  {
    return 1U;
  } else
  {
    unsigned int mbp = gt_maxbasepower(numofchars);
    if (mbp >= 1U)
    {
      return MIN(mbp,prefixlength);
    } else
    {
      return prefixlength;
    }
  }
}

unsigned int gt_whatisthemaximalprefixlength(unsigned int numofchars,
                                             unsigned long totallength,
                                             unsigned int prefixlenbits,
                                             bool withspecialsuffixes)
{
  unsigned int maxprefixlen, mbp;

  maxprefixlen = prefixlengthwithmaxspace(numofchars,totallength,
                                          GT_MAXMULTIPLIEROFTOTALLENGTH,
                                          totallength+1,
                                          withspecialsuffixes);
  mbp = gt_maxbasepower(numofchars);
  maxprefixlen = MIN(mbp,maxprefixlen);
  if (prefixlenbits > 0)
  {
    unsigned int tmplength;
    tmplength
      = prefixlengthwithmaxspace(numofchars,
                                 (unsigned long)
                                 MAXREMAININGAFTERPREFIXLEN(prefixlenbits),
                                 GT_RECOMMENDED_MULTIPLIER_DEFAULT,
                                 totallength+1,
                                 withspecialsuffixes);
    if (tmplength > 0 && maxprefixlen > tmplength)
    {
      maxprefixlen = tmplength;
    }
    tmplength = GT_MAXVALUEWITHBITS(prefixlenbits);
    if (tmplength > 0 && maxprefixlen > tmplength)
    {
      maxprefixlen = tmplength;
    }
  }
  if (maxprefixlen == 0)
  {
    return 1U;
  }
  return maxprefixlen;
}

int gt_checkprefixlength(unsigned int maxprefixlen,
                      unsigned int prefixlength,
                      GtError *err)
{
  gt_error_check(err);
  if (maxprefixlen < prefixlength)
  {
    gt_error_set(err,"prefix length %u is too large, maximal prefix length "
                      "for this input size and alphabet size is %u",
                      prefixlength,
                      maxprefixlen);
    return -1;
  }
  return 0;
}

void gt_showmaximalprefixlength(GtLogger *logger,
                             unsigned int maxprefixlen,
                             unsigned int recommended)
{
  gt_logger_log(logger,
              "for this input size and alphabet size, "
              "the maximal prefixlength");
  gt_logger_log(logger,"(argument of option -pl) is %u,",maxprefixlen);
  gt_logger_log(logger,"the recommended prefixlength is %u",recommended);
}
