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
#include "libgtcore/error.h"
#include "libgtcore/minmax.h"
#include "intcode-def.h"
#include "seqpos-def.h"
#include "verbose-def.h"

#include "bcktab.h"
#include "initbasepower.pr"

/*
  We need \texttt{prefixlenbits} bits to store the length of
  a matching prefix. So we can store the following maximal value
  in the remaining bits.
*/

#define MAXREMAININGAFTERPREFIXLEN(PFXLENBITS)\
        ((((Codetype) 1) << (32-(PFXLENBITS))) - 1)

/*
  We allow to choose the prefixlength \(l\) in such a way that the size of
  table bcktab (which is $8\cdot k^{l}$ never exceeds the
  $\texttt{MAXMULTIPLIEROFTOTALLENGTH} \cdot n$, where \(k\) is the size
  of the alphabet and \(n\) is the total length of the input sequence.
*/

#define MAXMULTIPLIEROFTOTALLENGTH 4

#define MAXVALUEWITHBITS(BITNUM)    ((1U << (BITNUM)) - 1)

static unsigned int prefixlengthwithmaxspace(unsigned int numofchars,
                                             Seqpos maxbytes,
                                             unsigned int factor)
{
  unsigned int prefixlength;
  unsigned long sizeofrep;

  for (prefixlength = 1U; /* Nothing */; prefixlength++)
  {
    sizeofrep = sizeofbuckettable(numofchars,prefixlength);
    if ((Seqpos) (sizeofrep / factor) > maxbytes)
    {
      return prefixlength-1;
    }
  }
  /*@ignore@*/
  return 1U;
  /*@end@*/
}

unsigned int recommendedprefixlength(unsigned int numofchars,
                                     Seqpos totallength)
{
  unsigned int prefixlength;

  prefixlength = prefixlengthwithmaxspace(numofchars,totallength,1U);
  if (prefixlength == 0)
  {
    return 1U;
  } else
  {
    unsigned int mbp = maxbasepower(numofchars);
    if (mbp >= 1U)
    {
      return MIN(mbp,prefixlength);
    } else
    {
      return prefixlength;
    }
  }
}

unsigned int whatisthemaximalprefixlength(unsigned int numofchars,
                                          Seqpos totallength,
                                          unsigned int prefixlenbits)
{
  unsigned int maxprefixlen, mbp;

  maxprefixlen = prefixlengthwithmaxspace(numofchars,totallength,
                                          (unsigned int)
                                          MAXMULTIPLIEROFTOTALLENGTH);
  mbp = maxbasepower(numofchars);
  maxprefixlen = MIN(mbp,maxprefixlen);
  if (prefixlenbits > 0)
  {
    unsigned int tmplength;
    tmplength
      = prefixlengthwithmaxspace(numofchars,
                                 (Seqpos)
                                 MAXREMAININGAFTERPREFIXLEN(prefixlenbits),
                                 1U);
    if (tmplength > 0 && maxprefixlen > tmplength)
    {
      maxprefixlen = tmplength;
    }
    tmplength = MAXVALUEWITHBITS(prefixlenbits);
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

int checkprefixlength(unsigned int maxprefixlen,
                      unsigned int prefixlength,
                      Error *err)
{
  error_check(err);
  if (maxprefixlen < prefixlength)
  {
    error_set(err,"prefix length %u is too large, maximal prefix length "
                      "for this input size and alphabet size is %u",
                      prefixlength,
                      maxprefixlen);
    return -1;
  }
  return 0;
}

void showmaximalprefixlength(Verboseinfo *verboseinfo,
                             unsigned int maxprefixlen,
                             unsigned int recommended)
{
  /* XXX: use one call of showverbose to display this */
  showverbose(verboseinfo,
              "for this input size and alphabet size, "
              "the maximal prefixlength");
  showverbose(verboseinfo,"(argument of option -pl) is %u,",maxprefixlen);
  showverbose(verboseinfo,"the recommended prefixlength is %u",recommended);
}
