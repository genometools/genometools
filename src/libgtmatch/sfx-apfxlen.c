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
#include "libgtcore/env.h"
#include "intcode-def.h"
#include "seqpos-def.h"

#define SIZEOFBCKENTRY (2 * sizeof (Seqpos))

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

#define MAXVALUEWITHBITS(BITNUM)    ((((unsigned int) 1) << (BITNUM)) - 1)

static unsigned int logalphasize(unsigned int numofchars,double value)
{
  unsigned int retval;
  double logtmp1, logtmp2;

  if (value <= (double) numofchars)
  {
    return (unsigned int) 1;
  }
  logtmp1 = log(value);
  logtmp2 = log((double) numofchars);
  retval = (unsigned int) floor(logtmp1/logtmp2);
  return retval;
}

unsigned int recommendedprefixlength(unsigned int numofchars,
                                     Seqpos totallength)
{
  unsigned int prefixlength;

  prefixlength = logalphasize(numofchars,
                              (double) totallength/SIZEOFBCKENTRY);
  if (prefixlength == 0)
  {
    return (unsigned int) 1;
  } else
  {
    return prefixlength;
  }
}

unsigned int whatisthemaximalprefixlength(unsigned int numofchars,
                                          Seqpos totallength,
                                          unsigned int prefixlenbits)
{
  unsigned int maxprefixlen;

  maxprefixlen = logalphasize(numofchars,
                           (double) totallength/
                                (SIZEOFBCKENTRY/MAXMULTIPLIEROFTOTALLENGTH));
  if (prefixlenbits > 0)
  {
    unsigned int tmplength;
    tmplength = logalphasize(numofchars,
                             (double)
                             MAXREMAININGAFTERPREFIXLEN(prefixlenbits));
    if (maxprefixlen > tmplength)
    {
      maxprefixlen = tmplength;
    }
    tmplength = MAXVALUEWITHBITS(prefixlenbits);
    if (maxprefixlen > tmplength)
    {
      maxprefixlen = tmplength;
    }
  }
  return maxprefixlen;
}

int checkprefixlength(unsigned int maxprefixlen,
                      unsigned int prefixlength,Env *env)
{
  env_error_check(env);
  if (maxprefixlen < prefixlength)
  {
    env_error_set(env,"prefix length %u is too large, maximal prefix length "
                      "for this input size and alphabet size is %u",
                      prefixlength,
                      maxprefixlen);
    return -1;
  }
  return 0;
}

void showmaximalprefixlength(unsigned int maxprefixlen,
                             unsigned int recommended)
{
  printf("# for this input size and alphabet size, the maximal prefixlength\n"
         "# (argument of option -pl) is %u,\n"
         "# the recommended prefixlength is %u\n",
         maxprefixlen,
         recommended);
}
