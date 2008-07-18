/*
  Copyright (c) 2008 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
  Copyright (c) 2008 Center for Bioinformatics, University of Hamburg

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

#include <string.h>
#include "libgtcore/symboldef.h"
#include "libgtcore/unused.h"
#include "libgtcore/ma.h"
#include "seqpos-def.h"
#include "encseq-def.h"
#include "defined-types.h"
#include "initeqsvec.h"

typedef struct
{
  unsigned long prefixofsuffix;
} Parallelmstats;

typedef struct
{
  unsigned long patternlength,
                *eqsvector;
} Matchtaskinfo;

#ifdef SKDEBUG

static void apm_showParallelmstats(const DECLAREPTRDFSSTATE(aliascol),
                                   unsigned long depth,
                                   const void *dfsconstinfo)
{
  const Matchtaskinfo *mti = (const Matchtaskinfo *) dfsconstinfo;
  const Myerscolumn *col = (const Myerscolumn *) aliascol;
  bool first = true;

  unsigned long idx, backmask;

  printf("[");
  for (idx=0, backmask = 1UL; idx<mti->patternlength; idx++, backmask <<= 1)
  {
    if (col->prefixofsuffix & backmask)
    {
      if (first)
      {
        printf("%lu",idx);
        first = false;
      } else
      {
        printf(",%lu",idx);
      }
    }
  }
  printf("\n");
}

#endif
