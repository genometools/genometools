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

#include <inttypes.h>
#include "libgtcore/strarray.h"
#include "libgtcore/seqiterator.h"
#include "libgtcore/chardef.h"
#include "spacedef.h"
#include "intcode-def.h"
#include "substriter.h"
#include "qgram2code.h"

#include "initbasepower.pr"

struct Substriter
{
  const Uchar *currentptr, *start;
  unsigned long remaining;
  Codetype currentcode;
  unsigned int qvalue,
               numofchars;
  Codetype **multimappower;
};

Substriter *substriter_new(const Alphabet *alphabet,unsigned int qvalue)
{
  Substriter *substriter;
  ALLOCASSIGNSPACE(substriter,NULL,Substriter,1);
  substriter->qvalue = qvalue;
  substriter->numofchars = getnumofcharsAlphabet(alphabet);
  substriter->multimappower = initmultimappower(substriter->numofchars,qvalue);
  return substriter;
}

void substriter_init(Substriter *substriter,const Uchar *start,
                    unsigned long len)
{
  substriter->start = substriter->currentptr = start;
  substriter->remaining = len;
  assert(substriter->remaining > 0);
}

int substriter_next(Substriter *substriter)
{
  unsigned int firstspecial;

  while (true)
  {
    if (substriter->remaining >= (unsigned long) substriter->qvalue)
    {
      firstspecial = qgram2code(&substriter->currentcode,
                                (const Codetype **) substriter->multimappower,
                                substriter->qvalue,
                                substriter->currentptr);
      if (firstspecial == substriter->qvalue)
      {
        substriter->remaining--;
        substriter->currentptr++;
        break;
      }
      substriter->remaining -= firstspecial + 1;
      substriter->currentptr += firstspecial + 1;
    } else
    {
      return 0;
    }
  }
  return 1;
}

void substriter_delete(Substriter **substriter)
{
  multimappowerfree(&(*substriter)->multimappower);
  FREESPACE(*substriter);
}
