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
#include "substriter.h"
#include "spacedef.h"

#include "initbasepower.pr"

struct Substriter
{
  uint64_t unitnum;
  SeqIterator *seqit;
  unsigned int qvalue, numofchars, *mappower;
  bool newseq;
};

Substriter *substriter_new(const StrArray *queryfilenames,
                           const Alphabet *alphabet,
                           unsigned int qvalue)
{
  Substriter *substriter;
  ALLOCASSIGNSPACE(substriter,NULL,Substriter,1);
  substriter->unitnum = 0;
  substriter->seqit
    = seqiterator_new(queryfilenames,getsymbolmapAlphabet(alphabet),true);
  substriter->mappower = initmappower(getnumofcharsAlphabet(alphabet),qvalue);
  substriter->newseq = true;
  return substriter;
}

static unsigned int qgram2code(unsigned int *code,unsigned int numofchars,
                               unsigned int qvalue,const Uchar *qgram)
{
  unsigned int i, tmpcode;
  Uchar a;

  a = qgram[0];
  if (ISSPECIAL(a))
  {
    return 0;
  }
  tmpcode = (unsigned int) a;
  for (i=1U; i < qvalue; i++)
  {
    a = qgram[i];
    if (ISSPECIAL(a))
    {
      return i;
    }
    tmpcode *= numofchars;
    tmpcode += (unsigned int) a;
  }
  *code = tmpcode;
  return qvalue;
}

int substriter_next(Substring *substring,Substriter *substriter,Error *err)
{
  while (true)
  {
    if (substriter->newseq)
    {
      int retval;
      unsigned int firstspecial;

      retval = seqiterator_next(substriter->seqit,
                                &substring->querystart,
                                &substring->remaining,
                                &substring->desc,
                                err);
      if (retval <= 0)
      {
        FREESPACE(substring->desc);
        return retval;
      }
      substring->queryptr = substring->querystart;
      assert(substring->remaining > 0);
      substriter->newseq = false;
      if (substring->remaining >= (unsigned long) substriter->qvalue)
      {
        firstspecial 
          = qgram2code(&substring->currentcode,substriter->numofchars,
                       substriter->qvalue,substring->querystart);
      }
      break;
    }
    assert(substring->remaining > 0);
    substring->remaining--;
    substring->queryptr++;
    if (substring->remaining > 0)
    {
      break;
    }
    substriter->newseq = true;
    substriter->unitnum++;
    FREESPACE(substring->desc);
  }
  return 1;
}

uint64_t substriter_unitnum(const Substriter *substriter)
{
  return substriter->unitnum;
}

void substriter_delete(Substriter **substriter)
{
  seqiterator_delete((*substriter)->seqit);
  FREESPACE((*substriter)->mappower);
  FREESPACE(*substriter);
}
