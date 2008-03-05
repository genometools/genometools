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
#include <assert.h>
#include "spacedef.h"
#include "intcode-def.h"
#include "encseq-def.h"
#include "sfx-enumcodes.h"

#include "initbasepower.pr"

struct Enumcodeatposition
{
  Sequencerange previousrange;
  Specialrangeiterator *sri;
  bool moveforward;
  Seqpos totallength;
  bool exhausted;
  const Encodedsequence *encseq;
  Readmode readmode;
  unsigned int prefixlength;
  Codetype **multimappower, *filltable;
};

Enumcodeatposition *newEnumcodeatposition(const Encodedsequence *encseq,
                                          Readmode readmode,
                                          unsigned int prefixlength,
                                          unsigned int numofchars)
{
  Enumcodeatposition *ecp;

  ALLOCASSIGNSPACE(ecp,NULL,Enumcodeatposition,1);
  ecp->encseq = encseq;
  ecp->readmode = readmode;
  ecp->multimappower = initmultimappower(numofchars,prefixlength);
  ecp->filltable = initfilltable(numofchars,prefixlength);
  ecp->prefixlength = prefixlength;
  ecp->moveforward = ISDIRREVERSE(readmode) ? true : false;
  ecp->totallength = getencseqtotallength(encseq);
  if (ecp->moveforward)
  {
    ecp->previousrange.leftpos = ecp->previousrange.rightpos = 0;
  } else
  {
    ecp->previousrange.leftpos = ecp->previousrange.rightpos = ecp->totallength;
  }
  ecp->exhausted = false;
  if (hasspecialranges(encseq))
  {
    ecp->sri = newspecialrangeiterator(encseq,ecp->moveforward);
  } else
  {
    ecp->sri = NULL;
  }
  return ecp;
}

static bool newcodelistelem(Specialcontext *specialcontext,
                            Seqpos smallerval,
                            Seqpos largerval,
                            const Enumcodeatposition *ecp)
{
  if (smallerval < largerval)
  {
    Seqpos distance = largerval - smallerval;

    if (distance > (Seqpos) (ecp->prefixlength-1))
    {
      distance = (Seqpos) (ecp->prefixlength-1);
    }
    specialcontext->maxprefixindex = (unsigned int) distance;
    if (ecp->moveforward)
    {
      specialcontext->position = ecp->totallength - smallerval;
    } else
    {
      specialcontext->position = largerval;
    }
    assert(specialcontext->position >= (Seqpos) specialcontext->maxprefixindex);
    return true;
  }
  return false;
}

bool nextEnumcodeatposition(Specialcontext *specialcontext,
                            Enumcodeatposition *ecp)
{
  Sequencerange currentrange;
  bool done;

  if (ecp->exhausted)
  {
    return false;
  }
  while (ecp->sri != NULL)
  {
    if (!nextspecialrangeiterator(&currentrange,ecp->sri))
    {
      freespecialrangeiterator(&ecp->sri);
      ecp->sri = NULL;
      break;
    }
    if (ecp->moveforward)
    {
      if (newcodelistelem(specialcontext,
                          ecp->previousrange.rightpos,
                          currentrange.leftpos,
                          ecp))
      {
        ecp->previousrange = currentrange;
        return true;
      }
    } else
    {
      if (newcodelistelem(specialcontext,
                          currentrange.rightpos,
                          ecp->previousrange.leftpos,
                          ecp))
      {
        ecp->previousrange = currentrange;
        return true;
      }
    }
    ecp->previousrange = currentrange;
  }
  ecp->exhausted = true;
  if (ecp->moveforward)
  {
    done = newcodelistelem(specialcontext,
                           ecp->previousrange.rightpos,
                           ecp->totallength,
                           ecp);
  } else
  {
    done = newcodelistelem(specialcontext,
                           0,
                           ecp->previousrange.leftpos,
                           ecp);
  }
  return done;
}

void freeEnumcodeatposition(Enumcodeatposition **ecp)
{
  FREESPACE((*ecp)->filltable);
  multimappowerfree(&(*ecp)->multimappower);
  FREESPACE(*ecp);
}

Codetype computefilledqgramcode(const Enumcodeatposition *ecp,
                                unsigned int prefixindex,
                                Seqpos pos)
{
  Codetype code;
  unsigned int idx;
  Uchar cc;

  assert(prefixindex > 0 && prefixindex < ecp->prefixlength);
  code = ecp->filltable[prefixindex];
  for (idx=0; idx<prefixindex; idx++)
  {
    assert((Seqpos) (pos + idx) < ecp->totallength);
    cc = getencodedcharnospecial(ecp->encseq,pos + idx, ecp->readmode);
    assert(ISNOTSPECIAL(cc));
    code += ecp->multimappower[idx][cc];
  }
  return code;
}

bool computefilledqgramcodestopatmax(Codetype *code,
                                     const Enumcodeatposition *ecp,
                                     unsigned int prefixindex,
                                     Seqpos pos,
                                     Codetype stopcode)
{
  Codetype tmpcode;
  unsigned int idx;
  Uchar cc;

  assert(prefixindex > 0 && prefixindex < ecp->prefixlength);
  tmpcode = ecp->filltable[prefixindex];
  if (tmpcode > stopcode)
  {
    return false;
  }
  for (idx=0; idx<prefixindex; idx++)
  {
    assert((Seqpos) (pos + idx) < ecp->totallength);
    cc = getencodedcharnospecial(ecp->encseq,pos + idx, ecp->readmode);
    assert(ISNOTSPECIAL(cc));
    tmpcode += ecp->multimappower[idx][cc];
    if (tmpcode > stopcode)
    {
      return false;
    }
  }
  *code = tmpcode;
  return true;
}
