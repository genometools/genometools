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

#include "core/chardef.h"
#include "core/minmax.h"
#include "encseq-def.h"
#include "lcpoverflow.h"

int comparetwosuffixes(const Encodedsequence *encseq,
                       Readmode readmode,
                       Seqpos *maxlcp,
                       bool specialsareequal,
                       bool specialsareequalatdepth0,
                       Seqpos maxdepth,
                       Seqpos start1,
                       Seqpos start2,
                       Encodedsequencescanstate *esr1,
                       Encodedsequencescanstate *esr2)
{
  GtUchar cc1, cc2;
  Seqpos pos1, pos2, end1, end2;
  int retval;

  end1 = end2 = getencseqtotallength(encseq);
  if (maxdepth > 0)
  {
    if (end1 > start1 + maxdepth)
    {
      end1 = start1 + maxdepth;
    }
    if (end2 > start2 + maxdepth)
    {
      end2 = start2 + maxdepth;
    }
  }
  if (esr1 != NULL && esr2 != NULL)
  {
    initEncodedsequencescanstate(esr1,encseq,readmode,start1);
    initEncodedsequencescanstate(esr2,encseq,readmode,start2);
  } else
  {
    gt_assert(esr1 == NULL && esr2 == NULL);
  }
  for (pos1=start1, pos2=start2; /* Nothing */; pos1++, pos2++)
  {
    if (pos1 >= end1 || pos2 >= end2)
    {
      *maxlcp = pos1 - start1;
      retval = 0;
      break;
    }
    if (esr1 != NULL)
    {
      cc1 = sequentialgetencodedchar(encseq,esr1,pos1,readmode);
      CHECKENCCHAR(cc1,encseq,pos1,readmode);
    } else
    {
      cc1 = getencodedchar(encseq,pos1,readmode);
    }
    if (esr2 != NULL)
    {
      cc2 = sequentialgetencodedchar(encseq,esr2,pos2,readmode);
      CHECKENCCHAR(cc2,encseq,pos2,readmode);
    } else
    {
      cc2 = getencodedchar(encseq,pos2,readmode);
    }
    if (ISSPECIAL(cc1))
    {
      if (ISSPECIAL(cc2))
      {
        if (specialsareequal || (pos1 == start1 && specialsareequalatdepth0))
        {
          *maxlcp = pos1 - start1 + 1;
          retval = 0;
          break;
        }
        if (pos1 < pos2)
        {
          *maxlcp = pos1  - start1;
          retval = -1; /* a < b */
          break;
        }
        if (pos1 > pos2)
        {
          *maxlcp = pos1 - start1;
          retval = 1; /* a > b */
          break;
        }
        *maxlcp = pos1 - start1 + 1;
        retval = 0; /* a = b */
        break;
      }
      *maxlcp = pos1 - start1;
      retval = 1; /* a > b */
      break;
    } else
    {
      if (ISSPECIAL(cc2))
      {
        *maxlcp = pos1 - start1;
        retval = -1; /* a < b */
        break;
      }
      if (cc1 < cc2)
      {
        *maxlcp = pos1 - start1;
        retval = -1; /* a < b */
        break;
      }
      if (cc1 > cc2)
      {
        *maxlcp = pos1 - start1;
        retval = 1; /* a > b */
        break;
      }
    }
  }
  return retval;
}

static Seqpos derefcharboundaries(const Encodedsequence *encseq,
                                  bool fwd,
                                  bool complement,
                                  Seqpos start,
                                  Seqpos maxoffset,
                                  Seqpos currentoffset,
                                  Seqpos totallength)
{
  if (fwd)
  {
    if (start + currentoffset == totallength)
    {
      return totallength + GT_COMPAREOFFSET;
    }
    start += currentoffset;
  } else
  {
    if (start < currentoffset)
    {
      return currentoffset - start + (Seqpos) GT_COMPAREOFFSET;
    }
    start -= currentoffset;
  }
  if (currentoffset <= maxoffset)
  {
    GtUchar cc;
    cc = getencodedchar(encseq,start,Forwardmode);
    if (ISSPECIAL(cc))
    {
      return start + GT_COMPAREOFFSET;
    }
    if (complement)
    {
      cc = COMPLEMENTBASE(cc);
    }
    return cc;
  }
  return  start + GT_COMPAREOFFSET;
}

int comparetwostrings(const Encodedsequence *encseq,
                      bool fwd,
                      bool complement,
                      Seqpos *maxcommon,
                      Seqpos pos1,
                      Seqpos pos2,
                      Seqpos maxdepth)
{
  Seqpos currentoffset, maxoffset, cc1, cc2,
         totallength = getencseqtotallength(encseq);

  if (fwd)
  {
    gt_assert(pos1 < totallength);
    gt_assert(pos2 < totallength);
    maxoffset = MIN(totallength - pos1,totallength - pos2);
  } else
  {
    maxoffset = MIN(pos1+1,pos2+1);
  }
  if (*maxcommon > 0)
  {
    maxoffset = MIN(*maxcommon,maxoffset);
  }
  if (maxdepth > 0)
  {
    maxoffset = MIN(maxoffset,maxdepth);
  }
  for (currentoffset = 0; currentoffset <= maxoffset; currentoffset++)
  {
    cc1 = derefcharboundaries(encseq,fwd,complement,
                              pos1,maxoffset,currentoffset,totallength);
    cc2 = derefcharboundaries(encseq,fwd,complement,
                              pos2,maxoffset,currentoffset,totallength);
    *maxcommon = currentoffset;
    if (cc1 != cc2)
    {
      if (!fwd && cc1 >= (Seqpos) GT_COMPAREOFFSET
               && cc2 >= (Seqpos) GT_COMPAREOFFSET)
      {
        return cc1 > cc2 ? -1 : 1;
      }
      return cc1 < cc2 ? -1 : 1;
    }
    if (pos1 == pos2 && cc1 >= (Seqpos) GT_COMPAREOFFSET)
    {
      return 0;
    }
  }
  *maxcommon = maxoffset;
  return 0;
}

int comparetwostringsgeneric(const Encodedsequence *encseq,
                             bool fwd,
                             bool complement,
                             Seqpos *maxcommon,
                             Seqpos pos1,
                             Seqpos pos2,
                             Seqpos depth,
                             Seqpos maxdepth)
{
  Seqpos totallength = getencseqtotallength(encseq);
  int retval;
  bool leftspecial, rightspecial;

  if (fwd)
  {
    Seqpos endpos1, endpos2;

    if (maxdepth == 0)
    {
      endpos1 = endpos2 = totallength;
    } else
    {
      gt_assert(maxdepth >= depth);
      endpos1 = MIN(pos1 + maxdepth,totallength);
      endpos2 = MIN(pos2 + maxdepth,totallength);
    }
    if (pos1 + depth < endpos1 && pos2 + depth < endpos2)
    {
      retval = comparetwostrings(encseq,
                                 fwd,
                                 complement,
                                 maxcommon,
                                 pos1+depth,
                                 pos2+depth,
                                 maxdepth > 0 ? (maxdepth - depth) : 0);
    } else
    {
      retval = comparewithonespecial(&leftspecial,
                                     &rightspecial,
                                     encseq,
                                     fwd,
                                     complement,
                                     pos1,
                                     pos2,
                                     depth,
                                     maxdepth);
    }
  } else
  {
    if (maxdepth > 0)
    {
      gt_assert(false);
    }
    if (pos1 >= depth && pos2 >= depth)
    {
      retval = comparetwostrings(encseq,
                                 fwd,
                                 complement,
                                 maxcommon,
                                 pos1-depth,
                                 pos2-depth,
                                 maxdepth > 0 ? (maxdepth - depth) : 0);
    } else
    {
      retval = comparewithonespecial(&leftspecial,
                                     &rightspecial,
                                     encseq,
                                     fwd,
                                     complement,
                                     pos1,
                                     pos2,
                                     depth,
                                     maxdepth);
    }
  }
  *maxcommon += depth;
  return retval;
}
