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

#include "libgtcore/chardef.h"
#include "libgtcore/minmax.h"
#include "encseq-def.h"

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
  Uchar cc1, cc2;
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
    assert(esr1 == NULL && esr2 == NULL);
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

#define COMPAREOFFSET   (UCHAR_MAX + 1)

static Seqpos extractsinglecharacter(const Encodedsequence *encseq,
                                     bool fwd,
                                     bool complement,
                                     Seqpos pos,
                                     Seqpos depth,
                                     Seqpos totallength)
{
  Seqpos cc;

  if (fwd)
  {
    if (pos + depth >= totallength)
    {
      cc = pos + depth + COMPAREOFFSET;
    } else
    {
      cc = getencodedchar(encseq,pos + depth,Forwardmode);
      if (ISSPECIAL(cc))
      {
        cc = pos + depth + COMPAREOFFSET;
      } else
      {
        if (complement)
        {
          cc = COMPLEMENTBASE(cc);
        }
      }
    }
  } else
  {
    if (pos < depth)
    {
      cc = depth - pos + COMPAREOFFSET;
    } else
    {
      cc = getencodedchar(encseq,pos - depth,Forwardmode);
      if (ISSPECIAL(cc))
      {
        cc = pos - depth + COMPAREOFFSET;
      } else
      {
        if (complement)
        {
          cc = COMPLEMENTBASE(cc);
        }
      }
    }
  }
  return cc;
}

int comparewithonespecial(const Encodedsequence *encseq,
                          bool fwd,
                          bool complement,
                          Seqpos pos1,
                          Seqpos pos2,
                          Seqpos depth)
{
  Seqpos cc1, cc2, totallength = getencseqtotallength(encseq);

  cc1 = extractsinglecharacter(encseq,
                               fwd,
                               complement,
                               pos1,
                               depth,
                               totallength);
  cc2 = extractsinglecharacter(encseq,
                               fwd,
                               complement,
                               pos2,
                               depth,
                               totallength);
  assert(cc1 != cc2);
  if (!fwd && cc1 >= (Seqpos) COMPAREOFFSET && cc2 >= (Seqpos) COMPAREOFFSET)
  {
    return cc1 > cc2 ? -1 : 1;
  }
  return cc1 < cc2 ? -1 : 1;
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
      return totallength + COMPAREOFFSET;
    }
    start += currentoffset;
  } else
  {
    if (start < currentoffset)
    {
      return currentoffset - start + (Seqpos) COMPAREOFFSET;
    }
    start -= currentoffset;
  }
  if (currentoffset <= maxoffset)
  {
    Uchar cc;
    cc = getencodedchar(encseq,start,Forwardmode);
    if (ISSPECIAL(cc))
    {
      return start + COMPAREOFFSET;
    }
    if (complement)
    {
      cc = COMPLEMENTBASE(cc);
    }
    return cc;
  }
  return  start + COMPAREOFFSET;
}

int comparetwostrings(const Encodedsequence *encseq,
                      bool fwd,
                      bool complement,
                      Seqpos *maxcommon,
                      Seqpos pos1,
                      Seqpos pos2)
{
  Seqpos currentoffset, maxoffset, cc1, cc2,
         totallength = getencseqtotallength(encseq);

  if (fwd)
  {
    assert(pos1 < totallength);
    assert(pos2 < totallength);
    maxoffset = MIN(totallength - pos1,totallength - pos2);
    if (*maxcommon > 0)
    {
      maxoffset = MIN(*maxcommon,maxoffset);
    }
  } else
  {
    maxoffset = MIN(pos1+1,pos2+1);
    if (*maxcommon > 0)
    {
      maxoffset = MIN(*maxcommon,maxoffset);
    }
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
      if (!fwd && cc1 >= (Seqpos) COMPAREOFFSET
               && cc2 >= (Seqpos) COMPAREOFFSET)
      {
        return cc1 > cc2 ? -1 : 1;
      }
      return cc1 < cc2 ? -1 : 1;
    }
    if (pos1 == pos2 && cc1 >= (Seqpos) COMPAREOFFSET)
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
                             Seqpos depth)
{
  Seqpos totallength = getencseqtotallength(encseq);
  int retval;

  if (fwd)
  {
    if (pos1 + depth < totallength && pos2 + depth < totallength)
    {
      retval = comparetwostrings(encseq,
                                 fwd,
                                 complement,
                                 maxcommon,
                                 pos1+depth,
                                 pos2+depth);
    } else
    {
      retval = comparewithonespecial(encseq,
                                     fwd,
                                     complement,
                                     pos1,
                                     pos2,
                                     depth);
    }
  } else
  {
    if (pos1 >= depth && pos2 >= depth)
    {
      retval = comparetwostrings(encseq,
                                 fwd,
                                 complement,
                                 maxcommon,
                                 pos1-depth,
                                 pos2-depth);
    } else
    {
      retval = comparewithonespecial(encseq,
                                     fwd,
                                     complement,
                                     pos1,
                                     pos2,
                                     depth);
    }
  }
  *maxcommon += depth;
  return retval;
}
