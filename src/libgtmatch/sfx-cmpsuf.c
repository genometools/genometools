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
#include "encseq-def.h"

int comparetwosuffixes(const Encodedsequence *encseq,
                       Readmode readmode,
                       Seqpos *maxlcp,
                       bool specialsareequal,
                       bool specialsareequalatdepth0,
                       Seqpos depth,
                       Seqpos start1,
                       Seqpos start2)
{
  Uchar cc1, cc2;
  Seqpos pos1, pos2, end1, end2;

  end1 = end2 = getencseqtotallength(encseq);
  if (depth > 0)
  {
    if (end1 > start1 + depth)
    {
      end1 = start1 + depth;
    }
    if (end2 > start2 + depth)
    {
      end2 = start2 + depth;
    }
  }
  for (pos1=start1, pos2=start2; pos1 < end1 && pos2 < end2; pos1++, pos2++)
  {
    cc1 = getencodedchar(encseq,pos1,readmode);
    /* printf("pos1=" FormatSeqpos "cc1=%u\n",pos1,cc1); */
    cc2 = getencodedchar(encseq,pos2,readmode);
    /* printf("pos2=" FormatSeqpos "cc2=%u\n",pos2,cc2); */
    if (ISSPECIAL(cc1))
    {
      if (ISSPECIAL(cc2))
      {
        if (specialsareequal || (pos1 == start1 && specialsareequalatdepth0))
        {
          *maxlcp = pos1 - start1 + 1;
          return 0;
        }
        if (pos1 < pos2)
        {
          *maxlcp = pos1  - start1;
          return -1; /* a < b */
        }
        if (pos1 > pos2)
        {
          *maxlcp = pos1 - start1;
          return 1; /* a > b */
        }
        *maxlcp = pos1 - start1 + 1;
        return 0; /* a = b */
      }
      *maxlcp = pos1 - start1;
      return 1; /* a > b */
    } else
    {
      if (ISSPECIAL(cc2))
      {
        *maxlcp = pos1 - start1;
        return -1; /* a < b */
      }
      if (cc1 < cc2)
      {
        *maxlcp = pos1 - start1;
        return -1;
      }
      if (cc1 > cc2)
      {
        *maxlcp = pos1 - start1;
        return 1;
      }
    }
  }
  *maxlcp = pos1 - start1;
  return 0;
}
