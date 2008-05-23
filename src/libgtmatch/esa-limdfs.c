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
#include "libgtcore/symboldef.h"
#include "libgtcore/unused.h"
#include "libgtcore/arraydef.h"
#include "sarr-def.h"
#include "seqpos-def.h"
#include "esa-splititv.h"
#include "spacedef.h"

typedef struct
{
  Lcpinterval lcpitv;
  void *info;
} Lcpintervalwithinfo;

DECLAREARRAYSTRUCT(Lcpintervalwithinfo);

void limiteddfs(const Encodedsequence *encseq,
                const Seqpos *suftab,
                unsigned int mapsize)
{
  ArrayLcpintervalwithinfo stack;
  Lcpintervalwithinfo *stackptr;
  Rightboundwithchar *rbwc;
  unsigned long idx, rboundscount;
  Uchar extendchar, alphasize = (Uchar) (mapsize-1);
  Seqpos lbound, rbound, offset, totallength = getencseqtotallength(encseq);

  INITARRAY(&stack,Lcpintervalwithinfo);
  ALLOCASSIGNSPACE(rbwc,NULL,Rightboundwithchar,alphasize);
  rboundscount = lcpintervalsplitwithoutspecial(rbwc,
                                                alphasize,
                                                encseq,
                                                suftab,
                                                0,
                                                0,
                                                totallength);
  for (idx=0; idx < rboundscount; idx++)
  {
    lbound = rbwc[idx].bound;
    rbound = rbwc[idx+1].bound-1;
    assert(lbound <= rbound);
    if (lbound < rbound)
    {
      assert(stack.spaceLcpintervalwithinfo != NULL);
      GETNEXTFREEINARRAY(stackptr,&stack,Lcpintervalwithinfo,128);
      stackptr->lcpitv.left = lbound;
      stackptr->lcpitv.right = rbound;
      stackptr->lcpitv.offset = (Seqpos) 1;
    } else
    {
      printf("singleton " FormatSeqpos "\n",PRINTSeqposcast(lbound));
    }
  }
  while (stack.nextfreeLcpintervalwithinfo > 0)
  {
    assert(stack.spaceLcpintervalwithinfo != NULL);
    stackptr = stack.spaceLcpintervalwithinfo +
               stack.nextfreeLcpintervalwithinfo - 1;
    extendchar = lcpintervalextendlcp(encseq,
                                      suftab,
                                      &stackptr->lcpitv,
                                      alphasize);
    if (extendchar < alphasize)
    {
      stackptr->lcpitv.offset++;
    } else
    {
      rboundscount = lcpintervalsplitwithoutspecial(rbwc,
                                                    alphasize,
                                                    encseq,
                                                    suftab,
                                                    stackptr->lcpitv.offset,
                                                    stackptr->lcpitv.left,
                                                    stackptr->lcpitv.right);
      for (idx=0; idx < rboundscount; idx++)
      {
        lbound = rbwc[idx].bound;
        rbound = rbwc[idx+1].bound-1;
        offset = stackptr->lcpitv.offset + 1;
        assert(lbound <= rbound);
        if (lbound < rbound)
        {
          GETNEXTFREEINARRAY(stackptr,&stack,Lcpintervalwithinfo,128);
          stackptr->lcpitv.left = lbound;
          stackptr->lcpitv.right = rbound;
          stackptr->lcpitv.offset = offset;
        } else
        {
          printf("singleton " FormatSeqpos "\n",PRINTSeqposcast(lbound));
        }
      }
    }
  }
  FREESPACE(rbwc);
  FREEARRAY(&stack,Lcpintervalwithinfo);
}
