/*
  Copyright (c) 2007 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
  Copyright (c) 2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#include "libgtcore/minmax.h"
#include "sarr-def.h"
#include "seqpos-def.h"
#include "divmodmul.h"
#include "symboldef.h"
#include "chardef.h"

#define COMPARE(OFFSET)\
        for(sidx = (OFFSET) + lcplen;\
            /* Nothing */; sidx++, lcplen++)\
        {\
          if(lcplen >= (Seqpos) querylen)\
          {\
            retcode = 0;\
            break;\
          }\
          if(sidx >= totallength)\
          {\
            retcode = -1;\
            break;\
          }\
          currentchar = getencodedchar(encseq,sidx,readmode);\
          retcode = (int) (query[lcplen] - currentchar);\
          if(retcode == 0)\
          {\
            if(ISSPECIAL(currentchar) && ISSPECIAL(query[lcplen]))\
            {\
              retcode = (int) -1;\
              break;\
            }\
          } else\
          {\
            break;\
          }\
        }

typedef struct
{
  Seqpos offset, 
         left, 
         right;
} Lcpinterval;

bool mmsearch(const Encodedsequence *encseq,
              const Seqpos *suftab,
              Readmode readmode,
              Lcpinterval *vnode,
              Uchar *query,
              unsigned long querylen)
{
  Seqpos sidx, left, leftsave, mid, right, lpref, rpref, totallength;
  int retcode = 0;
  Uchar currentchar;
  Seqpos lcplen;

  totallength = getencseqtotallength(encseq);
  leftsave = left = vnode->left;
  right = vnode->right;
  lcplen = vnode->offset;
  COMPARE(suftab[left]);
  if(retcode > 0)
  {
    lpref = lcplen;
    lcplen = vnode->offset;
    COMPARE(suftab[right]);
    if(retcode > 0)
    {
      return false;
    } else
    {
      rpref = lcplen;
      while(right > left + 1)
      {
        mid = DIV2(left+right);
        lcplen = MIN(lpref,rpref);
        COMPARE(suftab[mid]);
        if(retcode <= 0)
        {
          right = mid;
          rpref = lcplen;
        } else
        {
          left = mid;
          lpref = lcplen;
        }
      }
      vnode->left = right;
    }
  }

  left = leftsave;
  right = vnode->right;
  lcplen = vnode->offset;
  COMPARE(suftab[left]);
  if(retcode < 0)
  {
    return false;
  } else
  {
    lpref = lcplen;
    lcplen = vnode->offset;
    COMPARE(suftab[right]);
    if(retcode >= 0)
    {
      vnode->right = right;
    } else
    {
      rpref = lcplen;
      while(right > left + 1)
      {
        mid = DIV2(left+right);
        lcplen = MIN(lpref,rpref);
        COMPARE(suftab[mid]);
        if(retcode >= 0)
        {
          left = mid;
          lpref = lcplen;
        } else
        {
          right = mid;
          rpref = lcplen;
        }
      }
      vnode->right = left;
    }
  }
  return true;
}
