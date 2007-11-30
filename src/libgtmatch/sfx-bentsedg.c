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

#include <limits.h>
#include <assert.h>
#include "libgtcore/arraydef.h"
#include "libgtcore/chardef.h"
#include "libgtcore/env.h"
#include "libgtcore/minmax.h"
#include "divmodmul.h"
#include "sfx-codespec.h"
#include "intcode-def.h"
#include "lcpinterval.h"
#include "sfx-lcpsub.h"
#include "spacedef.h"
#include "encseq-def.h"

#define COMPAREOFFSET   (UCHAR_MAX + 1)
#define UNIQUEINT(P)    ((Seqpos) ((P) + COMPAREOFFSET))
#define ACCESSCHAR(POS) getencodedchar(encseq,POS,readmode) /* XXX */
#define ISNOTEND(POS)   ((POS) < totallength && ISNOTSPECIAL(ACCESSCHAR(POS)))

#define DECLARETMPC Uchar tmpsvar, tmptvar
#define DEREF(VAR,A,S)\
        (((A) < totallength && ISNOTSPECIAL(VAR = ACCESSCHAR(S))) ?\
        ((Seqpos) VAR) : UNIQUEINT(S))

#define PTR2INT(VAR,I) DEREF(VAR,cptr = *(I)+depth,cptr)

#define STRINGCOMPARE(S,T,OFFSET)\
        for (sptr = (S)+(OFFSET), tptr = (T)+(OFFSET); /* Nothing */;\
             sptr++, tptr++)\
        {\
          ccs = DEREF(tmpsvar,sptr,sptr);\
          cct = DEREF(tmptvar,tptr,tptr);\
          if (ccs != cct)\
          {\
            break;\
          }\
        }

#define SWAP(A,B)\
        if ((A) != (B))\
        {\
          temp = *(A);\
          *(A) = *(B);\
          *(B) = temp;\
        }

#define VECSWAP(A,B,N)\
        aptr = A;\
        bptr = B;\
        while ((N)-- > 0)\
        {\
          temp = *aptr;\
          *aptr++ = *bptr;\
          *bptr++ = temp;\
        }

#define SMALLSIZE (Seqpos) 6

#define SUBSORT(WIDTH,BORDER,LEFT,RIGHT,DEPTH)\
        if ((WIDTH) <= (BORDER))\
        {\
          insertionsort(encseq,readmode,totallength,DEPTH,LEFT,RIGHT);\
        } else\
        {\
          PUSHMKVSTACK(LEFT,RIGHT,DEPTH);\
        }

#define PUSHMKVSTACK(L,R,D)\
        CHECKARRAYSPACE(mkvauxstack,MKVstack,1024);\
        mkvauxstack->spaceMKVstack[mkvauxstack->nextfreeMKVstack].left = L;\
        mkvauxstack->spaceMKVstack[mkvauxstack->nextfreeMKVstack].right = R;\
        mkvauxstack->spaceMKVstack[mkvauxstack->nextfreeMKVstack++].depth = D

#define POPMKVstack(L,R,D)\
        L = mkvauxstack->spaceMKVstack[--mkvauxstack->nextfreeMKVstack].left;\
        R = mkvauxstack->spaceMKVstack[mkvauxstack->nextfreeMKVstack].right;\
        D = mkvauxstack->spaceMKVstack[mkvauxstack->nextfreeMKVstack].depth;\
        width = (Seqpos) ((R) - (L) + 1)

typedef Seqpos Suffixptr;

static Suffixptr *medianof3(const Encodedsequence *encseq,
                            Readmode readmode,
                            Seqpos totallength,
                            Seqpos depth,
                            Suffixptr *a,
                            Suffixptr *b,
                            Suffixptr *c)
{
  Suffixptr cptr;
  Seqpos vala, valb, valc;
  Uchar tmpsvar, tmptvar;

  vala = PTR2INT(tmpsvar,a);
  valb = PTR2INT(tmptvar,b);
  if (vala == valb)
  {
    return a;
  }
  if ((valc = PTR2INT(tmpsvar,c)) == vala || valc == valb)
  {
    return c;
  }
  return vala < valb ?
        (valb < valc ? b : (vala < valc ? c : a))
      : (valb > valc ? b : (vala < valc ? a : c));
}

static void insertionsort(const Encodedsequence *encseq,
                          Readmode readmode,
                          Seqpos totallength,
                          Seqpos depth,
                          Suffixptr *left,
                          Suffixptr *right)
{
  Suffixptr sptr, tptr;
  Suffixptr *pi, *pj, temp;
  Seqpos ccs, cct;
  Uchar tmpsvar, tmptvar;

  for (pi = left + 1; pi <= right; pi++)
  {
    for (pj = pi; pj > left; pj--)
    {
      STRINGCOMPARE(*(pj-1),*pj,depth);
      if (ccs < cct)
      {
        break;
      }
      SWAP(pj,pj-1);
    }
  }
}

typedef struct
{
  Suffixptr *left,
            *right;
  Seqpos depth;
} MKVstack;

DECLAREARRAYSTRUCT(MKVstack);

static void bentleysedgewick(const Encodedsequence *encseq,
                             Readmode readmode,
                             Seqpos totallength,
                             ArrayMKVstack *mkvauxstack,
                             Suffixptr *l,Suffixptr *r,Seqpos d,
                             Lcpsubtab *lcpsubtab,
                             Env *env)
{
  Suffixptr *left, *right, *leftplusw;
  Seqpos w, val, partval, depth, offset, doubleoffset, width;
  Suffixptr *pa, *pb, *pc, *pd, *pl, *pm, *pr, *aptr, *bptr, cptr, temp;
  Uchar tmpsvar;

  env_error_check(env);
  width = (Seqpos) (r - l + 1);
  if (width <= SMALLSIZE)
  {
    insertionsort(encseq,readmode,totallength,d,l,r);
    return;
  }
  left = l;
  right = r;
  depth = d;
  mkvauxstack->nextfreeMKVstack = 0;

  for (;;)
  {
    pl = left;
    pm = left + DIV2(width);
    pr = right;
    if (width > (Seqpos) 30)
    { /* On big arrays, pseudomedian of 9 */
      offset = DIV8(width);
      doubleoffset = MULT2(offset);
      pl = medianof3(encseq,readmode,totallength,depth,
                     pl,pl+offset,pl+doubleoffset);
      pm = medianof3(encseq,readmode,totallength,depth,
                     pm-offset,pm,pm+offset);
      pr = medianof3(encseq,readmode,totallength,depth,
                     pr-doubleoffset,pr-offset,pr);
    }
    pm = medianof3(encseq,readmode,totallength,depth,pl,pm,pr);
    SWAP(left, pm);
    partval = PTR2INT(tmpsvar,left);
    pa = pb = left + 1;
    pc = pd = right;
    for (;;)
    {
      while (pb <= pc)
      {
        if ((val = PTR2INT(tmpsvar,pb)) > partval)
        {
          break;
        }
        if (val == partval)
        {
          SWAP(pa, pb);
          pa++;
        }
        pb++;
      }
      while (pb <= pc)
      {
        if ((val = PTR2INT(tmpsvar,pc)) < partval)
        {
          break;
        }
        if (val == partval)
        {
          SWAP(pc, pd);
          pd--;
        }
        pc--;
      }
      if (pb > pc)
      {
        break;
      }
      SWAP(pb, pc);
      pb++;
      pc--;
    }

    assert(pa >= left);
    assert(pb >= pa);
    w = MIN((Seqpos) (pa-left),(Seqpos) (pb-pa));
    VECSWAP(left,  pb-w, w);
    pr = right + 1;
    assert(pd >= pc);
    assert(pr > pd);
    w = MIN((Seqpos) (pd-pc), (Seqpos) (pr-pd-1));
    VECSWAP(pb, pr-w, w);

    assert(pd >= pc);
    if ((w = (Seqpos) (pd-pc)) > 0)
    {
      SUBSORT(w,SMALLSIZE,right-w+1,right,depth);
    }
    assert(pb >= pa);
    w = (Seqpos) (pb-pa);
    leftplusw = left + w;
    cptr = *leftplusw + depth;
    if (ISNOTEND(cptr))
    {
      right -= (pd-pb);
      width = (Seqpos) (right-leftplusw);
      SUBSORT(width,SMALLSIZE,leftplusw,right-1,depth+1);
    }
    if (w > 0)
    {
      SUBSORT(w,SMALLSIZE,left,leftplusw-1,depth);
    }
    if (mkvauxstack->nextfreeMKVstack == 0)
    {
      break;
    }
    POPMKVstack(left,right,depth);
  }
}

static unsigned int bucketboundaries(Lcpinterval *lcpitv,
                                     const Seqpos *leftborder,
                                     const Seqpos *countspecialcodes,
                                     Codetype code,
                                     Codetype maxcode,
                                     Seqpos totalwidth,
                                     unsigned int rightchar,
                                     unsigned int numofchars)
{
  lcpitv->left = leftborder[code];
  if (code == maxcode)
  {
    assert(totalwidth > 0);
    lcpitv->right = totalwidth - 1;
  } else
  {
    if (leftborder[code+1] > 0)
    {
      lcpitv->right = leftborder[code+1] - 1;
    } else
    {
      lcpitv->right = 0;
    }
  }
  assert(rightchar == code % numofchars);
  if (rightchar == numofchars - 1)
  {
    Seqpos specialcodes
      = countspecialcodes[FROMCODE2SPECIALCODE(code,numofchars)];
    if (lcpitv->right >= specialcodes)
    {
      lcpitv->right -= specialcodes;
    } else
    {
      lcpitv->right = 0;
    }
     rightchar = 0;
   } else
   {
     rightchar++;
   }
   return rightchar;
}

static Seqpos determinemaxbucketsize(const Seqpos *leftborder,
                                     const Seqpos *countspecialcodes,
                                     const Codetype mincode,
                                     const Codetype maxcode,
                                     Seqpos totalwidth,
                                     unsigned int numofchars)
{
  Seqpos maxbucketsize = (Seqpos) 1, bsize;
  unsigned int rightchar = mincode % numofchars;
  Lcpinterval lcpitv;
  Codetype code;

  for (code = mincode; code <= maxcode; code++)
  {
    rightchar = bucketboundaries(&lcpitv,
                                 leftborder,
                                 countspecialcodes,
                                 code,
                                 maxcode,
                                 totalwidth,
                                 rightchar,
                                 numofchars);
    if (lcpitv.left < lcpitv.right)
    {
      bsize = lcpitv.right - lcpitv.left + 1;
      if (bsize > maxbucketsize)
      {
        maxbucketsize = bsize;
      }
    }
  }
  return maxbucketsize;
}

 struct Lcpsubtab
{
  Seqpos *spaceSeqpos;
  unsigned long nextfreeSeqpos, allocatedSeqpos;
};

void freelcpsubtab(Lcpsubtab **lcpsubtab,Env *env)
{
  FREEARRAY(*lcpsubtab,Seqpos);
  FREESPACE(*lcpsubtab);
  return;
}

Lcpsubtab *newlcpsubtab(Env *env)
{
  Lcpsubtab *lcpsubtab;

  ALLOCASSIGNSPACE(lcpsubtab,NULL,Lcpsubtab,1);
  INITARRAY(lcpsubtab,Seqpos);
  return lcpsubtab;
}

void sortallbuckets(Seqpos *suftabptr,
                    const Encodedsequence *encseq,
                    Readmode readmode,
                    const Seqpos *leftborder,
                    const Seqpos *countspecialcodes,
                    unsigned int numofchars,
                    unsigned int prefixlength,
                    Codetype mincode,
                    Codetype maxcode,
                    Seqpos totalwidth,
                    Lcpsubtab *lcpsubtab,
                    Env *env)
{
  Codetype code;
  unsigned int rightchar = mincode % numofchars;
  Seqpos totallength = getencseqtotallength(encseq);
  ArrayMKVstack mkvauxstack;
  Lcpinterval lcpitv;
  Seqpos maxbucketsize;

  env_error_check(env);
  maxbucketsize = determinemaxbucketsize(leftborder,
                                         countspecialcodes,
                                         mincode,
                                         maxcode,
                                         totalwidth,
                                         numofchars);
  if (maxbucketsize > lcpsubtab->allocatedSeqpos)
  {
    lcpsubtab->allocatedSeqpos = maxbucketsize;
    ALLOCASSIGNSPACE(lcpsubtab->spaceSeqpos,lcpsubtab->spaceSeqpos,Seqpos,
                     lcpsubtab->allocatedSeqpos);
  }
  INITARRAY(&mkvauxstack,MKVstack);
  for (code = mincode; code <= maxcode; code++)
  {
    rightchar = bucketboundaries(&lcpitv,
                                 leftborder,
                                 countspecialcodes,
                                 code,
                                 maxcode,
                                 totalwidth,
                                 rightchar,
                                 numofchars);
    if (lcpitv.left < lcpitv.right)
    {
      bentleysedgewick(encseq,
                       readmode,
                       totallength,
                       &mkvauxstack,
                       suftabptr + lcpitv.left,
                       suftabptr + lcpitv.right,
                       (Seqpos) prefixlength,
                       lcpsubtab,
                       env);
    }
  }
  FREEARRAY(&mkvauxstack,MKVstack);
}
