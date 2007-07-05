/*
  Copyright (c) 2007 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
  Copyright (c) 2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#include <limits.h>
#include <assert.h>
#include "libgtcore/env.h"
#include "types.h"
#include "chardef.h"
#include "arraydef.h"
#include "divmodmul.h"
#include "minmax.h"
#include "codespec.h"
#include "intcode-def.h"
#include "encseq-def.h"

#define COMPAREOFFSET   (UCHAR_MAX + 1)
#define UNIQUEINT(P)    ((Seqpos) ((P) + COMPAREOFFSET))
#define ACCESSCHAR(T)   getencodedchar(encseq,T)
#define ISNOTEND(P)     ((P) < totallength && ISNOTSPECIAL(ACCESSCHAR(P)))

#define DECLARETMPC Uchar tmpsvar, tmptvar
#define GENDEREF(VAR,A,S)\
        (((A) < totallength && ISNOTSPECIAL(VAR = ACCESSCHAR(S))) ?\
        ((Seqpos) VAR) : UNIQUEINT(S))
#define DEREF(VAR,S)   GENDEREF(VAR,S,S)

#define PTR2INT(VAR,I) GENDEREF(VAR,cptr = *(I)+depth,cptr)

#define STRINGCOMPARE(S,T,OFFSET)\
        for (sptr = (S)+(OFFSET), tptr = (T)+(OFFSET); /* Nothing */;\
             sptr++, tptr++)\
        {\
          ccs = DEREF(tmpsvar,sptr);\
          cct = DEREF(tmptvar,tptr);\
          if(ccs != cct)\
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

#define SMALLSIZE 6

#define SUBSORT(WIDTH,BORDER,LEFT,RIGHT,DEPTH)\
        if ((WIDTH) <= (BORDER))\
        {\
          insertionsort(encseq,totallength,DEPTH,LEFT,RIGHT);\
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
                             Seqpos totallength,
                             ArrayMKVstack *mkvauxstack,
                             Suffixptr *l,Suffixptr *r,Seqpos d,
                             Env *env)
{
  Suffixptr *left, *right, *leftplusw;
  Seqpos w, val, partval, depth, offset, doubleoffset, width;
  Suffixptr *pa, *pb, *pc, *pd, *pl, *pm, *pr, *aptr, *bptr, cptr, temp;
  Uchar tmpsvar;

  width = (Seqpos) (r - l + 1);
  if (width <= (Seqpos) (SMALLSIZE))
  {
    insertionsort(encseq,totallength,d,l,r);
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
      pl = medianof3(encseq,totallength,depth,pl,pl+offset,pl+doubleoffset);
      pm = medianof3(encseq,totallength,depth,pm-offset,pm,pm+offset);
      pr = medianof3(encseq,totallength,depth,pr-doubleoffset,pr-offset,pr);
    }
    pm = medianof3(encseq,totallength,depth,pl,pm,pr);
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
      SUBSORT(w,(Seqpos) (SMALLSIZE),right-w+1,right,depth);
    }
    assert(pb >= pa);
    w = (Seqpos) (pb-pa);
    leftplusw = left + w;
    cptr = *leftplusw + depth;
    if (ISNOTEND(cptr))
    {
      right -= (pd-pb);
      width = (Seqpos) (right-leftplusw);
      SUBSORT(width,(Seqpos) (SMALLSIZE),leftplusw,right-1,depth+1);
    }
    if (w > 0)
    {
      SUBSORT(w,(Seqpos) (SMALLSIZE),left,leftplusw-1,depth);
    }
    if (mkvauxstack->nextfreeMKVstack == 0)
    {
      break;
    }
    POPMKVstack(left,right,depth);
  }
}

void sortallbuckets(Seqpos *suftabptr,
                    const Encodedsequence *encseq,
                    const Seqpos *leftborder,
                    const Seqpos *countspecialcodes,
                    uint32_t numofchars,
                    uint32_t prefixlength,
                    Codetype mincode,
                    Codetype maxcode,
                    uint64_t widthofpart,
                    Env *env)
{
  Codetype code;
  uint32_t rightchar = mincode % numofchars;
  Seqpos left, right, specialcodes;
  ArrayMKVstack mkvauxstack;
  Seqpos totallength = getencseqtotallength(encseq);

  INITARRAY(&mkvauxstack,MKVstack);
  for (code=mincode; code<=maxcode; code++)
  {
    left = leftborder[code];
    if (code == maxcode)
    {
      assert(widthofpart > 0);
      right = widthofpart - 1;
    } else
    {
      if (leftborder[code+1] > 0)
      {
        right = leftborder[code+1] - 1;
      } else
      {
        right = 0;
      }
    }
    if (rightchar == numofchars - 1)
    {
      specialcodes = countspecialcodes[FROMCODE2SPECIALCODE(code,numofchars)];
      if (right >= specialcodes)
      {
        right -= specialcodes;
      } else
      {
        right = 0;
      }
      rightchar = 0;
    } else
    {
      rightchar++;
    }
    if (left < right)
    {
      bentleysedgewick(encseq,
                       totallength,
                       &mkvauxstack,
                       suftabptr + left,
                       suftabptr + right,
                       (Seqpos) prefixlength,
                       env);
    }
  }
  FREEARRAY(&mkvauxstack,MKVstack);
}
