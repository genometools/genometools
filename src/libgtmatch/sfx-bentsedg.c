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
#include "libgtcore/chardef.h"
#include "libgtcore/minmax.h"
#include "libgtcore/xansi.h"
#include "libgtcore/fa.h"
#include "libgtcore/arraydef.h"
#include "divmodmul.h"
#include "spacedef.h"
#include "encseq-def.h"
#include "turnwheels.h"
#include "esafileend.h"
#include "sfx-outlcp.h"
#include "bcktab.h"

#include "sfx-cmpsuf.pr"
#include "opensfxfile.pr"

#define COMPAREOFFSET   (UCHAR_MAX + 1)
#define UNIQUEINT(P)    ((Seqpos) ((P) + COMPAREOFFSET))
#define ACCESSCHAR(POS) getencodedchar(encseq,POS,readmode) /* XXX */
#define ISNOTEND(POS)   ((POS) < totallength && ISNOTSPECIAL(ACCESSCHAR(POS)))

#define DECLARETMPC Uchar tmpsvar, tmptvar
#define DEREF(VAR,A,S)\
        (((A) < totallength && ISNOTSPECIAL(VAR = ACCESSCHAR(S))) ?\
        ((Seqpos) VAR) : UNIQUEINT(S))

#define PTR2INT(VAR,I) DEREF(VAR,cptr = *(I)+depth,cptr)

#define LCPINDEX(I)       (Seqpos) ((I) - lcpsubtab->suftabbase)
#define SETLCP(I,V)       lcpsubtab->spaceSeqpos[I] = V

#define STRINGCOMPARE(S,T,OFFSET,LL)\
        for (sptr = (S)+(OFFSET), tptr = (T)+(OFFSET); /* Nothing */;\
             sptr++, tptr++)\
        {\
          ccs = DEREF(tmpsvar,sptr,sptr);\
          cct = DEREF(tmptvar,tptr,tptr);\
          if (ccs != cct)\
          {\
            LL = (Seqpos) (tptr - (T));\
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
          insertionsort(encseq,lcpsubtab,readmode,totallength,\
                        DEPTH,LEFT,RIGHT);\
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

DECLAREARRAYSTRUCT(Largelcpvalue);

typedef struct
{
  Seqpos *spaceSeqpos;
  Uchar *smalllcpvalues;
  ArrayLargelcpvalue largelcpvalues;
  unsigned long nextfreeSeqpos, allocatedSeqpos;
  const Seqpos *suftabbase;
} Lcpsubtab;

 struct Outlcpinfo
{
  FILE *outfplcptab,
       *outfpllvtab;
  Seqpos totallength,
         countoutputlcpvalues,
         numoflargelcpvalues,
         maxbranchdepth;
  Turningwheel *tw;
  Lcpsubtab lcpsubtab;
};

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
                          Lcpsubtab *lcpsubtab,
                          Readmode readmode,
                          Seqpos totallength,
                          Seqpos depth,
                          Suffixptr *leftptr,
                          Suffixptr *rightptr)
{
  Suffixptr sptr, tptr, *pi, *pj, temp;
  Seqpos ccs, cct, lcpindex, lcplen;
  Uchar tmpsvar, tmptvar;

  for (pi = leftptr + 1; pi <= rightptr; pi++)
  {
    for (pj = pi; pj > leftptr; pj--)
    {
      STRINGCOMPARE(*(pj-1),*pj,depth,lcplen);
      if (lcpsubtab != NULL)
      {
        lcpindex = LCPINDEX(pj);
        if (ccs > cct && pj < pi)
        {
          SETLCP(lcpindex+1,lcpsubtab->spaceSeqpos[lcpindex]);
        }
        SETLCP(lcpindex,lcplen);
      }
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
                             Suffixptr *l,
                             Suffixptr *r,
                             Seqpos d,
                             Lcpsubtab *lcpsubtab)
{
  Suffixptr *left, *right, *leftplusw;
  Seqpos w, val, partval, depth, offset, doubleoffset, width;
  Suffixptr *pa, *pb, *pc, *pd, *pl, *pm, *pr, *aptr, *bptr, cptr, temp;
  Uchar tmpsvar;

  width = (Seqpos) (r - l + 1);
  if (width <= SMALLSIZE)
  {
    insertionsort(encseq,lcpsubtab,readmode,totallength,d,l,r);
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
      if (lcpsubtab != NULL)
      {
        SETLCP(LCPINDEX(right-w+1),depth);
      }
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
      if (lcpsubtab != NULL)
      {
        SETLCP(LCPINDEX(leftplusw),depth);
      }
      SUBSORT(w,SMALLSIZE,left,leftplusw-1,depth);
    }
    if (mkvauxstack->nextfreeMKVstack == 0)
    {
      break;
    }
    POPMKVstack(left,right,depth);
  }
}

#define NUMBEROFZEROS 1024

static void outmany0lcpvalues(Seqpos many,Outlcpinfo *outlcpinfo)
{
  Seqpos i, countout;
  Uchar outvalues[NUMBEROFZEROS] = {0};

  countout = many/NUMBEROFZEROS;
  for (i=0; i<countout; i++)
  {
    xfwrite(outvalues,sizeof (Uchar),(size_t) NUMBEROFZEROS,
            outlcpinfo->outfplcptab);
  }
  xfwrite(outvalues,sizeof (Uchar),(size_t) many % NUMBEROFZEROS,
          outlcpinfo->outfplcptab);
  outlcpinfo->countoutputlcpvalues += many;
}

static unsigned long determinemaxbucketsize(const Bcktab *bcktab,
                                            const Codetype mincode,
                                            const Codetype maxcode,
                                            Seqpos totalwidth,
                                            unsigned int numofchars)
{
  unsigned long maxbucketsize = 1UL;
  unsigned int rightchar = mincode % numofchars;
  Bucketspecification bucketspec;
  Codetype code;

  for (code = mincode; code <= maxcode; code++)
  {
    rightchar = calcbucketboundsparts(&bucketspec,
                                      bcktab,
                                      code,
                                      maxcode,
                                      totalwidth,
                                      rightchar,
                                      numofchars);
    if (bucketspec.nonspecialsinbucket > maxbucketsize)
    {
      maxbucketsize = bucketspec.nonspecialsinbucket;
    }
    if (bucketspec.specialsinbucket > maxbucketsize)
    {
      maxbucketsize = bucketspec.specialsinbucket;
    }
  }
  return maxbucketsize;
}

static void multilcpvalue(Outlcpinfo *outlcpinfo,
                          unsigned long bucketsize,
                          Seqpos posoffset)
{
  unsigned long i;
  Seqpos lcpvalue;
  Largelcpvalue *largelcpvalueptr;

  outlcpinfo->lcpsubtab.largelcpvalues.nextfreeLargelcpvalue = 0;
  for (i=0; i<bucketsize; i++)
  {
    lcpvalue = outlcpinfo->lcpsubtab.spaceSeqpos[i];
    if (outlcpinfo->maxbranchdepth < lcpvalue)
    {
      outlcpinfo->maxbranchdepth = lcpvalue;
    }
    if (lcpvalue >= (Seqpos) UCHAR_MAX)
    {
      outlcpinfo->numoflargelcpvalues++;
      GETNEXTFREEINARRAY(largelcpvalueptr,&outlcpinfo->lcpsubtab.largelcpvalues,
                         Largelcpvalue,32);
      largelcpvalueptr->position = posoffset+i;
      largelcpvalueptr->value = lcpvalue;
      outlcpinfo->lcpsubtab.smalllcpvalues[i] = (Uchar) UCHAR_MAX;
    } else
    {
      outlcpinfo->lcpsubtab.smalllcpvalues[i] = (Uchar) lcpvalue;
    }
  }
  outlcpinfo->countoutputlcpvalues += bucketsize;
  xfwrite(outlcpinfo->lcpsubtab.smalllcpvalues,
          sizeof (Uchar),(size_t) bucketsize,outlcpinfo->outfplcptab);
  xfwrite(outlcpinfo->lcpsubtab.largelcpvalues.spaceLargelcpvalue,
          sizeof (Largelcpvalue),
          (size_t) outlcpinfo->lcpsubtab.largelcpvalues.nextfreeLargelcpvalue,
          outlcpinfo->outfpllvtab);
}

static Seqpos computelocallcpvalue(const Encodedsequence *encseq,
                                   Readmode readmode,
                                   const Suffixwithcode *previoussuffix,
                                   Seqpos suffixpos2,
                                   Encodedsequencescanstate *esr1,
                                   Encodedsequencescanstate *esr2)
{
  Seqpos lcpvalue;
  int cmp;

  cmp = comparetwosuffixes(encseq,
                           readmode,
                           &lcpvalue,
                           false,
                           false,
                           0,
                           previoussuffix->startpos,
                           suffixpos2,
                           esr1,
                           esr2);
  if (cmp > 0)
  {
    fprintf(stderr,"cmp " FormatSeqpos
            " " FormatSeqpos " = %d, lcpval=" FormatSeqpos "\n",
            PRINTSeqposcast(previoussuffix->startpos),
            PRINTSeqposcast(suffixpos2),
            cmp,
            PRINTSeqposcast(lcpvalue));
    exit(EXIT_FAILURE);
  }
  return lcpvalue;
}

static unsigned int bucketends(const Encodedsequence *encseq,
                               Readmode readmode,
                               Encodedsequencescanstate *esr1,
                               Encodedsequencescanstate *esr2,
                               Outlcpinfo *outlcpinfo,
                               Suffixwithcode *previoussuffix,
                               Seqpos firstspecialsuffix,
                               unsigned long specialsinbucket,
                               Codetype code,
                               unsigned int numofchars,
                               const Bcktab *bcktab)
{
  Seqpos lcpvalue;
  unsigned int maxprefixindex, minprefixindex;
  Suffixwithcode firstspecialsuffixwithcode;

  /*
     there is at least one element in the bucket. if there is more than
     one element in the bucket, then we insert them using the
     information from the bcktab
  */
  if (specialsinbucket > 1UL)
  {
    maxprefixindex = pfxidx2lcpvalues(&minprefixindex,
                                      outlcpinfo->lcpsubtab.smalllcpvalues,
                                      specialsinbucket,
                                      bcktab,
                                      code);
    if (outlcpinfo->maxbranchdepth < (Seqpos) maxprefixindex)
    {
      outlcpinfo->maxbranchdepth = (Seqpos) maxprefixindex;
    }
  } else
  {
    minprefixindex = maxprefixindex = singletonmaxprefixindex(bcktab,code);
  }
  firstspecialsuffixwithcode.code = code;
  firstspecialsuffixwithcode.startpos = firstspecialsuffix;
  firstspecialsuffixwithcode.prefixindex = maxprefixindex;
  consistencyofsuffix(__LINE__,
                      encseq,readmode,bcktab,numofchars,
                      &firstspecialsuffixwithcode);
  lcpvalue = computelocallcpvalue(encseq,
                                  readmode,
                                  previoussuffix,
                                  firstspecialsuffix,
                                  esr1,
                                  esr2);
  if (outlcpinfo->maxbranchdepth < lcpvalue)
  {
    outlcpinfo->maxbranchdepth = lcpvalue;
  }
  outlcpinfo->lcpsubtab.smalllcpvalues[0] = (Uchar) lcpvalue;
  outlcpinfo->countoutputlcpvalues += specialsinbucket;
  xfwrite(outlcpinfo->lcpsubtab.smalllcpvalues,
          sizeof (Uchar),(size_t) specialsinbucket,outlcpinfo->outfplcptab);
  return minprefixindex;
}

Outlcpinfo *newlcpoutfileinfo(const Str *indexname,
                              unsigned int prefixlength,
                              unsigned int numofchars,
                              Seqpos totallength,
                              Error *err)
{
  bool haserr = false;
  Outlcpinfo *outlcpinfo;

  ALLOCASSIGNSPACE(outlcpinfo,NULL,Outlcpinfo,1);
  if (indexname == NULL)
  {
    outlcpinfo->outfplcptab = NULL;
    outlcpinfo->outfpllvtab = NULL;
  } else
  {
    outlcpinfo->outfplcptab = opensfxfile(indexname,LCPTABSUFFIX,"wb",err);
    if (outlcpinfo->outfplcptab == NULL)
    {
      haserr = true;
    }
    if (!haserr)
    {
      outlcpinfo->outfpllvtab
        = opensfxfile(indexname,LARGELCPTABSUFFIX,"wb",err);
      if (outlcpinfo->outfpllvtab == NULL)
      {
        haserr = true;
      }
    }
  }
  outlcpinfo->numoflargelcpvalues = 0;
  outlcpinfo->maxbranchdepth = 0;
  outlcpinfo->countoutputlcpvalues = 0;
  outlcpinfo->totallength = totallength;
  INITARRAY(&outlcpinfo->lcpsubtab,Seqpos);
  INITARRAY(&outlcpinfo->lcpsubtab.largelcpvalues,Largelcpvalue);
  outlcpinfo->lcpsubtab.smalllcpvalues = NULL;
  outlcpinfo->tw = newTurningwheel(prefixlength,numofchars);
  if (haserr)
  {
    FREESPACE(outlcpinfo);
    return NULL;
  }
  return outlcpinfo;
}

void freeoutlcptab(Outlcpinfo **outlcpinfo)
{
  if ((*outlcpinfo)->countoutputlcpvalues < (*outlcpinfo)->totallength + 1)
  {
    outmany0lcpvalues((*outlcpinfo)->totallength + 1 -
                      (*outlcpinfo)->countoutputlcpvalues,
                      *outlcpinfo);
  }
  assert((*outlcpinfo)->countoutputlcpvalues == (*outlcpinfo)->totallength + 1);
  fa_fclose((*outlcpinfo)->outfplcptab);
  fa_fclose((*outlcpinfo)->outfpllvtab);
  FREEARRAY(&(*outlcpinfo)->lcpsubtab,Seqpos);
  freeTurningwheel(&(*outlcpinfo)->tw);
  FREEARRAY(&(*outlcpinfo)->lcpsubtab.largelcpvalues,Largelcpvalue);
  FREESPACE(*outlcpinfo);
}

Seqpos getnumoflargelcpvalues(const Outlcpinfo *outlcpinfo)
{
  return outlcpinfo->numoflargelcpvalues;
}

Seqpos getmaxbranchdepth(const Outlcpinfo *outlcpinfo)
{
  return outlcpinfo->maxbranchdepth;
}

void sortallbuckets(Seqpos *suftabptr,
                    const Encodedsequence *encseq,
                    Readmode readmode,
                    Codetype mincode,
                    Codetype maxcode,
                    Seqpos totalwidth,
                    Suffixwithcode *previoussuffix,
                    const Bcktab *bcktab,
                    unsigned int numofchars,
                    unsigned int prefixlength,
                    Outlcpinfo *outlcpinfo)
{
  Codetype code;
  unsigned int rightchar = mincode % numofchars, minprefixindex;
  Seqpos totallength = getencseqtotallength(encseq);
  ArrayMKVstack mkvauxstack;
  Bucketspecification bucketspec;
  unsigned long maxbucketsize;
  Seqpos lcpvalue;
  Encodedsequencescanstate *esr1, *esr2;
  Lcpsubtab *lcpsubtab;

  if (outlcpinfo == NULL)
  {
    lcpsubtab = NULL;
  } else
  {
    lcpsubtab = &outlcpinfo->lcpsubtab;
  }
  esr1 = newEncodedsequencescanstate();
  esr2 = newEncodedsequencescanstate();
  maxbucketsize = determinemaxbucketsize(bcktab,
                                         mincode,
                                         maxcode,
                                         totalwidth,
                                         numofchars);
  if (lcpsubtab != NULL && maxbucketsize > lcpsubtab->allocatedSeqpos)
  {
    lcpsubtab->allocatedSeqpos = maxbucketsize;
    ALLOCASSIGNSPACE(lcpsubtab->spaceSeqpos,
                     lcpsubtab->spaceSeqpos,Seqpos,
                     lcpsubtab->allocatedSeqpos);
    lcpsubtab->smalllcpvalues = (Uchar *) lcpsubtab->spaceSeqpos;
  }
  INITARRAY(&mkvauxstack,MKVstack);
  for (code = mincode; code <= maxcode; code++)
  {
    rightchar = calcbucketboundsparts(&bucketspec,
                                      bcktab,
                                      code,
                                      maxcode,
                                      totalwidth,
                                      rightchar,
                                      numofchars);
    if (outlcpinfo != NULL && code > 0)
    {
      (void) nextTurningwheel(outlcpinfo->tw);
    }
    if (bucketspec.nonspecialsinbucket > 0)
    {
      if (bucketspec.nonspecialsinbucket > 1UL)
      {
        if (lcpsubtab != NULL)
        {
          lcpsubtab->suftabbase = suftabptr + bucketspec.left;
        }
        bentleysedgewick(encseq,
                         readmode,
                         totallength,
                         &mkvauxstack,
                         suftabptr + bucketspec.left,
                         suftabptr + bucketspec.left +
                                     bucketspec.nonspecialsinbucket - 1,
                         (Seqpos) prefixlength,
                         lcpsubtab);
      }
      if (outlcpinfo != NULL)
      {
        if (previoussuffix->defined)
        {
          /* compute lcpvalue of first element of bucket with
             last element of previous bucket */
          lcpvalue = computelocallcpvalue(encseq,
                                          readmode,
                                          previoussuffix,
                                          suftabptr[bucketspec.left],
                                          esr1,
                                          esr2);
        } else
        {
          /* first part first code */
          lcpvalue = 0;
        }
        assert(lcpsubtab != NULL);
        SETLCP(0,lcpvalue);
        /* all other lcp-values are computed and they can be output */
        multilcpvalue(outlcpinfo,
                      bucketspec.nonspecialsinbucket,
                      bucketspec.left);
        /* previoussuffix becomes last nonspecial element in current bucket */
        previoussuffix->startpos
          = suftabptr[bucketspec.left + bucketspec.nonspecialsinbucket - 1];
        previoussuffix->code = code;
        previoussuffix->prefixindex = prefixlength;
        consistencyofsuffix(__LINE__,
                            encseq,readmode,bcktab,numofchars,previoussuffix);
      }
    }
    if (outlcpinfo != NULL)
    {
      if (bucketspec.specialsinbucket > 0)
      {
        minprefixindex = bucketends(encseq,
                                    readmode,
                                    esr1,
                                    esr2,
                                    outlcpinfo,
                                    previoussuffix,
                                    /* first special element in bucket */
                                    suftabptr[bucketspec.left +
                                              bucketspec.nonspecialsinbucket],
                                    bucketspec.specialsinbucket,
                                    code,
                                    numofchars,
                                    bcktab);
        /* there is at least one special element: this is the last element
           in the bucket, and thus the previoussuffix for the next round */
        previoussuffix->startpos = suftabptr[bucketspec.left +
                                             bucketspec.nonspecialsinbucket +
                                             bucketspec.specialsinbucket - 1];
        previoussuffix->defined = true;
        previoussuffix->code = code;
        previoussuffix->prefixindex = minprefixindex;
        consistencyofsuffix(__LINE__,
                            encseq,readmode,bcktab,numofchars,previoussuffix);
      } else
      {
        if (bucketspec.nonspecialsinbucket > 0)
        {
          /* if there is at least one element in the bucket, then the last
             one becomes the next previous suffix */
          previoussuffix->startpos
             = suftabptr[bucketspec.left + bucketspec.nonspecialsinbucket - 1];
          previoussuffix->defined = true;
          previoussuffix->code = code;
          previoussuffix->prefixindex = prefixlength;
          consistencyofsuffix(__LINE__,
                              encseq,readmode,bcktab,numofchars,previoussuffix);
        }
      }
    }
  }
  freeEncodedsequencescanstate(&esr1);
  freeEncodedsequencescanstate(&esr2);
  FREEARRAY(&mkvauxstack,MKVstack);
}
