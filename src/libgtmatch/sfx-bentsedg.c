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
#include "libgtcore/unused.h"
#include "divmodmul.h"
#include "spacedef.h"
#include "encseq-def.h"
#include "turnwheels.h"
#include "esafileend.h"
#include "sfx-outlcp.h"
#include "bcktab.h"

#include "sfx-cmpsuf.pr"
#include "opensfxfile.pr"
#include "kmer2string.pr"

#define COMPAREOFFSET   (UCHAR_MAX + 1)
#define UNIQUEINT(P)    ((Seqpos) ((P) + COMPAREOFFSET))
#define ACCESSCHAR(POS) getencodedchar(encseq,POS,readmode) /* XXX */
#define ISNOTEND(POS)   ((POS) < totallength && ISNOTSPECIAL(ACCESSCHAR(POS)))

#define DEREF(VAR,PTR,STOPPOS)\
        (((PTR) < (STOPPOS) && ISNOTSPECIAL(VAR = ACCESSCHAR(PTR))) ?\
        ((Seqpos) VAR) : UNIQUEINT(PTR))

#define LCPINDEX(I)       (Seqpos) ((I) - lcpsubtab->suftabbase)
#define SETLCP(I,V)       lcpsubtab->spaceSeqpos[I] = V

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

#define SMALLSIZE (Seqpos) 5

#define SUBSORT(WIDTH,BORDER,LEFT,RIGHT,DEPTH)\
        /*checksuffixrange(encseq,\
                         fwd,\
                         complement,\
                         baseptr,\
                         lcpsubtab,\
                         LEFT,\
                         RIGHT,\
                         DEPTH,\
                         __LINE__);*/\
        if (!maxdepth->defined ||\
            (DEPTH) < (Seqpos) maxdepth->valueunsignedint)\
        {\
          if ((WIDTH) <= (BORDER))\
          {\
            if ((LEFT) < (RIGHT))\
            {\
              insertionsort(encseq,esr1,esr2,\
                            lcpsubtab,readmode,totallength,\
                            LEFT,RIGHT,DEPTH,maxdepth,cmpcharbychar);\
            }\
          } else\
          {\
            PUSHMKVSTACK(LEFT,RIGHT,DEPTH);\
          }\
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

#define UPDATELCP(MINVAL,MAXVAL)\
        if ((MINVAL) > commonunits)\
        {\
          MINVAL = commonunits;\
        }\
        if ((MAXVAL) < commonunits)\
        {\
          MAXVAL = commonunits;\
        }

DECLAREARRAYSTRUCT(Largelcpvalue);

typedef struct
{
  Seqpos *spaceSeqpos;
  Uchar *smalllcpvalues;
  ArrayLargelcpvalue largelcpvalues;
  unsigned long nextfreeSeqpos, allocatedSeqpos;
  const Seqpos *suftabbase;
} Lcpsubtab;

typedef struct
{
  bool defined;
  Codetype code;
  unsigned int prefixindex;
#ifdef SKDEBUG
  Seqpos startpos;
#endif
} Suffixwithcode;

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
  Suffixwithcode previoussuffix;
  bool previousbucketwasempty;
};

typedef Seqpos Suffixptr;

#define CMPCHARBYCHARPTR2INT(VAR,TMPVAR,I)\
        VAR = (((cptr = *(I)+depth) < totallength &&\
              ISNOTSPECIAL(TMPVAR = ACCESSCHAR(cptr))) ? ((Seqpos) TMPVAR)\
                                                       : UNIQUEINT(cptr))

typedef EndofTwobitencoding Sfxcmp;

#define PTR2INT(VAR,IDXPTR)\
        {\
          Seqpos pos = *(IDXPTR);\
          if (fwd)\
          {\
            if (pos + depth < totallength)\
            {\
              pos += depth;\
              initEncodedsequencescanstategeneric(esr1,encseq,true,pos);\
              extract2bitenc(true,&VAR,encseq,esr1,pos);\
            } else\
            {\
              VAR.tbe = 0;\
              VAR.unitsnotspecial = 0;\
              VAR.position = pos;\
            }\
          } else\
          {\
            pos = REVERSEPOS(totallength,pos);\
            if (pos >= depth)\
            {\
              pos -= depth;\
              initEncodedsequencescanstategeneric(esr1,encseq,false,pos);\
              extract2bitenc(false,&VAR,encseq,esr1,pos);\
            } else\
            {\
              VAR.tbe = 0;\
              VAR.unitsnotspecial = 0;\
              VAR.position = pos;\
            }\
          }\
        }

#define Sfxdocompare(COMMONUNITS,X,Y)\
        ret##X##Y = compareTwobitencodings(fwd,complement,COMMONUNITS,&X,&Y)

#define SfxcmpEQUAL(X,Y)      (ret##X##Y == 0)
#define SfxcmpSMALLER(X,Y)    (ret##X##Y < 0)
#define SfxcmpGREATER(X,Y)    (ret##X##Y > 0)


#ifdef SKDEBUG
static void showsuffixrange(const Encodedsequence *encseq,
                            bool fwd,
                            bool complement,
                            Seqpos baseptr,
                            const Lcpsubtab *lcpsubtab,
                            const Suffixptr *leftptr,
                            const Suffixptr *rightptr,
                            Seqpos depth)
{
  const Suffixptr *pi;

  printf("of %d suffixes [%d,%d] at depth %d:\n",
         (int) ((rightptr) - (leftptr) + 1),
         (int) baseptr + LCPINDEX(leftptr),
         (int) baseptr + LCPINDEX(rightptr),
         (int) depth);
  for (pi = leftptr; pi <= rightptr; pi++)
  {
    printf("suffix %d:",*pi);
    showsequenceatstartpos(stdout,
                           fwd,
                           complement,
                           encseq,
                           *pi);
  }
}
#endif

#undef CHECKSUFFIXRANGE
#ifdef CHECKSUFFIXRANGE
static void checksuffixrange(const Encodedsequence *encseq,
                             bool fwd,
                             bool complement,
                             Seqpos *left,
                             Seqpos *right,
                             Seqpos depth,
                             int line)
{
  Seqpos *sufptr, newdepth = depth, pos1, pos2;

#ifdef SKDEBUG
  printf("checksuffixrange ");
  showsuffixrange(encseq,
                  fwd,
                  complement,
                  baseptr,
                  lcpsubtab,
                  left,
                  right,
                  depth);
#endif
  for (sufptr=left; sufptr<right; sufptr++)
  {
    if (fwd)
    {
      pos1 = *sufptr;
      pos2 = *(sufptr+1);
    } else
    {
      pos1 = REVERSEPOS(getencseqtotallength(encseq),*sufptr);
      pos2 = REVERSEPOS(getencseqtotallength(encseq),*(sufptr+1));
    }
    (void) comparetwostrings(encseq,
                             fwd,
                             complement,
                             &newdepth,
                             pos1,
                             pos2);
    if (depth > newdepth)
    {
      fprintf(stderr,"line %d: "
                     "depth=" FormatSeqpos " > " FormatSeqpos "=newdepth\n",
                     line,
                     PRINTSeqposcast(depth),
                     PRINTSeqposcast(newdepth));
      fprintf(stderr,"suffix " FormatSeqpos " vs " FormatSeqpos "\n",
                     PRINTSeqposcast(*sufptr),
                     PRINTSeqposcast(*(sufptr+1)));
      fprintf(stderr,"in range of length " FormatSeqpos "\n",
                     PRINTSeqposcast(right - left + 1));
      exit(EXIT_FAILURE);
    }
  }
}
#endif

static Suffixptr *medianof3cmpcharbychar(const Encodedsequence *encseq,
                                         Readmode readmode,
                                         Seqpos totallength,
                                         Seqpos depth,
                                         Suffixptr *a,
                                         Suffixptr *b,
                                         Suffixptr *c)
{
  Seqpos vala, valb, valc;
  Suffixptr cptr;
  Uchar tmpavar, tmpbvar;

  CMPCHARBYCHARPTR2INT(vala,tmpavar,a);
  CMPCHARBYCHARPTR2INT(valb,tmpbvar,b);
  if (vala == valb)
  {
    return a;
  }
  CMPCHARBYCHARPTR2INT(valc,tmpavar,c);
  if (vala == valc || valb == valc)
  {
    return c;
  }
  return vala < valb ?
        (valb < valc ? b : (vala < valc ? c : a))
      : (valb > valc ? b : (vala < valc ? a : c));
}

static Suffixptr *medianof3(const Encodedsequence *encseq,
                            Encodedsequencescanstate *esr1,
                            bool fwd,
                            bool complement,
                            Seqpos totallength,
                            Seqpos depth,
                            Suffixptr *a,
                            Suffixptr *b,
                            Suffixptr *c)
{
  Sfxcmp vala, valb, valc;
  int retvalavalb, retvalavalc, retvalbvalc;

  PTR2INT(vala,a);
  PTR2INT(valb,b);
  Sfxdocompare(NULL,vala,valb);
  if (SfxcmpEQUAL(vala,valb))
  {
    return a;
  }
  PTR2INT(valc,c);
  Sfxdocompare(NULL,vala,valc);
  if (SfxcmpEQUAL(vala,valc))
  {
    return c;
  }
  Sfxdocompare(NULL,valb,valc);
  if (SfxcmpEQUAL(valb,valc))
  {
    return c;
  }
  return SfxcmpSMALLER(vala,valb) ?
        (SfxcmpSMALLER(valb,valc) ? b : (SfxcmpSMALLER(vala,valc) ? c : a))
      : (SfxcmpGREATER(valb,valc) ? b : (SfxcmpSMALLER(vala,valc) ? a : c));
}

static void insertionsort(const Encodedsequence *encseq,
                          Encodedsequencescanstate *esr1,
                          Encodedsequencescanstate *esr2,
                          Lcpsubtab *lcpsubtab,
                          Readmode readmode,
                          Seqpos totallength,
                          Suffixptr *leftptr,
                          Suffixptr *rightptr,
                          Seqpos depth,
                          const Definedunsignedint *maxdepth,
                          bool cmpcharbychar)
{
  Suffixptr *pi, *pj;
  Seqpos lcpindex, lcplen = 0;
  int retval;
  Suffixptr sptr, tptr, temp;
  bool fwd = ISDIRREVERSE(readmode) ? false : true,
       complement = ISDIRCOMPLEMENT(readmode) ? true : false;

#ifdef SKDEBUG
  printf("insertion sort ");
  showsuffixrange(encseq,fwd,complement,baseptr,lcpsubtab,leftptr,rightptr,
                  depth);
#endif
  for (pi = leftptr + 1; pi <= rightptr; pi++)
  {
    for (pj = pi; pj > leftptr; pj--)
    {
      if (cmpcharbychar)
      {
        for (sptr = (*(pj-1))+depth, tptr = (*pj)+depth; /* Nothing */;
             sptr++, tptr++)
        {
          Seqpos ccs, cct;
          Uchar tmpsvar, tmptvar;
          if (maxdepth->defined)
          {
            ccs = DEREF(tmpsvar,sptr,
                        MIN(totallength,*(pj-1)+maxdepth->valueunsignedint));
            cct = DEREF(tmptvar,tptr,
                        MIN(totallength,*pj+maxdepth->valueunsignedint));
          } else
          {
            ccs = DEREF(tmpsvar,sptr,totallength);
            cct = DEREF(tmptvar,tptr,totallength);
          }
          if (ccs != cct)
          {
            lcplen = (Seqpos) (tptr - *pj);
            retval = (ccs < cct) ? -1 : 1;
            break;
          }
        }
      } else
      {
#ifdef SKDEBUG
        printf("compareEncseqsequences[%d,%d] at depth %d\n",
                       (int) *(pj-1),(int) *pj,(int) depth);
        showsequenceatstartpos(stdout,
                            fwd,
                            complement,
                            encseq,
                            *(pj-1));
        showsequenceatstartpos(stdout,
                            fwd,
                            complement,
                            encseq,
                            *pj);
#endif
        retval = compareEncseqsequences(&lcplen,encseq,fwd,complement,
                                        esr1,esr2,*(pj-1),*pj,depth);
      }
      if (lcpsubtab != NULL)
      {
        lcpindex = LCPINDEX(pj);
        if (retval > 0 && pj < pi)
        {
          SETLCP(lcpindex+1,lcpsubtab->spaceSeqpos[lcpindex]);
        }
        SETLCP(lcpindex,lcplen);
      }
      if (retval < 0)
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

static unsigned long quicksortsteps = 0;
static unsigned long quicksortdiff = 0;

static void bentleysedgewick(const Encodedsequence *encseq,
                             Encodedsequencescanstate *esr1,
                             Encodedsequencescanstate *esr2,
                             Readmode readmode,
                             Seqpos totallength,
                             ArrayMKVstack *mkvauxstack,
                             Suffixptr *l,
                             Suffixptr *r,
                             Seqpos d,
                             Lcpsubtab *lcpsubtab,
                             const Definedunsignedint *maxdepth,
                             bool cmpcharbychar)
{
  Suffixptr *left, *right, *leftplusw;
  Seqpos pivotcmpcharbychar = 0, valcmpcharbychar;
  Sfxcmp pivot, val;
  Seqpos w, depth, offset, doubleoffset, width;
  Suffixptr *pa, *pb, *pc, *pd, *pl, *pm, *pr, *aptr, *bptr, cptr, temp;
  bool fwd = ISDIRREVERSE(readmode) ? false : true,
       complement = ISDIRCOMPLEMENT(readmode) ? true : false;
  int retvalpivot;
  Uchar tmpvar;
  unsigned int commonunits, smallermaxlcp, greatermaxlcp,
               smallerminlcp, greaterminlcp;
  const int commonunitsequal = cmpcharbychar ? 1 : UNITSIN2BITENC;

  width = (Seqpos) (r - l + 1);
  if (width <= SMALLSIZE)
  {
    if (l < r)
    {
      insertionsort(encseq,esr1,esr2,lcpsubtab,readmode,totallength,
                    l,r,d,maxdepth,cmpcharbychar);
    }
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
      if (cmpcharbychar)
      {
        pl = medianof3cmpcharbychar(encseq,readmode,totallength,depth,
                       pl,pl+offset,pl+doubleoffset);
        pm = medianof3cmpcharbychar(encseq,readmode,totallength,depth,
                       pm-offset,pm,pm+offset);
        pr = medianof3cmpcharbychar(encseq,readmode,totallength,depth,
                       pr-doubleoffset,pr-offset,pr);
      } else
      {
        pl = medianof3(encseq,esr1,fwd,complement,totallength,depth,
                       pl,pl+offset,pl+doubleoffset);
        pm = medianof3(encseq,esr1,fwd,complement,totallength,depth,
                       pm-offset,pm,pm+offset);
        pr = medianof3(encseq,esr1,fwd,complement,totallength,depth,
                       pr-doubleoffset,pr-offset,pr);
      }
    }
    if (cmpcharbychar)
    {
      pm = medianof3cmpcharbychar(encseq,readmode,totallength,depth,pl,pm,pr);
      SWAP(left, pm);
      CMPCHARBYCHARPTR2INT(pivotcmpcharbychar,tmpvar,left);
    } else
    {
      pm = medianof3(encseq,esr1,fwd,complement,totallength,depth,pl,pm,pr);
      SWAP(left, pm);
      PTR2INT(pivot,left);
    }
    /* now pivot element is at index left */
    /* all elements to be compared are between pb and pc */
    /* pa is the position at which the next element smaller than the
       pivot element is inserted at */
    /* pd is the position at which the next element greater than the
       pivot element is inserted at */
    pa = pb = left + 1;
    pc = pd = right;
    if (cmpcharbychar)
    {
      smallerminlcp = greaterminlcp = smallermaxlcp = greatermaxlcp = 0;
      for (;;)
      {
        while (pb <= pc)
        {
          CMPCHARBYCHARPTR2INT(valcmpcharbychar,tmpvar,pb);
          if (valcmpcharbychar > pivotcmpcharbychar)
          {
            break;
          }
          if (valcmpcharbychar == pivotcmpcharbychar)
          {
            SWAP(pa, pb);
            pa++;
          }
          pb++;
        }
        while (pb <= pc)
        {
          CMPCHARBYCHARPTR2INT(valcmpcharbychar,tmpvar,pc);
          if (valcmpcharbychar < pivotcmpcharbychar)
          { /* stop for elements < pivot */
            break;
          }
          if (valcmpcharbychar == pivotcmpcharbychar)
          {
            SWAP(pc, pd);  /* exchange equal element and element at index pd */
            pd--;
          }
          pc--;
        }
        if (pb > pc)
        { /* no elements to compare to pivot */
          break;
        }
        SWAP(pb, pc);
        pb++;
        pc--;
      }
    } else
    {
      smallermaxlcp = greatermaxlcp = 0;
      smallerminlcp = greaterminlcp = (unsigned int) UNITSIN2BITENC;
      for (;;)
      {
        /* look for elements identical or smaller than pivot from left */
        while (pb <= pc)
        {
          PTR2INT(val,pb);
          Sfxdocompare(&commonunits,val,pivot);
          if (SfxcmpGREATER(val,pivot))
          { /* stop for elements val > pivot */
            UPDATELCP(greaterminlcp,greatermaxlcp);
            break;
          }
          if (SfxcmpEQUAL(val,pivot))
          {
            SWAP(pa, pb); /* exchange equal element and element at index pa */
            pa++;
          } else /* smaller */
          {
            UPDATELCP(smallerminlcp,smallermaxlcp);
          }
          pb++;
        }
        /* look for elements identical or greater than pivot from right */
        while (pb <= pc)
        {
          PTR2INT(val,pc);
          Sfxdocompare(&commonunits,val,pivot);
          if (SfxcmpSMALLER(val,pivot))
          { /* stop for elements val < pivot */
            UPDATELCP(smallerminlcp,smallermaxlcp);
            break;
          }
          if (SfxcmpEQUAL(val,pivot))
          {
            SWAP(pc, pd); /* exchange equal element and element at index pa */
            pd--;
          } else /* greater */
          {
            UPDATELCP(greaterminlcp,greatermaxlcp);
          }
          pc--;
        }
        if (pb > pc)
        { /* interval is empty small */
          break;
        }
        SWAP(pb, pc);
        pb++;
        pc--;
      }
    }
    assert(pa >= left);
    assert(pb >= pa);
    w = MIN((Seqpos) (pa-left),(Seqpos) (pb-pa));
    /* move w elements at the left to the middle */
    VECSWAP(left,  pb-w, w);

    assert(pd >= pc);
    assert(right >= pd);
    w = MIN((Seqpos) (pd-pc), (Seqpos) (right-pd));
    /* move w elements at the right to the middle */
    VECSWAP(pb, right+1-w, w);

    /* all elements equal to the pivot are now in the middle namely in the
       range [left + (pb-pa) and right - (pd-pc)] */
    /* hence we have to sort the elements in the intervals
       [left..left+(pb-pa)-1] and
       [right-(pd-pc)+1..right] */

    assert(pb >= pa);
    if ((w = (Seqpos) (pb-pa)) > 0)
    {
      leftplusw = left + w;
      if (lcpsubtab != NULL)
      {
        SETLCP(LCPINDEX(leftplusw),depth + smallermaxlcp);
      }
      SUBSORT(w,SMALLSIZE,left,leftplusw-1,depth + smallerminlcp);
    } else
    {
      leftplusw = left;
    }

    cptr = *leftplusw + depth;
    if (ISNOTEND(cptr))
    {
      width = (Seqpos) (right-(pd-pb)-leftplusw);
      SUBSORT(width,SMALLSIZE,leftplusw,right-(pd-pb)-1,depth+commonunitsequal);
    }

    assert(pd >= pc);
    if ((w = (Seqpos) (pd-pc)) > 0)
    {
      if (lcpsubtab != NULL)
      {
        SETLCP(LCPINDEX(right-w+1),depth + greatermaxlcp);
      }
      SUBSORT(w,SMALLSIZE,right-w+1,right,depth + greaterminlcp);
    }
    quicksortsteps++;
    if ((unsigned long) (pb-pa) < (unsigned long) (pd-pc))
    {
      quicksortdiff += (unsigned long) (pd-pc) - (unsigned long) (pb-pa);
    } else
    {
      quicksortdiff += (unsigned long) (pb-pa) - (unsigned long) (pd-pc);
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
  unsigned int rightchar = (unsigned int) (mincode % numofchars);
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

#ifdef SKDEBUG
static void showSuffixwithcode(FILE *fp,const Suffixwithcode *suffix)
{
  char buffer[18+1];

  kmercode2string(buffer,
                  suffix->code,
                  4,
                  8,
                  "acgt");
  fprintf(fp,"(startpos=%lu,code=%u,prefixindex=%u,\"%s\")",
              (unsigned long) suffix->startpos,
              suffix->code,suffix->prefixindex,
              buffer);
}
#endif

#ifdef SKDEBUG
static Seqpos computelocallcpvalue(const Encodedsequence *encseq,
                                   Readmode readmode,
                                   const Suffixwithcode *previoussuffix,
                                   const Suffixwithcode *currentsuffix,
                                   unsigned int minchanged,
                                   Encodedsequencescanstate *esr1,
                                   Encodedsequencescanstate *esr2)
{
  Seqpos lcpvalue;
  unsigned int lcpvalue2;
  int cmp;

  cmp = comparetwosuffixes(encseq,
                           readmode,
                           &lcpvalue,
                           false,
                           false,
                           0,
                           previoussuffix->startpos,
                           currentsuffix->startpos,
                           esr1,
                           esr2);
  if (cmp > 0)
  {
    fprintf(stderr,"cmp " FormatSeqpos
            " " FormatSeqpos " = %d, lcpval=" FormatSeqpos "\n",
            PRINTSeqposcast(previoussuffix->startpos),
            PRINTSeqposcast(currentsuffix->startpos),
            cmp,
            PRINTSeqposcast(lcpvalue));
    exit(EXIT_FAILURE);
  }
  if (previoussuffix->code == currentsuffix->code)
  {
    assert(lcpvalue == MIN(previoussuffix->prefixindex,
                           currentsuffix->prefixindex));
  } else
  {
    assert(previoussuffix->code < currentsuffix->code);
    lcpvalue2 = MIN(minchanged,MIN(previoussuffix->prefixindex,
                                   currentsuffix->prefixindex));
    if (lcpvalue != lcpvalue2)
    {
      fprintf(stderr,"lcpvalue = %lu != %u = lcpvalue2\n",
              (unsigned long) lcpvalue,
              lcpvalue2);
      fprintf(stderr,"previoussuffix=");
      showSuffixwithcode(stderr,previoussuffix);
      fprintf(stderr,"\ncurrentsuffix=");
      showSuffixwithcode(stderr,currentsuffix);
      fprintf(stderr,"\n");
      exit(EXIT_FAILURE);
    }
  }
  return lcpvalue;
}
#endif

static Seqpos computelocallcpvalue(const Suffixwithcode *previoussuffix,
                                   const Suffixwithcode *currentsuffix,
                                   unsigned int minchanged)
{
  unsigned int lcpvalue;

  if (previoussuffix->code == currentsuffix->code)
  {
    lcpvalue = MIN(previoussuffix->prefixindex,
                   currentsuffix->prefixindex);
  } else
  {
    assert(previoussuffix->code < currentsuffix->code);
    lcpvalue = MIN(minchanged,MIN(previoussuffix->prefixindex,
                                  currentsuffix->prefixindex));
  }
  return (Seqpos) lcpvalue;
}

static unsigned int bucketends(Outlcpinfo *outlcpinfo,
                               Suffixwithcode *previoussuffix,
                               UNUSED Seqpos firstspecialsuffix,
                               unsigned int minchanged,
                               unsigned long specialsinbucket,
                               Codetype code,
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
  firstspecialsuffixwithcode.prefixindex = maxprefixindex;
#ifdef SKDEBUG
  firstspecialsuffixwithcode.startpos = firstspecialsuffix;
  consistencyofsuffix(__LINE__,
                      encseq,readmode,bcktab,numofchars,
                      &firstspecialsuffixwithcode);
#endif
  lcpvalue = computelocallcpvalue(previoussuffix,
                                  &firstspecialsuffixwithcode,
                                  minchanged);
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

Outlcpinfo *newlcpoutinfo(const Str *indexname,
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
#ifdef SKDEBUG
  outlcpinfo->previoussuffix.startpos = 0;
#endif
  outlcpinfo->previoussuffix.code = 0;
  outlcpinfo->previoussuffix.prefixindex = 0;
  outlcpinfo->previoussuffix.defined = false;
  outlcpinfo->previousbucketwasempty = false;
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
                    const Bcktab *bcktab,
                    unsigned int numofchars,
                    unsigned int prefixlength,
                    Outlcpinfo *outlcpinfo,
                    const Definedunsignedint *maxdepth,
                    bool cmpcharbychar,
                    unsigned long long *bucketiterstep)
{
  Codetype code;
  unsigned int rightchar = (unsigned int) (mincode % numofchars),
               minprefixindex,
               minchanged = 0;
  Seqpos totallength = getencseqtotallength(encseq);
  ArrayMKVstack mkvauxstack;
  Bucketspecification bucketspec;
  unsigned long maxbucketsize;
  Seqpos lcpvalue;
  Lcpsubtab *lcpsubtab;
  Suffixwithcode firstsuffixofbucket;
  Seqpos baseptr;
  Encodedsequencescanstate *esr1 = NULL,
                           *esr2 = NULL;

  if (outlcpinfo == NULL)
  {
    lcpsubtab = NULL;
  } else
  {
    lcpsubtab = &outlcpinfo->lcpsubtab;
  }
  if (!cmpcharbychar && hasspecialranges(encseq))
  {
    esr1 = newEncodedsequencescanstate();
    esr2 = newEncodedsequencescanstate();
  }
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
    (*bucketiterstep)++;
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
      if (outlcpinfo->previousbucketwasempty)
      {
        minchanged = MIN(minchanged,minchangedTurningwheel(outlcpinfo->tw));
      } else
      {
        minchanged = minchangedTurningwheel(outlcpinfo->tw);
      }
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
                         esr1,
                         esr2,
                         readmode,
                         totallength,
                         &mkvauxstack,
                         suftabptr + bucketspec.left,
                         suftabptr + bucketspec.left +
                                     bucketspec.nonspecialsinbucket - 1,
                         (Seqpos) prefixlength,
                         lcpsubtab,
                         maxdepth,
                         cmpcharbychar);
      }
      if (outlcpinfo != NULL)
      {
        if (outlcpinfo->previoussuffix.defined)
        {
          /* compute lcpvalue of first element of bucket with
             last element of previous bucket */
          firstsuffixofbucket.code = code;
          firstsuffixofbucket.prefixindex = prefixlength;
#ifdef SKDEBUG
          firstsuffixofbucket.startpos = suftabptr[bucketspec.left];
          consistencyofsuffix(__LINE__,
                              encseq,readmode,bcktab,numofchars,
                              &firstsuffixofbucket);
#endif
          lcpvalue = computelocallcpvalue(&outlcpinfo->previoussuffix,
                                          &firstsuffixofbucket,
                                          minchanged);
        } else
        {
          /* first part first code */
          lcpvalue = 0;
        }
        assert(lcpsubtab != NULL);
        baseptr = bucketspec.left;
        SETLCP(0,lcpvalue);
        /* all other lcp-values are computed and they can be output */
        multilcpvalue(outlcpinfo,
                      bucketspec.nonspecialsinbucket,
                      bucketspec.left);
        /* previoussuffix becomes last nonspecial element in current bucket */
        outlcpinfo->previoussuffix.code = code;
        outlcpinfo->previoussuffix.prefixindex = prefixlength;
#ifdef SKDEBUG
        outlcpinfo->previoussuffix.startpos
          = suftabptr[bucketspec.left + bucketspec.nonspecialsinbucket - 1];
        consistencyofsuffix(__LINE__,
                            encseq,readmode,bcktab,numofchars,
                            &outlcpinfo->previoussuffix);
#endif
      }
    }
    if (outlcpinfo != NULL)
    {
      if (bucketspec.specialsinbucket > 0)
      {
        minprefixindex = bucketends(outlcpinfo,
                                    &outlcpinfo->previoussuffix,
                                    /* first special element in bucket */
                                    suftabptr[bucketspec.left +
                                              bucketspec.nonspecialsinbucket],
                                    minchanged,
                                    bucketspec.specialsinbucket,
                                    code,
                                    bcktab);
        /* there is at least one special element: this is the last element
           in the bucket, and thus the previoussuffix for the next round */
        outlcpinfo->previoussuffix.defined = true;
        outlcpinfo->previoussuffix.code = code;
        outlcpinfo->previoussuffix.prefixindex = minprefixindex;
#ifdef SKDEBUG
        outlcpinfo->previoussuffix.startpos
           = suftabptr[bucketspec.left + bucketspec.nonspecialsinbucket +
                                         bucketspec.specialsinbucket - 1];
        consistencyofsuffix(__LINE__,
                            encseq,readmode,bcktab,numofchars,
                            &outlcpinfo->previoussuffix);
#endif
      } else
      {
        if (bucketspec.nonspecialsinbucket > 0)
        {
          /* if there is at least one element in the bucket, then the last
             one becomes the next previous suffix */
          outlcpinfo->previoussuffix.defined = true;
          outlcpinfo->previoussuffix.code = code;
          outlcpinfo->previoussuffix.prefixindex = prefixlength;
#ifdef SKDEBUG
          outlcpinfo->previoussuffix.startpos
            = suftabptr[bucketspec.left + bucketspec.nonspecialsinbucket - 1];
          consistencyofsuffix(__LINE__,
                              encseq,readmode,bcktab,numofchars,
                              &outlcpinfo->previoussuffix);
#endif
        }
      }
    }
    if (outlcpinfo != NULL)
    {
      if (bucketspec.nonspecialsinbucket + bucketspec.specialsinbucket == 0)
      {
        outlcpinfo->previousbucketwasempty = true;
      } else
      {
        outlcpinfo->previousbucketwasempty = false;
      }
    }
  }
  if (!cmpcharbychar && hasspecialranges(encseq))
  {
    assert(esr1 != NULL);
    freeEncodedsequencescanstate(&esr1);
    assert(esr2 != NULL);
    freeEncodedsequencescanstate(&esr2);
  }
  FREEARRAY(&mkvauxstack,MKVstack);
  printf("# quicksortsteps: %lu, avg diff %.2f\n",
          quicksortsteps,(double) quicksortdiff/quicksortsteps);

}
