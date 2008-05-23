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
#include "libgtcore/unused.h"
#include "libgtcore/strarray.h"
#include "libgtcore/ma.h"
#include "libgtcore/error.h"
#include "libgtcore/fileutils.h"
#include "libgtcore/seqiterator.h"
#include "libgtmatch/tagerator.h"
#include "libgtmatch/sarr-def.h"
#include "libgtmatch/esa-mmsearch-def.h"

#include "libgtmatch/esa-map.pr"

#define MAXTAGSIZE 64
#define UNDEFINDEX      ((short) (-1))

static void exactpatternmatching(const Encodedsequence *encseq,
                                 const Seqpos *suftab,
                                 Readmode readmode,
                                 Seqpos totallength,
                                 const Uchar *pattern,
                                 unsigned long patternlength)
{
  MMsearchiterator *mmsi;
  Seqpos dbstartpos;

  mmsi = newmmsearchiterator(encseq,
                             suftab,
                             0,  /* leftbound */
                             totallength, /* rightbound */
                             0, /* offset */
                             readmode,
                             pattern,
                             patternlength);
  while (nextmmsearchiterator(&dbstartpos,mmsi))
  {
    printf(" " FormatSeqpos,PRINTSeqposcast(dbstartpos));
  }
  printf("\n");
  freemmsearchiterator(&mmsi);
}

typedef struct
{
  unsigned long *Eqs,         /* bit vector for reverse order match */
                maxdistance,  /* distance threshold */
                patternlength;   /* pattern length */
} MyersBVinfo;

#define SKDEBUG

typedef struct
{
  unsigned long Pv,    /* the plus-vector for Myers Algorithm */
                Mv;    /* the minus-vector for Myers Algorithm */
  short maxleqk;       /* \(\max\{i\in[0,m]\mid D(i)\leq k\}\) where
                          \(m\) is the length of the pattern, \(k\) is the
                          distance threshold, and \(D\) is
                          the current distance column */
#ifdef SKDEBUG
  unsigned long scorevalue;    /* the score for the given depth */
#endif
} MyersColumn;

#ifdef APPROX

#ifdef SKDEBUG
static void verifycolumnvalues(unsigned long maxdistance,
                               unsigned long patternlength,
                               const MyersColumn *col,
                               unsigned long startscore)
{
  unsigned long idx, score = startscore, minscore, mask;
  short maxleqkindex;

  if (score <= maxdistance)
  {
    maxleqkindex = 0;
    minscore = score;
  } else
  {
    maxleqkindex = UNDEFINDEX;
    minscore = 0;
  }
  assert(patternlength <= (unsigned long) SHRT_MAX);
  for (idx=1UL, mask = 1UL; idx <= patternlength; idx++, mask <<= 1)
  {
    if (col->Pv & mask)
    {
      score++;
    } else
    {
      if (col->Mv & mask)
      {
        score--;
      }
    }
    if (score <= maxdistance)
    {
      maxleqkindex = (short) idx;
      minscore = score;
    }
  }
  if (maxleqkindex != col->maxleqk)
  {
    fprintf(stderr,"correct maxleqkindex = %hd != %hd = col->maxleqk\n",
                   maxleqkindex,
                   col->maxleqk);
    exit(EXIT_FAILURE);
  }
  if (maxleqkindex != UNDEFINDEX)
  {
    if (minscore != col->scorevalue)
    {
      fprintf(stderr,"correct score = %lu != %lu = col->score\n",
                   minscore,
                   col->scorevalue);
      exit(EXIT_FAILURE);
    }
  }
}
#endif

static void nextEDcolumn(MyersBVinfo *apminfo,
                         MyersColumn *outcol,
                         Uchar currentchar,
                         MyersColumn *incol)
{
  unsigned long Eq, Xv, Xh, Ph, Mh, /* as in Myers Paper */
                backmask;           /* only one bit is on */
  short idx;                        /* a counter */
  unsigned long score;              /* current score */

  Eq = apminfo->Eqs[(unsigned long) currentchar];
  Xv = Eq | incol->Mv;
  Xh = (((Eq & incol->Pv) + incol->Pv) ^ incol->Pv) | Eq;

  Ph = incol->Mv | ~ (Xh | incol->Pv);
  Mh = incol->Pv & Xh;

  Ph = (Ph << 1) | 1UL;
  outcol->Pv = (Mh << 1) | ~ (Xv | Ph);
  outcol->Mv = Ph & Xv;
  /* printf("incol->maxleqk %ld\n",(Showsint) incol->maxleqk); */
#ifdef SKDEBUG
  if ((unsigned long) incol->maxleqk == apminfo->patternlength)
  {
    fprintf(stderr,"incol->maxleqk = %lu = patternlength not allowed\n",
            apminfo->patternlength);
    exit(EXIT_FAILURE);
  }
  if (incol->maxleqk == UNDEFINDEX)
  {
    fprintf(stderr,"incol->maxleqk = UNDEFINDEX not allowed\n");
    exit(EXIT_FAILURE);
  }
#endif
  backmask = 1UL << incol->maxleqk;
  if (Eq & backmask || Mh & backmask)
  {
    outcol->maxleqk = incol->maxleqk + (short) 1;
#ifdef SKDEBUG
    outcol->scorevalue = incol->scorevalue;
#endif
  } else
  {
    if (Ph & backmask)
    {
      score = apminfo->maxdistance+1;
      outcol->maxleqk = UNDEFINDEX;
      for (idx = incol->maxleqk - (short) 1, backmask >>= 1;
           idx >= 0;
           idx--, backmask >>= 1)
      {
        if (outcol->Pv & backmask)
        {
          score--;
          if (score <= apminfo->maxdistance)
          {
            outcol->maxleqk = idx;
#ifdef SKDEBUG
            outcol->scorevalue = score;
#endif
            break;
          }
        } else
        {
          if (outcol->Mv & backmask)
          {
            score++;
          }
        }
      }
    } else
    {
      outcol->maxleqk = incol->maxleqk;
#ifdef SKDEBUG
      outcol->scorevalue = incol->scorevalue;
#endif
    }
  }
}

typedef struct
{
  unsigned long offset;
  Seqpos left, right;
  MyersColumn col;
} Lcpintervalwithcolumn;

static void approxpatternmatch(const Encodedsequence *encseq,
                               const Seqpos *suftab,
                               Readmode readmode,
                               Seqpos totallength,
                               const Uchar *pattern,
                               unsigned long patternlength,
                               unsigned long maxdistance)
{
  Genericstack gstack;
  Vbound *vbounds, *vboundsptr;
  Scoredvnode scvnode, scvnodenew;
  Uint l, r, vboundmaxsize, vboundscount, countprocessed = 0;

  emptygenericStack(&gstack,UintConst(1024));
  scvnode.score = (ProfScore) 0;
  scvnode.offset = 0;
  scvnode.left = 0;
  scvnode.right = virtualtree->multiseq.totallength;
  pushGenericstack(&gstack,scvnode);
  vboundmaxsize = virtualtree->alpha.mapsize-1+1;
  ALLOCASSIGNSPACE(vbounds,NULL,Vbound,vboundmaxsize);
  while (!stackisempty(&gstack))
  {
    countprocessed++;
    scvnode = popGenericstack(&gstack);
    vboundscount = splitnodewithcharbinwithoutspecial(&virtualtree->multiseq,
                                                      virtualtree->suftab,
                                                      vbounds,
                                                      vboundmaxsize,
                                                      scvnode.offset,
                                                      scvnode.left,
                                                      scvnode.right);
    for (vboundsptr = vbounds + vboundscount -1;
         vboundsptr >= vbounds;
         vboundsptr--)
    {
      l = vboundsptr->bound;
      r = (vboundsptr+1)->bound-1;
      if (ISSPECIAL(vboundsptr->inchar))
      {
        NOTSUPPOSED;
      }
      DEBUG4(2,"%lu-(%lu,%lu)%c",
                (Showuint) vboundsptr->inchar,(Showuint) l,
                (Showuint) r,(vboundsptr == vbounds + vboundscount - 1)
                              ? '\n'
                              : ' ');
      if (l == r)
      {
        (void) evaluateedge(virtualtree,
                            prof,
                            results,
                            NULL,
                            l,r,
                            scvnode.offset,
                            scvnode.score);
      } else
      {
        if (evaluateedge(virtualtree,
                        prof,
                        results,
                        &scvnodenew,
                        l,r,
                        scvnode.offset,
                        scvnode.score))
        {
          pushGenericstack(&gstack,scvnodenew);
        }
      }
    }
  }
  wrapGenericstack(&gstack);
  FREESPACE(vbounds);
}
#endif

int runtagerator(const TageratorOptions *tageratoroptions,Error *err)
{
  Suffixarray suffixarray;
  Seqpos totallength;
  SeqIterator *seqit;
  Uchar charcode;
  bool haserr = false;
  char *desc;
  int retval;
  unsigned long idx, taglen, tagnumber;
  unsigned int demand = SARR_SUFTAB | SARR_ESQTAB;
  const Uchar *symbolmap, *currenttag;
  Uchar transformedtag[MAXTAGSIZE];

  if (mapsuffixarray(&suffixarray,
                     &totallength,
                     demand,
                     tageratoroptions->indexname,
                     NULL,
                     err) != 0)
  {
    haserr = true;
  }
  symbolmap = getsymbolmapAlphabet(suffixarray.alpha);
  seqit = seqiterator_new(tageratoroptions->tagfiles, NULL, true);
  for (tagnumber = 0; /* Nothing */; tagnumber++)
  {
    retval = seqiterator_next(seqit, &currenttag, &taglen, &desc, err);
    if (retval != 1)
    {
      break;
    }
    if (taglen > (unsigned long) MAXTAGSIZE)
    {
      error_set(err,"tag of length %lu; tags must not be longer than %d",
                     taglen,MAXTAGSIZE);
      haserr = true;
      break;
    }
    for (idx = 0; idx < taglen; idx++)
    {
      charcode = symbolmap[currenttag[idx]];
      if (charcode == (Uchar) UNDEFCHAR)
      {
        error_set(err,"undefed character '%c' in tag number %lu",
                  currenttag[idx],
                  tagnumber);
        haserr = true;
        break;
      }
      transformedtag[idx] = charcode;
    }
    if (tageratoroptions->maxdistance == 0)
    {
      printf("tag %lu:",tagnumber);
      exactpatternmatching(suffixarray.encseq,
                           suffixarray.suftab,
                           suffixarray.readmode,
                           totallength,
                           transformedtag,
                           taglen);
    }
    ma_free(desc);
  }
  seqiterator_delete(seqit);
  freesuffixarray(&suffixarray);
  return haserr ? -1 : 0;
}
