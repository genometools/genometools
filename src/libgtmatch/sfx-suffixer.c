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

#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <limits.h>
#include <errno.h>
#include "libgtcore/arraydef.h"
#include "libgtcore/error.h"
#include "spacedef.h"
#include "intbits.h"
#include "divmodmul.h"
#include "measure-time-if.h"
#include "intcode-def.h"
#include "encseq-def.h"
#include "safecast-gen.h"
#include "sfx-codespec.h"
#include "sfx-partssuf-def.h"
#include "sfx-suffixer.h"
#include "sfx-outlcp.h"
#include "stamp.h"

#include "sfx-mappedstr.pr"

#define CODEBITS        (32-PREFIXLENBITS)
#define MAXPREFIXLENGTH ((1U << PREFIXLENBITS) - 1)
#define MAXCODEVALUE    ((1U << CODEBITS) - 1)

typedef struct
{
  unsigned int maxprefixlen:PREFIXLENBITS;
  unsigned int code:CODEBITS;
  Seqpos position; /* get rid of this by using information from encseq */
} Codeatposition;

DECLAREARRAYSTRUCT(Codeatposition);

#define LONGOUTPUT
#undef LONGOUTPUT

DECLAREARRAYSTRUCT(Seqpos);

 struct Sfxiterator
{
  bool storespecials;
  Codetype currentmincode,
           currentmaxcode;
  unsigned int *filltable,
               *basepower;
  Seqpos specialcharacters,
         *leftborder, /* export this for bcktab */
         *countspecialcodes, /* export this for bcktab */
         *suftab,
         *suftabptr;
  unsigned long nextfreeCodeatposition;
  Codeatposition *spaceCodeatposition;
  Suftabparts *suftabparts;
  const Encodedsequence *encseq;
  Readmode readmode;
  Seqpos widthofpart,
         totallength;
  Outlcpinfo *outlcpinfo;
  unsigned int part,
               numofchars,
               prefixlength;
  ArraySeqpos fusp;
  Specialrangeiterator *sri;
  Sequencerange overhang;
  Seqpos previoussuffix;
  bool exhausted;
};

static int initbasepower(unsigned int **basepower,
                         unsigned int **filltable,
                         unsigned int base,
                         unsigned int len,
                         Error *err)
{
  unsigned int thepower = 1U, i, minfailure;
  bool haserr = false;

  error_check(err);
  ALLOCASSIGNSPACE(*basepower,NULL,unsigned int,len+1);
  ALLOCASSIGNSPACE(*filltable,NULL,unsigned int,len);
  minfailure = UINT_MAX/base;
  for (i=0; /* Nothing */; i++)
  {
    (*basepower)[i] = thepower;
    if (i == len)
    {
      break;
    }
    if (thepower >= minfailure)
    {
      error_set(err,"overflow when computing %u * %u",thepower,base);
      haserr = true;
      break;
    }
    thepower *= base;
  }
  if (!haserr)
  {
    for (i=0; i<len; i++)
    {
      (*filltable)[i] = (*basepower)[len-i]-1;
    }
  }
  if (haserr)
  {
    FREESPACE(*basepower);
    FREESPACE(*filltable);
    return -1;
  }
  return 0;
}

static void updatekmercount(void *processinfo,
                            Codetype code,
                            Seqpos position,
                            const Firstspecialpos *firstspecial)
{
  Sfxiterator *sfi = (Sfxiterator *) processinfo;

  if (firstspecial->defined)
  {
    if (sfi->storespecials)
    {
      if (firstspecial->specialpos > 0)
      {
        Codeatposition *cp;

        cp = sfi->spaceCodeatposition + sfi->nextfreeCodeatposition++;
        cp->code = code;
        cp->maxprefixlen = firstspecial->specialpos;
        cp->position = position + firstspecial->specialpos;
        sfi->storespecials = false;
        sfi->leftborder[code]++;
      }
    } else
    {
      if (firstspecial->specialpos > 0)
      {
        sfi->leftborder[code]++;
      } else
      {
        sfi->storespecials = true;
      }
    }
  } else
  {
    sfi->leftborder[code]++;
  }
}

static void insertwithoutspecial(void *processinfo,
                                 Codetype code,
                                 Seqpos position,
                                 const Firstspecialpos *firstspecial)
{
  if (!firstspecial->defined)
  {
    Sfxiterator *sfi = (Sfxiterator *) processinfo;

    if (code >= sfi->currentmincode && code <= sfi->currentmaxcode)
    {
      Seqpos stidx;

      stidx = --sfi->leftborder[code];
#ifdef LONGOUTPUT
      printf("insert suffix " FormatSeqpos " at location " FormatSeqpos "\n",
              PRINTSeqposcast(position),
              PRINTSeqposcast(stidx));
#endif
      sfi->suftabptr[stidx] = position;
    }
  }
}

static void reversespecialcodes(Codeatposition *spaceCodeatposition,
                                unsigned long nextfreeCodeatposition)
{
  Codeatposition *front, *back, tmp;

  for (front = spaceCodeatposition,
       back = spaceCodeatposition + nextfreeCodeatposition - 1;
       front < back; front++, back--)
  {
    tmp = *front;
    *front = *back;
    *back = tmp;
  }
}

static Codetype codedownscale(const unsigned int *filltable,
                              const unsigned int *basepower,
                              Codetype code,
                              unsigned int prefixindex,
                              unsigned int maxprefixlen)
{
  unsigned int remain;

  code -= filltable[maxprefixlen];
  remain = maxprefixlen-prefixindex;
  code %= (filltable[remain]+1);
  code *= basepower[remain];
  code += filltable[prefixindex];
  return code;
}

static void derivespecialcodes(Sfxiterator *sfi,bool deletevalues)
{
  Codetype code;
  unsigned int prefixindex;
  unsigned long insertindex, j;
  Seqpos stidx;

  for (prefixindex=0; prefixindex < sfi->prefixlength; prefixindex++)
  {
    for (j=0, insertindex = 0; j < sfi->nextfreeCodeatposition; j++)
    {
      if (prefixindex <= sfi->spaceCodeatposition[j].maxprefixlen)
      {
        code = codedownscale(sfi->filltable,
                             sfi->basepower,
                             sfi->spaceCodeatposition[j].code,
                             prefixindex,
                             sfi->spaceCodeatposition[j].maxprefixlen);
        if (code >= sfi->currentmincode && code <= sfi->currentmaxcode &&
            (prefixindex > 0 || code != sfi->filltable[0]))
        {
          sfi->countspecialcodes[FROMCODE2SPECIALCODE(code,sfi->numofchars)]++;
          stidx = --sfi->leftborder[code];
#ifdef LONGOUTPUT
          printf("insert special_suffix " FormatSeqpos
                 " (code %u) at location " FormatSeqpos "\n",
                 PRINTSeqposcast(sfi->spaceCodeatposition[j].position -
                                 prefixindex),
                 (unsigned int) code,
                 PRINTSeqposcast(stidx));
#endif
          sfi->suftabptr[stidx] = sfi->spaceCodeatposition[j].position -
                                  prefixindex;
        }
      }
      if (deletevalues)
      {
        if (prefixindex < sfi->prefixlength - 1 &&
            prefixindex < sfi->spaceCodeatposition[j].maxprefixlen)
        {
          if (insertindex < j)
          {
            sfi->spaceCodeatposition[insertindex] =
              sfi->spaceCodeatposition[j];
          }
          insertindex++;
        }
      }
    }
    if (deletevalues)
    {
      sfi->nextfreeCodeatposition = insertindex;
    }
  }
}

void freeSfxiterator(Sfxiterator **sfi)
{
  Codetype specialcode;

  specialcode = FROMCODE2SPECIALCODE((*sfi)->filltable[0],
                                     (*sfi)->numofchars);
  (*sfi)->countspecialcodes[specialcode] += ((*sfi)->specialcharacters + 1);
  if ((*sfi)->sri != NULL)
  {
    freespecialrangeiterator(&(*sfi)->sri);
  }
  FREESPACE((*sfi)->spaceCodeatposition);
  FREESPACE((*sfi)->filltable);
  FREESPACE((*sfi)->basepower);
  FREESPACE((*sfi)->leftborder);
  FREESPACE((*sfi)->countspecialcodes);
  FREESPACE((*sfi)->suftab);
  freesuftabparts((*sfi)->suftabparts);
  FREESPACE(*sfi);
}

 DECLARESAFECASTFUNCTION(Seqpos,Seqpos,unsigned long,unsigned_long)

Sfxiterator *newSfxiterator(Seqpos specialcharacters,
                            Seqpos specialranges,
                            const Encodedsequence *encseq,
                            Readmode readmode,
                            unsigned int numofchars,
                            unsigned int prefixlength,
                            unsigned int numofparts,
                            Outlcpinfo *outlcpinfo,
                            Measuretime *mtime,
                            Verboseinfo *verboseinfo,
                            Error *err)
{
  Sfxiterator *sfi = NULL;
  unsigned int numofallcodes = 0, numofspecialcodes;
  Seqpos *optr;
  bool haserr = false;

  error_check(err);
  if (prefixlength == 0 || prefixlength > MAXPREFIXLENGTH)
  {
    error_set(err,"argument for option -pl must be in the range [1,%u]",
                  MAXPREFIXLENGTH);
    haserr = true;
  } else
  {
    ALLOCASSIGNSPACE(sfi,NULL,Sfxiterator,1);
    ALLOCASSIGNSPACE(sfi->spaceCodeatposition,NULL,
                     Codeatposition,specialranges+1);
    sfi->nextfreeCodeatposition = 0;
    sfi->filltable = NULL;
    sfi->basepower = NULL;
    sfi->leftborder = NULL;
    sfi->countspecialcodes = NULL;
    sfi->suftab = NULL;
    sfi->suftabptr = NULL;
    sfi->suftabparts = NULL;
    sfi->encseq = encseq;
    sfi->readmode = readmode;
    sfi->numofchars = numofchars;
    sfi->prefixlength = prefixlength;
    sfi->totallength = getencseqtotallength(encseq);
    sfi->specialcharacters = specialcharacters;
    sfi->previoussuffix = 0;
    sfi->outlcpinfo = outlcpinfo;
    sfi->sri = NULL;
    sfi->part = 0;
    sfi->exhausted = false;
  }
  if (!haserr)
  {
    assert(sfi != NULL);
    if (initbasepower(&sfi->basepower,
                      &sfi->filltable,
                      numofchars,
                      prefixlength,
                      err) != 0)
    {
      haserr = true;
    }
  }
  if (!haserr)
  {
    assert(sfi != NULL);
    assert(sfi->basepower != NULL);
    numofallcodes = sfi->basepower[prefixlength];
    if (numofallcodes-1 > MAXCODEVALUE)
    {
      error_set(err,"alphasize^prefixlength-1 = %u does not fit into "
                    " %u bits: choose smaller value for prefixlength",
                    numofallcodes-1,
                    (unsigned int) CODEBITS);
      haserr = true;
    }
  }
  if (!haserr)
  {
    assert(sfi != NULL);
    assert(sfi->basepower != NULL);
    numofspecialcodes = sfi->basepower[prefixlength-1];
    ALLOCASSIGNSPACE(sfi->leftborder,NULL,Seqpos,numofallcodes+1);
    memset(sfi->leftborder,0,
           sizeof (*sfi->leftborder) * (size_t) numofallcodes);
    ALLOCASSIGNSPACE(sfi->countspecialcodes,NULL,Seqpos,numofspecialcodes);
    memset(sfi->countspecialcodes,0,
           sizeof (*sfi->countspecialcodes) *
                  (size_t) numofspecialcodes);
    sfi->storespecials = true;
    if (mtime != NULL)
    {
      deliverthetime(stdout,mtime,"counting prefix distribution");
    }
    getencseqkmers(encseq,
                   readmode,
                   updatekmercount,
                   sfi,
                   numofchars,
                   prefixlength,
                   err);
    assert(specialranges+1 >= (Seqpos) sfi->nextfreeCodeatposition);
    assert(sfi->filltable != NULL);
    assert(sfi->leftborder != NULL);
    /* printf("leftborder[0]=%u\n",sfi.leftborder[0]); */
    for (optr = sfi->leftborder + 1;
         optr < sfi->leftborder + numofallcodes; optr++)
    {
      /*/ printf("leftborder[%u]=%u\n",(unsigned int) (optr - sfi->leftborder),
                                   *optr); */
      *optr += *(optr-1);
    }
    sfi->leftborder[numofallcodes] = sfi->totallength - specialcharacters;
    sfi->suftabparts = newsuftabparts(numofparts,
                                      sfi->leftborder,
                                      numofallcodes,
                                      sfi->totallength - specialcharacters,
                                      specialcharacters + 1,
                                      verboseinfo);
    assert(sfi->suftabparts != NULL);
    ALLOCASSIGNSPACE(sfi->suftab,NULL,Seqpos,
                     stpgetlargestwidth(sfi->suftabparts));
    reversespecialcodes(sfi->spaceCodeatposition,sfi->nextfreeCodeatposition);
    if (hasspecialranges(sfi->encseq))
    {
      sfi->sri = newspecialrangeiterator(sfi->encseq,
                                         ISDIRREVERSE(sfi->readmode)
                                           ? false : true);
    } else
    {
      sfi->sri = NULL;
    }
    sfi->fusp.spaceSeqpos = sfi->suftab;
    sfi->fusp.allocatedSeqpos
      = CALLCASTFUNC(Seqpos,unsigned_long,
                     stpgetlargestwidth(sfi->suftabparts));
    sfi->overhang.leftpos = sfi->overhang.rightpos = 0;
  }
  if (haserr)
  {
    if (sfi != NULL)
    {
      freeSfxiterator(&sfi);
    }
    return NULL;
  }
  return sfi;
}

static void preparethispart(Sfxiterator *sfi,
                            Measuretime *mtime,
                            Error *err)
{
  Seqpos totalwidth;

  sfi->currentmincode = stpgetcurrentmincode(sfi->part,sfi->suftabparts);
  sfi->currentmaxcode = stpgetcurrentmaxcode(sfi->part,sfi->suftabparts);
  sfi->widthofpart = stpgetcurrentwidthofpart(sfi->part,sfi->suftabparts);
  sfi->suftabptr = sfi->suftab -
                   stpgetcurrentsuftaboffset(sfi->part,sfi->suftabparts);
  derivespecialcodes(sfi,
                     (stpgetnumofparts(sfi->suftabparts) == 1U)
                       ? true : false);
  if (mtime != NULL)
  {
    deliverthetime(stdout,mtime,"inserting suffixes into buckets");
  }
  getencseqkmers(sfi->encseq,
                 sfi->readmode,
                 insertwithoutspecial,
                 sfi,
                 sfi->numofchars,
                 sfi->prefixlength,
                 err);
  if (mtime != NULL)
  {
    deliverthetime(stdout,mtime,"sorting the buckets");
  }
  totalwidth = stpgetcurrentsumofwdith(sfi->part,sfi->suftabparts);
  sortallbuckets(sfi->suftabptr,
                 sfi->encseq,
                 sfi->readmode,
                 sfi->leftborder,
                 sfi->countspecialcodes,
                 sfi->numofchars,
                 sfi->prefixlength,
                 sfi->currentmincode,
                 sfi->currentmaxcode,
                 totalwidth,
                 sfi->previoussuffix,
                 sfi->outlcpinfo);
  assert(totalwidth > 0);
  sfi->previoussuffix = sfi->suftab[sfi->widthofpart-1];
  sfi->part++;
}

static void insertfullspecialrange(Sfxiterator *sfi,
                                   Seqpos leftpos,
                                   Seqpos rightpos)
{
  Seqpos pos;

  assert(leftpos < rightpos);
  if (ISDIRREVERSE(sfi->readmode))
  {
    pos = rightpos - 1;
  } else
  {
    pos = leftpos;
  }
  while (true)
  {
    if (ISDIRREVERSE(sfi->readmode))
    {
      sfi->fusp.spaceSeqpos[sfi->fusp.nextfreeSeqpos++]
        = REVERSEPOS(sfi->totallength,pos);
      if (pos == leftpos)
      {
        break;
      }
      pos--;
    } else
    {
      sfi->fusp.spaceSeqpos[sfi->fusp.nextfreeSeqpos++] = pos;
      if (pos == rightpos-1)
      {
        break;
      }
      pos++;
    }
  }
}

static void fillspecialnextpage(Sfxiterator *sfi)
{
  Sequencerange range;
  Seqpos width;

  while (true)
  {
    if (sfi->overhang.leftpos < sfi->overhang.rightpos)
    {
      width = sfi->overhang.rightpos - sfi->overhang.leftpos;
      if (sfi->fusp.nextfreeSeqpos + width > sfi->fusp.allocatedSeqpos)
      {
        /* does not fit into the buffer, so only output a part */
        unsigned long rest = sfi->fusp.nextfreeSeqpos +
                             width - sfi->fusp.allocatedSeqpos;
        assert(rest > 0);
        if (ISDIRREVERSE(sfi->readmode))
        {
          insertfullspecialrange(sfi,sfi->overhang.leftpos + rest,
                                 sfi->overhang.rightpos);
          sfi->overhang.rightpos = sfi->overhang.leftpos + rest;
        } else
        {
          insertfullspecialrange(sfi,sfi->overhang.leftpos,
                                     sfi->overhang.rightpos - rest);
          sfi->overhang.leftpos = sfi->overhang.rightpos - rest;
        }
        break;
      }
      if (sfi->fusp.nextfreeSeqpos + width == sfi->fusp.allocatedSeqpos)
      { /* overhang fits into the buffer and buffer is full */
        insertfullspecialrange(sfi,sfi->overhang.leftpos,
                               sfi->overhang.rightpos);
        sfi->overhang.leftpos = sfi->overhang.rightpos = 0;
        break;
      }
      /* overhang fits into the buffer and buffer is not full */
      insertfullspecialrange(sfi,sfi->overhang.leftpos,
                             sfi->overhang.rightpos);
      sfi->overhang.leftpos = sfi->overhang.rightpos = 0;
    } else
    {
      if (sfi->sri != NULL && nextspecialrangeiterator(&range,sfi->sri))
      {
        width = range.rightpos - range.leftpos;
        assert(width > 0);
        if (sfi->fusp.nextfreeSeqpos + width > sfi->fusp.allocatedSeqpos)
        { /* does not fit into the buffer, so only output a part */
          unsigned long rest = sfi->fusp.nextfreeSeqpos +
                               width - sfi->fusp.allocatedSeqpos;
          if (ISDIRREVERSE(sfi->readmode))
          {
            insertfullspecialrange(sfi,range.leftpos + rest,
                                   range.rightpos);
            sfi->overhang.leftpos = range.leftpos;
            sfi->overhang.rightpos = range.leftpos + rest;
          } else
          {
            insertfullspecialrange(sfi,range.leftpos,range.rightpos - rest);
            sfi->overhang.leftpos = range.rightpos - rest;
            sfi->overhang.rightpos = range.rightpos;
          }
          break;
        }
        if (sfi->fusp.nextfreeSeqpos + width == sfi->fusp.allocatedSeqpos)
        { /* overhang fits into the buffer and buffer is full */
          insertfullspecialrange(sfi,range.leftpos,range.rightpos);
          sfi->overhang.leftpos = sfi->overhang.rightpos = 0;
          break;
        }
        insertfullspecialrange(sfi,range.leftpos,range.rightpos);
        sfi->overhang.leftpos = sfi->overhang.rightpos = 0;
      } else
      {
        if (sfi->fusp.nextfreeSeqpos < sfi->fusp.allocatedSeqpos)
        {
          sfi->fusp.spaceSeqpos[sfi->fusp.nextfreeSeqpos++] = sfi->totallength;
          sfi->exhausted = true;
        }
        break;
      }
    }
  }
}

const Seqpos *nextSfxiterator(Seqpos *numberofsuffixes,bool *specialsuffixes,
                              Measuretime *mtime,Sfxiterator *sfi,Error *err)
{
  error_check(err);
  if (sfi->part < stpgetnumofparts(sfi->suftabparts))
  {
    preparethispart(sfi,mtime,err);
    *numberofsuffixes = sfi->widthofpart;
    *specialsuffixes = false;
    return sfi->suftab;
  }
  if (sfi->exhausted)
  {
    return NULL;
  }
  sfi->fusp.nextfreeSeqpos = 0;
  fillspecialnextpage(sfi);
  assert(sfi->fusp.nextfreeSeqpos > 0);
  *numberofsuffixes = (Seqpos) sfi->fusp.nextfreeSeqpos;
  *specialsuffixes = true;
  return sfi->suftab;
}

int bcktab2file(FILE *fp,
                const Sfxiterator *sfi,
                unsigned int prefixlength,
                Error *err)
{
  unsigned int numofallcodes = sfi->basepower[prefixlength],
               numofspecialcodes = sfi->basepower[prefixlength-1];

  if (fwrite(sfi->leftborder,
             sizeof (*sfi->leftborder),
             (size_t) (numofallcodes+1),
             fp)
             != (size_t) (numofallcodes+1))
  {
    error_set(err,"cannot write %u items of size %u: errormsg=\"%s\"",
              numofallcodes+1,
              (unsigned int) sizeof (*sfi->leftborder),
              strerror(errno));
    return -1;
  }
  if (fwrite(sfi->countspecialcodes,
             sizeof (*sfi->countspecialcodes),
             (size_t) numofspecialcodes,
             fp)
             != (size_t) numofspecialcodes)
  {
    error_set(err,"cannot write %u items of size %u: errormsg=\"%s\"",
              numofspecialcodes,
              (unsigned int) sizeof (*sfi->countspecialcodes),
              strerror(errno));
    return -2;
  }
  return 0;
}
