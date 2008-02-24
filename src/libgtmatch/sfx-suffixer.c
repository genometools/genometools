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
#include "libgtcore/unused.h"
#include "spacedef.h"
#include "measure-time-if.h"
#include "intcode-def.h"
#include "encseq-def.h"
#include "safecast-gen.h"
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
  unsigned int maxprefixindex:PREFIXLENBITS;
  unsigned int code:CODEBITS;
  Seqpos position; /* get rid of this by using information from encseq */
} Codeatposition;

DECLAREARRAYSTRUCT(Seqpos);

struct Sfxiterator
{
  bool storespecials;
  Codetype currentmincode,
           currentmaxcode;
  Seqpos specialcharacters,
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
  const Uchar *characters;
  ArraySeqpos fusp;
  Specialrangeiterator *sri;
  Sequencerange overhang;
  bool exhausted;
  Bcktab *bcktab;
  Codetype numofallcodes;
  Seqpos *leftborder; /* points to bcktab->leftborder */
  const Definedunsignedint *maxdepth;
};

static void newcodelistelem(Codeatposition *cp,Seqpos prevright,Seqpos currleft,
                            unsigned int prefixlength)
{
 Seqpos distance = currleft - prevright;

  if (distance > (Seqpos) (prefixlength-1))
  {
    distance = (Seqpos) (prefixlength-1);
  }
  cp->maxprefixindex = (unsigned int) distance;
  cp->position = currleft;
}

static unsigned long produceCodeatpositionlist(Codeatposition *codelist,
                                               const  Encodedsequence *encseq,
                                               UNUSED bool moveforward,
                                               Seqpos realspecialranges,
                                               unsigned int prefixlength)
{
  Sequencerange previousrange;
  unsigned long insertindex = 0;
  Seqpos totallength;

  previousrange.rightpos = 0;
  if (hasspecialranges(encseq))
  {
    Specialrangeiterator *sri;
    Sequencerange currentrange;
    bool firstrange = true;

    sri = newspecialrangeiterator(encseq,true);
    while (nextspecialrangeiterator(&currentrange,sri))
    {
      if (firstrange)
      {
        firstrange = false;
      }
      if (currentrange.leftpos > previousrange.rightpos)
      {
        assert(insertindex < (unsigned long) (realspecialranges + 1));
        newcodelistelem(codelist + insertindex,previousrange.rightpos,
                        currentrange.leftpos,prefixlength);
        insertindex++;
      }
      previousrange = currentrange;
    }
    freespecialrangeiterator(&sri);
  }
  totallength = getencseqtotallength(encseq);
  if (totallength > previousrange.rightpos)
  {
    assert(insertindex < (unsigned long) (realspecialranges + 1));
    newcodelistelem(codelist + insertindex,previousrange.rightpos,totallength,
                    prefixlength);
    insertindex++;
  }
  return insertindex;
}

static void compareCodeatpositionlists(UNUSED const Codeatposition *codelist1,
                                       unsigned long len1,
                                       UNUSED const Codeatposition *codelist2,
                                       unsigned long len2)
{
  // unsigned long idx;

  if (len1 != len2)
  {
    fprintf(stderr,"len1 = %lu != %lu = len2\n",len1,len2);
    exit(EXIT_FAILURE); /* program error */
  }
  /*
  for (idx=0; idx<len1; idx++)
  {
    if (codelist1[idx].position != codelist2[idx].position)
    {
      fprintf(stderr,"idx %lu, codelist1.position = " FormatSeqpos " != "
                      FormatSeqpos " = codelist2.position\n",idx,
                      PRINTSeqposcast(codelist1[idx].position),
                      PRINTSeqposcast(codelist2[idx].position));
      exit(EXIT_FAILURE); 
    }
    if (codelist1[idx].maxprefixindex != codelist2[idx].maxprefixindex)
    {
      fprintf(stderr,"idx %lu, codelist1.maxprefixindex = %u != %u = "
                     "codelist2.maxprefixindex\n",idx,
                      codelist1[idx].maxprefixindex,
                      codelist2[idx].maxprefixindex);
      exit(EXIT_FAILURE);
    }
  }
  */
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
        cp->maxprefixindex = firstspecial->specialpos;
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
      sfi->suftabptr[--sfi->leftborder[code]] = position;
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

static void derivespecialcodes(Sfxiterator *sfi,bool deletevalues)
{
  Codetype code;
  unsigned int prefixindex;
  unsigned long insertindex, j;
  Seqpos stidx;

  for (prefixindex=1U; prefixindex < sfi->prefixlength; prefixindex++)
  {
    for (j=0, insertindex = 0; j < sfi->nextfreeCodeatposition; j++)
    {
      if (prefixindex <= sfi->spaceCodeatposition[j].maxprefixindex)
      {
        code = codedownscale(sfi->bcktab,
                             sfi->spaceCodeatposition[j].code,
                             prefixindex,
                             sfi->spaceCodeatposition[j].maxprefixindex);
        if (code >= sfi->currentmincode && code <= sfi->currentmaxcode)
        {
          updatebckspecials(sfi->bcktab,code,sfi->numofchars,prefixindex);
          stidx = --sfi->leftborder[code];
          sfi->suftabptr[stidx] = sfi->spaceCodeatposition[j].position -
                                  prefixindex;
        }
      }
      if (deletevalues)
      {
        if (prefixindex < sfi->prefixlength - 1 &&
            prefixindex < sfi->spaceCodeatposition[j].maxprefixindex)
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
#ifndef NDEBUG
  checkcountspecialcodes((*sfi)->bcktab);
#endif
  addfinalbckspecials((*sfi)->bcktab,(*sfi)->numofchars,
                      (*sfi)->specialcharacters);
  if ((*sfi)->sri != NULL)
  {
    freespecialrangeiterator(&(*sfi)->sri);
  }
  FREESPACE((*sfi)->spaceCodeatposition);
  FREESPACE((*sfi)->suftab);
  freesuftabparts((*sfi)->suftabparts);
  if ((*sfi)->bcktab != NULL)
  {
    freebcktab(&(*sfi)->bcktab);
  }
  FREESPACE(*sfi);
}

 DECLARESAFECASTFUNCTION(Seqpos,Seqpos,unsigned long,unsigned_long)

Sfxiterator *newSfxiterator(Seqpos specialcharacters,
                            Seqpos realspecialranges,
                            const Encodedsequence *encseq,
                            Readmode readmode,
                            unsigned int numofchars,
                            const Uchar *characters,
                            unsigned int prefixlength,
                            const Definedunsignedint *maxdepth,
                            unsigned int numofparts,
                            Outlcpinfo *outlcpinfo,
                            Measuretime *mtime,
                            Verboseinfo *verboseinfo,
                            Error *err)
{
  Sfxiterator *sfi = NULL;
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
                     Codeatposition,realspecialranges+1);
    printf("# sizeof (spaceCodeatposition)=%lu\n",
              (unsigned long ) sizeof (Codeatposition) * (realspecialranges+1));
    sfi->nextfreeCodeatposition = 0;
    sfi->suftab = NULL;
    sfi->suftabptr = NULL;
    sfi->suftabparts = NULL;
    sfi->encseq = encseq;
    sfi->readmode = readmode;
    sfi->numofchars = numofchars;
    sfi->characters = characters;
    sfi->prefixlength = prefixlength;
    assert(maxdepth != NULL);
    sfi->maxdepth = maxdepth;
    sfi->totallength = getencseqtotallength(encseq);
    sfi->specialcharacters = specialcharacters;
    sfi->outlcpinfo = outlcpinfo;
    sfi->sri = NULL;
    sfi->part = 0;
    sfi->exhausted = false;
    sfi->bcktab = allocBcktab(sfi->totallength,
                              numofchars,
                              prefixlength,
                              (unsigned int) CODEBITS,
                              (unsigned int) MAXCODEVALUE,
                              err);
    if (sfi->bcktab == NULL)
    {
      haserr = true;
      sfi->leftborder = NULL;
      sfi->numofallcodes = 0;
    } else
    {
      sfi->leftborder = bcktab_leftborder(sfi->bcktab);
      sfi->numofallcodes = bcktab_numofallcodes(sfi->bcktab);
    }
  }
  if (!haserr)
  {
    unsigned long nextfreeCodeatposition2;
    Codeatposition *spaceCodeatposition2;

    assert(sfi != NULL);
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
                   prefixlength);
    assert(realspecialranges+1 >= (Seqpos) sfi->nextfreeCodeatposition);
    ALLOCASSIGNSPACE(spaceCodeatposition2,NULL,Codeatposition,
                     realspecialranges+1);
    nextfreeCodeatposition2
      = produceCodeatpositionlist(spaceCodeatposition2,
                                  encseq,
                                  ISDIRREVERSE(readmode) ? false : true,
                                  realspecialranges,
                                  prefixlength);
    compareCodeatpositionlists(sfi->spaceCodeatposition,
                               sfi->nextfreeCodeatposition,
                               spaceCodeatposition2,
                               nextfreeCodeatposition2);
    FREESPACE(spaceCodeatposition2);
    assert(sfi->leftborder != NULL);
    for (optr = sfi->leftborder + 1;
         optr < sfi->leftborder + sfi->numofallcodes; optr++)
    {
      *optr += *(optr-1);
    }
    sfi->leftborder[sfi->numofallcodes]
      = sfi->totallength - specialcharacters;
    sfi->suftabparts = newsuftabparts(numofparts,
                                      sfi->leftborder,
                                      sfi->numofallcodes,
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
                            Measuretime *mtime)
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
                 sfi->prefixlength);
  if (mtime != NULL)
  {
    deliverthetime(stdout,mtime,"sorting the buckets");
  }
  totalwidth = stpgetcurrentsumofwdith(sfi->part,sfi->suftabparts);
  if (!sfi->maxdepth->defined ||
      sfi->prefixlength < sfi->maxdepth->valueunsignedint)
  {
    sortallbuckets(sfi->suftabptr,
                   sfi->encseq,
                   sfi->readmode,
                   sfi->currentmincode,
                   sfi->currentmaxcode,
                   totalwidth,
                   sfi->bcktab,
                   sfi->numofchars,
                   sfi->prefixlength,
                   sfi->maxdepth,
                   sfi->outlcpinfo);
  }
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
                              Measuretime *mtime,Sfxiterator *sfi)
{
  if (sfi->part < stpgetnumofparts(sfi->suftabparts))
  {
    preparethispart(sfi,mtime);
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

int sfibcktab2file(FILE *fp,
                   const Sfxiterator *sfi,
                   Error *err)
{
  error_check(err);
  return bcktab2file(fp,sfi->bcktab,err);
}
