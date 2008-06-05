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
#include "tagerator.h"
#include "sarr-def.h"
#include "esa-mmsearch-def.h"
#include "intbits.h"
#include "alphadef.h"
#include "esa-myersapm.h"
#include "esa-limdfs.h"

#include "libgtmatch/esa-map.pr"

#define MAXTAGSIZE INTWORDSIZE

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

static void storematch(UNUSED void *processinfo,Seqpos startpos,Seqpos len)
{
  printf("match " FormatSeqpos " " FormatSeqpos "\n",
          PRINTSeqposcast(startpos),
          PRINTSeqposcast(len));
}

int runtagerator(const TageratorOptions *tageratoroptions,Error *err)
{
  Suffixarray suffixarray;
  Seqpos totallength;
  SeqIterator *seqit;
  Uchar charcode;
  bool haserr = false;
  char *desc = NULL;
  int retval;
  unsigned long idx, taglen, tagnumber;
  unsigned int mapsize, demand = SARR_SUFTAB | SARR_ESQTAB;
  const Uchar *symbolmap, *currenttag;
  Limdfsresources *limdfsresources = NULL;
  Myersonlineresources *mor = NULL;
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
  mapsize = getmapsizeAlphabet(suffixarray.alpha);
  if (tageratoroptions->online)
  {
    mor = newMyersonlineresources(mapsize,suffixarray.encseq,
                                  storematch,
                                  NULL);
  }
  if (tageratoroptions->maxdistance > 0)
  {
    limdfsresources = newLimdfsresources(suffixarray.encseq,
                                         suffixarray.readmode,
                                         mapsize,
                                         suffixarray.suftab,
                                         storematch,
                                         NULL);
  }
  seqit = seqiterator_new(tageratoroptions->tagfiles, NULL, true);
  for (tagnumber = 0; /* Nothing */; tagnumber++)
  {
    retval = seqiterator_next(seqit, &currenttag, &taglen, &desc, err);
    if (retval != 1)
    {
      if (retval < 0)
      {
        ma_free(desc);
      }
      break;
    }
    if (taglen > (unsigned long) MAXTAGSIZE)
    {
      error_set(err,"tag of length %lu; tags must not be longer than %lu",
                     taglen, (unsigned long) MAXTAGSIZE);
      haserr = true;
      ma_free(desc);
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
        ma_free(desc);
        break;
      }
      transformedtag[idx] = charcode;
    }
    printf("# patternlength=%lu\n",taglen);
    printf("# maxdistance=%lu\n",tageratoroptions->maxdistance);
    printf("# tag=");
    showsymbolstringgeneric(stdout,suffixarray.alpha,transformedtag,taglen);
    printf("\n");
    if (tageratoroptions->online)
    {
      edistmyersbitvectorAPM(mor,
                             transformedtag,
                             taglen,
                             tageratoroptions->maxdistance);
    } else
    {
      if (tageratoroptions->maxdistance == 0)
      {
        printf("tag %lu:",tagnumber);
        exactpatternmatching(suffixarray.encseq,
                             suffixarray.suftab,
                             suffixarray.readmode,
                             totallength,
                             transformedtag,
                             taglen);
      } else
      {
        esalimiteddfs(limdfsresources,
                      transformedtag,
                      taglen,
                      tageratoroptions->maxdistance);
      }
    }
    ma_free(desc);
  }
  if (limdfsresources != NULL)
  {
    freeLimdfsresources(&limdfsresources);
  }
  if (mor != NULL)
  {
    freeMyersonlineresources(&mor);
  }
  seqiterator_delete(seqit);
  freesuffixarray(&suffixarray);
  return haserr ? -1 : 0;
}
