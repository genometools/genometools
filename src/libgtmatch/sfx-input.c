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
#include <stdbool.h>
#include <errno.h>
#include "libgtcore/fa.h"
#include "libgtcore/filelengthvalues.h"
#include "spacedef.h"
#include "alphadef.h"
#include "sfx-optdef.h"
#include "encseq-def.h"
#include "measure-time-if.h"
#include "esafileend.h"
#include "verbose-def.h"
#include "sarr-def.h"
#include "stamp.h"
#include "sfx-input.h"

#include "opensfxfile.pr"
#include "esa-map.pr"
#include "fillsci.pr"

static int outal1file(const Str *indexname,const Alphabet *alpha,Error *err)
{
  FILE *al1fp;
  bool haserr = false;

  error_check(err);
  al1fp = opensfxfile(indexname,ALPHABETFILESUFFIX,"wb",err);
  if (al1fp == NULL)
  {
    haserr = true;
  }
  if (!haserr)
  {
    outputalphabet(al1fp,alpha);
    fa_xfclose(al1fp);
  }
  return haserr ? -1 : 0;
}

static unsigned long *initcharacterdistribution(const Alphabet *alpha)
{
  unsigned long *characterdistribution;
  unsigned int mapsize, idx;

  mapsize = getmapsizeAlphabet(alpha);
  ALLOCASSIGNSPACE(characterdistribution,NULL,unsigned long,mapsize-1);
  for (idx=0; idx<mapsize-1; idx++)
  {
    characterdistribution[idx] = 0;
  }
  return characterdistribution;
}

int fromfiles2Sfxseqinfo(Sfxseqinfo *sfxseqinfo,
                         Measuretime *mtime,
                         const Suffixeratoroptions *so,
                         Verboseinfo *verboseinfo,
                         Error *err)
{
  Seqpos totallength;
  bool haserr = false;

  error_check(err);
  sfxseqinfo->voidptr2suffixarray = NULL;
  sfxseqinfo->filelengthtab = NULL;
  sfxseqinfo->encseq = NULL;
  sfxseqinfo->characterdistribution = NULL;
  sfxseqinfo->readmode = so->readmode;
  sfxseqinfo->filenametab = so->filenametab;
  sfxseqinfo->alpha = assigninputalphabet(so->isdna,
                                          so->isprotein,
                                          so->str_smap,
                                          so->filenametab,
                                          err);
  if (sfxseqinfo->alpha == NULL)
  {
    haserr = true;
  }
  if (!haserr)
  {
    if (!so->isplain)
    {
      sfxseqinfo->characterdistribution
        = initcharacterdistribution(sfxseqinfo->alpha);
    }
    if (fasta2sequencekeyvalues(so->str_indexname,
                                &sfxseqinfo->numofsequences,
                                &totallength,
                                &sfxseqinfo->specialcharinfo,
                                so->filenametab,
                                &sfxseqinfo->filelengthtab,
                                getsymbolmapAlphabet(sfxseqinfo->alpha),
                                so->isplain,
                                so->outdestab,
                                sfxseqinfo->characterdistribution,
                                verboseinfo,
                                err) != 0)
    {
      haserr = true;
    }
  }
  if (!haserr)
  {
    if (outal1file(so->str_indexname,sfxseqinfo->alpha,err) != 0)
    {
      haserr = true;
    }
  }
  if (!haserr)
  {
    if (mtime != NULL)
    {
      deliverthetime(stdout,mtime,"computing sequence encoding");
    }
    sfxseqinfo->encseq
      = files2encodedsequence(true,
                              so->filenametab,
                              so->isplain,
                              totallength,
                              sfxseqinfo->specialcharinfo.specialranges,
                              sfxseqinfo->alpha,
                              str_length(so->str_sat) > 0
                                ? str_get(so->str_sat)
                                : NULL,
                              verboseinfo,
                              err);
    if (sfxseqinfo->encseq == NULL)
    {
      haserr = true;
    } else
    {
      if (so->outtistab)
      {
        if (flushencseqfile(so->str_indexname,sfxseqinfo->encseq,err) != 0)
        {
          haserr = true;
        }
      }
    }
  }
  return haserr ? -1 : 0;
}

int fromsarr2Sfxseqinfo(Sfxseqinfo *sfxseqinfo,
                        const Str *indexname,
                        Verboseinfo *verboseinfo,
                        Error *err)
{
  Seqpos totallength;
  bool haserr = false;
  Suffixarray *suffixarray;

  ALLOCASSIGNSPACE(suffixarray,NULL,Suffixarray,1);
  sfxseqinfo->characterdistribution = NULL;
  sfxseqinfo->voidptr2suffixarray = NULL;
  if (mapsuffixarray(suffixarray,
                     &totallength,
                     SARR_ESQTAB,
                     indexname,
                     verboseinfo,
                     err) != 0)
  {
    haserr = true;
    FREESPACE(suffixarray);
  } else
  {
    sfxseqinfo->voidptr2suffixarray = suffixarray; /* for freeing it later */
    sfxseqinfo->numofsequences = suffixarray->numofdbsequences;
    sfxseqinfo->alpha = suffixarray->alpha;
    sfxseqinfo->specialcharinfo = suffixarray->specialcharinfo;
    sfxseqinfo->filelengthtab = suffixarray->filelengthtab;
    sfxseqinfo->readmode = suffixarray->readmode;
    sfxseqinfo->filenametab = suffixarray->filenametab;
    assert(sfxseqinfo->filelengthtab != NULL);
    sfxseqinfo->encseq = suffixarray->encseq;
  }
  return haserr ? -1 : 0;
}

void freeSfxseqinfo(Sfxseqinfo *sfxseqinfo,bool mapped)
{
  FREESPACE(sfxseqinfo->characterdistribution);
  if (mapped)
  {
    if (sfxseqinfo->voidptr2suffixarray != NULL)
    {
      freesuffixarray((Suffixarray *) sfxseqinfo->voidptr2suffixarray);
    }
    FREESPACE(sfxseqinfo->voidptr2suffixarray);
  } else
  {
    FREESPACE(sfxseqinfo->filelengthtab);
    if (sfxseqinfo->alpha != NULL)
    {
      freeAlphabet(&sfxseqinfo->alpha);
    }
    freeEncodedsequence(&sfxseqinfo->encseq);
  }
}
