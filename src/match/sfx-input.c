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
#include "core/alphabet.h"
#include "core/filelengthvalues.h"
#include "spacedef.h"
#include "sfx-optdef.h"
#include "encseq-def.h"
#include "sfx-progress.h"
#include "verbose-def.h"
#include "sfx-input.h"

static unsigned long *initcharacterdistribution(const GtAlphabet *alpha)
{
  unsigned long *characterdistribution;
  unsigned int numofchars, idx;

  numofchars = gt_alphabet_num_of_chars(alpha);
  ALLOCASSIGNSPACE(characterdistribution,NULL,unsigned long,numofchars);
  for (idx=0; idx<numofchars; idx++)
  {
    characterdistribution[idx] = 0;
  }
  return characterdistribution;
}

int fromfiles2Sfxseqinfo(Sfxseqinfo *sfxseqinfo,
                         Sfxprogress *sfxprogress,
                         const Suffixeratoroptions *so,
                         Verboseinfo *verboseinfo,
                         GtError *err)
{
  Seqpos totallength;
  bool haserr = false;
  unsigned int forcetable;
  Specialcharinfo specialcharinfo;
  const GtAlphabet *alpha = NULL;
  bool alphaisbound = false;
  Filelengthvalues *filelengthtab = NULL;
  Seqpos specialrangestab[3];
  unsigned long *characterdistribution = NULL;

  gt_error_check(err);
  sfxseqinfo->encseq = NULL;
  sfxseqinfo->readmode = so->readmode;
  GT_INITARRAY(&sfxseqinfo->sequenceseppos,Seqpos);
  if (gt_str_length(so->str_sat) > 0)
  {
    int retval = getsatforcevalue(gt_str_get(so->str_sat),err);
    if (retval < 0)
    {
      haserr = true;
    } else
    {
      forcetable = (unsigned int) retval;
    }
  } else
  {
    forcetable = 3U;
  }
  if (!haserr)
  {
    alpha = gt_alphabet_new(so->isdna, so->isprotein, so->str_smap,
                            so->filenametab, err);
    if (alpha == NULL)
    {
      haserr = true;
    }
  }
  if (!haserr)
  {
    if (gt_outal1file(so->str_indexname,alpha,err) != 0)
    {
      haserr = true;
    }
  }
  if (!haserr)
  {
    characterdistribution = initcharacterdistribution(alpha);
    if (gt_inputfiles2sequencekeyvalues(
                                so->str_indexname,
                                &totallength,
                                &specialcharinfo,
                                forcetable,
                                specialrangestab,
                                so->filenametab,
                                &filelengthtab,
                                alpha,
                                so->isplain,
                                so->outdestab,
                                so->outsdstab,
                                characterdistribution,
                                so->outssptab,
                                &sfxseqinfo->sequenceseppos,
                                verboseinfo,
                                err) != 0)
    {
      haserr = true;
    }
  }
  if (!haserr)
  {
    if (sfxprogress != NULL)
    {
      sfxprogress_deliverthetime(stdout,
                                 sfxprogress,"computing sequence encoding");
    }
    sfxseqinfo->encseq
      = files2encodedsequence(true,
                              so->filenametab,
                              filelengthtab,
                              so->isplain,
                              totallength,
                              sfxseqinfo->sequenceseppos.nextfreeSeqpos+1,
                              specialrangestab,
                              alpha,
                              gt_str_length(so->str_sat) > 0
                                ? gt_str_get(so->str_sat)
                                : NULL,
                              characterdistribution,
                              &specialcharinfo,
                              verboseinfo,
                              err);
    if (sfxseqinfo->encseq == NULL)
    {
      haserr = true;
      FREESPACE(characterdistribution);
    } else
    {
      alphaisbound = true;
      if (so->outtistab)
      {
        if (flushencseqfile(so->str_indexname,sfxseqinfo->encseq,err) != 0)
        {
          haserr = true;
        }
      }
    }
  }
  if (haserr)
  {
    FREESPACE(characterdistribution);
    gt_free(filelengthtab);
    filelengthtab = NULL;
    if (alpha != NULL && !alphaisbound)
    {
      gt_alphabet_delete((GtAlphabet*) alpha);
    }
  }
  return haserr ? -1 : 0;
}

int fromsarr2Sfxseqinfo(Sfxseqinfo *sfxseqinfo,
                        const GtStr *indexname,
                        Readmode readmodeoption,
                        Verboseinfo *verboseinfo,
                        GtError *err)
{
  sfxseqinfo->readmode = readmodeoption;
  GT_INITARRAY(&sfxseqinfo->sequenceseppos,Seqpos);
  sfxseqinfo->encseq = mapencodedsequence(true,
                                          indexname,
                                          true,
                                          false,
                                          false,
                                          false,
                                          verboseinfo,
                                          err);
  if (sfxseqinfo->encseq == NULL)
  {
    return -1;
  }
  return 0;
}

void freeSfxseqinfo(Sfxseqinfo *sfxseqinfo)
{
  if (sfxseqinfo->encseq != NULL)
  {
    encodedsequence_free(&sfxseqinfo->encseq);
  }
  GT_FREEARRAY(&sfxseqinfo->sequenceseppos,Seqpos);
}
