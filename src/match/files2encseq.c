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
#include "core/ma_api.h"
#include "core/filelengthvalues.h"
#include "encseq-def.h"
#include "sfx-progress.h"
#include "verbose-def.h"
#include "files2encseq.h"

static unsigned long *initcharacterdistribution(const GtAlphabet *alpha)
{
  unsigned long *characterdistribution;
  unsigned int numofchars, idx;

  numofchars = gt_alphabet_num_of_chars(alpha);
  characterdistribution = gt_malloc(sizeof (*characterdistribution) *
                                    numofchars);
  for (idx=0; idx<numofchars; idx++)
  {
    characterdistribution[idx] = 0;
  }
  return characterdistribution;
}

Encodedsequence *fromfiles2encseq(ArraySeqpos *sequenceseppos,
                                  Sfxprogress *sfxprogress,
                                  const GtStr *str_indexname,
                                  const GtStr *str_smap,
                                  const GtStr *str_sat,
                                  const GtStrArray *filenametab,
                                  bool isdna,
                                  bool isprotein,
                                  bool isplain,
                                  bool outtistab,
                                  bool outdestab,
                                  bool outsdstab,
                                  bool outssptab,
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
  Encodedsequence *encseq = NULL;

  gt_error_check(err);
  encseq = NULL;
  if (gt_str_length(str_sat) > 0)
  {
    int retval = getsatforcevalue(gt_str_get(str_sat),err);
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
    alpha = gt_alphabet_new(isdna, isprotein,str_smap, filenametab, err);
    if (alpha == NULL)
    {
      haserr = true;
    }
  }
  if (!haserr)
  {
    if (gt_outal1file(str_indexname,alpha,err) != 0)
    {
      haserr = true;
    }
  }
  if (!haserr)
  {
    characterdistribution = initcharacterdistribution(alpha);
    if (gt_inputfiles2sequencekeyvalues(str_indexname,
                                        &totallength,
                                        &specialcharinfo,
                                        forcetable,
                                        specialrangestab,
                                        filenametab,
                                        &filelengthtab,
                                        alpha,
                                        isplain,
                                        outdestab,
                                        outsdstab,
                                        characterdistribution,
                                        outssptab,
                                        sequenceseppos,
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
    encseq = files2encodedsequence(true,
                                   filenametab,
                                   filelengthtab,
                                   isplain,
                                   totallength,
                                   sequenceseppos->nextfreeSeqpos+1,
                                   specialrangestab,
                                   alpha,
                                   gt_str_length(str_sat) > 0
                                     ? gt_str_get(str_sat)
                                     : NULL,
                                   characterdistribution,
                                   &specialcharinfo,
                                   verboseinfo,
                                   err);
    if (encseq == NULL)
    {
      haserr = true;
    } else
    {
      alphaisbound = true;
      if (outtistab)
      {
        if (flushencseqfile(str_indexname,encseq,err) != 0)
        {
          haserr = true;
        }
      }
    }
  }
  if (haserr)
  {
    gt_free(characterdistribution);
    gt_free(filelengthtab);
    filelengthtab = NULL;
    if (alpha != NULL && !alphaisbound)
    {
      gt_alphabet_delete((GtAlphabet*) alpha);
    }
  }
  return haserr ? NULL : encseq;
}
