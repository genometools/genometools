/*
  Copyright (c) 2009 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
  Copyright (c) 2009 Center for Bioinformatics, University of Hamburg

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

#include "core/unused_api.h"
#include "core/ma_api.h"
#include "core/seqiterator.h"
#include "sarr-def.h"
#include "intbits.h"
#include "format64.h"
#include "idx-limdfs.h"
#include "idxlocali.h"
#include "idxlocalidp.h"
#include "absdfstrans-def.h"
#include "esa-map.h"
#include "stamp.h"

static void showmatch(GT_UNUSED void *processinfo,
                      Seqpos dbstartpos,
                      Seqpos dblen,
                      GT_UNUSED const Uchar *dbsubstring,
                      unsigned long pprefixlen,
                      unsigned long distance)
{
  printf(FormatSeqpos "\t",PRINTSeqposcast(dblen));
  printf(FormatSeqpos "\t",PRINTSeqposcast(dbstartpos));
  printf("%lu\t%lu\n",pprefixlen,distance);
}

int runidxlocali(const IdxlocaliOptions *arguments,GtError *err)
{
  Genericindex *genericindex = NULL;
  bool haserr = false;
  Verboseinfo *verboseinfo;

  verboseinfo = newverboseinfo(arguments->verbose);

  genericindex = genericindex_new(arguments->indexname,
                                  arguments->withesa,
                                  arguments->withesa,0,
                                  verboseinfo,err);
  if (genericindex == NULL)
  {
    haserr = true;
  }
  if (!haserr)
  {
    GtSeqIterator *seqit;
    const Uchar *query;
    unsigned long querylen;
    char *desc = NULL;
    uint64_t unitnum;
    int retval;
    Limdfsresources *limdfsresources = NULL;
    const AbstractDfstransformer *dfst;
    const Encodedsequence *encseq;

    dfst = locali_AbstractDfstransformer();
    gt_assert(genericindex != NULL);
    limdfsresources = newLimdfsresources(genericindex,
                                         true,
                                         0,
                                         0,
                                         showmatch,
                                         NULL, /* processmatch info */
                                         NULL, /* processresult */
                                         NULL, /* processresult info */
                                         dfst);
    encseq = genericindex_getencseq(genericindex);
    seqit = gt_seqiterator_new(arguments->queryfiles,
                               getencseqAlphabetsymbolmap(encseq),
                               true);
    for (unitnum = 0; /* Nothing */; unitnum++)
    {
      retval = gt_seqiterator_next(seqit,
                                   &query,
                                   &querylen,
                                   &desc,
                                   err);
      if (retval < 0)
      {
        haserr = true;
        break;
      }
      if (retval == 0)
      {
        break;
      }
      printf("process sequence " Formatuint64_t "\n",
              PRINTuint64_tcast(unitnum));
      indexbasedlocali(limdfsresources,
                       arguments->matchscore,
                       arguments->mismatchscore,
                       arguments->gapstart,
                       arguments->gapextend,
                       arguments->threshold,
                       query,
                       querylen,
                       dfst);
      gt_free(desc);
    }
    if (limdfsresources != NULL)
    {
      freeLimdfsresources(&limdfsresources,dfst);
    }
    gt_seqiterator_delete(seqit);
  }
  genericindex_delete(genericindex);
  freeverboseinfo(&verboseinfo);
  return haserr ? -1 : 0;
}
