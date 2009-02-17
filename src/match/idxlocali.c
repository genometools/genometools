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
#include "extended/alignment.h"
#include "sarr-def.h"
#include "intbits.h"
#include "format64.h"
#include "idx-limdfs.h"
#include "idxlocali.h"
#include "idxlocalidp.h"
#include "absdfstrans-def.h"
#include "esa-map.h"
#include "stamp.h"

typedef struct
{
  const Uchar *characters;
  Uchar wildcardshow;
  bool showalignment;
  uint64_t unitnum;
  const Encodedsequence *encseq;
} Showmatchinfo;

static void showmatch(void *processinfo,const GtMatch *match)
{
  Showmatchinfo *showmatchinfo = (Showmatchinfo *) processinfo;
  Seqinfo seqinfo;
  unsigned long seqnum;

  seqnum = getencseqfrompos2seqnum(showmatchinfo->encseq,match->dbstartpos);
  getencseqSeqinfo(&seqinfo,showmatchinfo->encseq,seqnum);
  gt_assert(seqinfo.seqstartpos <= match->dbstartpos);
  printf("%lu\t" FormatSeqpos "\t",
         seqnum,
         PRINTSeqposcast(match->dbstartpos - seqinfo.seqstartpos));
  printf(FormatSeqpos "\t",PRINTSeqposcast(match->dblen));
  printf("\t" Formatuint64_t "\t%lu\t%lu\t%lu\n",
              PRINTuint64_tcast(showmatchinfo->unitnum),
              match->querystartpos,
              match->querylen,
              match->distance);
  if (showmatchinfo->showalignment)
  {
    gt_alignment_showwithmappedcharacters(
                (const GtAlignment *) match->alignment,
                showmatchinfo->characters,
                showmatchinfo->wildcardshow,
                stdout);
  }
}

int runidxlocali(const IdxlocaliOptions *arguments,GtError *err)
{
  Genericindex *genericindex = NULL;
  bool haserr = false;
  Verboseinfo *verboseinfo;

  verboseinfo = newverboseinfo(arguments->verbose);

  genericindex = genericindex_new(arguments->indexname,
                                  arguments->withesa,
                                  arguments->withesa,
                                  false,
                                  true,
                                  0,
                                  verboseinfo,
                                  err);
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
    int retval;
    Limdfsresources *limdfsresources = NULL;
    const AbstractDfstransformer *dfst;
    Showmatchinfo showmatchinfo;

    showmatchinfo.encseq = genericindex_getencseq(genericindex);
    showmatchinfo.characters 
      = getencseqAlphabetcharacters(showmatchinfo.encseq);
    showmatchinfo.wildcardshow 
      = getencseqAlphabetwildcardshow(showmatchinfo.encseq);
    showmatchinfo.showalignment = arguments->showalignment;
    dfst = locali_AbstractDfstransformer();
    gt_assert(genericindex != NULL);
    limdfsresources = newLimdfsresources(genericindex,
                                         true,
                                         0,
                                         0,    /* maxpathlength */
                                         true, /* keepexpandedonstack */
                                         showmatch,
                                         &showmatchinfo,
                                         NULL, /* processresult */
                                         NULL, /* processresult info */
                                         dfst);
    seqit = gt_seqiterator_new(arguments->queryfiles,
                               getencseqAlphabetsymbolmap(showmatchinfo.encseq),
                               true);
    for (showmatchinfo.unitnum = 0; /* Nothing */; showmatchinfo.unitnum++)
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
      printf("process sequence " Formatuint64_t " of length %lu\n",
              PRINTuint64_tcast(showmatchinfo.unitnum),querylen);
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
