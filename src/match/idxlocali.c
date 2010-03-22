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
#include "core/seqiterator_sequence_buffer.h"
#include "extended/alignment.h"
#include "sarr-def.h"
#include "core/intbits.h"
#include "core/format64.h"
#include "idx-limdfs.h"
#include "idxlocali.h"
#include "idxlocalidp.h"
#include "idxlocalisw.h"
#include "absdfstrans-def.h"
#include "esa-map.h"
#include "stamp.h"

typedef struct
{
  const GtUchar *characters;
  uint64_t queryunit;
  GtUchar wildcardshow;
  bool showalignment;
  const GtEncodedsequence *encseq;
} Showmatchinfo;

static void showmatch(void *processinfo,const GtMatch *match)
{
  Showmatchinfo *showmatchinfo = (Showmatchinfo *) processinfo;
  unsigned long seqnum;
  unsigned long relpos;

  if (match->dbabsolute)
  {
    GtSeqinfo seqinfo;

    seqnum = getencseqfrompos2seqnum(showmatchinfo->encseq,match->dbstartpos);
    gt_encodedsequence_seqinfo(showmatchinfo->encseq,&seqinfo,seqnum);
    gt_assert(seqinfo.seqstartpos <= match->dbstartpos);
    relpos = match->dbstartpos - seqinfo.seqstartpos;
  } else
  {
    relpos = match->dbstartpos;
    seqnum = match->dbseqnum;
  }
  printf("%lu\t%lu\t",seqnum,relpos);
  printf("%lu\t",match->dblen);
  printf("\t" Formatuint64_t "\t%lu\t%lu\t%lu\n",
              PRINTuint64_tcast(showmatchinfo->queryunit),
              match->querystartpos,
              match->querylen,
              match->distance);
  if (showmatchinfo->showalignment)
  {
    gt_alignment_show_with_mapped_chars(
                (const GtAlignment *) match->alignment,
                showmatchinfo->characters,
                showmatchinfo->wildcardshow,
                stdout);
  }
}

typedef struct
{
  const GtEncodedsequence *encseq;
  GtBitsequence *hasmatch;
} Storematchinfo;

void initstorematch(Storematchinfo *storematch,
                    const GtEncodedsequence *encseq)
{
  unsigned long numofdbsequences = gt_encodedsequence_num_of_sequences(encseq);

  storematch->encseq = encseq;
  GT_INITBITTAB(storematch->hasmatch,numofdbsequences);
}

static void storematch(void *info,const GtMatch *match)
{
  Storematchinfo *storematch = (Storematchinfo *) info;
  unsigned long seqnum;

  if (match->dbabsolute)
  {
    seqnum = getencseqfrompos2seqnum(storematch->encseq,match->dbstartpos);
  } else
  {
    seqnum = match->dbseqnum;
  }
  if (!GT_ISIBITSET(storematch->hasmatch,seqnum))
  {
    GT_SETIBIT(storematch->hasmatch,seqnum);
  }
}

void checkandresetstorematch(GT_UNUSED uint64_t queryunit,
                             Storematchinfo *storeonline,
                             Storematchinfo *storeoffline)
{
  unsigned long seqnum, countmatchseq = 0,
    numofdbsequences = gt_encodedsequence_num_of_sequences(storeonline->encseq);

  for (seqnum = 0; seqnum < numofdbsequences; seqnum++)
  {
#ifndef NDEBUG
    if (GT_ISIBITSET(storeonline->hasmatch,seqnum) &&
        !GT_ISIBITSET(storeoffline->hasmatch,seqnum))
    {
      fprintf(stderr,"query " Formatuint64_t " refseq %lu: "
                     "online has match but offline not\n",
                     PRINTuint64_tcast(queryunit),seqnum);
      exit(GT_EXIT_PROGRAMMING_ERROR);
    }
    if (!GT_ISIBITSET(storeonline->hasmatch,seqnum) &&
        GT_ISIBITSET(storeoffline->hasmatch,seqnum))
    {
      fprintf(stderr,"query " Formatuint64_t " refseq %lu: "
                     "offline has match but online not\n",
                     PRINTuint64_tcast(queryunit),seqnum);
      exit(GT_EXIT_PROGRAMMING_ERROR);
    }
#endif
    if (GT_ISIBITSET(storeonline->hasmatch,seqnum))
    {
      countmatchseq++;
    }
  }
  GT_CLEARBITTAB(storeonline->hasmatch,numofdbsequences);
  GT_CLEARBITTAB(storeoffline->hasmatch,numofdbsequences);
  printf("matching sequences: %lu\n",countmatchseq);
}

void freestorematch(Storematchinfo *storematch)
{
  gt_free(storematch->hasmatch);
}

int runidxlocali(const IdxlocaliOptions *idxlocalioptions,GtError *err)
{
  Genericindex *genericindex = NULL;
  bool haserr = false;
  GtLogger *logger;
  const GtEncodedsequence *encseq = NULL;

  logger = gt_logger_new(idxlocalioptions->verbose,
                         GT_LOGGER_DEFLT_PREFIX, stdout);

  if (idxlocalioptions->doonline)
  {
    encseq = gt_encodedsequence_new_from_index (true,
                                 idxlocalioptions->indexname,
                                 true,
                                 false,
                                 false,
                                 true,
                                 logger,
                                 err);
    if (encseq == NULL)
    {
      haserr = true;
    }
  } else
  {
    genericindex = genericindex_new(idxlocalioptions->indexname,
                                    idxlocalioptions->withesa,
                                    idxlocalioptions->withesa ||
                                    idxlocalioptions->docompare,
                                    false,
                                    true,
                                    0,
                                    logger,
                                    err);
    if (genericindex == NULL)
    {
      haserr = true;
    } else
    {
      encseq = genericindex_getencseq(genericindex);
    }
  }
  if (!haserr)
  {
    GtSeqIterator *seqit;
    const GtUchar *query;
    unsigned long querylen;
    char *desc = NULL;
    int retval;
    Limdfsresources *limdfsresources = NULL;
    const AbstractDfstransformer *dfst;
    SWdpresource *swdpresource = NULL;
    Showmatchinfo showmatchinfo;
    Processmatch processmatch;
    void *processmatchinfoonline, *processmatchinfooffline;
    Storematchinfo storeonline, storeoffline;

    if (idxlocalioptions->docompare)
    {
      processmatch = storematch;
      initstorematch(&storeonline,encseq);
      initstorematch(&storeoffline,encseq);
      processmatchinfoonline = &storeonline;
      processmatchinfooffline = &storeoffline;
    } else
    {
      processmatch = showmatch;
      showmatchinfo.encseq = encseq;
      showmatchinfo.characters = gt_encodedsequence_alphabetcharacters(encseq);
      showmatchinfo.wildcardshow =
                               gt_encodedsequence_alphabetwildcardshow(encseq);
      showmatchinfo.showalignment = idxlocalioptions->showalignment;
      processmatchinfoonline = processmatchinfooffline = &showmatchinfo;
    }
    if (idxlocalioptions->doonline || idxlocalioptions->docompare)
    {
      swdpresource = newSWdpresource(idxlocalioptions->matchscore,
                                     idxlocalioptions->mismatchscore,
                                     idxlocalioptions->gapextend,
                                     idxlocalioptions->threshold,
                                     idxlocalioptions->showalignment,
                                     processmatch,
                                     processmatchinfoonline);
    }
    dfst = locali_AbstractDfstransformer();
    if (!idxlocalioptions->doonline || idxlocalioptions->docompare)
    {
      gt_assert(genericindex != NULL);
      limdfsresources = newLimdfsresources(genericindex,
                                           true,
                                           0,
                                           0,    /* maxpathlength */
                                           true, /* keepexpandedonstack */
                                           processmatch,
                                           processmatchinfooffline,
                                           NULL, /* processresult */
                                           NULL, /* processresult info */
                                           dfst);
    }
    seqit = gt_seqiterator_sequence_buffer_new(idxlocalioptions->queryfiles,
                                               err);
    if (!seqit)
      haserr = true;
    if (!haserr)
    {
      gt_seqiterator_set_symbolmap(seqit,
                                  gt_encodedsequence_alphabetsymbolmap(encseq));
      for (showmatchinfo.queryunit = 0; /* Nothing */;
           showmatchinfo.queryunit++)
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
                PRINTuint64_tcast(showmatchinfo.queryunit),querylen);
        if (idxlocalioptions->doonline || idxlocalioptions->docompare)
        {
          multiapplysmithwaterman(swdpresource,encseq,query,querylen);
        }
        if (!idxlocalioptions->doonline || idxlocalioptions->docompare)
        {
          indexbasedlocali(limdfsresources,
                           idxlocalioptions->matchscore,
                           idxlocalioptions->mismatchscore,
                           idxlocalioptions->gapstart,
                           idxlocalioptions->gapextend,
                           idxlocalioptions->threshold,
                           query,
                           querylen,
                           dfst);
        }
        if (idxlocalioptions->docompare)
        {
          checkandresetstorematch(showmatchinfo.queryunit,
                                  &storeonline,&storeoffline);
        }
        gt_free(desc);
      }
      if (limdfsresources != NULL)
      {
        freeLimdfsresources(&limdfsresources,dfst);
      }
      if (swdpresource != NULL)
      {
        freeSWdpresource(swdpresource);
        swdpresource = NULL;
      }
      gt_seqiterator_delete(seqit);
    }
    if (idxlocalioptions->docompare)
    {
      freestorematch(&storeonline);
      freestorematch(&storeoffline);
    }
  }
  if (genericindex == NULL)
  {
    gt_assert(encseq != NULL);
    gt_encodedsequence_delete((GtEncodedsequence *) encseq);
    encseq = NULL;
  } else
  {
    genericindex_delete(genericindex);
  }
  gt_logger_delete(logger);
  logger = NULL;
  return haserr ? -1 : 0;
}
