/*
  Copyright (c) 2008 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
  Copyright (c) 2008 Center for Bioinformatics, University of Hamburg

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

#include <errno.h>
#include "core/alphabet.h"
#include "core/divmodmul.h"
#include "core/fa.h"
#include "core/format64.h"
#include "core/logger.h"
#include "core/spacecalc.h"
#include "core/str.h"
#include "core/unused_api.h"
#include "core/xansi_api.h"
#include "core/ma_api.h"
#include "esa-seqread.h"
#include "esa-mmsearch.h"
#include "tyr-basic.h"
#include "tyr-mkindex.h"
#include "echoseq.h"

typedef struct /* information stored for each node of the lcp interval tree */
{
  unsigned long leftmostleaf,
                rightmostleaf,
                suftabrightmostleaf,
                lcptabrightmostleafplus1;
} TyrDfsinfo;

typedef int (*Processoccurrencecount)(unsigned long,
                                      unsigned long,
                                      void *,
                                      GtError *);

typedef struct listUlong
{
  unsigned long position;
  struct listUlong *nextptr;
} ListUlong;

typedef struct
{
  unsigned long occcount;
  ListUlong *positionlist;
} Countwithpositions;

GT_DECLAREARRAYSTRUCT(Countwithpositions);

typedef struct /* global information */
{
  unsigned long mersize,
                totallength,
                minocc,
                maxocc;
  const GtEncseq *encseq;
  GtReadmode readmode;
  Processoccurrencecount processoccurrencecount;
  GtArrayCountwithpositions occdistribution;
  FILE *merindexfpout,
       *countsfilefpout;
  GtEncseqReader *esrspace;
  bool performtest;
  bool storecounts;
  GtUchar *bytebuffer;
  unsigned long sizeofbuffer;
  GtArrayLargecount largecounts;
  unsigned long countoutputmers;
  const ESASuffixptr *suftab; /* only necessary for performtest */
  GtUchar *currentmer;    /* only necessary for performtest */
} TyrDfsstate;

#include "esa-dfs.h"

static uint64_t bruteforcecountnumofmers(const TyrDfsstate *state)
{
  unsigned long idx;
  uint64_t numofmers = 0;

  for (idx=0; idx <= state->totallength - state->mersize; idx++)
  {
    if (!gt_encseq_contains_special(state->encseq,
                                    state->readmode,
                                    state->esrspace,
                                    idx,
                                    state->mersize))
    {
      numofmers++;
    }
  }
  return numofmers;
}

static void checknumofmers(const TyrDfsstate *state,
                           uint64_t dnumofmers)
{
  uint64_t bfnumofmers = bruteforcecountnumofmers(state);

  if (dnumofmers != bfnumofmers)
  {
    fprintf(stderr,"numofmers(distribution) = " Formatuint64_t " != "
                   Formatuint64_t,
                   PRINTuint64_tcast(dnumofmers),
                   PRINTuint64_tcast(bfnumofmers));
    exit(GT_EXIT_PROGRAMMING_ERROR);
  }
}

static void checknumberofoccurrences(const TyrDfsstate *dfsstate,
                                     unsigned long countocc,
                                     unsigned long position)
{
  GtMMsearchiterator *mmsi;
  unsigned long idx, bfcount;

  for (idx = 0; idx < dfsstate->mersize; idx++)
  {
    dfsstate->currentmer[idx] =
              gt_encseq_get_encoded_char(dfsstate->encseq,position+idx,
                                                dfsstate->readmode);
  }
  mmsi = gt_mmsearchiterator_new_complete_olain(dfsstate->encseq,
                                              dfsstate->suftab,
                                              0,
                                              dfsstate->totallength,
                                              0,
                                              dfsstate->readmode,
                                              dfsstate->currentmer,
                                              dfsstate->mersize);
  bfcount = gt_mmsearchiterator_count(mmsi);
  if (bfcount != countocc)
  {
    fprintf(stderr,"bfcount = %lu != %lu = countocc\n",bfcount,countocc);
    exit(GT_EXIT_PROGRAMMING_ERROR);
  }
  gt_mmsearchiterator_delete(mmsi);
}

static /*@null@*/ ListUlong *insertListUlong(ListUlong *liststart,
                                             unsigned long position)
{
  ListUlong *newnode;

  newnode = gt_malloc(sizeof *newnode);
  newnode->position = position;
  newnode->nextptr = liststart;
  return newnode;
}

static void wrapListUlong(ListUlong *node)
{
  ListUlong *tmpnext;

  if (node != NULL)
  {
    while (true)
    {
      if (node->nextptr == NULL)
      {
        gt_free(node);
        break;
      }
      tmpnext = node->nextptr;
      gt_free(node);
      node = tmpnext;
    }
  }
}

static void showListUlong(const GtEncseq *encseq,
                          unsigned long mersize,
                          const ListUlong *node)
{
  const ListUlong *tmp;

  for (tmp = node; tmp != NULL; tmp = tmp->nextptr)
  {
    gt_fprintfencseq(stdout,encseq,tmp->position,mersize);
    (void) putchar((int) '\n');
  }
}

static bool decideifocc(const TyrDfsstate *state,unsigned long countocc)
{
  if (state->minocc > 0)
  {
    if (state->maxocc > 0)
    {
      if (countocc >= state->minocc && countocc <= state->maxocc)
      {
        return true;
      }
    } else
    {
      if (countocc >= state->minocc)
      {
        return true;
      }
    }
  } else
  {
    if (state->maxocc > 0)
    {
      if (countocc <= state->maxocc)
      {
        return true;
      }
    }
  }
  return false;
}

static uint64_t addupdistribution(const GtArrayCountwithpositions *distribution)
{
  unsigned long idx;
  uint64_t addcount = 0;

  for (idx=0; idx < distribution->nextfreeCountwithpositions; idx++)
  {
    if (distribution->spaceCountwithpositions[idx].occcount > 0)
    {
      addcount += (uint64_t)
                  (idx * distribution->spaceCountwithpositions[idx].occcount);
    }
  }
  return addcount;
}

static void showmerdistribution(const TyrDfsstate *state)
{
  unsigned long countocc;

  for (countocc = 0;
      countocc < state->occdistribution.nextfreeCountwithpositions;
      countocc++)
  {
    if (state->occdistribution.spaceCountwithpositions[countocc].occcount > 0)
    {
      printf("%lu %lu\n",countocc,
                         state->occdistribution.
                                spaceCountwithpositions[countocc].occcount);
      if (decideifocc(state,countocc))
      {
        showListUlong(state->encseq,
                      state->mersize,
                      state->occdistribution.spaceCountwithpositions[countocc].
                                             positionlist);
        wrapListUlong(state->occdistribution.spaceCountwithpositions[countocc].
                                             positionlist);
      }
    }
  }
}

static void showfinalstatistics(const TyrDfsstate *state,
                                const char *inputindex,
                                GtLogger *logger)
{
  uint64_t dnumofmers = addupdistribution(&state->occdistribution);

  if (state->performtest)
  {
    checknumofmers(state,dnumofmers);
  }
  gt_logger_log(logger,
              "the following output refers to the set of all sequences");
  gt_logger_log(logger,
              "represented by the index \"%s\"",inputindex);
  gt_logger_log(logger,
              "number of %lu-mers in the sequences not containing a "
              "wildcard: " Formatuint64_t,
              (unsigned long) state->mersize,
              PRINTuint64_tcast(dnumofmers));
  gt_logger_log(logger,
              "show the distribution of the number of occurrences of %lu-mers",
               (unsigned long) state->mersize);
  gt_logger_log(logger,"not containing a wildcard as rows of the form "
              "i d where");
  gt_logger_log(logger,
              "d is the number of events that a %lu-mer occurs exactly i times",
              (unsigned long) state->mersize);
  showmerdistribution(state);
}

static void incrementdistribcounts(GtArrayCountwithpositions *occdistribution,
                                   unsigned long countocc,unsigned long value)
{
  if (countocc >= occdistribution->allocatedCountwithpositions)
  {
    const unsigned long addamount = 128UL;
    unsigned long idx;

    occdistribution->spaceCountwithpositions
      = gt_realloc(occdistribution->spaceCountwithpositions,
                   sizeof *occdistribution->spaceCountwithpositions
                   * (countocc+addamount));
    for (idx=occdistribution->allocatedCountwithpositions;
         idx<countocc+addamount; idx++)
    {
      occdistribution->spaceCountwithpositions[idx].occcount = 0;
      occdistribution->spaceCountwithpositions[idx].positionlist = NULL;
    }
    occdistribution->allocatedCountwithpositions = countocc+addamount;
  }
  if (countocc + 1 > occdistribution->nextfreeCountwithpositions)
  {
    occdistribution->nextfreeCountwithpositions = countocc+1;
  }
  occdistribution->spaceCountwithpositions[countocc].occcount += value;
}

static int adddistpos2distribution(unsigned long countocc,
                                   unsigned long position,
                                   void *adddistposinfo,
                                   GT_UNUSED GtError *err)
{
  TyrDfsstate *state = (TyrDfsstate *) adddistposinfo;

  incrementdistribcounts(&state->occdistribution,countocc,1UL);
  if (decideifocc(state,countocc))
  {
    state->occdistribution.spaceCountwithpositions[countocc].positionlist
      = insertListUlong(state->occdistribution.
                               spaceCountwithpositions[countocc].positionlist,
                         position);
  }
  if (state->performtest)
  {
    checknumberofoccurrences(state,countocc,position);
  }
  return 0;
}

#define MAXSMALLMERCOUNT UCHAR_MAX

static int outputsortedstring2indexviafileptr(const GtEncseq *encseq,
                                              unsigned long mersize,
                                              GtUchar *bytebuffer,
                                              unsigned long sizeofbuffer,
                                              FILE *merindexfpout,
                                              FILE *countsfilefpout,
                                              unsigned long position,
                                              unsigned long countocc,
                                              GtArrayLargecount *largecounts,
                                              unsigned long countoutputmers,
                                              GT_UNUSED GtError *err)
{
  gt_encseq_sequence2bytecode(bytebuffer,encseq,position,mersize);
  gt_xfwrite(bytebuffer, sizeof (*bytebuffer), (size_t) sizeofbuffer,
             merindexfpout);
  if (countsfilefpout != NULL)
  {
    GtUchar smallcount;

    if (countocc <= MAXSMALLMERCOUNT)
    {
      smallcount = (GtUchar) countocc;
    } else
    {
      Largecount *lc;

      GT_GETNEXTFREEINARRAY(lc,largecounts,Largecount,32);
      lc->idx = countoutputmers;
      lc->value = countocc;
      smallcount = 0;
    }
    gt_xfwrite(&smallcount, sizeof (smallcount),(size_t) 1,countsfilefpout);
  }
  return 0;
}

static int outputsortedstring2index(unsigned long countocc,
                                    unsigned long position,
                                    void *adddistposinfo,
                                    GtError *err)
{
  TyrDfsstate *state = (TyrDfsstate *) adddistposinfo;

  if (decideifocc(state,countocc))
  {
    if (outputsortedstring2indexviafileptr(state->encseq,
                                           state->mersize,
                                           state->bytebuffer,
                                           state->sizeofbuffer,
                                           state->merindexfpout,
                                           state->countsfilefpout,
                                           position,
                                           countocc,
                                           &state->largecounts,
                                           state->countoutputmers,
                                           err) != 0)
    {
      return -1;
    }
    state->countoutputmers++;
  }
  return 0;
}

static Dfsinfo* tyr_allocateDfsinfo(GT_UNUSED Dfsstate *state)
{
  TyrDfsinfo *dfsinfo;

  dfsinfo = gt_malloc(sizeof *dfsinfo);
  return (Dfsinfo*) dfsinfo;
}

static void tyr_freeDfsinfo(Dfsinfo *adfsinfo, GT_UNUSED Dfsstate *state)
{
  TyrDfsinfo *dfsinfo = (TyrDfsinfo*) adfsinfo;
  gt_free(dfsinfo);
}

static int tyr_processleafedge(GT_UNUSED bool firstsucc,
                           unsigned long fatherdepth,
                           GT_UNUSED Dfsinfo *father,
                           unsigned long leafnumber,
                           Dfsstate *astate,
                           GtError *err)
{
  TyrDfsstate *state = (TyrDfsstate*) astate;
  gt_error_check(err);
  if (fatherdepth < state->mersize &&
      leafnumber + state->mersize <= state->totallength &&
      !gt_encseq_contains_special(state->encseq,
                                  state->readmode,
                                  state->esrspace,
                                  leafnumber + fatherdepth,
                                  state->mersize - fatherdepth))
  {
    if (state->processoccurrencecount(1UL,leafnumber,state,err) != 0)
    {
      return -1;
    }
  }
  return 0;
}

static int tyr_processcompletenode(unsigned long nodeptrdepth,
                               Dfsinfo *anodeptr,
                               unsigned long nodeptrminusonedepth,
                               Dfsstate *astate,
                               GtError *err)
{
  TyrDfsinfo *nodeptr = (TyrDfsinfo*) anodeptr;
  TyrDfsstate *state = (TyrDfsstate*) astate;
  gt_error_check(err);
  if (state->mersize <= nodeptrdepth)
  {
    unsigned long fatherdepth;

    fatherdepth = nodeptr->lcptabrightmostleafplus1;
    if (fatherdepth < nodeptrminusonedepth)
    {
      fatherdepth = nodeptrminusonedepth;
    }
    if (fatherdepth < state->mersize)
    {
      if (state->processoccurrencecount((unsigned long)
                                        (nodeptr->rightmostleaf -
                                         nodeptr->leftmostleaf + 1),
                                        nodeptr->suftabrightmostleaf,
                                        state,
                                        err) != 0)
      {
        return -1;
      }
    }
  }
  return 0;
}

static void tyr_assignleftmostleaf(Dfsinfo *adfsinfo,unsigned long leftmostleaf,
                                   GT_UNUSED Dfsstate *dfsstate)
{
  TyrDfsinfo *dfsinfo = (TyrDfsinfo*) adfsinfo;
  dfsinfo->leftmostleaf = leftmostleaf;
}

static void tyr_assignrightmostleaf(Dfsinfo *adfsinfo,
                                    unsigned long currentindex,
                                    unsigned long previoussuffix,
                                    unsigned long currentlcp,
                                    GT_UNUSED Dfsstate *dfsstate)
{
  TyrDfsinfo *dfsinfo = (TyrDfsinfo*) adfsinfo;
  dfsinfo->rightmostleaf = currentindex;
  dfsinfo->suftabrightmostleaf = previoussuffix;
  dfsinfo->lcptabrightmostleafplus1 = currentlcp;
}

static void outputbytewiseUlongvalue(FILE *fpout,unsigned long value)
{
  size_t i;

  for (i=0; i < sizeof (value); i++)
  {
    (void) putc((int) (value & UCHAR_MAX),fpout);
    value >>= 8;
  }
}

static int enumeratelcpintervals(const char *inputindex,
                                 Sequentialsuffixarrayreader *ssar,
                                 const char *storeindex,
                                 bool storecounts,
                                 unsigned long mersize,
                                 unsigned long minocc,
                                 unsigned long maxocc,
                                 bool performtest,
                                 GtLogger *logger,
                                 GtError *err)
{
  TyrDfsstate *state;
  bool haserr = false;
  unsigned int alphasize;

  gt_error_check(err);
  state = gt_malloc(sizeof (*state));
  GT_INITARRAY(&state->occdistribution,Countwithpositions);
  state->esrspace = gt_encseq_create_reader_with_readmode(
                                   gt_encseqSequentialsuffixarrayreader(ssar),
                                   gt_readmodeSequentialsuffixarrayreader(ssar),
                                   0);
  state->mersize = (unsigned long) mersize;
  state->encseq = gt_encseqSequentialsuffixarrayreader(ssar);
  alphasize = gt_alphabet_num_of_chars(gt_encseq_alphabet(state->encseq));
  state->readmode = gt_readmodeSequentialsuffixarrayreader(ssar);
  state->storecounts = storecounts;
  state->minocc = minocc;
  state->maxocc = maxocc;
  state->totallength = gt_encseq_total_length(state->encseq);
  state->performtest = performtest;
  state->countoutputmers = 0;
  state->merindexfpout = NULL;
  state->countsfilefpout = NULL;
  GT_INITARRAY(&state->largecounts,Largecount);
  if (strlen(storeindex) == 0)
  {
    state->sizeofbuffer = 0;
    state->bytebuffer = NULL;
  } else
  {
    state->sizeofbuffer = MERBYTES(mersize);
    state->bytebuffer = gt_malloc(sizeof *state->bytebuffer
                                  * state->sizeofbuffer);
  }
  if (performtest)
  {
    state->currentmer = gt_malloc(sizeof *state->currentmer
                                  * state->mersize);
    state->suftab = gt_suftabSequentialsuffixarrayreader(ssar);
  } else
  {
    state->currentmer = NULL;
    state->suftab = NULL;
  }
  if (state->mersize > state->totallength)
  {
    gt_error_set(err,"mersize %lu > %lu = totallength not allowed",
                 state->mersize,
                 state->totallength);
    haserr = true;
  } else
  {
    if (strlen(storeindex) == 0)
    {
      state->processoccurrencecount = adddistpos2distribution;
    } else
    {
      state->merindexfpout = gt_fa_fopen_with_suffix(storeindex,MERSUFFIX,
                                                    "wb",err);
      if (state->merindexfpout == NULL)
      {
        haserr = true;
      } else
      {
        if (state->storecounts)
        {
          state->countsfilefpout
            = gt_fa_fopen_with_suffix(storeindex,COUNTSSUFFIX,"wb",err);
          if (state->countsfilefpout == NULL)
          {
            haserr = true;
          }
        }
      }
      state->processoccurrencecount = outputsortedstring2index;
    }
    if (!haserr)
    {
      if (gt_depthfirstesa(ssar,
                          tyr_allocateDfsinfo,
                          tyr_freeDfsinfo,
                          tyr_processleafedge,
                          NULL,
                          tyr_processcompletenode,
                          tyr_assignleftmostleaf,
                          tyr_assignrightmostleaf,
                          (Dfsstate*) state,
                          logger,
                          err) != 0)
      {
        haserr = true;
      }
      if (strlen(storeindex) == 0)
      {
        showfinalstatistics(state,inputindex,logger);
      }
    }
    if (!haserr)
    {
      if (state->countsfilefpout != NULL)
      {
        gt_logger_log(logger,"write %lu mercounts > %lu to file \"%s%s\"",
                    state->largecounts.nextfreeLargecount,
                    (unsigned long) MAXSMALLMERCOUNT,
                    storeindex,
                    COUNTSSUFFIX);
        gt_xfwrite(state->largecounts.spaceLargecount, sizeof (Largecount),
                  (size_t) state->largecounts.nextfreeLargecount,
                  state->countsfilefpout);
      }
    }
    if (!haserr)
    {
      gt_logger_log(logger,"number of %lu-mers in index: %lu",
                  mersize,
                  state->countoutputmers);
      gt_logger_log(logger,"index size: %.2f megabytes\n",
                  GT_MEGABYTES(state->countoutputmers * state->sizeofbuffer +
                               sizeof (unsigned long) * EXTRAINTEGERS));
    }
  }
  /* now out EXTRAINTEGERS integer values */
  if (!haserr && state->merindexfpout != NULL)
  {
    outputbytewiseUlongvalue(state->merindexfpout,
                             (unsigned long) state->mersize);
    outputbytewiseUlongvalue(state->merindexfpout,(unsigned long) alphasize);
  }
  gt_fa_xfclose(state->merindexfpout);
  gt_fa_xfclose(state->countsfilefpout);
  GT_FREEARRAY(&state->occdistribution,Countwithpositions);
  gt_free(state->currentmer);
  gt_free(state->bytebuffer);
  GT_FREEARRAY(&state->largecounts,Largecount);
  gt_encseq_reader_delete(state->esrspace);
  gt_free(state);
  return haserr ? -1 : 0;
}

int gt_merstatistics(const char *inputindex,
                  unsigned long mersize,
                  unsigned long minocc,
                  unsigned long maxocc,
                  const char *storeindex,
                  bool storecounts,
                  bool scanfile,
                  bool performtest,
                  GtLogger *logger,
                  GtError *err)
{
  bool haserr = false;
  Sequentialsuffixarrayreader *ssar;

  gt_error_check(err);
  ssar = gt_newSequentialsuffixarrayreaderfromfile(inputindex,
                                                SARR_LCPTAB |
                                                SARR_SUFTAB |
                                                SARR_ESQTAB,
                                                (scanfile && !performtest)
                                                  ? SEQ_scan : SEQ_mappedboth,
                                                logger,
                                                err);
  if (ssar == NULL)
  {
    haserr = true;
  }
  if (!haserr)
  {
    if (enumeratelcpintervals(inputindex,
                              ssar,
                              storeindex,
                              storecounts,
                              mersize,
                              minocc,
                              maxocc,
                              performtest,
                              logger,
                              err) != 0)
    {
      haserr = true;
    }
  }
  if (ssar != NULL)
  {
    gt_freeSequentialsuffixarrayreader(&ssar);
  }
  return haserr ? -1 : 0;
}
