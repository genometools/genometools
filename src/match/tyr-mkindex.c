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
#include "core/divmodmul.h"
#include "core/str.h"
#include "core/unused_api.h"
#include "core/fa.h"
#include "esa-seqread.h"
#include "alphadef.h"
#include "verbose-def.h"
#include "spacedef.h"
#include "format64.h"
#include "esa-mmsearch-def.h"
#include "opensfxfile.h"
#include "tyr-basic.h"
#include "tyr-mkindex.h"
#include "echoseq.h"

struct Dfsinfo /* information stored for each node of the lcp interval tree */
{
  Seqpos leftmostleaf,
         rightmostleaf,
         suftabrightmostleaf,
         lcptabrightmostleafplus1;
};

typedef int (*Processoccurrencecount)(unsigned long,Seqpos,void *,GtError *);

typedef struct _listSeqpos
{
  Seqpos position;
  struct _listSeqpos *nextptr;
} ListSeqpos;

typedef struct
{
  unsigned long occcount;
  ListSeqpos *positionlist;
} Countwithpositions;

GT_DECLAREARRAYSTRUCT(Countwithpositions);

struct Dfsstate /* global information */
{
  Seqpos mersize,
         totallength;
  unsigned long minocc,
                maxocc;
  const Encodedsequence *encseq;
  Readmode readmode;
  Processoccurrencecount processoccurrencecount;
  GtArrayCountwithpositions occdistribution;
  FILE *merindexfpout,
       *countsfilefpout;
  bool moveforward;
  Encodedsequencescanstate *esrspace;
  bool performtest;
  bool storecounts;
  GtUchar *bytebuffer;
  unsigned long sizeofbuffer;
  GtArrayLargecount largecounts;
  unsigned long countoutputmers;
  const Seqpos *suftab; /* only necessary for performtest */
  GtUchar *currentmer;    /* only necessary for performtest */
};

#include "esa-dfs.h"

static uint64_t bruteforcecountnumofmers(const Dfsstate *state)
{
  Seqpos idx;
  uint64_t numofmers = 0;

  for (idx=0; idx <= state->totallength - state->mersize; idx++)
  {
    if (!containsspecial(state->encseq,
                         state->moveforward,
                         state->esrspace,
                         idx,
                         state->mersize))
    {
      numofmers++;
    }
  }
  return numofmers;
}

static void checknumofmers(const Dfsstate *state,
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

static void checknumberofoccurrences(const Dfsstate *dfsstate,
                                     unsigned long countocc,
                                     Seqpos position)
{
  MMsearchiterator *mmsi;
  Seqpos idx, bfcount;

  for (idx = 0; idx < dfsstate->mersize; idx++)
  {
    dfsstate->currentmer[idx] = getencodedchar(dfsstate->encseq,position+idx,
                                               dfsstate->readmode);
  }
  mmsi = newmmsearchiterator(dfsstate->encseq,
                             dfsstate->suftab,
                             0,
                             dfsstate->totallength,
                             0,
                             dfsstate->readmode,
                             dfsstate->currentmer,
                             (unsigned long) dfsstate->mersize);
  bfcount = countmmsearchiterator(mmsi);
  if (bfcount != (Seqpos) countocc)
  {
    fprintf(stderr,"bfcount = " FormatSeqpos " != %lu = countocc\n",
                   PRINTSeqposcast(bfcount),
                   countocc);
    exit(GT_EXIT_PROGRAMMING_ERROR);
  }
  freemmsearchiterator(&mmsi);
}

static /*@null@*/ ListSeqpos *insertListSeqpos(ListSeqpos *liststart,
                                               Seqpos position)
{
  ListSeqpos *newnode;

  ALLOCASSIGNSPACE(newnode,NULL,ListSeqpos,1);
  newnode->position = position;
  newnode->nextptr = liststart;
  return newnode;
}

static void wrapListSeqpos(ListSeqpos *node)
{
  ListSeqpos *tmpnext;

  if (node != NULL)
  {
    while (true)
    {
      if (node->nextptr == NULL)
      {
        FREESPACE(node);
        break;
      }
      tmpnext = node->nextptr;
      FREESPACE(node);
      node = tmpnext;
    }
  }
}

static void showListSeqpos(const Encodedsequence *encseq,
                           Seqpos mersize,
                           const ListSeqpos *node)
{
  const ListSeqpos *tmp;

  for (tmp = node; tmp != NULL; tmp = tmp->nextptr)
  {
    fprintfencseq(stdout,encseq,tmp->position,mersize);
    (void) putchar((int) '\n');
  }
}

static bool decideifocc(const Dfsstate *state,unsigned long countocc)
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

static void showmerdistribution(const Dfsstate *state)
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
        showListSeqpos(state->encseq,
                       state->mersize,
                       state->occdistribution.spaceCountwithpositions[countocc].
                                              positionlist);
        wrapListSeqpos(state->occdistribution.spaceCountwithpositions[countocc].
                                              positionlist);
      }
    }
  }
}

static void showfinalstatistics(const Dfsstate *state,
                                const GtStr *inputindex,
                                Verboseinfo *verboseinfo)
{
  uint64_t dnumofmers = addupdistribution(&state->occdistribution);

  if (state->performtest)
  {
    checknumofmers(state,dnumofmers);
  }
  showverbose(verboseinfo,
              "the following output refers to the set of all sequences");
  showverbose(verboseinfo,
              "represented by the index \"%s\"",gt_str_get(inputindex));
  showverbose(verboseinfo,
              "number of %lu-mers in the sequences not containing a "
              "wildcard: " Formatuint64_t,
              (unsigned long) state->mersize,
              PRINTuint64_tcast(dnumofmers));
  showverbose(verboseinfo,
              "show the distribution of the number of occurrences of %lu-mers",
               (unsigned long) state->mersize);
  showverbose(verboseinfo,"not containing a wildcard as rows of the form "
              "i d where");
  showverbose(verboseinfo,
              "d is the number of events that a %lu-mer occurs exactly i times",
              (unsigned long) state->mersize);
  showmerdistribution(state);
}

static void incrementdistribcounts(GtArrayCountwithpositions *occdistribution,
                                   unsigned long countocc,unsigned long value)
{
  if (countocc >= occdistribution->allocatedCountwithpositions)
  {
    const Seqpos addamount = (Seqpos) 128;
    unsigned long idx;

    ALLOCASSIGNSPACE(occdistribution->spaceCountwithpositions,
                     occdistribution->spaceCountwithpositions,
                     Countwithpositions,countocc+addamount);
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
                                   Seqpos position,
                                   void *adddistposinfo,
                                   GT_UNUSED GtError *err)
{
  Dfsstate *state = (Dfsstate *) adddistposinfo;

  incrementdistribcounts(&state->occdistribution,countocc,1UL);
  if (decideifocc(state,countocc))
  {
    state->occdistribution.spaceCountwithpositions[countocc].positionlist
      = insertListSeqpos(state->occdistribution.
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

static int outputsortedstring2indexviafileptr(const Encodedsequence *encseq,
                                              Seqpos mersize,
                                              GtUchar *bytebuffer,
                                              unsigned long sizeofbuffer,
                                              FILE *merindexfpout,
                                              FILE *countsfilefpout,
                                              Seqpos position,
                                              unsigned long countocc,
                                              GtArrayLargecount *largecounts,
                                              unsigned long countoutputmers,
                                              GtError *err)
{
  sequence2bytecode(bytebuffer,encseq,position,(Seqpos) mersize);
  if (fwrite(bytebuffer,sizeof(*bytebuffer),(size_t) sizeofbuffer,merindexfpout)
            != (size_t) sizeofbuffer)
  {
    gt_error_set(err,"cannot write %lu items of size %u: errormsg=\"%s\"",
                  (unsigned long) sizeofbuffer,
                  (unsigned int) sizeof (*bytebuffer),
                  strerror(errno));
    return -1;
  }
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
    if (fwrite(&smallcount, sizeof (smallcount),(size_t) 1,countsfilefpout)
              != (size_t) 1)
    {
      gt_error_set(err,"cannot write 1 item of size %u: errormsg=\"%s\"",
                  (unsigned int) sizeof (smallcount),
                  strerror(errno));
      return -1;
    }
  }
  return 0;
}

static int outputsortedstring2index(unsigned long countocc,
                                    Seqpos position,
                                    void *adddistposinfo,
                                    GtError *err)
{
  Dfsstate *state = (Dfsstate *) adddistposinfo;

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

static Dfsinfo *allocateDfsinfo(GT_UNUSED Dfsstate *state)
{
  Dfsinfo *dfsinfo;

  ALLOCASSIGNSPACE(dfsinfo,NULL,Dfsinfo,1);
  return dfsinfo;
}

static void freeDfsinfo(Dfsinfo *dfsinfo, GT_UNUSED Dfsstate *state)
{
  FREESPACE(dfsinfo);
}

#ifdef WITHcontainsspecial2
static bool containsspecial2(const Encodedsequence *encseq,
                     GT_UNUSED bool moveforward,
                     GT_UNUSED Encodedsequencescanstate *esrspace,
                     Seqpos startpos,
                     Seqpos len)
{
  Seqpos pos;
  bool result = false, result2;

  for (pos=startpos; pos<startpos+len; pos++)
  {
    if (ISSPECIAL(getencodedchar(encseq,pos,Forwardmode)))
    {
      result = true;
      break;
    }
  }
  result2 = containsspecial(encseq,
                            moveforward,
                            esrspace,
                            startpos,len);
  if ((result && !result2) || (!result && result2))
  {
    fprintf(stderr,"pos = %lu, len = %lu: result = %s != %s = result2\n",
                    (unsigned long) startpos,
                    (unsigned long) len,
                    result ? "true" : "false",
                    result2 ? "true" : "false");
    exit(GT_EXIT_PROGRAMMING_ERROR);
  }
  return result;
}
#endif

static int processleafedge(GT_UNUSED bool firstsucc,
                           Seqpos fatherdepth,
                           GT_UNUSED Dfsinfo *father,
                           Seqpos leafnumber,
                           Dfsstate *state,
                           GtError *err)
{
  gt_error_check(err);
  if (fatherdepth < state->mersize &&
      leafnumber + state->mersize <= state->totallength &&
      !containsspecial(state->encseq,
                       state->moveforward,
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

static int processcompletenode(Seqpos nodeptrdepth,
                               Dfsinfo *nodeptr,
                               Seqpos nodeptrminusonedepth,
                               Dfsstate *state,
                               GtError *err)
{
  gt_error_check(err);
  if (state->mersize <= nodeptrdepth)
  {
    Seqpos fatherdepth;

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

static void assignleftmostleaf(Dfsinfo *dfsinfo,Seqpos leftmostleaf,
                               GT_UNUSED Dfsstate *dfsstate)
{
  dfsinfo->leftmostleaf = leftmostleaf;
}

static void assignrightmostleaf(Dfsinfo *dfsinfo,Seqpos currentindex,
                                Seqpos previoussuffix,Seqpos currentlcp,
                                GT_UNUSED Dfsstate *dfsstate)
{
  dfsinfo->rightmostleaf = currentindex;
  dfsinfo->suftabrightmostleaf = previoussuffix;
  dfsinfo->lcptabrightmostleafplus1 = currentlcp;
}

#define MEGABYTES(V)  ((double) (V)/((((unsigned long) 1) << 20) - 1))

static void outputbytewiseUlongvalue(FILE *fpout,unsigned long value)
{
  size_t i;

  for (i=0; i < sizeof (value); i++)
  {
    (void) putc((int) (value & UCHAR_MAX),fpout);
    value >>= 8;
  }
}

static int enumeratelcpintervals(const GtStr *str_inputindex,
                                 Sequentialsuffixarrayreader *ssar,
                                 const GtStr *str_storeindex,
                                 bool storecounts,
                                 unsigned long mersize,
                                 unsigned long minocc,
                                 unsigned long maxocc,
                                 bool performtest,
                                 Verboseinfo *verboseinfo,
                                 GtError *err)
{
  Dfsstate state;
  bool haserr = false;
  unsigned int alphasize;

  gt_error_check(err);
  GT_INITARRAY(&state.occdistribution,Countwithpositions);
  state.esrspace = newEncodedsequencescanstate();
  state.mersize = (Seqpos) mersize;
  state.encseq = encseqSequentialsuffixarrayreader(ssar);
  alphasize = getencseqAlphabetnumofchars(state.encseq);
  state.readmode = readmodeSequentialsuffixarrayreader(ssar);
  state.storecounts = storecounts;
  state.minocc = minocc;
  state.maxocc = maxocc;
  state.moveforward = ISDIRREVERSE(state.readmode) ? false : true;
  state.totallength = getencseqtotallength(state.encseq);
  state.performtest = performtest;
  state.countoutputmers = 0;
  state.merindexfpout = NULL;
  state.countsfilefpout = NULL;
  GT_INITARRAY(&state.largecounts,Largecount);
  if (gt_str_length(str_storeindex) == 0)
  {
    state.sizeofbuffer = 0;
    state.bytebuffer = NULL;
  } else
  {
    state.sizeofbuffer = MERBYTES(mersize);
    ALLOCASSIGNSPACE(state.bytebuffer,NULL,GtUchar,state.sizeofbuffer);
  }
  if (performtest)
  {
    ALLOCASSIGNSPACE(state.currentmer,NULL,GtUchar,state.mersize);
    state.suftab = suftabSequentialsuffixarrayreader(ssar);
  } else
  {
    state.currentmer = NULL;
    state.suftab = NULL;
  }
  if (state.mersize > state.totallength)
  {
    gt_error_set(err,"mersize " FormatSeqpos " > " FormatSeqpos
                     " = totallength not allowed",
                 PRINTSeqposcast(state.mersize),
                 PRINTSeqposcast(state.totallength));
    haserr = true;
  } else
  {
    if (gt_str_length(str_storeindex) == 0)
    {
      state.processoccurrencecount = adddistpos2distribution;
    } else
    {
      state.merindexfpout = opensfxfile(str_storeindex,
                                        MERSUFFIX,
                                        "wb",
                                        err);
      if (state.merindexfpout == NULL)
      {
        haserr = true;
      } else
      {
        if (state.storecounts)
        {
          state.countsfilefpout = opensfxfile(str_storeindex,
                                              COUNTSSUFFIX,"wb",err);
          if (state.countsfilefpout == NULL)
          {
            haserr = true;
          }
        }
      }
      state.processoccurrencecount = outputsortedstring2index;
    }
    if (!haserr)
    {
      if (depthfirstesa(ssar,
                        allocateDfsinfo,
                        freeDfsinfo,
                        processleafedge,
                        NULL,
                        processcompletenode,
                        assignleftmostleaf,
                        assignrightmostleaf,
                        &state,
                        verboseinfo,
                        err) != 0)
      {
        haserr = true;
      }
      if (gt_str_length(str_storeindex) == 0)
      {
        showfinalstatistics(&state,str_inputindex,verboseinfo);
      }
    }
    if (!haserr)
    {
      if (state.countsfilefpout != NULL)
      {
        showverbose(verboseinfo,"write %lu mercounts > %lu to file \"%s%s\"",
                    state.largecounts.nextfreeLargecount,
                    (unsigned long) MAXSMALLMERCOUNT,
                    gt_str_get(str_storeindex),
                    COUNTSSUFFIX);
        if (fwrite(state.largecounts.spaceLargecount,
                  sizeof (Largecount),
                  (size_t) state.largecounts.nextfreeLargecount,
                  state.countsfilefpout) !=
                  (size_t) state.largecounts.nextfreeLargecount)
        {
          gt_error_set(err,
                       "cannot write %lu items of size %u: errormsg=\"%s\"",
                       (unsigned long) state.largecounts.nextfreeLargecount,
                       (unsigned int) sizeof (Largecount),
                       strerror(errno));
          haserr = true;
        }
      }
    }
    if (!haserr)
    {
      showverbose(verboseinfo,"number of %lu-mers in index: %lu",
                  mersize,
                  state.countoutputmers);
      showverbose(verboseinfo,"index size: %.2f megabytes\n",
                  MEGABYTES(state.countoutputmers * state.sizeofbuffer +
                            sizeof (unsigned long) * EXTRAINTEGERS));
    }
  }
  /* now out EXTRAINTEGERS integer values */
  if (!haserr && state.merindexfpout != NULL)
  {
    outputbytewiseUlongvalue(state.merindexfpout,
                             (unsigned long) state.mersize);
    outputbytewiseUlongvalue(state.merindexfpout,(unsigned long) alphasize);
  }
  gt_fa_xfclose(state.merindexfpout);
  gt_fa_xfclose(state.countsfilefpout);
  GT_FREEARRAY(&state.occdistribution,Countwithpositions);
  FREESPACE(state.currentmer);
  FREESPACE(state.bytebuffer);
  GT_FREEARRAY(&state.largecounts,Largecount);
  freeEncodedsequencescanstate(&state.esrspace);
  return haserr ? -1 : 0;
}

int merstatistics(const GtStr *str_inputindex,
                  unsigned long mersize,
                  unsigned long minocc,
                  unsigned long maxocc,
                  const GtStr *str_storeindex,
                  bool storecounts,
                  bool scanfile,
                  bool performtest,
                  Verboseinfo *verboseinfo,
                  GtError *err)
{
  bool haserr = false;
  Sequentialsuffixarrayreader *ssar;

  gt_error_check(err);
  ssar = newSequentialsuffixarrayreaderfromfile(str_inputindex,
                                                SARR_LCPTAB |
                                                SARR_SUFTAB |
                                                SARR_ESQTAB,
                                                (scanfile && !performtest)
                                                  ? SEQ_scan : SEQ_mappedboth,
                                                err);
  if (ssar == NULL)
  {
    haserr = true;
  }
  if (!haserr)
  {
    if (enumeratelcpintervals(str_inputindex,
                              ssar,
                              str_storeindex,
                              storecounts,
                              mersize,
                              minocc,
                              maxocc,
                              performtest,
                              verboseinfo,
                              err) != 0)
    {
      haserr = true;
    }
  }
  if (ssar != NULL)
  {
    freeSequentialsuffixarrayreader(&ssar);
  }
  return haserr ? -1 : 0;
}
