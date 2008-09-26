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

#include "core/str.h"
#include "core/unused_api.h"
#include "core/arraydef.h"
#include "esa-seqread.h"
#include "verbose-def.h"
#include "spacedef.h"
#include "tyr-mkindex.h"

struct Dfsinfo /* information stored for each node of the lcp interval tree */
{
  Seqpos depth,
         leftmostleaf,
         rightmostleaf,
         suftabrightmostleaf,
         lcptabrightmostleafplus1;
};

typedef int (*Processoccurrencecount)(Seqpos,Seqpos,void *);

typedef struct _listSeqpos
{
  Seqpos position; 
  struct _listSeqpos *nextptr;
} ListSeqpos;

typedef struct
{
  Seqpos occcount;
  ListSeqpos *positionlist;
} Countwithpositions;

DECLAREARRAYSTRUCT(Countwithpositions);

struct Dfsstate /* global information */
{
  Seqpos searchlength,
         minocc,
         maxocc,
         totallength;
  const Encodedsequence *encseq;
  Readmode readmode;
  Processoccurrencecount processoccurrencecount;
  ArrayCountwithpositions occdistribution;
  FILE *merindexfpout;
};

#include "esa-dfs.h"

static /*@null@*/ ListSeqpos *insertListSeqpos(ListSeqpos *liststart,
                                               Seqpos position)
{
  ListSeqpos *newnode;

  ALLOCASSIGNSPACE(newnode,NULL,ListSeqpos,1);
  newnode->position = position;
  newnode->nextptr = liststart;
  return newnode;
}

static bool decideifocc(const Dfsstate *state,Seqpos countocc)
{
  if(state->minocc > 0)
  {
    if(state->maxocc > 0)
    {
      if(countocc >= state->minocc && countocc <= state->maxocc)
      {
        return true;
      }
    } else
    {
      if(countocc >= state->minocc)
      {
        return true;
      }
    }
  } else
  {
    if(state->maxocc > 0)
    {
      if(countocc <= state->maxocc)
      {
        return true;
      }
    } 
  }
  return false;
}

static void incrementdistribcounts(ArrayCountwithpositions *occdistribution,
                                   Seqpos countocc,Seqpos value)
{
  if(countocc >= occdistribution->allocatedCountwithpositions)
  {
    const Seqpos addamount = (Seqpos) 128;
    unsigned long idx;

    ALLOCASSIGNSPACE(occdistribution->spaceCountwithpositions,
                     occdistribution->spaceCountwithpositions,
                     Countwithpositions,countocc+addamount);
    for(idx=occdistribution->allocatedCountwithpositions; 
        idx<countocc+addamount; idx++)
    {
      occdistribution->spaceCountwithpositions[idx].occcount = 0; 
      occdistribution->spaceCountwithpositions[idx].positionlist = NULL;
    }
    occdistribution->allocatedCountwithpositions = countocc+addamount;
  } 
  if(countocc + 1 > occdistribution->nextfreeCountwithpositions)
  {
    occdistribution->nextfreeCountwithpositions = countocc+1;
  }
  occdistribution->spaceCountwithpositions[countocc].occcount += value;
}

static int adddistpos2distribution(Seqpos countocc,Seqpos position,
                                   void *adddistposinfo)
{
  Dfsstate *state = (Dfsstate *) adddistposinfo;

  incrementdistribcounts(&state->occdistribution,countocc,(Seqpos) 1);
  if(decideifocc(state,countocc))
  {
    state->occdistribution.spaceCountwithpositions[countocc].positionlist
      = insertListSeqpos(state->occdistribution.
                              spaceCountwithpositions[countocc].positionlist,
                       position);
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

static bool containsspecial(GT_UNUSED const Encodedsequence *encseq,
                            GT_UNUSED Seqpos startindex,
                            GT_UNUSED Seqpos len)
{
  return false;
}

static int processleafedge(GT_UNUSED bool firstsucc,
                           GT_UNUSED Seqpos fatherdepth,
                           Dfsinfo *father,
                           Seqpos leafnumber,
                           GT_UNUSED Dfsstate *state,
                           GT_UNUSED GtError *err)
{
  if (father->depth < state->searchlength &&
      leafnumber + state->searchlength <=
      state->totallength &&
      !containsspecial(state->encseq,
                      leafnumber + father->depth,
                      state->searchlength - father->depth))
  {
    if (state->processoccurrencecount((Seqpos) 1,leafnumber,state) != 0)
    {
      return -1;
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

static int enumeratelcpintervals(Sequentialsuffixarrayreader *ssar,
                                 const GtStr *str_storeindex,
                                 unsigned long searchlength,
                                 unsigned long minocc,
                                 unsigned long maxocc,
                                 Verboseinfo *verboseinfo,
                                 GtError *err)
{
  Dfsstate state;
  bool haserr = false;
  char *merindexoutfilename = NULL;

  state.searchlength = (Seqpos) searchlength;
  state.encseq = encseqSequentialsuffixarrayreader(ssar);
  state.totallength = getencseqtotallength(state.encseq);
  state.readmode = readmodeSequentialsuffixarrayreader(ssar);
  state.minocc = minocc;
  state.maxocc = maxocc;
  INITARRAY(&state.occdistribution,Countwithpositions);
  if(gt_str_length(str_storeindex) == 0)
  {
    state.merindexfpout = NULL;
    state.processoccurrencecount = adddistpos2distribution;
  } else
  {  
    state.merindexfpout = NULL;
    assert(false);
  }
  if (depthfirstesa(ssar,
                    allocateDfsinfo,
                    freeDfsinfo,
                    processleafedge,
                    NULL,
                    NULL,
                    assignleftmostleaf,
                    assignrightmostleaf,
                    &state,
                    verboseinfo,
                    err) != 0)
  {
    haserr = true;
  }
  FREESPACE(merindexoutfilename);
  FREEARRAY(&state.occdistribution,Countwithpositions);
  return haserr ? -1 : 0;
}

int merstatistics(const GtStr *str_inputindex,
                  unsigned long searchlength,
                  unsigned long minocc,
                  unsigned long maxocc,
                  const GtStr *str_storeindex,
                  GT_UNUSED bool storecounts,
                  bool scanfile,
                  bool verbose,
                  GtError *err)
{
  bool haserr = false;
  Sequentialsuffixarrayreader *ssar;

  gt_error_check(err);
  ssar = newSequentialsuffixarrayreaderfromfile(str_inputindex,
                                                SARR_LCPTAB |
                                                SARR_SUFTAB |
                                                SARR_ESQTAB,
                                                scanfile
                                                  ? SEQ_scan : SEQ_mappedboth,
                                                err);
  if (ssar == NULL)
  {
    haserr = true;
  }
  if (!haserr)
  {
    Verboseinfo *verboseinfo = newverboseinfo(verbose);
    if (enumeratelcpintervals(ssar,
                              str_storeindex,
                              searchlength,
                              minocc,
                              maxocc,
                              verboseinfo,
                              err) != 0)
    {
      haserr = true;
    }
    freeverboseinfo(&verboseinfo);
  }
  if (ssar != NULL)
  {
    freeSequentialsuffixarrayreader(&ssar);
  }
  return haserr ? -1 : 0;
}
