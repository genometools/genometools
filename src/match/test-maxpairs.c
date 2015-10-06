/*
  Copyright (c) 2007-2013 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
  Copyright (c) 2007-2013 Center for Bioinformatics, University of Hamburg

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
#include "core/alphabet.h"
#include "core/array.h"
#include "core/arraydef.h"
#include "core/divmodmul.h"
#include "core/encseq.h"
#include "core/format64.h"
#include "core/logger.h"
#include "core/unused_api.h"
#include "core/ma_api.h"
#include "core/yarandom.h"
#include "esa-mmsearch.h"
#include "echoseq.h"
#include "esa-maxpairs.h"
#include "sfx-sain.h"
#include "test-maxpairs.h"

static int gt_constructsarrandrunmaxpairs(GtUchar *dbseqquery,
                                          GtUword dbseqquerylen,
                                          GtReadmode readmode,
                                          unsigned int searchlength,
                                          unsigned int numofchars,
                                          GtProcessmaxpairs processmaxpairs,
                                          void *processmaxpairsinfo,
                                          GtError *err)
{
  bool haserr = false;
  GtSainSufLcpIterator *suflcpiterator;

  suflcpiterator = gt_sain_suf_lcp_iterator_new(true,
                                                dbseqquery,
                                                dbseqquerylen,
                                                readmode,
                                                (GtUword) numofchars,
                                                err);
  if (suflcpiterator == NULL)
  {
    haserr = true;
  } else
  {
    if (gt_enumeratemaxpairs_sain(suflcpiterator,
                                  searchlength,
                                  processmaxpairs,
                                  processmaxpairsinfo,
                                  err) != 0)
    {
      haserr = true;
    }
  }
  gt_sain_suf_lcp_iterator_delete(suflcpiterator);
  return haserr ? -1 : 0;
}

static GtUword gt_samplesubstring(bool replacespecialchars,
                                  GtUchar *seqspace,
                                  const GtEncseq *encseq,
                                  GtUword substringlength)
{
  GtUword start, totallength = gt_encseq_total_length(encseq);

  start = (GtUword) (random() % totallength);
  if (start + substringlength > totallength)
  {
    substringlength = totallength - start;
  }
  gt_assert(substringlength > 0);
  gt_encseq_extract_encoded(encseq,seqspace,start,start+substringlength-1);
  if (replacespecialchars)
  {
    GtUword idx;
    const unsigned int numofchars = gt_encseq_alphabetnumofchars(encseq);

    for (idx = 0; idx < substringlength; idx++)
    {
      if (ISSPECIAL(seqspace[idx]))
      {
        seqspace[idx] = (GtUchar) (random() % numofchars);
      }
    }
  }
  return substringlength;
}

typedef struct
{
  GtUword len,
          dbstart,
          querystart;
  uint64_t queryseqnum;
} Substringmatch;

static void gt_storemaxmatchquery(void *info, const GtQuerymatch *querymatch)
{
  GtArray *tab = (GtArray *) info;
  Substringmatch subm;

  subm.len = gt_querymatch_querylen(querymatch);
  subm.dbstart = gt_querymatch_dbstart(querymatch);
  subm.querystart = gt_querymatch_querystart(querymatch);
  subm.queryseqnum = gt_querymatch_queryseqnum(querymatch);
  gt_array_add(tab,subm);
}

typedef struct
{
  GtArray *results;
  GtUword dblen, *querymarkpos, querylen, numofquerysequences;
} GtMaxmatchselfinfo;

static int gt_storemaxmatchself(void *info,
                                GT_UNUSED const GtGenericEncseq *genericencseq,
                                GtUword len,
                                GtUword pos1,
                                GtUword pos2,
                                GT_UNUSED GtError *err)
{
  GtMaxmatchselfinfo *maxmatchselfinfo = (GtMaxmatchselfinfo *) info;
  GtUword dbstart, querystart;

  if (pos1 < pos2)
  {
    dbstart = pos1;
    querystart = pos2;
  } else
  {
    dbstart = pos2;
    querystart = pos1;
  }
  if (dbstart < maxmatchselfinfo->dblen &&
      maxmatchselfinfo->dblen < querystart)
  {
    Substringmatch subm;
    GtUword pos;

    subm.len = len;
    subm.dbstart = dbstart;
    pos = querystart - (maxmatchselfinfo->dblen + 1);
    if (maxmatchselfinfo->querymarkpos == NULL)
    {
      subm.queryseqnum = 0;
      subm.querystart = pos;
    } else
    {
      GtUword queryseqnum
        = gt_encseq_sep2seqnum(maxmatchselfinfo->querymarkpos,
                               maxmatchselfinfo->numofquerysequences,
                               maxmatchselfinfo->querylen,
                               pos);
      if (queryseqnum == maxmatchselfinfo->numofquerysequences)
      {
        return -1;
      }
      if (queryseqnum == 0)
      {
        subm.querystart = pos;
      } else
      {
        subm.querystart = pos -
                          (maxmatchselfinfo->querymarkpos[queryseqnum-1] + 1);
      }
      subm.queryseqnum = (uint64_t) queryseqnum;
    }
    gt_array_add(maxmatchselfinfo->results,subm);
  }
  return 0;
}

static int gt_orderSubstringmatch(const void *a,const void *b)
{
  Substringmatch *m1 = (Substringmatch *) a,
                 *m2 = (Substringmatch *) b;

  if (m1->queryseqnum < m2->queryseqnum)
  {
    return -1;
  }
  if (m1->queryseqnum > m2->queryseqnum)
  {
    return 1;
  }
  if (m1->querystart < m2->querystart)
  {
    return -1;
  }
  if (m1->querystart > m2->querystart)
  {
    return 1;
  }
  if (m1->dbstart < m2->dbstart)
  {
    return -1;
  }
  if (m1->dbstart > m2->dbstart)
  {
    return 1;
  }
  if (m1->len < m2->len)
  {
    return -1;
  }
  if (m1->len > m2->len)
  {
    return 1;
  }
  return 0;
}

static int gt_showSubstringmatch(void *a, GT_UNUSED void *info,
                                 GT_UNUSED GtError *err)
{
  Substringmatch *m = (Substringmatch *) a;

  printf(""GT_WU" "GT_WU" " Formatuint64_t " "GT_WU"\n",
           m->len,
           m->dbstart,
           PRINTuint64_tcast(m->queryseqnum),
           m->querystart);
  return 0;
}

static GtUword *gt_sequence2markpositions(GtUword *numofsequences,
                                          const GtUchar *seq,
                                          GtUword seqlen)
{
  GtUword *spacemarkpos, idx, allocatedmarkpos, nextfreemarkpos;

  *numofsequences = 1UL;
  for (idx=0; idx<seqlen; idx++)
  {
    if (seq[idx] == (GtUchar) SEPARATOR)
    {
      (*numofsequences)++;
    }
  }
  if (*numofsequences == 1UL)
  {
    return NULL;
  }
  allocatedmarkpos = (*numofsequences)-1;
  spacemarkpos = gt_malloc(sizeof *spacemarkpos * allocatedmarkpos);
  for (idx=0, nextfreemarkpos = 0; idx<seqlen; idx++)
  {
    if (seq[idx] == (GtUchar) SEPARATOR)
    {
      spacemarkpos[nextfreemarkpos++] = idx;
    }
  }
  return spacemarkpos;
}

int gt_testmaxpairs(const char *indexname,
                    GtUword samples,
                    unsigned int minlength,
                    GtUword substringlength,
                    GtLogger *logger,
                    GtError *err)
{
  bool haserr = false;
  GtEncseq *encseq;
  GtUword idx, totallength = 0, dblen, querylen;
  GtUchar *dbseq = NULL, *query = NULL;
  GtAlphabet *alphabet = NULL;
  GtArray *tabmaxquerymatches;
  GtMaxmatchselfinfo maxmatchselfinfo;
  GtEncseqLoader *el;

  gt_logger_log(logger,"draw "GT_WU" samples",samples);
  el = gt_encseq_loader_new();
  gt_encseq_loader_do_not_require_des_tab(el);
  gt_encseq_loader_do_not_require_ssp_tab(el);
  gt_encseq_loader_do_not_require_sds_tab(el);
  gt_encseq_loader_set_logger(el, logger);
  encseq = gt_encseq_loader_load(el, indexname, err);
  gt_encseq_loader_delete(el);

  if (encseq == NULL)
  {
    haserr = true;
  } else
  {
    totallength = gt_encseq_total_length(encseq);
  }
  if (!haserr)
  {
    if (substringlength > totallength/2)
    {
      substringlength = totallength/2;
    }
    /* use one memory area for dbseq and query, since in one method
       we concatenate both sequences */
    dbseq = gt_malloc(sizeof *dbseq * (1 + GT_MULT2(substringlength)));
    alphabet = gt_encseq_alphabet(encseq);
  }
  for (idx = 0; idx < samples && !haserr; idx++)
  {
    dblen = gt_samplesubstring(false,dbseq,encseq,substringlength);
    gt_assert(dbseq != NULL);
    query = dbseq + dblen + 1;
    gt_assert(query != NULL);
    dbseq[dblen] = SEPARATOR;
    querylen = gt_samplesubstring(true,query,encseq,substringlength);
    if (querylen < (GtUword) minlength || dblen < (GtUword) minlength ||
        dbseq[0] == SEPARATOR || query[0] == SEPARATOR ||
        dbseq[substringlength-1] == SEPARATOR ||
        query[substringlength-1] == SEPARATOR)
    {
      continue;
    }
    gt_logger_log(logger,"run query match for dblen="GT_WU""
                         ",querylen= "GT_WU", minlength=%u",
                         dblen,
                         querylen,
                         minlength);
    tabmaxquerymatches = gt_array_new(sizeof (Substringmatch));
    if (gt_sarrquerysubstringmatch(dbseq,
                                   dblen,
                                   query,
                                   querylen,
                                   minlength,
                                   alphabet,
                                   gt_storemaxmatchquery,
                                   tabmaxquerymatches,
                                   logger,
                                   err) != 0)
    {
      haserr = true;
      break;
    }
    gt_logger_log(logger,"run self match for dblen="GT_WU""
                         ",querylen= "GT_WU", minlength=%u",
                         dblen,
                         querylen,
                         minlength);
    maxmatchselfinfo.results = gt_array_new(sizeof (Substringmatch));
    maxmatchselfinfo.dblen = dblen;
    maxmatchselfinfo.querylen = querylen;
    maxmatchselfinfo.querymarkpos
      = gt_sequence2markpositions(&maxmatchselfinfo.numofquerysequences,
                                  query,querylen);
    if (gt_constructsarrandrunmaxpairs(
                 dbseq,
                 dblen + 1 + querylen,
                 GT_READMODE_FORWARD,
                 minlength,
                 gt_encseq_alphabetnumofchars(encseq),
                 gt_storemaxmatchself,
                 &maxmatchselfinfo,
                 err) != 0)
    {
      haserr = true;
      break;
    }
    gt_array_sort(tabmaxquerymatches,gt_orderSubstringmatch);
    gt_array_sort(maxmatchselfinfo.results,gt_orderSubstringmatch);
    if (!gt_array_equal(tabmaxquerymatches,maxmatchselfinfo.results,
                        gt_orderSubstringmatch))
    {
      const GtUword width = 60UL;
      fprintf(stderr,"failure for query of length "GT_WU"\n",querylen);
      fprintf(stderr,"querymatches\n");
      (void) gt_array_iterate(tabmaxquerymatches,gt_showSubstringmatch,NULL,
                              err);
      fprintf(stderr,"dbmatches\n");
      (void) gt_array_iterate(maxmatchselfinfo.results,gt_showSubstringmatch,
                              NULL,err);
      gt_symbolstring2fasta(stderr,"dbseq", alphabet, dbseq, dblen, width);
      gt_symbolstring2fasta(stderr,"queryseq", alphabet, query, querylen,width);
      exit(GT_EXIT_PROGRAMMING_ERROR);
    }
    gt_free(maxmatchselfinfo.querymarkpos);
    gt_array_delete(tabmaxquerymatches);
    gt_array_delete(maxmatchselfinfo.results);
  }
  gt_free(dbseq);
  gt_encseq_delete(encseq);
  return haserr ? -1 : 0;
}
