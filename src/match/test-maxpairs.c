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

#include <limits.h>
#include "core/alphabet.h"
#include "core/array.h"
#include "core/arraydef.h"
#include "core/divmodmul.h"
#include "core/encseq.h"
#include "core/format64.h"
#include "core/logger.h"
#include "core/timer_api.h"
#include "core/unused_api.h"
#include "core/ma_api.h"
#include "esa-mmsearch.h"
#include "echoseq.h"
#include "sfx-suffixer.h"
#include "sfx-apfxlen.h"
#include "esa-maxpairs.h"
#include "esa-seqread.h"

typedef struct
{
  unsigned int minlength;
  GtEncseq *encseq;
  Processmaxpairs processmaxpairs;
  void *processmaxpairsinfo;
} Substringmatchinfo;

static int constructsarrandrunmaxpairs(
                 Substringmatchinfo *ssi,
                 GtReadmode readmode,
                 unsigned int prefixlength,
                 unsigned int numofparts,
                 unsigned long maximumspace,
                 GtTimer *sfxprogress,
                 bool withprogressbar,
                 GT_UNUSED GtLogger *logger,
                 GtError *err)
{
  const GtSuffixsortspace *suffixsortspace;
  unsigned long numberofsuffixes;
  bool haserr = false;
  Sfxiterator *sfi;
  bool specialsuffixes = false;

  Sfxstrategy sfxstrategy;

  defaultsfxstrategy(&sfxstrategy,
                     gt_encseq_bitwise_cmp_ok(ssi->encseq) ? false : true);
  sfi = gt_Sfxiterator_new(ssi->encseq,
                           readmode,
                           prefixlength,
                           numofparts,
                           maximumspace,
                           &sfxstrategy,
                           sfxprogress,
                           withprogressbar,
                           NULL, /* logger */
                           err);
  if (sfi == NULL)
  {
    haserr = true;
  } else
  {
    Sequentialsuffixarrayreader *ssar = NULL;
    bool firstpage = true;

    ssar = gt_newSequentialsuffixarrayreaderfromRAM(ssi->encseq,
                                                    readmode);
    while (true)
    {
      suffixsortspace = gt_Sfxiterator_next(&numberofsuffixes,&specialsuffixes,
                                            sfi);
      if (suffixsortspace == NULL || specialsuffixes)
      {
        break;
      }
      gt_updateSequentialsuffixarrayreaderfromRAM(
               ssar,
               (const ESASuffixptr *)
               gt_suffixsortspace_ulong_get(suffixsortspace),
               firstpage,
               numberofsuffixes);
      firstpage = false;
      if (gt_enumeratemaxpairs(ssar,
                               ssi->encseq,
                               readmode,
                               ssi->minlength,
                               ssi->processmaxpairs,
                               ssi->processmaxpairsinfo,
                               err) != 0)
      {
        haserr = true;
      }
    }
    if (ssar != NULL)
    {
      gt_freeSequentialsuffixarrayreader(&ssar);
    }
  }
  if (gt_Sfxiterator_delete(sfi,err) != 0)
  {
    haserr = true;
  }
  return haserr ? -1 : 0;
}

static int sarrselfsubstringmatch(const GtUchar *dbseq,
                                  unsigned long dblen,
                                  const GtUchar *query,
                                  unsigned long querylen,
                                  unsigned int minlength,
                                  GtAlphabet *alpha,
                                  Processmaxpairs processmaxpairs,
                                  void *processmaxpairsinfo,
                                  GtLogger *logger,
                                  GtError *err)
{
  Substringmatchinfo ssi;
  unsigned int numofchars, recommendedprefixlength;
  GtEncseqBuilder *eb;
  bool haserr = false;

  eb = gt_encseq_builder_new(alpha);
  gt_encseq_builder_disable_multiseq_support(eb);
  gt_encseq_builder_disable_description_support(eb);
  gt_encseq_builder_set_logger(eb, logger);
  gt_encseq_builder_add_encoded(eb, dbseq, dblen, NULL);
  gt_encseq_builder_add_encoded(eb, query, querylen, NULL);
  ssi.encseq = gt_encseq_builder_build(eb, err);
  gt_encseq_builder_delete(eb);

  ssi.minlength = minlength;
  ssi.processmaxpairs = processmaxpairs;
  ssi.processmaxpairsinfo = processmaxpairsinfo;
  numofchars = gt_alphabet_num_of_chars(alpha);
  recommendedprefixlength
    = gt_recommendedprefixlength(numofchars,
                                 dblen+querylen+1,
                                 GT_RECOMMENDED_MULTIPLIER_DEFAULT,
                                 true);
  if (constructsarrandrunmaxpairs(&ssi,
                                  GT_READMODE_FORWARD,
                                  recommendedprefixlength,
                                  1U, /* parts */
                                  0, /* maximumspace */
                                  NULL,
                                  false,
                                  logger,
                                  err) != 0)
  {
    haserr = true;
  }
  gt_encseq_delete(ssi.encseq);
  ssi.encseq = NULL;
  return haserr ? -1 : 0;
}

static unsigned long samplesubstring(GtUchar *seqspace,
                              const GtEncseq *encseq,
                              unsigned long substringlength)
{
  unsigned long start, totallength;

  totallength = gt_encseq_total_length(encseq);
  start = (unsigned long) (random() % totallength);
  if (start + substringlength > totallength)
  {
    substringlength = totallength - start;
  }
  gt_assert(substringlength > 0);
  gt_encseq_extract_encoded(encseq,seqspace,start,
                                       start+substringlength-1);
  return substringlength;
}

typedef struct
{
  unsigned long len,
         dbstart,
         querystart;
  uint64_t queryseqnum;
} Substringmatch;

static int storemaxmatchquery(void *info,
                              GT_UNUSED const GtEncseq *encseq,
                              const GtQuerymatch *querymatch,
                              GT_UNUSED const GtUchar *query,
                              GT_UNUSED unsigned long query_totallength,
                              GT_UNUSED GtError *err)
{
  GtArray *tab = (GtArray *) info;
  Substringmatch subm;

  subm.len = gt_querymatch_querylen(querymatch);
  subm.dbstart = gt_querymatch_dbstart(querymatch);
  subm.querystart = gt_querymatch_querystart(querymatch);
  subm.queryseqnum = gt_querymatch_queryseqnum(querymatch);
  gt_array_add(tab,subm);
  return 0;
}

typedef struct
{
  GtArray *results;
  unsigned long dblen, *querymarkpos, querylen;
  unsigned long numofquerysequences;
} Maxmatchselfinfo;

static int storemaxmatchself(void *info,
                             GT_UNUSED const GtEncseq *encseq,
                             unsigned long len,
                             unsigned long pos1,
                             unsigned long pos2,
                             GT_UNUSED GtError *err)
{
  Maxmatchselfinfo *maxmatchselfinfo = (Maxmatchselfinfo *) info;
  unsigned long dbstart, querystart;

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
    unsigned long pos;

    subm.len = len;
    subm.dbstart = dbstart;
    pos = querystart - (maxmatchselfinfo->dblen + 1);
    if (maxmatchselfinfo->querymarkpos == NULL)
    {
      subm.queryseqnum = 0;
      subm.querystart = pos;
    } else
    {
      unsigned long queryseqnum
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

static int orderSubstringmatch(const void *a,const void *b)
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

static int showSubstringmatch(void *a, GT_UNUSED void *info,
                              GT_UNUSED GtError *err)
{
  Substringmatch *m = (Substringmatch *) a;

  printf("%lu %lu " Formatuint64_t " %lu\n",
           m->len,
           m->dbstart,
           PRINTuint64_tcast(m->queryseqnum),
           m->querystart);
  return 0;
}

static unsigned long *sequence2markpositions(unsigned long *numofsequences,
                                      const GtUchar *seq,
                                      unsigned long seqlen)
{
  unsigned long *spacemarkpos, idx;
  unsigned long allocatedmarkpos, nextfreemarkpos;

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
                    unsigned long samples,
                    unsigned int minlength,
                    unsigned long substringlength,
                    GtLogger *logger,
                    GtError *err)
{
  GtEncseq *encseq;
  unsigned long totallength = 0, dblen, querylen;
  GtUchar *dbseq = NULL, *query = NULL;
  bool haserr = false;
  unsigned long s;
  GtArray *tabmaxquerymatches;
  Maxmatchselfinfo maxmatchselfinfo;
  GtEncseqLoader *el;

  gt_logger_log(logger,"draw %lu samples",samples);

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
    dbseq = gt_malloc(sizeof *dbseq * substringlength);
    query = gt_malloc(sizeof *query * substringlength);
  }
  for (s=0; s<samples && !haserr; s++)
  {
    dblen = samplesubstring(dbseq,encseq,substringlength);
    querylen = samplesubstring(query,encseq,substringlength);
    gt_logger_log(logger,"run query match for dblen=%lu"
                         ",querylen= %lu, minlength=%u",
                         dblen,
                         querylen,
                         minlength);
    tabmaxquerymatches = gt_array_new(sizeof (Substringmatch));
    if (gt_sarrquerysubstringmatch(dbseq,
                                   dblen,
                                   query,
                                   (unsigned long) querylen,
                                   minlength,
                                   gt_encseq_alphabet(encseq),
                                   storemaxmatchquery,
                                   tabmaxquerymatches,
                                   logger,
                                   err) != 0)
    {
      haserr = true;
      break;
    }
    gt_logger_log(logger,"run self match for dblen=%lu"
                         ",querylen= %lu, minlength=%u",
                         dblen,
                         querylen,
                         minlength);
    maxmatchselfinfo.results = gt_array_new(sizeof (Substringmatch));
    maxmatchselfinfo.dblen = dblen;
    maxmatchselfinfo.querylen = querylen;
    maxmatchselfinfo.querymarkpos
      = sequence2markpositions(&maxmatchselfinfo.numofquerysequences,
                               query,querylen);
    if (sarrselfsubstringmatch(dbseq,
                               dblen,
                               query,
                               (unsigned long) querylen,
                               minlength,
                               gt_encseq_alphabet(encseq),
                               storemaxmatchself,
                               &maxmatchselfinfo,
                               logger,
                               err) != 0)
    {
      haserr = true;
      break;
    }
    gt_array_sort(tabmaxquerymatches,orderSubstringmatch);
    gt_array_sort(maxmatchselfinfo.results,orderSubstringmatch);
    if (!gt_array_equal(tabmaxquerymatches,maxmatchselfinfo.results,
                        orderSubstringmatch))
    {
      const unsigned long width = 60UL;
      printf("failure for query of length %lu\n",(unsigned long) querylen);
      printf("querymatches\n");
      (void) gt_array_iterate(tabmaxquerymatches,showSubstringmatch,NULL,
                           err);
      printf("dbmatches\n");
      (void) gt_array_iterate(maxmatchselfinfo.results,showSubstringmatch,
                           NULL,err);
      gt_symbolstring2fasta(stdout,"dbseq",
                         gt_encseq_alphabet(encseq),
                         dbseq,
                         (unsigned long) dblen,
                         width);
      gt_symbolstring2fasta(stdout,"queryseq",
                         gt_encseq_alphabet(encseq),
                         query,
                         (unsigned long) querylen,
                         width);
      exit(GT_EXIT_PROGRAMMING_ERROR);
    }
    gt_free(maxmatchselfinfo.querymarkpos);
    printf("# numberofmatches=%lu\n",gt_array_size(tabmaxquerymatches));
    gt_array_delete(tabmaxquerymatches);
    gt_array_delete(maxmatchselfinfo.results);
  }
  gt_free(dbseq);
  gt_free(query);
  gt_encseq_delete(encseq);
  encseq = NULL;
  return haserr ? -1 : 0;
}
