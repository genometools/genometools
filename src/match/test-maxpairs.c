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

#ifdef INLINEDSequentialsuffixarrayreader

#include "core/error.h"
#include "core/str.h"
#include "core/logger.h"
#include "seqpos-def.h"

int testmaxpairs(GT_UNUSED const GtStr *indexname,
                 GT_UNUSED unsigned long samples,
                 GT_UNUSED unsigned int minlength,
                 GT_UNUSED Seqpos substringlength,
                 GT_UNUSED GtLogger *logger,
                 GT_UNUSED GtError *err)
{
  return 0;
}

#else

#include <limits.h>
#include "core/alphabet.h"
#include "core/array.h"
#include "core/arraydef.h"
#include "core/divmodmul.h"
#include "core/unused_api.h"
#include "spacedef.h"
#include "esa-mmsearch.h"
#include "format64.h"
#include "echoseq.h"
#include "encodedsequence.h"
#include "seqpos-def.h"
#include "core/logger.h"
#include "sfx-suffixer.h"
#include "sfx-apfxlen.h"
#include "core/progress_timer.h"
#include "esa-maxpairs.h"
#include "esa-seqread.h"

#include "arrcmp.pr"

typedef struct
{
  unsigned int minlength;
  GtEncodedsequence *encseq;
  Processmaxpairs processmaxpairs;
  void *processmaxpairsinfo;
} Substringmatchinfo;

static int constructsarrandrunmaxpairs(
                 Substringmatchinfo *ssi,
                 GtReadmode readmode,
                 unsigned int prefixlength,
                 unsigned int numofparts,
                 GtProgressTimer *sfxprogress,
                 GtLogger *logger,
                 GtError *err)
{
  const Seqpos *suftabptr;
  Seqpos numberofsuffixes;
  bool haserr = false;
  Sfxiterator *sfi;
  bool specialsuffixes = false;

  sfi = newSfxiterator(ssi->encseq,
                       readmode,
                       prefixlength,
                       numofparts,
                       NULL, /* oulcpinfo */
                       NULL, /* sfxstrategy */
                       sfxprogress,
                       NULL, /* verbosinfo */
                       err);
  if (sfi == NULL)
  {
    haserr = true;
  } else
  {
    Sequentialsuffixarrayreader *ssar = NULL;
    bool firstpage = true;

    ssar = newSequentialsuffixarrayreaderfromRAM(ssi->encseq,
                                                 readmode);
    while (true)
    {
      suftabptr = nextSfxiterator(&numberofsuffixes,&specialsuffixes,sfi);
      if (suftabptr == NULL || specialsuffixes)
      {
        break;
      }
      updateSequentialsuffixarrayreaderfromRAM(ssar,
                                               suftabptr,
                                               firstpage,
                                               numberofsuffixes);
      firstpage = false;
      if (enumeratemaxpairs(ssar,
                            ssi->encseq,
                            readmode,
                            ssi->minlength,
                            ssi->processmaxpairs,
                            ssi->processmaxpairsinfo,
                            logger,
                            err) != 0)
      {
        haserr = true;
      }
    }
    if (ssar != NULL)
    {
      freeSequentialsuffixarrayreader(&ssar);
    }
  }
  if (sfi != NULL)
  {
    freeSfxiterator(&sfi);
  }
  return haserr ? -1 : 0;
}

static int sarrselfsubstringmatch(const GtUchar *dbseq,
                                  Seqpos dblen,
                                  const GtUchar *query,
                                  unsigned long querylen,
                                  unsigned int minlength,
                                  const GtAlphabet *alpha,
                                  Processmaxpairs processmaxpairs,
                                  void *processmaxpairsinfo,
                                  GtLogger *logger,
                                  GtError *err)
{
  Substringmatchinfo ssi;
  unsigned int numofchars;
  bool haserr = false;

  ssi.encseq = plain2encodedsequence(true,
                                     dbseq,
                                     dblen,
                                     query,
                                     querylen,
                                     alpha,
                                     logger);
  ssi.minlength = minlength;
  ssi.processmaxpairs = processmaxpairs;
  ssi.processmaxpairsinfo = processmaxpairsinfo;
  numofchars = gt_alphabet_num_of_chars(alpha);
  if (constructsarrandrunmaxpairs(&ssi,
                                  GT_READMODE_FORWARD,
                                  recommendedprefixlength(numofchars,
                                                          dblen+querylen+1),
                                  1U, /* parts */
                                  NULL,
                                  logger,
                                  err) != 0)
  {
    haserr = true;
  }
  removealpharef(ssi.encseq);
  gt_encodedsequence_delete(ssi.encseq);
  ssi.encseq = NULL;
  return haserr ? -1 : 0;
}

static Seqpos samplesubstring(GtUchar *seqspace,
                              const GtEncodedsequence *encseq,
                              Seqpos substringlength)
{
  Seqpos start, totallength;

  totallength = gt_encodedsequence_total_length(encseq);
  start = (Seqpos) (drand48() * (double) totallength);
  if (start + substringlength > totallength)
  {
    substringlength = totallength - start;
  }
  gt_assert(substringlength > 0);
  gt_encodedsequence_extract_substring(encseq,seqspace,start,
                                       start+substringlength-1);
  return substringlength;
}

typedef struct
{
  Seqpos len,
         dbstart,
         querystart;
  uint64_t queryseqnum;
} Substringmatch;

static int storemaxmatchquery(void *info,
                              GT_UNUSED const GtEncodedsequence *encseq,
                              const Querymatch *querymatch,
                              GT_UNUSED GtError *err)
{
  GtArray *tab = (GtArray *) info;
  Substringmatch subm;

  subm.len = querymatch_len(querymatch);
  subm.dbstart = querymatch_dbstart(querymatch);
  subm.querystart = querymatch_querystart(querymatch);
  subm.queryseqnum = querymatch_queryseqnum(querymatch);
  gt_array_add(tab,subm);
  return 0;
}

typedef struct
{
  GtArray *results;
  Seqpos dblen, *querymarkpos, querylen;
  unsigned long numofquerysequences;
} Maxmatchselfinfo;

static int storemaxmatchself(void *info,
                             GT_UNUSED const GtEncodedsequence *encseq,
                             Seqpos len,
                             Seqpos pos1,
                             Seqpos pos2,
                             GT_UNUSED GtError *err)
{
  Maxmatchselfinfo *maxmatchselfinfo = (Maxmatchselfinfo *) info;
  Seqpos dbstart, querystart;

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
    Seqpos pos;

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
        = getrecordnumSeqpos(maxmatchselfinfo->querymarkpos,
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

  printf(FormatSeqpos " " FormatSeqpos " " Formatuint64_t " " FormatSeqpos "\n",
           PRINTSeqposcast(m->len),
           PRINTSeqposcast(m->dbstart),
           PRINTuint64_tcast(m->queryseqnum),
           PRINTSeqposcast(m->querystart));
  return 0;
}

static Seqpos *sequence2markpositions(unsigned long *numofsequences,
                                      const GtUchar *seq,
                                      Seqpos seqlen)
{
  Seqpos *spacemarkpos, idx;
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
  ALLOCASSIGNSPACE(spacemarkpos,NULL,Seqpos,allocatedmarkpos);
  for (idx=0, nextfreemarkpos = 0; idx<seqlen; idx++)
  {
    if (seq[idx] == (GtUchar) SEPARATOR)
    {
      spacemarkpos[nextfreemarkpos++] = idx;
    }
  }
  return spacemarkpos;
}

int testmaxpairs(const GtStr *indexname,
                 unsigned long samples,
                 unsigned int minlength,
                 Seqpos substringlength,
                 GtLogger *logger,
                 GtError *err)
{
  GtEncodedsequence *encseq;
  Seqpos totallength = 0, dblen, querylen;
  GtUchar *dbseq = NULL, *query = NULL;
  bool haserr = false;
  unsigned long s;
  GtArray *tabmaxquerymatches;
  Maxmatchselfinfo maxmatchselfinfo;

  gt_logger_log(logger,"draw %lu samples",samples);
  encseq = gt_encodedsequence_new_from_index(true,
                              indexname,
                              true,
                              false,
                              false,
                              false,
                              logger,
                              err);
  if (encseq == NULL)
  {
    haserr = true;
  } else
  {
    totallength = gt_encodedsequence_total_length(encseq);
  }
  if (!haserr)
  {
    srand48(42349421);
    if (substringlength > totallength/2)
    {
      substringlength = totallength/2;
    }
    ALLOCASSIGNSPACE(dbseq,NULL,GtUchar,substringlength);
    ALLOCASSIGNSPACE(query,NULL,GtUchar,substringlength);
  }
  for (s=0; s<samples && !haserr; s++)
  {
    dblen = samplesubstring(dbseq,encseq,substringlength);
    querylen = samplesubstring(query,encseq,substringlength);
    gt_logger_log(logger,"run query match for dblen=" FormatSeqpos
                            ",querylen= " FormatSeqpos ", minlength=%u",
           PRINTSeqposcast(dblen),PRINTSeqposcast(querylen),minlength);
    tabmaxquerymatches = gt_array_new(sizeof (Substringmatch));
    if (sarrquerysubstringmatch(dbseq,
                                dblen,
                                query,
                                (unsigned long) querylen,
                                minlength,
                                gt_encodedsequence_alphabet(encseq),
                                storemaxmatchquery,
                                tabmaxquerymatches,
                                logger,
                                err) != 0)
    {
      haserr = true;
      break;
    }
    gt_logger_log(logger,"run self match for dblen=" FormatSeqpos
                            ",querylen= " FormatSeqpos ", minlength=%u",
           PRINTSeqposcast(dblen),PRINTSeqposcast(querylen),minlength);
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
                               gt_encodedsequence_alphabet(encseq),
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
    if (array_compare(tabmaxquerymatches,maxmatchselfinfo.results,
                      orderSubstringmatch) != 0)
    {
      const unsigned long width = 60UL;
      printf("failure for query of length %lu\n",(unsigned long) querylen);
      printf("querymatches\n");
      (void) gt_array_iterate(tabmaxquerymatches,showSubstringmatch,NULL,
                           err);
      printf("dbmatches\n");
      (void) gt_array_iterate(maxmatchselfinfo.results,showSubstringmatch,
                           NULL,err);
      symbolstring2fasta(stdout,"dbseq",
                         gt_encodedsequence_alphabet(encseq),
                         dbseq,
                         (unsigned long) dblen,
                         width);
      symbolstring2fasta(stdout,"queryseq",
                         gt_encodedsequence_alphabet(encseq),
                         query,
                         (unsigned long) querylen,
                         width);
      exit(GT_EXIT_PROGRAMMING_ERROR);
    }
    FREESPACE(maxmatchselfinfo.querymarkpos);
    printf("# numberofmatches=%lu\n",gt_array_size(tabmaxquerymatches));
    gt_array_delete(tabmaxquerymatches);
    gt_array_delete(maxmatchselfinfo.results);
  }
  FREESPACE(dbseq);
  FREESPACE(query);
  gt_encodedsequence_delete(encseq);
  encseq = NULL;
  return haserr ? -1 : 0;
}

#endif /* INLINEDSequentialsuffixarrayreader */
