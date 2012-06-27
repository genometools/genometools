/*
  Copyright (c) 2010 Dirk Willrodt <willrodt@zbh.uni-hamburg.de>
  Copyright (c) 2010 Center for Bioinformatics, University of Hamburg

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

#include <stdio.h>
#include <float.h>

#include "core/array2dim_api.h"
#include "core/encseq_api.h"
#include "core/log_api.h"
#include "core/logger.h"
#include "core/ma.h"
#include "core/seqiterator.h"
#include "core/seqiterator_sequence_buffer.h"
#include "match/eis-voiditf.h"
#include "match/idx-limdfs.h"
#include "match/shu-dfs.h"
#include "match/shu-divergence.h"
#include "match/shu-genomediff-pck-simple.h"

int gt_genomediff_pck_shu_simple(GtLogger *logger,
                                 const GtGenomediffArguments *arguments,
                                 GtError *err)
{
  int had_err = 0;
  int retval;
  GtSeqIterator *queries = NULL;
  const GtUchar *symbolmap, *currentQuery;
  const GtAlphabet *alphabet;
  GtUchar c_sym = 0,
          g_sym = 0;
  uint64_t queryNo;
  char *description = NULL;
  unsigned long queryLength,
                subjectLength = 0,
                currentSuffix;
  double avgShuLength,
         currentShuLength = 0.0,
         /*gc_subject,*/
         gc_query /*, gc*/;
  const FMindex *subjectindex = NULL;
  Genericindex *genericindexSubject;
  const GtEncseq *encseq = NULL;
  double *ln_n_fac;

  /* get the precalculation of ln(n!) for 0<n<max_ln_n_fac */
  ln_n_fac = gt_get_ln_n_fac(arguments->max_ln_n_fac);
  gt_log_log("ln(max_ln_n_fac!) = %f\n",
             ln_n_fac[arguments->max_ln_n_fac]);

  genericindexSubject = genericindex_new(gt_str_get(
                                           arguments->indexname),
                                         arguments->with_esa,
                                         true,
                                         false,
                                         true,
                                         arguments->user_max_depth,
                                         logger,
                                         err);
  if (genericindexSubject == NULL)
  {
    had_err = 1;
  }
  else
  {
    encseq = genericindex_getencseq(genericindexSubject);
  }

  if (!had_err)
  {
    subjectLength = genericindex_get_totallength(genericindexSubject) - 1;
    /*subjectLength /= 2;*/
    /*gt_log_log("subject length: %lu", subjectLength);*/
    subjectindex = genericindex_get_packedindex(genericindexSubject);

    queries = gt_seqiterator_sequence_buffer_new(
                                          arguments->queryname,
                                          err);
    gt_assert(queries);
    alphabet = gt_encseq_alphabet(encseq);
    /* makes assumption that alphabet is dna, it has to calculate the gc! */
    if (!gt_alphabet_is_dna(alphabet))
    {
      fprintf(stderr, "error: Sequences need to be dna");
      had_err = 1;
    }
    else
    {
      symbolmap = gt_alphabet_symbolmap(alphabet);
      gt_seqiterator_set_symbolmap(queries, symbolmap);
      c_sym = gt_alphabet_encode(alphabet, 'c');
      g_sym = gt_alphabet_encode(alphabet, 'g');
    }
  }

  for (queryNo = 0; !had_err; queryNo++)
  {
    retval = gt_seqiterator_next(queries,
                                 &currentQuery,
                                 &queryLength,
                                 &description,
                                 err);
    if ( retval != 1)
    {
      if (retval < 0)
      {
        gt_free(description);
      }
      break;
    }
    gt_logger_log(logger,
                  "found query of length: %lu",
                  queryLength);
    avgShuLength = 0.0;
    gc_query = 0.0;
    for (currentSuffix = 0; currentSuffix < queryLength; currentSuffix++)
    {
      currentShuLength = (double) gt_pck_getShuStringLength(
                    subjectindex,
                    &currentQuery[currentSuffix],
                    queryLength - currentSuffix);
      avgShuLength += currentShuLength;
      if (currentQuery[currentSuffix] == c_sym ||
          currentQuery[currentSuffix] == g_sym)
      {
        gc_query++;
      }
    }
    if (arguments->shulen_only)
    {
      printf("# Query %d sum of shulen:\n %.0f\n",
             (int) queryNo, avgShuLength);
    }
    else
    {
      avgShuLength /= (double) queryLength;
      gc_query /= (double) queryLength;

      gt_logger_log(logger, "Query %d has an average SHUstring length "
                            "of\n# shulength: %f",
                            (int) queryNo, avgShuLength);
      gt_logger_log(logger, "Query description: %s", description);
      gt_log_log("Query (i): %s", description);

  /* XXX Fehlerabfragen einbauen */

      if ( !had_err )
      {
        double div, kr;

        gt_logger_log(logger, "shulen:\n%f", avgShuLength);
        gt_log_log("shu: %f, gc: %f, len: %lu",
            avgShuLength, gc_query, subjectLength);
        div =  gt_divergence(arguments->divergence_rel_err,
                             arguments->divergence_abs_err,
                             arguments->divergence_m,
                             arguments->divergence_threshold,
                             avgShuLength,
                             subjectLength,
                             gc_query,
                             ln_n_fac,
                             arguments->max_ln_n_fac);
        gt_logger_log(logger, "divergence:\n%f", div);

        kr = gt_calculateKr(div);

        printf("# Kr:\n%f\n", kr);
      }
    }
  }
  gt_free(ln_n_fac);
  gt_seqiterator_delete(queries);
  genericindex_delete(genericindexSubject);
  return had_err;
}
