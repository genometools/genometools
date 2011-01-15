/*
  Copyright (c) 2010 Dirk Willrodt <dwillrodt@zbh.uni-hamburg.de>
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

#include "core/array2dim_api.h"
#include "core/encseq_api.h"
#include "core/format64.h"
#include "core/log_api.h"
#include "core/logger.h"
#include "core/ma.h"
#include "core/mathsupport.h"

#include "match/eis-voiditf.h"
#include "match/esa-seqread.h"
#include "match/esa-shulen.h"
#include "match/idx-limdfs.h"
#include "match/shu-dfs.h"
#include "match/shu-divergence.h"
#include "match/shu-encseq-gc.h"
#include "match/sarr-def.h"

#include "match/shu-genomediff.h"

int gt_genomediff_shu(GtLogger *logger,
                      const GtGenomediffArguments *arguments,
                      GtProgressTimer *timer,
                      GtError *err)
{
  int had_err = 0;
  double **div = NULL,
         *gc_contents = NULL;
  unsigned long numoffiles = 0,
                *filelength = NULL,
                i_idx, j_idx;
  uint64_t **shulen = NULL;
  const GtStrArray *filenames = NULL;
  const GtEncseq *encseq = NULL;
  Genericindex *genericindexSubject = NULL;
  Sequentialsuffixarrayreader *ssar = NULL;

  if (arguments->withesa)
  {
    gt_error_check(err);
    if (timer != NULL)
    {
      gt_progress_timer_start_new_state(timer,
                                        "load sequential sa reader",
                                        stdout);
    }
    ssar =
      gt_newSequentialsuffixarrayreaderfromfile(gt_str_get(
                                                  arguments->indexname),
                                                SARR_LCPTAB |
                                                SARR_SUFTAB |
                                                SARR_ESQTAB |
                                                SARR_SSPTAB,
                                                arguments->scan
                                                ? SEQ_scan : SEQ_mappedboth,
                                                err);
    if (ssar == NULL)
    {
      had_err = -1;
    }
    if (!had_err)
    {
      encseq = gt_encseqSequentialsuffixarrayreader(ssar);
      numoffiles = gt_encseq_num_of_files(encseq);
      gt_array2dim_calloc(shulen, numoffiles, numoffiles);
      filelength = gt_calloc((size_t) numoffiles, sizeof (unsigned long));
      if (timer != NULL)
      {
        gt_progress_timer_start_new_state(timer,
                                          "dfs esa index",
                                          stdout);
      }
      if (gt_get_multiesashulengthdist(ssar,
                                       encseq,
                                       shulen,
                                       logger,
                                       err) != 0)
      {
        had_err = -1;
      }
    }
  } else
  {
    unsigned long numofchars = 0,
                  totallength = 0;
    const GtAlphabet *alphabet;
    const FMindex *subjectindex = NULL;

    if (timer != NULL)
    {
      gt_progress_timer_start_new_state(timer,
                                        "map generic index",
                                        stdout);
    }
    gt_assert(!arguments->withesa);
    genericindexSubject = genericindex_new(gt_str_get(arguments->indexname),
                                           arguments->withesa,
                                           true,
                                           false,
                                           true,
                                           arguments->user_max_depth,
                                           logger,
                                           err);
    if (genericindexSubject == NULL)
    {
      gt_assert(gt_error_is_set(err));
      had_err = -1;
    } else
    {
      encseq = genericindex_getencseq(genericindexSubject);
      numoffiles = gt_encseq_num_of_files(encseq);
      gt_array2dim_calloc(shulen, numoffiles, numoffiles);
      filelength = gt_calloc((size_t) numoffiles, sizeof (unsigned long));
    }

    if (!had_err)
    {
      alphabet = gt_encseq_alphabet(encseq);
      if (!gt_alphabet_is_dna(alphabet) &&
          (!arguments->traverse_only ||
           !arguments->shulen_only))
      {
        gt_error_set(err,"error: sequences need to be dna to calculate gc!");
        had_err = -1;
      } else
      {
        numofchars = (unsigned long) gt_alphabet_num_of_chars(alphabet);
        totallength = genericindex_get_totallength(genericindexSubject);
        subjectindex = genericindex_get_packedindex(genericindexSubject);
      }
    }

    if (!had_err)
    {
      gt_assert(shulen);
      gt_assert(subjectindex);
      had_err = gt_pck_calculate_shulen(subjectindex,
                                        encseq,
                                        shulen,
                                        numofchars,
                                        totallength,
                                        !arguments->traverse_only,
                                        timer,
                                        logger,
                                        err);
    }
    if (!had_err && arguments->traverse_only)
    {
      gt_logger_log(logger, "stopping after traversal");
      genericindex_delete(genericindexSubject);
      gt_free(filelength);
      if (shulen != NULL)
        gt_array2dim_delete(shulen);
      gt_assert(timer != NULL);
      return had_err;
    }
  }

  if (!had_err)
  {
    gt_assert(filelength);
    filenames = gt_encseq_filenames(encseq);
    for (i_idx = 0UL; i_idx < numoffiles && !had_err; i_idx++)
    {
      filelength[i_idx] =
        (unsigned long) gt_encseq_effective_filelength(encseq, i_idx) - 1;
    }
  }
  if (!had_err)
  {
    gt_assert(shulen);
    if (timer != NULL)
    {
      gt_progress_timer_start_new_state(timer,
                                        "calculate avg",
                                        stdout);
    }
    if (arguments->shulen_only)
    {
      printf("# sum of shulen\n%lu\n", numoffiles);
    }
    gt_array2dim_calloc(div, numoffiles, numoffiles);
    gt_assert(filelength);
    for (i_idx = 0; i_idx < numoffiles; i_idx++)
    {
      unsigned long length_i;
      length_i = filelength[i_idx];
      if (arguments->shulen_only)
        printf("%s\t", gt_str_array_get(filenames, i_idx));
      for (j_idx = 0; j_idx < numoffiles; j_idx++)
      {
        if (j_idx == i_idx)
        {
          if (arguments->shulen_only)
            printf("0\t");
        } else
        {
          if (arguments->shulen_only)
            printf(Formatuint64_t"\t", arguments->withesa ?
                              PRINTuint64_tcast(shulen[j_idx][i_idx]) :
                              PRINTuint64_tcast(shulen[i_idx][j_idx]));
          else
            div[i_idx][j_idx] = arguments->withesa ?
                                  ((double) shulen[j_idx][i_idx]) / length_i :
                                  ((double) shulen[i_idx][j_idx]) / length_i;
        }
      }
      if (arguments->shulen_only)
        printf("\n");
    }
  }
  if (!arguments->shulen_only && !had_err)
  {
    gt_assert(div);
    if (timer != NULL)
    {
      gt_progress_timer_start_new_state(timer,
                                        "calculate gc",
                                        stdout);
    }
    gc_contents = gt_encseq_get_gc(encseq,
                                  true,
                                  false,
                                  err);
    if (gc_contents == NULL)
      had_err = -1;

    gt_logger_log(logger, "table of avg shulens");
    if (!had_err && gt_logger_enabled(logger))
    {
      for (i_idx = 0; i_idx < numoffiles; i_idx++)
      {
        printf("# %s\t", gt_str_array_get(filenames, i_idx));
        for (j_idx = 0; j_idx < numoffiles; j_idx++)
        {
          if (i_idx == j_idx)
            printf("0\t\t");
          else
            printf("%f\t", div[i_idx][j_idx]);
        }
        printf("\n");
      }
    }
    if (!had_err)
    {
      double *ln_n_fac;

      gt_assert(div);
      gt_assert(filelength);
      gt_assert(gc_contents);

      if (timer != NULL)
      {
        gt_progress_timer_start_new_state(timer,
                                          "precalculate ln_n_fac",
                                          stdout);
      }
      ln_n_fac = gt_get_ln_n_fac(arguments->max_ln_n_fac);
      if (timer != NULL)
      {
        gt_progress_timer_start_new_state(timer,
                                          "calculate divergence",
                                          stdout);
      }
      for (i_idx = 0; i_idx < numoffiles; i_idx++)
      {
        for (j_idx = i_idx+1; j_idx < numoffiles; j_idx++)
        {
          double query_gc, query_shulen;
          unsigned long subject_len;
          if (gt_double_smaller_double(div[i_idx][j_idx],
                                       div[j_idx][i_idx]))
          { /* S=j_idx Q=i_idx */
            query_gc = gc_contents[i_idx];
            query_shulen = div[i_idx][j_idx];
            subject_len = filelength[j_idx];
          } else
          {
            if (gt_double_smaller_double(div[j_idx][i_idx],
                                         div[i_idx][j_idx]))
            { /* S=i_idx Q=j_idx */
              query_gc = gc_contents[j_idx];
              query_shulen = div[j_idx][i_idx];
              subject_len = filelength[i_idx];
            } else
            {
              if (gt_double_smaller_double(fabs(gc_contents[i_idx]-0.5),
                                           fabs(gc_contents[j_idx]-0.5)))
              { /* S=i_idx Q=j_idx XXX check this if right*/
                query_gc = gc_contents[j_idx];
                query_shulen = div[j_idx][i_idx];
                subject_len = filelength[i_idx];
              } else
                query_gc = gc_contents[i_idx];
                query_shulen = div[i_idx][j_idx];
                subject_len = filelength[j_idx];
              { /* S=j_idx Q=i_idx */
              }
            }
          }

          div[i_idx][j_idx] =
            gt_divergence(arguments->divergence_rel_err,
                          arguments->divergence_abs_err,
                          arguments->divergence_m,
                          arguments->divergence_threshold,
                          query_shulen,
                          subject_len,
                          query_gc,
                          ln_n_fac,
                          arguments->max_ln_n_fac);
          div[j_idx][i_idx] = div[i_idx][j_idx];
        }
      }
      gt_free(ln_n_fac);
    }
    gt_logger_log(logger, "table of divergences");
    if (!had_err && gt_logger_enabled(logger))
    {
      gt_assert(div);
      for (i_idx = 0; i_idx < numoffiles; i_idx++)
      {
        printf("# %s\t", gt_str_array_get(filenames, i_idx));
        for (j_idx = 0; j_idx < numoffiles; j_idx++)
        {
          if (i_idx == j_idx)
          {
            printf("0\t\t");
            continue;
          }
          printf("%f\t", div[i_idx][j_idx]);
        }
        printf("\n");
      }
    }
    if (!had_err)
    {
      gt_assert(div);
      if (timer != NULL)
      {
        gt_progress_timer_start_new_state(timer,
                                          "calculate kr",
                                          stdout);
      }
      printf("# Table of Kr\n%lu\n", numoffiles);
      for (i_idx = 0; i_idx < numoffiles; i_idx++)
      {
        printf("%s\t", gt_str_array_get(filenames, i_idx));
        for (j_idx = 0; j_idx < numoffiles; j_idx++)
        {
          if ( i_idx == j_idx )
            printf("0\t\t");
          else
            printf("%f\t", gt_calculateKr(div[i_idx][j_idx]));
        }
        printf("\n");
      }
    }
  }
  if (arguments->withesa)
  {
    if (ssar != NULL)
    {
      gt_freeSequentialsuffixarrayreader(&ssar);
    }
  } else
  {
    genericindex_delete(genericindexSubject);
  }

  gt_free(filelength);
  gt_free(gc_contents);
  if (shulen != NULL)
    gt_array2dim_delete(shulen);
  if (div != NULL)
    gt_array2dim_delete(div);
  return had_err;
}
