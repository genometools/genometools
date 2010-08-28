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
#include "core/log_api.h"
#include "core/logger.h"
#include "core/ma.h"
#include "core/mathsupport.h"

#include "match/eis-voiditf.h"
#include "match/idx-limdfs.h"
#include "match/shu-dfs.h"
#include "match/shu-divergence.h"
#include "match/shu_encseq_gc.h"

#include "match/shu-genomediff-kr2.h"

int gt_genomediff_run_kr2_search(Genericindex *genericindexSubject,
                                 const GtEncseq *encseq,
                                 GtLogger *logger,
                                 const GtGenomediffArguments *arguments,
                                 GtProgressTimer *timer,
                                 GtError *err)
{
  int had_err = 0;
  unsigned long numofchars,
                numoffiles,
                totallength,
                start = 0UL,
                end = 0UL,
                i_idx, j_idx,
                *filelength;
  double **shulen,
         *gc_contents = NULL;
  const GtAlphabet *alphabet;
  const GtStrArray *filenames;
  const FMindex *subjectindex;

  alphabet = gt_encseq_alphabet(encseq);
  numofchars = (unsigned long) gt_alphabet_num_of_chars(alphabet);
  totallength = genericindex_get_totallength(genericindexSubject);
  gt_logger_log(logger, "totallength=%lu", totallength);
  filenames = gt_encseq_filenames(encseq);
  subjectindex = genericindex_get_packedindex(genericindexSubject);
  numoffiles = gt_encseq_num_of_files(encseq);
  gt_logger_log(logger, "number of files=%lu", numoffiles);
  gt_array2dim_calloc(shulen, numoffiles, numoffiles);
  filelength = gt_calloc((size_t) numoffiles, sizeof (unsigned long));

  for (i_idx = 0UL; i_idx < numoffiles; i_idx++)
  {
    start = gt_encseq_filestartpos(encseq, i_idx);
    filelength[i_idx] =
      (unsigned long) gt_encseq_effective_filelength(encseq, i_idx) - 1;
    end = start + filelength[i_idx];
    gt_logger_log(logger,
           "File: %s (No: %lu)\tstart: %lu, end: %lu, sep: %lu",
           gt_str_array_get(filenames, i_idx),
           i_idx,
           start,
           end,
           end + 1);
  }
  had_err = gt_pck_calculate_shulen(subjectindex,
                                    encseq,
                                    shulen,
                                    (unsigned long) numofchars,
                                    totallength,
                                    !arguments->traverse_only,
                                    timer,
                                    logger,
                                    err);
  if (!had_err)
  {
    /*due to the nature of the fm-index calculated shulens are always too long
     * by 2, this adjusts for this */
    for (i_idx = 0; i_idx < numoffiles; i_idx++)
      for (j_idx = i_idx + 1; j_idx < numoffiles; j_idx++)
      {
        if (gt_double_smaller_double(2.0, shulen[i_idx][j_idx]))
          shulen[i_idx][j_idx] -= 2.0;
        if (gt_double_smaller_double(2.0, shulen[j_idx][i_idx]))
          shulen[j_idx][i_idx] -= 2.0;
      }
    if (arguments->traverse_only)
    {
      gt_free(filelength);
      gt_array2dim_delete(shulen);
      gt_assert(timer != NULL);
      return had_err;
    }
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
  }
  if (!had_err)
  {
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
          continue;
        }
        if (arguments->shulen_only)
          printf("%.0f\t", shulen[i_idx][j_idx]);
        shulen[i_idx][j_idx] = shulen[i_idx][j_idx] / length_i;
      }
      if (arguments->shulen_only)
        printf("\n");
    }
  }
  if (!arguments->shulen_only)
  {
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
            printf("%f\t", shulen[i_idx][j_idx]);
        }
        printf("\n");
      }
    }
    if (!had_err)
    {
      double *ln_n_fac;

      gt_assert(gc_contents != NULL);

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
          if (gt_double_smaller_double(shulen[i_idx][j_idx],
                                       shulen[j_idx][i_idx]))
          { /* S=j_idx Q=i_idx */
            query_gc = gc_contents[i_idx];
            query_shulen = shulen[i_idx][j_idx];
            subject_len = filelength[j_idx];
          } else
          {
            if (gt_double_smaller_double(shulen[j_idx][i_idx],
                                         shulen[i_idx][j_idx]))
            { /* S=i_idx Q=j_idx */
              query_gc = gc_contents[j_idx];
              query_shulen = shulen[j_idx][i_idx];
              subject_len = filelength[i_idx];
            } else
            {
              if (gt_double_smaller_double(fabs(gc_contents[i_idx]-0.5),
                                           fabs(gc_contents[j_idx]-0.5)))
              { /* S=i_idx Q=j_idx XXX check this if right*/
                query_gc = gc_contents[j_idx];
                query_shulen = shulen[j_idx][i_idx];
                subject_len = filelength[i_idx];
              } else
                query_gc = gc_contents[i_idx];
                query_shulen = shulen[i_idx][j_idx];
                subject_len = filelength[j_idx];
              { /* S=j_idx Q=i_idx */
              }
            }
          }

          shulen[i_idx][j_idx] =
            gt_divergence(arguments->divergence_rel_err,
                          arguments->divergence_abs_err,
                          arguments->divergence_m,
                          arguments->divergence_threshold,
                          query_shulen,
                          subject_len,
                          query_gc,
                          ln_n_fac,
                          arguments->max_ln_n_fac);
          shulen[j_idx][i_idx] = shulen[i_idx][j_idx];
        }
      }
      gt_free(ln_n_fac);
    }
    gt_logger_log(logger, "table of divergences");
    if (!had_err && gt_logger_enabled(logger))
    {
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
          printf("%f\t", shulen[i_idx][j_idx]);
        }
        printf("\n");
      }
    }
    if (!had_err)
    {
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
            printf("%f\t", gt_calculateKr(shulen[i_idx][j_idx]));
        }
        printf("\n");
      }
    }
  }
  gt_free(filelength);
  gt_free(gc_contents);
  gt_array2dim_delete(shulen);
  return had_err;
}
