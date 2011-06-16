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
#include "match/sarr-def.h"
#include "match/shu-dfs.h"
#include "match/shu-divergence.h"
#include "match/shu-encseq-gc.h"
#include "match/shu_unitfile.h"

#include "match/shu-genomediff.h"

static inline int parse_unit(const GtEncseq *encseq,
                             struct GtShuUnitFileInfo_tag *unit_info,
                             const GtGenomediffArguments *arguments,
                             GtTimer *timer,
                             GtLogger *logger,
                             GtError *err)
{
  int had_err = 0;
  unsigned long i_idx;
  unit_info->num_of_files = gt_encseq_num_of_files(encseq);
  unit_info->file_names = gt_encseq_filenames(encseq);
  if (arguments->with_units)
  {
    if (timer != NULL)
    {
      gt_timer_show_progress(timer, "parse unitfile", stdout);
    }
    had_err = gt_read_genomediff_unitfile(arguments->unitfile,
                                          unit_info,
                                          logger,
                                          err);
    if (!had_err)
    {
      gt_logger_log(logger, "successfully loaded unitfile");
    }
  }
  else
  {
    unit_info->num_of_genomes = unit_info->num_of_files;
    unit_info->genome_names = gt_str_array_new();
    for (i_idx = 0; i_idx < unit_info->num_of_files; i_idx++)
    {
      gt_str_array_add_cstr(unit_info->genome_names,
                            gt_str_array_get(unit_info->file_names, i_idx));
    }
  }
  return had_err;
}

int gt_genomediff_shu(GtLogger *logger,
                      const GtGenomediffArguments *arguments,
                      GtTimer *timer,
                      GtError *err)
{
  int had_err = 0;
  double **div = NULL,
         *seq_gc_content = NULL,
         *gc_content = NULL;
  unsigned long *genome_length = NULL,
                i_idx, j_idx;
  uint64_t **shulen = NULL;
  const GtEncseq *encseq = NULL;
  Genericindex *genericindexSubject = NULL;
  Sequentialsuffixarrayreader *ssar = NULL;
  struct GtShuUnitFileInfo_tag *unit_info;

  unit_info = gt_malloc(sizeof(*unit_info));
  unit_info->map_files = NULL;
  unit_info->genome_names = NULL;
  unit_info->num_of_genomes = 0;
  unit_info->num_of_files = 0;
  unit_info->file_names = NULL;

  /*XXX change esa-functions so that unitfile can be used*/
  if (arguments->with_esa)
  {
    gt_error_check(err);
    if (timer != NULL)
    {
      gt_timer_show_progress(timer, "load sequential sa reader", stdout);
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
      had_err = parse_unit(encseq,
                           unit_info,
                           arguments,
                           timer,
                           logger,
                           err);
      gt_array2dim_calloc(shulen,
                          unit_info->num_of_genomes,
                          unit_info->num_of_genomes);
      genome_length = gt_calloc((size_t) unit_info->num_of_genomes,
                                sizeof (unsigned long));
    }
    if (!had_err)
    {
      if (timer != NULL)
      {
        gt_timer_show_progress(timer, "dfs esa index", stdout);
      }
      if (gt_get_multiesashulengthdist(ssar,
                                       encseq,
                                       shulen,
                                       unit_info,
                                       err) != 0)
      {
        had_err = -1;
      }
    }
  }
  else
  {
    const GtAlphabet *alphabet;

    if (timer != NULL)
    {
      gt_timer_show_progress(timer, "map generic index", stdout);
    }
    genericindexSubject = genericindex_new(gt_str_get(arguments->indexname),
                                           arguments->with_esa,
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
    }
    /*get num of genomes from unitfile or encseq and allocate shulen-matrix*/
    if (!had_err)
    {
      encseq = genericindex_getencseq(genericindexSubject);
      had_err = parse_unit(encseq,
                           unit_info,
                           arguments,
                           timer,
                           logger,
                           err);
      gt_array2dim_calloc(shulen,
                          unit_info->num_of_genomes,
                          unit_info->num_of_genomes);
      genome_length = gt_calloc((size_t) unit_info->num_of_genomes,
                                sizeof (unsigned long));
    }
    /*check if alphabet is DNA, which is presumed for gc-calculation*/
    if (!had_err)
    {
      alphabet = gt_encseq_alphabet(encseq);
      if (!gt_alphabet_is_dna(alphabet) &&
          (!arguments->traverse_only ||
           !arguments->shulen_only))
      {
        gt_error_set(err,"error: sequences need to be dna to calculate gc!");
        had_err = -1;
      }
    }
    /*calculate shulens*/
    if (!had_err)
    {
      const FMindex *subjectindex = NULL;
      unsigned long num_of_chars,
                    total_length;
      num_of_chars = (unsigned long) gt_alphabet_num_of_chars(alphabet);
      total_length = genericindex_get_totallength(genericindexSubject);
      subjectindex = genericindex_get_packedindex(genericindexSubject);
      gt_assert(shulen);
      gt_assert(subjectindex);
      had_err = gt_pck_calculate_shulen(subjectindex,
                                        encseq,
                                        unit_info,
                                        shulen,
                                        num_of_chars,
                                        total_length,
                                        !arguments->traverse_only,
                                        timer,
                                        logger,
                                        err);
    }
    if (!had_err && arguments->traverse_only)
    {
      gt_assert(timer != NULL);
      gt_logger_log(logger, "stopping after traversal");
      genericindex_delete(genericindexSubject);
      gt_free(genome_length);
      gt_delete_unit_file_info(unit_info);
      if (shulen != NULL)
      {
        gt_array2dim_delete(shulen);
      }
      return had_err;
    }
  }
   /*get file respective genome length*/
  if (!had_err)
  {
    bool mirrored = gt_encseq_is_mirrored(encseq);
    unsigned long sep_per_file = 0,
                  seqs_passed = 0,
                  eff_file_length,
                  seqs_upto_fileend,
                  filestart;
    gt_assert(genome_length);
    for (i_idx = 0UL; i_idx < unit_info->num_of_files && !had_err; i_idx++)
    {
      eff_file_length = gt_encseq_effective_filelength(encseq, i_idx);
      filestart = gt_encseq_filestartpos(encseq,i_idx);
      if (i_idx == unit_info->num_of_files - 1) {
        if (mirrored) {
          sep_per_file = GT_DIV2(gt_encseq_num_of_sequences(encseq))
                         - seqs_passed - 1;
        } else {
          sep_per_file = gt_encseq_num_of_sequences(encseq) - seqs_passed - 1;
        }
      } else {
        seqs_upto_fileend = gt_encseq_seqnum(encseq,
            filestart +
            eff_file_length - 1 ) + 1;
        sep_per_file = seqs_upto_fileend - seqs_passed - 1;
        seqs_passed = seqs_upto_fileend;
      }
      if (unit_info->map_files != NULL)
      {
        genome_length[unit_info->map_files[i_idx]] +=
          eff_file_length - sep_per_file;
      }
      else
      {
        genome_length[i_idx] = eff_file_length - sep_per_file;
      }
    }
    for (i_idx = 0UL; i_idx < unit_info->num_of_files && mirrored; i_idx++) {
      genome_length[i_idx] += genome_length[i_idx];
    }
    if (gt_log_enabled()) {
      for  (i_idx = 0UL; i_idx < unit_info->num_of_files; i_idx++) {
        gt_log_log("file/genome %lu has length %lu",
            i_idx, genome_length[i_idx]);
      }
    }
    for (i_idx = 0UL; i_idx < unit_info->num_of_files && mirrored; i_idx++) {
      genome_length[i_idx] += genome_length[i_idx];
    }
    if (gt_log_enabled()) {
      for  (i_idx = 0UL; i_idx < unit_info->num_of_files; i_idx++) {
        gt_log_log("file/genome %lu has length %lu",
            i_idx, genome_length[i_idx]);
      }
    }
  }
  /*calculate avg shulen or print sum shulen*/
  if (!had_err)
  {
    gt_assert(shulen);
    if (timer != NULL)
    {
      if (arguments->shulen_only)
      {
        gt_timer_show_progress(timer, "print sums", stdout);
      }
      else
      {
        gt_timer_show_progress(timer, "calculate avg", stdout);
      }
    }
    if (arguments->shulen_only)
    {
      printf("# sum of shulen\n%lu\n", unit_info->num_of_genomes);
    }
    gt_array2dim_calloc(div,
                        unit_info->num_of_genomes,
                        unit_info->num_of_genomes);
    gt_assert(genome_length);
    for (i_idx = 0; i_idx < unit_info->num_of_genomes; i_idx++)
    {
      unsigned long length_i;
      length_i = genome_length[i_idx];
      if (arguments->shulen_only)
      {
        printf("%s\t", gt_str_array_get(unit_info->genome_names, i_idx));
      }
      for (j_idx = 0; j_idx < unit_info->num_of_genomes; j_idx++)
      {
        if (j_idx == i_idx)
        {
          if (arguments->shulen_only)
          {
            printf("0\t");
          }
        }
        else
        {
          if (arguments->shulen_only)
          {
            printf(Formatuint64_t"\t", arguments->with_esa ?
                              PRINTuint64_tcast(shulen[j_idx][i_idx]) :
                              PRINTuint64_tcast(shulen[i_idx][j_idx]));
          }
          else
          {
            div[i_idx][j_idx] = arguments->with_esa ?
                                  ((double) shulen[j_idx][i_idx]) / length_i :
                                  ((double) shulen[i_idx][j_idx]) / length_i;
          }
        }
      }
      if (arguments->shulen_only)
      {
        printf("\n");
      }
    }
  }
  if (!arguments->shulen_only && !had_err)
  {
    gt_assert(div);
    if (timer != NULL)
    {
      gt_timer_show_progress(timer, "calculate gc", stdout);
    }
    /*count gc, do not calculate ratio*/
    seq_gc_content = gt_encseq_get_gc(encseq,
                                  false,
                                  false,
                                  err);

    if (seq_gc_content != NULL)
    {
      /*uint64_t file_sep;*/
      unsigned long file_idx = 0UL,
                    seq_idx,
                    genome_idx,
                    num_of_seq;
      gc_content = gt_calloc((size_t)unit_info->num_of_genomes,
                             sizeof (*gc_content));
      num_of_seq = gt_encseq_num_of_sequences(encseq);
      /*file_sep = gt_encseq_effective_filelength(encseq, file_idx);*/
      for (seq_idx = 0; seq_idx < num_of_seq; seq_idx++) {
        file_idx = gt_encseq_filenum(encseq,
                                     gt_encseq_seqstartpos(encseq,
                                                           seq_idx));
        if (unit_info->map_files != NULL)
        {
          gc_content[unit_info->map_files[file_idx]] +=
                                                  seq_gc_content[seq_idx];
        }
        else
        {
          gc_content[file_idx] += seq_gc_content[seq_idx];
        }
      }
      for (genome_idx = 0; genome_idx < unit_info->num_of_genomes; genome_idx++)
      {
        gt_assert(genome_length != NULL);
        gt_log_log("gc content of genome %lu: %f",
                   genome_idx, gc_content[genome_idx]);
        gc_content[genome_idx] /= (double) genome_length[genome_idx];
      }
      gt_free(seq_gc_content);
    }
    else
    {
      had_err = -1;
    }
    if (!had_err && gt_logger_enabled(logger))
    {
      gt_logger_log(logger, "table of avg shulens");
      for (i_idx = 0; i_idx < unit_info->num_of_genomes; i_idx++)
      {
        printf("# %s\t", gt_str_array_get(unit_info->genome_names, i_idx));
        for (j_idx = 0; j_idx < unit_info->num_of_genomes; j_idx++)
        {
          if (i_idx == j_idx)
          {
            printf("0\t\t");
          }
          else
          {
            printf("%f\t", div[i_idx][j_idx]);
          }
        }
        printf("\n");
      }
    }
    /* calculation of divergence and ln_n_fac */
    if (!had_err)
    {
      double *ln_n_fac;

      gt_assert(div);
      gt_assert(genome_length);
      gt_assert(gc_content);

      if (timer != NULL)
      {
        gt_timer_show_progress(timer, "pre calculate ln_n_fac", stdout);
      }
      ln_n_fac = gt_get_ln_n_fac(arguments->max_ln_n_fac);
      if (timer != NULL)
      {
        gt_timer_show_progress(timer, "calculate divergence", stdout);
      }
      for (i_idx = 0; i_idx < unit_info->num_of_genomes; i_idx++)
      {
        for (j_idx = i_idx+1; j_idx < unit_info->num_of_genomes; j_idx++)
        {
          double query_gc, query_shulen;
          unsigned long subject_len;
          if (gt_double_smaller_double(div[i_idx][j_idx],
                                       div[j_idx][i_idx]))
          { /* S=j_idx Q=i_idx */
            query_gc = gc_content[i_idx];
            query_shulen = div[i_idx][j_idx];
            subject_len = genome_length[j_idx];
          }
          else
          {
            if (gt_double_smaller_double(div[j_idx][i_idx],
                                         div[i_idx][j_idx]))
            { /* S=i_idx Q=j_idx */
              query_gc = gc_content[j_idx];
              query_shulen = div[j_idx][i_idx];
              subject_len = genome_length[i_idx];
            }
            else
            {
              if (gt_double_smaller_double(fabs(gc_content[i_idx]-0.5),
                                           fabs(gc_content[j_idx]-0.5)))
              { /* S=i_idx Q=j_idx XXX check this if right*/
                query_gc = gc_content[j_idx];
                query_shulen = div[j_idx][i_idx];
                subject_len = genome_length[i_idx];
              }
              else
              { /* S=j_idx Q=i_idx */
                query_gc = gc_content[i_idx];
                query_shulen = div[i_idx][j_idx];
                subject_len = genome_length[j_idx];
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
    if (!had_err && gt_logger_enabled(logger))
    {
      gt_logger_log(logger, "table of divergences");
      gt_assert(div);
      for (i_idx = 0; i_idx < unit_info->num_of_genomes; i_idx++)
      {
        printf("# %s\t", gt_str_array_get(unit_info->genome_names, i_idx));
        for (j_idx = 0; j_idx < unit_info->num_of_genomes; j_idx++)
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
        gt_timer_show_progress(timer, "calculate kr", stdout);
      }
      printf("# Table of Kr\n%lu\n", unit_info->num_of_genomes);
      for (i_idx = 0; i_idx < unit_info->num_of_genomes; i_idx++)
      {
        printf("%s\t", gt_str_array_get(unit_info->genome_names, i_idx));
        for (j_idx = 0; j_idx < unit_info->num_of_genomes; j_idx++)
        {
          if ( i_idx == j_idx )
          {
            printf("0\t\t");
          }
          else
          {
            printf("%f\t", gt_calculateKr(div[i_idx][j_idx]));
          }
        }
        printf("\n");
      }
    }
  }
  if (arguments->with_esa)
  {
    if (ssar != NULL)
    {
      gt_freeSequentialsuffixarrayreader(&ssar);
    }
  }
  else
  {
    genericindex_delete(genericindexSubject);
  }

  gt_free(genome_length);
  gt_free(gc_content);
  gt_delete_unit_file_info(unit_info);
  if (shulen != NULL)
  {
    gt_array2dim_delete(shulen);
  }
  if (div != NULL)
  {
    gt_array2dim_delete(div);
  }
  return had_err;
}
