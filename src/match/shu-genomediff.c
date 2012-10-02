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
#include "core/safearith.h"
#include "match/eis-voiditf.h"
#include "match/esa-seqread.h"
#include "match/esa-shulen.h"
#include "match/idx-limdfs.h"
#include "match/sarr-def.h"
#include "match/shu-dfs.h"
#include "match/shu-divergence.h"
#include "match/shu-encseq-gc.h"
#include "match/shu-genomediff.h"
#include "match/shu_unitfile.h"

/* static inline int genomediff_parse_unit(GtShuUnitFileInfo *unit_info,
                                        const GtStr *filename,
                                        GtTimer *timer, GtLogger *logger,
                                        GtError *err)
{
  int had_err = 0;

  if (timer != NULL)
    gt_timer_show_progress(timer, "parse unitfile", stdout);
  had_err = gt_shu_unit_file_info_read(filename, unit_info, logger, err);
  if (!had_err)
    gt_logger_log(logger, "successfully loaded unitfile");
  return had_err;
} */

static unsigned long*
genomediff_calculate_genome_lengths(GtShuUnitFileInfo *unit_info)
{
  bool mirrored = gt_encseq_is_mirrored(unit_info->encseq);
  unsigned long eff_file_length,
                filestart,
                sep_per_file = 0,
                seqs_passed = 0,
                seqs_upto_fileend,
                i_idx;
  unsigned long *genome_lengths = gt_calloc((size_t) unit_info->num_of_genomes,
                                           sizeof (*genome_lengths));
  for (i_idx = 0UL; i_idx < unit_info->num_of_files; i_idx++) {
    eff_file_length = gt_safe_cast2ulong_64(
        gt_encseq_effective_filelength(unit_info->encseq, i_idx));
    filestart = gt_encseq_filestartpos(unit_info->encseq, i_idx);
    if (i_idx == unit_info->num_of_files - 1) {
      if (mirrored) {
        sep_per_file = GT_DIV2(gt_encseq_num_of_sequences(unit_info->encseq))
                       - seqs_passed - 1;
      }
      else {
        sep_per_file = gt_encseq_num_of_sequences(unit_info->encseq)
                       - seqs_passed - 1;
      }
    }
    else {
      seqs_upto_fileend = gt_encseq_seqnum(unit_info->encseq,
                                           filestart +
                                           eff_file_length - 1)
                          + 1;

      sep_per_file = seqs_upto_fileend - seqs_passed - 1;
      seqs_passed = seqs_upto_fileend;
    }
    if (unit_info->map_files != NULL) {
      genome_lengths[unit_info->map_files[i_idx]] +=
        eff_file_length - sep_per_file;
    }
    else {
      genome_lengths[i_idx] = eff_file_length - sep_per_file;
    }
  }
  if (mirrored) {
    for (i_idx = 0UL; i_idx < unit_info->num_of_genomes; i_idx++) {
      genome_lengths[i_idx] += genome_lengths[i_idx];
    }
  }
  if (gt_log_enabled()) {
    for  (i_idx = 0UL; i_idx < unit_info->num_of_genomes; i_idx++) {
      gt_log_log("file/genome %lu has length %lu",
          i_idx, genome_lengths[i_idx]);
    }
  }
  return genome_lengths;
}

static double *genomediff_calculate_gc(unsigned long *genome_lengths,
                                       GtShuUnitFileInfo *unit_info,
                                       GT_UNUSED GtError *err)
{
  unsigned long file_idx = 0UL,
                genome_idx,
                idx,
                *seq_gc_content = NULL,
                num_of_seq,
                *genome_gc_content;
  double *gc_content;

  gt_assert(genome_lengths != NULL);
  gt_assert(unit_info != NULL);

  gt_assert(genome_lengths != NULL);
  gt_assert(unit_info != NULL);

  /*count gc, do not calculate ratio*/
  seq_gc_content = gt_encseq_gc_count(unit_info->encseq);

  if (seq_gc_content == NULL)
    return NULL;

  gc_content = gt_calloc((size_t)unit_info->num_of_genomes,
                         sizeof (*gc_content));
  genome_gc_content = gt_calloc((size_t)unit_info->num_of_genomes,
                                sizeof (*genome_gc_content));
  num_of_seq = gt_encseq_num_of_sequences(unit_info->encseq);

  for (idx = 0; idx < num_of_seq; idx++) {
    file_idx = gt_encseq_filenum(unit_info->encseq,
                                 gt_encseq_seqstartpos(unit_info->encseq, idx));
    if (unit_info->map_files != NULL)
      genome_idx = unit_info->map_files[file_idx];
    else
      genome_idx = file_idx;
    gt_safe_add(genome_gc_content[genome_idx],
                genome_gc_content[genome_idx],
                seq_gc_content[idx]);
  }
  for (idx = 0; idx < unit_info->num_of_genomes; idx++) {
    gt_log_log("file/genome %lu has gc %lu", idx, genome_gc_content[idx]);
  }

  for (idx = 0; idx < unit_info->num_of_genomes; idx++) {
    gt_assert(genome_gc_content[idx] <= genome_lengths[idx]);
    gc_content[idx] = (double) genome_gc_content[idx] /
                      (double) genome_lengths[idx];
  }
  gt_free(seq_gc_content);
  gt_free(genome_gc_content);
  return gc_content;
}

static void genomediff_calculate_div(GtShuUnitFileInfo *unit_info,
                                     double **div, double *gc_content,
                                     unsigned long *genome_lengths,
                                     const GtGenomediffArguments *arguments,
                                     GtTimer *timer)
{
  double query_gc,
         query_shulen,
         *ln_n_fac;
  unsigned long subject_len,
                i_idx, j_idx,
                subject, query;

  if (timer != NULL)
    gt_timer_show_progress(timer, "pre calculate ln_n_fac", stdout);

  ln_n_fac = gt_get_ln_n_fac(arguments->max_ln_n_fac);
  if (timer != NULL)
    gt_timer_show_progress(timer, "calculate divergence", stdout);
  for (i_idx = 0; i_idx < unit_info->num_of_genomes; i_idx++) {
    for (j_idx = i_idx+1; j_idx < unit_info->num_of_genomes; j_idx++) {
      /* query is the one with smaller avg shulen */
      if (gt_double_smaller_double(div[i_idx][j_idx], div[j_idx][i_idx])) {
        subject = j_idx;
        query = i_idx;
      }
      else if (gt_double_smaller_double(div[j_idx][i_idx], div[i_idx][j_idx])) {
          subject = i_idx;
          query = j_idx;
      }
      /* if avg shulen is equal, choose query with gc content farther from .5 */
      else {
        if (gt_double_smaller_double(fabs(gc_content[i_idx]-0.5),
                                     fabs(gc_content[j_idx]-0.5))) {
          subject = i_idx;
          query = j_idx;
        }
        else {
          subject = j_idx;
          query = i_idx;
        }
      }
      query_gc = gc_content[query];
      query_shulen = div[query][subject];
      subject_len = genome_lengths[subject];

      div[i_idx][j_idx] = gt_divergence(arguments->divergence_rel_err,
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

uint64_t **gt_genomediff_shulen_sum(const GtGenomediffArguments *arguments,
                                    GtShuUnitFileInfo *unit_info,
                                    GtLogger *logger,
                                    GtTimer *timer,
                                    GtError *err)
{
  int had_err = 0;
  unsigned long i_idx, j_idx;
  uint64_t **shulensums = NULL;
  const GtEncseq *encseq = NULL;
  Genericindex *generic_index_subject = NULL;
  Sequentialsuffixarrayreader *ssar = NULL;

  gt_error_check(err);
  gt_assert(unit_info != NULL);

  /* prepare index */
  if (arguments->with_esa) {
    if (timer != NULL)
      gt_timer_show_progress(timer, "load sequential sa reader", stdout);

    ssar = gt_newSequentialsuffixarrayreaderfromfile(
                                    gt_str_array_get(arguments->filenames, 0),
                                    SARR_LCPTAB | SARR_SUFTAB | SARR_ESQTAB |
                                      SARR_SSPTAB,
                                    arguments->scan ? SEQ_scan : SEQ_mappedboth,
                                    logger, err);
    if (ssar == NULL)
      had_err = -1;

    if (!had_err) {
      encseq = gt_encseqSequentialsuffixarrayreader(ssar);
    }
  }
  else {
    gt_assert(arguments->with_pck);

    if (timer != NULL)
      gt_timer_show_progress(timer, "map generic index", stdout);

    generic_index_subject =
      genericindex_new(gt_str_array_get(arguments->filenames, 0),
                       false, true, false, true,
                       arguments->user_max_depth,
                       logger, err);
    if (generic_index_subject == NULL) {
      had_err = -1;
    }
    /*get num of genomes from unitfile or encseq and allocate shulen-matrix*/
    if (!had_err) {
      encseq = genericindex_getencseq(generic_index_subject);
    }
  }

  /*calculate sums of shulens*/
  if (!had_err) {
    gt_array2dim_calloc(shulensums,
                        unit_info->num_of_genomes,
                        unit_info->num_of_genomes);

    if (arguments->with_esa) {
      if (timer != NULL)
        gt_timer_show_progress(timer, "dfs esa index", stdout);

      had_err = gt_multiesa2shulengthdist(ssar, encseq, shulensums,
                                          unit_info, err);
    }
    else {
      const FMindex *subjectindex = NULL;
      unsigned long num_of_chars,
                    total_length;
      const GtAlphabet *alphabet = gt_encseq_alphabet(unit_info->encseq);

      num_of_chars = (unsigned long) gt_alphabet_num_of_chars(alphabet);
      total_length = genericindex_get_totallength(generic_index_subject);

      subjectindex = genericindex_get_packedindex(generic_index_subject);
      gt_assert(subjectindex != NULL);

      had_err = gt_pck_calculate_shulen(subjectindex,
                                        unit_info,
                                        shulensums,
                                        num_of_chars,
                                        total_length,
                                        timer,
                                        logger,
                                        err);
    }
  }

  if (arguments->with_esa) {
    if (ssar != NULL)
      gt_freeSequentialsuffixarrayreader(&ssar);
  }
  else
    genericindex_delete(generic_index_subject);

  if (had_err) {
    gt_log_log("foo");
    if (shulensums != NULL)
      gt_array2dim_delete(shulensums);
    shulensums = NULL;
  }
  else {
    gt_assert(shulensums != NULL);
    gt_log_log("sum of shulens");
    if (gt_log_enabled()) {
      for (i_idx = 0; i_idx < unit_info->num_of_genomes; i_idx++) {
        fprintf(stderr, "debug: %s\t",
                gt_str_array_get(unit_info->genome_names, i_idx));
        for (j_idx = 0; j_idx < unit_info->num_of_genomes; j_idx++) {
          if (i_idx == j_idx)
            fprintf(stderr, "0");
          else
            fprintf(stderr, Formatuint64_t"\t",
                    PRINTuint64_tcast(shulensums[i_idx][j_idx]));
        }
        fprintf(stderr, "\n");
      }
    }
  }
  return shulensums;
}

static void genomediff_print_table(double **tab,
                                   GtShuUnitFileInfo *unit_info)
{
  unsigned long i_idx, j_idx;

  for (i_idx = 0; i_idx < unit_info->num_of_genomes; i_idx++) {
    printf("# %s\t", gt_str_array_get(unit_info->genome_names, i_idx));
    for (j_idx = 0; j_idx < unit_info->num_of_genomes; j_idx++) {
      if (i_idx == j_idx)
        printf("%.6f\t",0.0);
      else
        printf("%f\t", tab[i_idx][j_idx]);
    }
    printf("\n");
  }
}

int gt_genomediff_calculate_div_from_avg(double **avgshu,
                                         const GtGenomediffArguments *arguments,
                                         GtShuUnitFileInfo *unit_info,
                                         GtLogger *logger,
                                         GtTimer *timer,
                                         GtError *err)
{
  int had_err = 0;
  const GtAlphabet *alphabet;
  double *gc_content = NULL;
  unsigned long i_idx, j_idx;
  unsigned long *genome_lengths = NULL;

  genome_lengths = genomediff_calculate_genome_lengths(unit_info);

  /*get alphabet and check if alphabet is dna, which is presumed for
    gc-calculation*/
  alphabet = gt_encseq_alphabet(unit_info->encseq);
  if (!gt_alphabet_is_dna(alphabet)) {
    gt_error_set(err,"error: sequences need to be dna to calculate gc!");
    had_err = -1;
  }

  if (!had_err && timer != NULL)
    gt_timer_show_progress(timer, "calculate gc", stdout);
  if (!had_err) {
    gc_content = genomediff_calculate_gc(genome_lengths, unit_info, err);
    if (gc_content == NULL)
      had_err = -1;
  }

    /* calculation of divergence */
  if (!had_err) {
    genomediff_calculate_div(unit_info, avgshu, gc_content, genome_lengths,
                             arguments, timer);

    if (gt_logger_enabled(logger)) {
      gt_logger_log(logger, "table of divergences");
      genomediff_print_table(avgshu, unit_info);
    }
  }

  if (!had_err && timer != NULL)
    gt_timer_show_progress(timer, "calculate kr", stdout);

  if (!had_err) {
    gt_logger_log(logger, "# Table of Kr");
    printf("%lu\n", unit_info->num_of_genomes);
    for (i_idx = 0; i_idx < unit_info->num_of_genomes; i_idx++) {
      printf("%s\t", gt_str_array_get(unit_info->genome_names, i_idx));
      for (j_idx = 0; j_idx < unit_info->num_of_genomes; j_idx++) {
        if ( i_idx == j_idx )
          printf("%.6f\t",0.0);
        else
          printf("%f\t", gt_calculateKr(avgshu[i_idx][j_idx]));
      }
      printf("\n");
    }
  }
  gt_free(gc_content);
  return had_err;
}

int gt_genomediff_kr_calc(uint64_t **shulensums,
                          const GtGenomediffArguments *arguments,
                          GtShuUnitFileInfo *unit_info,
                          bool pcktable,
                          GtLogger *logger,
                          GtTimer *timer,
                          GtError *err)
{
  int had_err = 0;
  double **div = NULL,
         *gc_content = NULL;
  unsigned long *genome_lengths = NULL,
                i_idx, j_idx;
  const GtAlphabet *alphabet;

  gt_assert(shulensums != NULL);
  gt_assert(unit_info != NULL);

  genome_lengths =
    genomediff_calculate_genome_lengths(unit_info);

  if (timer != NULL)
    gt_timer_show_progress(timer, "calculate avg shulen", stdout);

  gt_array2dim_calloc(div,
                      unit_info->num_of_genomes,
                      unit_info->num_of_genomes);

  for (i_idx = 0; i_idx < unit_info->num_of_genomes; i_idx++) {
    unsigned long length_i;
    length_i = genome_lengths[i_idx];

    for (j_idx = 0; j_idx < unit_info->num_of_genomes; j_idx++) {
      div[i_idx][j_idx] = pcktable ?
                            ((double) shulensums[i_idx][j_idx]) / length_i :
                            ((double) shulensums[j_idx][i_idx]) / length_i;
    }
  }
  if (gt_logger_enabled(logger)) {
    gt_logger_log(logger, "table of avg shulens");
    genomediff_print_table(div, unit_info);
  }

  /*get alphabet and check if alphabet is dna, which is presumed for
    gc-calculation*/
  alphabet = gt_encseq_alphabet(unit_info->encseq);
  if (!gt_alphabet_is_dna(alphabet)) {
    gt_error_set(err,"error: sequences need to be dna to calculate gc!");
    had_err = -1;
  }

  if (!had_err && timer != NULL)
    gt_timer_show_progress(timer, "calculate gc", stdout);
  if (!had_err) {
    gc_content = genomediff_calculate_gc(genome_lengths, unit_info, err);
    if (gc_content == NULL)
      had_err = -1;
  }

    /* calculation of divergence */
  if (!had_err) {
    genomediff_calculate_div(unit_info, div, gc_content, genome_lengths,
                             arguments, timer);

    if (gt_logger_enabled(logger)) {
      gt_logger_log(logger, "table of divergences");
      genomediff_print_table(div, unit_info);
    }
  }

  if (!had_err && timer != NULL)
    gt_timer_show_progress(timer, "calculate kr", stdout);

  if (!had_err) {
    gt_logger_log(logger, "# Table of Kr");
    printf("%lu\n", unit_info->num_of_genomes);
    for (i_idx = 0; i_idx < unit_info->num_of_genomes; i_idx++) {
      printf("%s\t", gt_str_array_get(unit_info->genome_names, i_idx));
      for (j_idx = 0; j_idx < unit_info->num_of_genomes; j_idx++) {
        if ( i_idx == j_idx )
          printf("%.6f\t",0.0);
        else
          printf("%f\t", gt_calculateKr(div[i_idx][j_idx]));
      }
      printf("\n");
    }
  }

  gt_free(genome_lengths);
  gt_free(gc_content);
  if (div != NULL)
    gt_array2dim_delete(div);
  return had_err;
}
