/*
  Copyright (c) 2007      Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
  Copyright (c)      2010 Sascha Steinbiss <steinbiss@zbh.uni-hamburg.de>
  Copyright (c) 2007-2010 Center for Bioinformatics, University of Hamburg

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

#include <string.h>
#include <inttypes.h>
#include "core/basename_api.h"
#include "core/error.h"
#include "core/logger.h"
#include "core/ma_api.h"
#include "core/option_api.h"
#include "core/readmode.h"
#include "core/spacecalc.h"
#include "core/str.h"
#include "core/unused_api.h"
#include "core/versionfunc.h"
#include "match/sfx-opt.h"
#include "match/sfx-strategy.h"
#include "match/stamp.h"

static GtOPrval parse_options(int *parsed_args,
                              bool doesa,
                              Suffixeratoroptions *so,
                              int argc, const char **argv, GtError *err)
{
  GtOptionParser *op;
  GtOption *option,
           *optionshowprogress,
           *optiongenomediff,
           *optionii;
  GtOPrval oprval;
  gt_error_check(err);

  op = gt_option_parser_new("[option ...] (-db file [...] | -ii index)",
                            doesa ? "Compute enhanced suffix array."
                                  : "Compute packed index.");
  gt_option_parser_set_mail_address(op, "<kurtz@zbh.uni-hamburg.de>");

  /* input info */
  so->indexname = gt_str_new();
  so->inputindex = gt_str_new();
  so->db = gt_str_array_new();

  /* register options for encoded sequence handling */
  so->encopts = gt_encseq_options_register_encoding(op, so->indexname, so->db);
  so->loadopts = gt_encseq_options_register_loading(op, so->indexname);

  /* register options for index handling */
  if (doesa)
    so->idxopts = gt_index_options_register_esa(op, so->encopts);
  else
    so->idxopts = gt_index_options_register_packedidx(op, so->indexname,
                                                      so->encopts);

  /* verbosity */
  option = gt_option_new_verbose(&so->beverbose);
  gt_option_parser_add_option(op, option);

  optionshowprogress = gt_option_new_bool("showprogress",
                                          "show a progress bar",
                                          &so->showprogress,
                                          false);
  gt_option_parser_add_option(op, optionshowprogress);

  optionii = gt_option_new_filename("ii", "specify existing encoded sequence",
                                    so->inputindex);
  gt_option_parser_add_option(op, optionii);
  gt_option_is_mandatory_either(gt_encseq_options_db_option(so->encopts),
                                optionii);
  gt_option_exclude(gt_encseq_options_db_option(so->encopts), optionii);
  gt_option_exclude(optionii, gt_encseq_options_smap_option(so->encopts));
  gt_option_exclude(optionii, gt_encseq_options_dna_option(so->encopts));
  gt_option_exclude(optionii, gt_encseq_options_protein_option(so->encopts));
  gt_option_exclude(optionii, gt_encseq_options_plain_option(so->encopts));
  gt_option_exclude(optionii, gt_encseq_options_sat_option(so->encopts));

  optiongenomediff = gt_option_new_bool("genomediff",
                                   "directly process the lcp intervals using "
                                   "the genomediff algorithm (suffix array and "
                                   "lcp-tables are not output)",
                                   &so->genomediff,
                                   false);
  gt_option_is_extended_option(optiongenomediff);
  if (gt_index_options_outsuftab_option(so->idxopts) != NULL) {
    gt_option_exclude(optiongenomediff,
                      gt_index_options_outsuftab_option(so->idxopts));
  }
  gt_option_parser_add_option(op, optiongenomediff);

  /* suffixerator and friends do not take arguments */
  gt_option_parser_set_min_max_args(op, 0U, 0U);

  oprval = gt_option_parser_parse(op, parsed_args, argc, argv, gt_versionfunc,
                                  err);

  if (gt_str_length(so->indexname) == 0UL) {
    /* we do not have an indexname yet, so there was none given in the
       -indexname option and it could not be derived from the input filenames.
       So it must be in the -ii parameter. */
    char *basenameptr;
    basenameptr = gt_basename(gt_str_get(so->inputindex));
    gt_str_set(so->indexname, basenameptr);
    gt_free(basenameptr);
  }

  gt_option_parser_delete(op);

  return oprval;
}

static void showoptions(const Suffixeratoroptions *so)
{
  unsigned long i;
  Sfxstrategy sfxtrategy;
  GtLogger *logger = gt_logger_new(true, GT_LOGGER_DEFLT_PREFIX, stdout);

  if (gt_str_length(gt_encseq_options_smap_value(so->encopts)) > 0)
  {
    gt_logger_log_force(logger, "smap=\"%s\"",
                        gt_str_get(gt_encseq_options_smap_value(so->encopts)));
  }
  if (gt_encseq_options_dna_value(so->encopts))
  {
    gt_logger_log_force(logger, "dna=yes");
  }
  if (gt_encseq_options_protein_value(so->encopts))
  {
    gt_logger_log_force(logger, "protein=yes");
  }
  if (gt_encseq_options_plain_value(so->encopts))
  {
    gt_logger_log_force(logger, "plain=yes");
  }
  gt_logger_log_force(logger, "indexname=\"%s\"",
                    gt_str_get(so->indexname));

  if (gt_index_options_prefixlength_value(so->idxopts)
                                                   == GT_PREFIXLENGTH_AUTOMATIC)
  {
    gt_logger_log_force(logger, "prefixlength=automatic");
  } else
  {
    gt_logger_log_force(logger, "prefixlength=%u",
                        gt_index_options_prefixlength_value(so->idxopts));
  }
  sfxtrategy = gt_index_options_sfxstrategy_value(so->idxopts);
  gt_logger_log_force(logger, "storespecialcodes=%s",
                        sfxtrategy.storespecialcodes ? "true" : "false");
  for (i=0; i<gt_str_array_size(so->db); i++)
  {
    gt_logger_log_force(logger, "inputfile[%lu]=%s", i,
                   gt_str_array_get(so->db, i));
  }
  if (gt_str_length(so->inputindex) > 0)
  {
    gt_logger_log_force(logger, "inputindex=%s",
                        gt_str_get(so->inputindex));
  }
  gt_assert(gt_str_length(so->indexname) > 0);
  gt_logger_log_force(logger, "indexname=%s",
                    gt_str_get(so->indexname));
  gt_logger_log_force(logger, "outtistab=%s,outsuftab=%s,outlcptab=%s,"
                              "outbwttab=%s,outbcktab=%s,outdestab=%s,"
                              "outsdstab=%s,outssptab=%s,outkystab=%s",
          gt_encseq_options_tis_value(so->encopts) ? "true" : "false",
          gt_index_options_outsuftab_value(so->idxopts) ? "true" : "false",
          gt_index_options_outlcptab_value(so->idxopts) ? "true" : "false",
          gt_index_options_outbwttab_value(so->idxopts) ? "true" : "false",
          gt_index_options_outbcktab_value(so->idxopts) ? "true" : "false",
          gt_encseq_options_des_value(so->encopts) ? "true" : "false",
          gt_encseq_options_sds_value(so->encopts) ? "true" : "false",
          gt_encseq_options_ssp_value(so->encopts) ? "true" : "false",
          gt_index_options_outkystab_value(so->idxopts) ?
             (gt_index_options_outkyssort_value(so->idxopts) ?
                              "true with sort" : "true") :
                              "false");

  if (gt_index_options_maximumspace_value(so->idxopts) > 0)
  {
    gt_assert(gt_index_options_numofparts_value(so->idxopts) == 1U);
    gt_logger_log_force(logger, "maximumspace=%.0f MB",
            GT_MEGABYTES(gt_index_options_maximumspace_value(so->idxopts)));
  } else
  {
    gt_logger_log_force(logger, "parts=%u",
                            gt_index_options_numofparts_value(so->idxopts));
  }
  gt_logger_log_force(logger, "maxinsertionsort=%lu",
                        sfxtrategy.maxinsertionsort);
  gt_logger_log_force(logger, "maxbltriesort=%lu",
                        sfxtrategy.maxbltriesort);
  gt_logger_log_force(logger, "maxcountingsort=%lu",
                        sfxtrategy.maxcountingsort);
  gt_logger_log_force(logger, "lcpdist=%s",
                        gt_index_options_lcpdist_value(so->idxopts)
                          ? "true"
                          : "false");
  gt_logger_delete(logger);
}

void gt_sfxoptions_delete(Suffixeratoroptions *so)
{
  /* no checking if error occurs, since errors have been output before */
  gt_index_options_delete(so->idxopts);
  gt_encseq_options_delete(so->encopts);
  gt_encseq_options_delete(so->loadopts);
  gt_str_delete(so->indexname);
  gt_str_delete(so->inputindex);
  gt_str_array_delete(so->db);
}

int gt_suffixeratoroptions(Suffixeratoroptions *so,
                        bool doesa,
                        int argc,
                        const char **argv,
                        GtError *err)
{
  int parsed_args, retval = 0;
  GtOPrval rval;

  gt_error_check(err);
  rval = parse_options(&parsed_args, doesa, so, argc, argv, err);
  if (rval == GT_OPTION_PARSER_ERROR)
  {
    retval = -1;
  } else
  {
    if (rval == GT_OPTION_PARSER_REQUESTS_EXIT)
    {
      retval = 2;
    } else
    {
      if (rval == GT_OPTION_PARSER_OK)
      {
        if (so->beverbose)
        {
          showoptions(so);
        }
      }
    }
  }
  return retval;
}
