/*
  Copyright (c) 2007-2015 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
  Copyright (c) 2007-2015 Center for Bioinformatics, University of Hamburg

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

#include <float.h>
#include "core/error_api.h"
#include "core/format64.h"
#include "core/log_api.h"
#include "core/logger.h"
#include "core/ma_api.h"
#include "core/option_api.h"
#include "core/str_api.h"
#include "core/tool_api.h"
#include "core/unused_api.h"
#include "core/versionfunc_api.h"
#include "core/minmax_api.h"
#include "core/encseq.h"
#include "core/showtime.h"
#include "core/timer_api.h"
#include "core/encseq_metadata.h"
#include "match/esa-maxpairs.h"
#include "match/esa-mmsearch.h"
#include "match/querymatch.h"
#include "match/test-maxpairs.h"
#include "match/seed-extend.h"
#include "match/esa-map.h"
#include "tools/gt_repfind.h"

typedef struct
{
  unsigned int userdefinedleastlength, /* minalilen */
               seedlength; /* like kmerlen, but here MEMs are used */
  GtXdropscore xdropbelowscore;
  GtUword samples,
          extendxdrop,
          maxfreq,
          historysize,   /* number of bits used for history of alignments */
          perc_mat_history, /* percent of matches in history */
          minidentity, /* We prefer to specify the minidentity. The use of the
           notion of error percentage may be misleading, as in Myers paper it
           refers to the percentage of errors in a read. If (for compatibility
           reasons, the option -err is used, then the minidentity contains the
           error rate, a value in the range from 1 to 30. */
          maxalignedlendifference, /* maxfrontdist */
          extendgreedy; /* determines which of the tables in
                           seed-extend-params.h is used */
  bool scanfile, beverbose, forward, reverse, reverse_complement, searchspm,
       check_extend_symmetry, trimstat_on, noxpolish, verify_alignment;
  GtStr *indexname, *query_indexname, *cam_string;
  GtStrArray *query_files;
  GtOption *refforwardoption,
           *refseedlengthoption,
           *refuserdefinedleastlengthoption,
           *refextendxdropoption,
           *refextendgreedyoption,
           *ref_op_evalue;
  GtStrArray *display_args;
  double evalue_threshold;
} GtMaxpairsoptions;

static int gt_exact_selfmatch_with_output(void *info,
                                          const GtGenericEncseq *genericencseq,
                                          GtUword len,
                                          GtUword pos1,
                                          GtUword pos2,
                                          GT_UNUSED GtError *err)
{
  GtUword queryseqnum, query_seqstart, query_seqlen, dbseqnum,
          db_seqstart, dbseqlen;
  const GtEncseq *encseq;
  GtProcessinfo_and_querymatchspaceptr *info_querymatch
    = (GtProcessinfo_and_querymatchspaceptr *) info;
  GtSeqorEncseq dbes, queryes;

  gt_assert(pos1 < pos2 && genericencseq != NULL && genericencseq->hasencseq);
  encseq = genericencseq->seqptr.encseq;
  if (gt_encseq_has_multiseq_support(encseq))
  {
    dbseqnum = gt_encseq_seqnum(encseq,pos1);
    db_seqstart = gt_encseq_seqstartpos(encseq,dbseqnum);
    dbseqlen = gt_encseq_seqlength(encseq,dbseqnum);
    queryseqnum = gt_encseq_seqnum(encseq,pos2);
    query_seqstart = gt_encseq_seqstartpos(encseq,queryseqnum);
    query_seqlen = gt_encseq_seqlength(encseq,queryseqnum);
  } else
  {
    dbseqnum = 0;
    db_seqstart = 0;
    dbseqlen = 0;
    queryseqnum = 0;
    query_seqstart = 0;
    query_seqlen = 0;
  }
  gt_assert(pos2 >= query_seqstart);
  GT_SEQORENCSEQ_INIT_ENCSEQ(&dbes,encseq);
  GT_SEQORENCSEQ_INIT_ENCSEQ(&queryes,encseq);
  if (gt_querymatch_complete(info_querymatch->querymatchspaceptr,
                             info_querymatch->out_display_flag,
                             len,
                             dbseqnum,
                             pos1 - db_seqstart,
                             db_seqstart,
                             dbseqlen,
                             0, /* score */
                             0, /* edist */
                             0, /* mismatches */
                             true,
                             (uint64_t) queryseqnum,
                             len,
                             pos2 - query_seqstart,
                             &dbes,
                             &queryes,
                             query_seqstart,
                             query_seqlen,
                             pos1 - db_seqstart,
                             pos2 - query_seqstart,
                             len,
                             false))
  {
    /* for exact matches we do not output evalues and bitscores */
    gt_querymatch_prettyprint(DBL_MAX,DBL_MAX,
                              info_querymatch->out_display_flag,
                              info_querymatch->querymatchspaceptr);
  }
  return 0;
}

static int gt_suffix_prefix_match_with_output(GT_UNUSED void *info,
                                              const GtGenericEncseq
                                                *genericencseq,
                                              GtUword matchlen,
                                              GtUword pos1,
                                              GtUword pos2,
                                              GT_UNUSED GtError *err)
{
  GtUword seqnum1, relpos1, seqnum2, relpos2, seqstartpos;
  const GtEncseq *encseq;

  gt_assert(pos1 < pos2);
  gt_assert(genericencseq != NULL && genericencseq->hasencseq);
  encseq = genericencseq->seqptr.encseq;
  seqnum1 = gt_encseq_seqnum(encseq,pos1);
  seqstartpos = gt_encseq_seqstartpos(encseq, seqnum1);
  gt_assert(seqstartpos <= pos1);
  relpos1 = pos1 - seqstartpos;
  seqnum2 = gt_encseq_seqnum(encseq,pos2);
  seqstartpos = gt_encseq_seqstartpos(encseq, seqnum2);
  gt_assert(seqstartpos <= pos2);
  relpos2 = pos2 - seqstartpos;
  if (relpos1 == 0)
  {
    GtUword seqlen2 = gt_encseq_seqlength(encseq,seqnum2);

    if (relpos2 + matchlen == seqlen2)
    {
      printf(GT_WU " " GT_WU " " GT_WU "\n",seqnum2,seqnum1,matchlen);
    }
  } else
  {
    if (relpos2 == 0)
    {
      GtUword seqlen1 = gt_encseq_seqlength(encseq,seqnum1);

      if (relpos1 + matchlen == seqlen1)
      {
        printf(GT_WU " " GT_WU " " GT_WU "\n",seqnum1,seqnum2,matchlen);
      }
    }
  }
  return 0;
}

static void *gt_repfind_arguments_new(void)
{
  GtMaxpairsoptions *arguments;

  arguments = gt_malloc(sizeof (*arguments));
  arguments->indexname = gt_str_new();
  arguments->query_indexname = gt_str_new();
  arguments->cam_string = gt_str_new();
  arguments->query_files = gt_str_array_new();
  arguments->display_args = gt_str_array_new();
  return arguments;
}

static void gt_repfind_arguments_delete(void *tool_arguments)
{
  GtMaxpairsoptions *arguments = tool_arguments;

  if (!arguments)
  {
    return;
  }
  gt_str_delete(arguments->indexname);
  gt_str_delete(arguments->query_indexname);
  gt_str_delete(arguments->cam_string);
  gt_str_array_delete(arguments->query_files);
  gt_str_array_delete(arguments->display_args);
  gt_option_delete(arguments->refforwardoption);
  gt_option_delete(arguments->refseedlengthoption);
  gt_option_delete(arguments->refuserdefinedleastlengthoption);
  gt_option_delete(arguments->refextendxdropoption);
  gt_option_delete(arguments->refextendgreedyoption);
  gt_option_delete(arguments->ref_op_evalue);
  gt_free(arguments);
}

static GtOptionParser *gt_repfind_option_parser_new(void *tool_arguments)
{
  const GtUword extension_sensitivity = 97;
  GtOptionParser *op;
  GtOption *option, *reverseoption, *reverse_complementoption,
           *option_query_files, *extendxdropoption,
           *extendgreedyoption, *scanoption, *sampleoption, *forwardoption,
           *spmoption, *seedlengthoption, *minidentityoption,
           *maxalilendiffoption, *leastlength_option, *char_access_mode_option,
           *check_extend_symmetry_option, *xdropbelowoption, *historyoption,
           *percmathistoryoption, *errorpercentageoption, *optiontrimstat,
           *optionnoxpolish, *verify_alignment_option, *option_query_indexname,
           *op_evalue;
  GtMaxpairsoptions *arguments = tool_arguments;

  op = gt_option_parser_new("[options] -ii indexname",
                                  "Compute maximal exact matches (and more).");
  gt_option_parser_set_mail_address(op,"<kurtz@zbh.uni-hamburg.de>");

  leastlength_option
    = gt_option_new_uint("l","Specify minimum length of matches",
                         &arguments->userdefinedleastlength,
                         0U);
  gt_option_parser_add_option(op, leastlength_option);
  arguments->refuserdefinedleastlengthoption
    = gt_option_ref(leastlength_option);

  forwardoption = gt_option_new_bool("f","Compute forward matches",
                                     &arguments->forward,
                                     true);
  gt_option_parser_add_option(op, forwardoption);
  arguments->refforwardoption = gt_option_ref(forwardoption);

  reverseoption = gt_option_new_bool("r","Compute reverse matches",
                                     &arguments->reverse,
                                     false);
  gt_option_parser_add_option(op, reverseoption);

  reverse_complementoption = gt_option_new_bool("p",
                                        "Compute matches on reverse strand",
                                        &arguments->reverse_complement,
                                        false);
  gt_option_parser_add_option(op, reverse_complementoption);

  seedlengthoption = gt_option_new_uint_min("seedlength",
                                             "Specify minimum length of seed",
                                             &arguments->seedlength,
                                             0,
                                             1UL);
  gt_option_parser_add_option(op, seedlengthoption);
  arguments->refseedlengthoption = gt_option_ref(seedlengthoption);

  option = gt_option_new_uword_min("maxfreq",
                                   "Specify maximal frequency of maximal exact "
                                   "matches in reference sequence",
                                   &arguments->maxfreq,
                                   0,
                                   2UL);
  gt_option_parser_add_option(op, option);

  extendxdropoption
    = gt_option_new_uword_min_max("extendxdrop",
                         "Extend seed to both sides using xdrop algorithm, "
                         "optional parameter specifies sensitivity",
                         &arguments->extendxdrop,
                         extension_sensitivity, 90, 100);
  gt_option_argument_is_optional(extendxdropoption);
  gt_option_parser_add_option(op, extendxdropoption);
  arguments->refextendxdropoption = gt_option_ref(extendxdropoption);

  xdropbelowoption = gt_option_new_word("xdropbelow",
                                        "Specify xdrop cutoff score "
                                        "(argument 0 means undefined). If "
                                        "undefined an optimal value is "
                                        "determined automatically depending "
                                        "on the error "
                                        "rate",
                                        &arguments->xdropbelowscore,
                                        0);
  gt_option_parser_add_option(op, xdropbelowoption);

  extendgreedyoption
    = gt_option_new_uword_min_max("extendgreedy",
                                  "Extend seed to both sides using "
                                  "greedy algorithm with trimming of waves, "
                                  "optional parameter specifies sensitivity",
                                  &arguments->extendgreedy,
                                  extension_sensitivity, 90, 100);
  gt_option_argument_is_optional(extendgreedyoption);
  gt_option_parser_add_option(op, extendgreedyoption);
  arguments->refextendgreedyoption = gt_option_ref(extendgreedyoption);

  errorpercentageoption
    = gt_option_new_uword_min_max("err","Specify error percentage of matches "
                         "as integer in the range from 1 to 30 "
                         "(for xdrop and greedy extension) [deprecated option, "
                         "kept for backwards compatibility]",
                         &arguments->minidentity,
                         10,
                         1,
                         100 - GT_EXTEND_MIN_IDENTITY_PERCENTAGE);
  gt_option_is_development_option(errorpercentageoption);
  gt_option_parser_add_option(op, errorpercentageoption);

  minidentityoption
    = gt_option_new_uword_min_max("minidentity",
                          "Specify minimum identity of matches\n"
                          "as integer in the range from 70 to 99 "
                          "(for xdrop and greedy extension)",
                          &arguments->minidentity,
                          80,
                          GT_EXTEND_MIN_IDENTITY_PERCENTAGE,99);
  gt_option_parser_add_option(op, minidentityoption);

  maxalilendiffoption
    = gt_option_new_uword("maxalilendiff","Specify maximum difference of "
                          "alignment length (trimming for greedy extension), "
                          "if option is not used or parameter 0 is specified, "
                          "then good value is automatically chosen",
                          &arguments->maxalignedlendifference,
                          0);
  gt_option_parser_add_option(op, maxalilendiffoption);
  gt_option_is_development_option(maxalilendiffoption);

  historyoption
    = gt_option_new_uword_min_max("history",
                                  "Specify size of history in range [1..64] "
                                  "(trimming for greedy extension)",
                                  &arguments->historysize,
                                  60,
                                  0,
                                  64);
  gt_option_parser_add_option(op, historyoption);
  gt_option_is_development_option(historyoption);

  percmathistoryoption
    = gt_option_new_uword_min_max("percmathistory",
                                  "percentage of matches required in history",
                                  &arguments->perc_mat_history,
                                  0,
                                  1,
                                  100);
  gt_option_parser_add_option(op, percmathistoryoption);
  gt_option_is_development_option(percmathistoryoption);

  option = gt_option_new_string("ii",
                                "Specify input index",
                                arguments->indexname, NULL);
  gt_option_parser_add_option(op, option);
  gt_option_is_mandatory(option);

  char_access_mode_option = gt_option_new_string("cam",
                                                 gt_cam_extendgreedy_comment(),
                                                 arguments->cam_string,
                                                 "any,any");
  gt_option_parser_add_option(op, char_access_mode_option);
  gt_option_is_development_option(char_access_mode_option);

  /* -outfmt */
  option = gt_option_new_string_array("outfmt",
                                      gt_querymatch_display_help(),
                                      arguments->display_args);
  gt_option_parser_add_option(op, option);

  /* -evalue */
  op_evalue = gt_option_new_double("evalue","switch on evalue filtering of "
                                            "matches (optional argument "
                                            "specifies evalue threshold)",
                                   &arguments->evalue_threshold,
                                   10.0);
  gt_option_parser_add_option(op, op_evalue);
  gt_option_argument_is_optional(op_evalue);
  arguments->ref_op_evalue = gt_option_ref(op_evalue);

  optionnoxpolish
    = gt_option_new_bool("noxpolish","do not polish X-drop extensions",
                         &arguments->noxpolish, false);
  gt_option_parser_add_option(op, optionnoxpolish);
  gt_option_is_development_option(optionnoxpolish);

  verify_alignment_option
    = gt_option_new_bool("verify-alignment","verify correctness of alignments",
                         &arguments->verify_alignment, false);
  gt_option_parser_add_option(op, verify_alignment_option);
  gt_option_is_development_option(verify_alignment_option);

  optiontrimstat = gt_option_new_bool("trimstat","show trimming statistics",
                                      &arguments->trimstat_on, false);
  gt_option_parser_add_option(op, optiontrimstat);
  gt_option_is_development_option(optiontrimstat);

  option_query_indexname = gt_option_new_string("qii",
                                                "Specify name of query index",
                                                arguments->query_indexname,
                                                NULL);
  gt_option_is_development_option(option_query_indexname);
  gt_option_parser_add_option(op, option_query_indexname);

  /* the following option are options special to repfind */

  option_query_files = gt_option_new_filename_array("q",
                                                    "Specify query files",
                                                    arguments->query_files);
  gt_option_is_development_option(option_query_files);
  gt_option_parser_add_option(op, option_query_files);

  sampleoption = gt_option_new_uword_min("samples","Specify number of samples",
                                         &arguments->samples,
                                         0,
                                         1UL);
  gt_option_is_development_option(sampleoption);
  gt_option_parser_add_option(op, sampleoption);

  spmoption = gt_option_new_bool("spm","Search for suffix prefix matches",
                                       &arguments->searchspm,
                                       false);
  gt_option_is_development_option(spmoption);
  gt_option_parser_add_option(op, spmoption);

  check_extend_symmetry_option
    = gt_option_new_bool("check_extend_symmetry",
                         "check that left/right greedy extension is symmetric "
                         "for sequences mirror around seed",
                         &arguments->check_extend_symmetry,
                         false);
  gt_option_parser_add_option(op, check_extend_symmetry_option);
  gt_option_is_development_option(check_extend_symmetry_option);

  scanoption = gt_option_new_bool("scan","scan index rather than map "
                                         "it to main memory",
                                  &arguments->scanfile,
                                  false);
  gt_option_parser_add_option(op, scanoption);

  option = gt_option_new_verbose(&arguments->beverbose);
  gt_option_parser_add_option(op, option);

  gt_option_exclude(option_query_files,sampleoption);
  gt_option_exclude(option_query_files,scanoption);
  gt_option_exclude(option_query_files,spmoption);
  gt_option_exclude(option_query_files,option_query_indexname);
  gt_option_exclude(sampleoption,spmoption);
  gt_option_exclude(reverseoption,spmoption);
  gt_option_exclude(reverse_complementoption,spmoption);
  gt_option_exclude(extendgreedyoption,extendxdropoption);
  gt_option_exclude(errorpercentageoption,minidentityoption);
  gt_option_imply(xdropbelowoption,extendxdropoption);
  gt_option_imply(historyoption,extendgreedyoption);
  gt_option_imply(maxalilendiffoption,extendgreedyoption);
  gt_option_imply(percmathistoryoption,extendgreedyoption);
  gt_option_imply(optiontrimstat,extendgreedyoption);
  gt_option_imply(optionnoxpolish,extendxdropoption);
  gt_option_imply_either_2(seedlengthoption,extendxdropoption,
                           extendgreedyoption);
  gt_option_imply_either_2(minidentityoption,extendxdropoption,
                           extendgreedyoption);
  gt_option_imply_either_2(errorpercentageoption,extendxdropoption,
                           extendgreedyoption);
  return op;
}

static int gt_repfind_arguments_check(GT_UNUSED int rest_argc,
                                      void *tool_arguments,
                                      GT_UNUSED GtError *err)
{
  GtMaxpairsoptions *arguments = tool_arguments;

  if (!gt_option_is_set(arguments->refforwardoption) &&
      (arguments->reverse || arguments->reverse_complement))
  {
    arguments->forward = false;
  }
  if (!gt_option_is_set(arguments->refuserdefinedleastlengthoption))
  {
    if (!gt_option_is_set(arguments->refseedlengthoption))
    {
      arguments->seedlength = arguments->userdefinedleastlength = 20U;
    } else
    {
      arguments->userdefinedleastlength = arguments->seedlength;
    }
  } else
  {
    if (!gt_option_is_set(arguments->refseedlengthoption))
    {
      arguments->seedlength = arguments->userdefinedleastlength;
    } else
    {
      if (arguments->seedlength > arguments->userdefinedleastlength)
      {
        arguments->seedlength = arguments->userdefinedleastlength;
      }
    }
  }
  if (!gt_option_is_set(arguments->ref_op_evalue))
  {
    arguments->evalue_threshold = DBL_MAX;
  }
  return 0;
}

static int gt_generic_extend_selfmatch_xdrop_with_output(
                                           void *info,
                                           const GtGenericEncseq *genericencseq,
                                           GtUword len,
                                           GtUword pos1,
                                           GtUword pos2,
                                           GtError *err)
{
  gt_assert(genericencseq != NULL && genericencseq->hasencseq);
  return gt_rf_xdrop_extend_selfmatch_with_output(info,
                                                  genericencseq->seqptr.encseq,
                                                  len,
                                                  pos1,
                                                  pos2,
                                                  err);
}

static int gt_generic_simplegreedyselfmatchoutput(
                                           void *processinfo,
                                           const GtGenericEncseq *genericencseq,
                                           GtUword len,
                                           GtUword pos1,
                                           GtUword pos2,
                                           GtError *err)
{
  gt_assert(genericencseq != NULL && genericencseq->hasencseq);
  return gt_rf_greedy_extend_selfmatch_with_output(processinfo,
                                                   genericencseq->seqptr.encseq,
                                                   len,
                                                   pos1,
                                                   pos2,
                                                   err);
}

typedef void (*Gt_extend_querymatch_func)(void *,
                                          const GtEncseq *,
                                          const GtQuerymatch *,
                                          const GtSeqorEncseq *,
                                          bool);

static int gt_callenumquerymatches(bool selfmatch,
                                   const char *indexname,
                                   const GtStrArray *query_files,
                                   const GtStr *query_indexname,
                                   GtReadmode query_readmode,
                                   unsigned int userdefinedleastlength,
                                   GtQuerymatchoutoptions *querymatchoutoptions,
                                   Gt_extend_querymatch_func eqmf,
                                   void *eqmf_data,
                                   const GtSeedExtendDisplayFlag
                                      *out_display_flag,
                                   GtLogger *logger,
                                   GtError *err)
{
  Suffixarray suffixarray;
  GtQuerysubstringmatchiterator *qsmi = NULL;
  bool haserr = false, query_encseq_own = false;
  GtEncseq *query_encseq = NULL;
  GtUword totallength = 0;

  if (gt_mapsuffixarray(&suffixarray,
                        SARR_ESQTAB | SARR_SUFTAB | SARR_SSPTAB,
                        indexname,
                        logger,
                        err) != 0)
  {
    haserr = true;
  }
  if (!haserr)
  {
    if (query_files == NULL || gt_str_array_size(query_files) == 0)
    {
      if (query_indexname == NULL || gt_str_length(query_indexname) == 0)
      {
        query_encseq = suffixarray.encseq;
      } else
      {
        GtEncseqLoader *el = gt_encseq_loader_new();

        gt_encseq_loader_require_ssp_tab(el);
        gt_encseq_loader_set_logger(el, logger);
        query_encseq = gt_encseq_loader_load(el, gt_str_get(query_indexname),
                                             err);
        gt_encseq_loader_delete(el);
        query_encseq_own = true;
        if (query_encseq == NULL)
        {
          haserr = true;
        }
      }
    } else
    {
      gt_assert(query_indexname == NULL || gt_str_length(query_indexname) == 0);
    }
  }
  if (!haserr)
  {
    totallength = gt_encseq_total_length(suffixarray.encseq);
    qsmi = gt_querysubstringmatchiterator_new(suffixarray.encseq,
                                              totallength,
                                              suffixarray.suftab,
                                              suffixarray.readmode,
                                              totallength + 1,
                                              query_files,
                                              query_encseq,
                                              query_readmode,
                                              userdefinedleastlength,
                                              err);
    if (qsmi == NULL)
    {
      haserr = true;
    }
  }
  if (!haserr)
  {
    int retval;
    GtSeqorEncseq query_seqorencseq;
    GtQuerymatch *exactseed = gt_querymatch_new();
    bool same_encseq;

    if (query_files == NULL || gt_str_array_size(query_files) == 0)
    {
      same_encseq = (suffixarray.encseq == query_encseq) ? true : false;
    } else
    {
      same_encseq = false;
    }
    if (querymatchoutoptions != NULL)
    {
      gt_querymatch_outoptions_set(exactseed,querymatchoutoptions);
    }
    gt_querymatch_query_readmode_set(exactseed,query_readmode);
    while (!haserr &&
           (retval = gt_querysubstringmatchiterator_next(qsmi, err)) == 0)
    {
      GtUword dbstart, dbseqnum, db_seqstart, dbseqlen,
              matchlength, query_seqlen, querystart, query_seqstart;
      uint64_t queryunitnum;

      dbstart = gt_querysubstringmatchiterator_dbstart(qsmi);
      if (gt_encseq_has_multiseq_support(suffixarray.encseq))
      {
        dbseqnum = gt_encseq_seqnum(suffixarray.encseq,dbstart);
        dbseqlen = gt_encseq_seqlength(suffixarray.encseq, dbseqnum);
        db_seqstart = gt_encseq_seqstartpos(suffixarray.encseq, dbseqnum);
      } else
      {
        dbseqnum = dbseqlen = db_seqstart = 0;
      }
      matchlength = gt_querysubstringmatchiterator_matchlength(qsmi);
      query_seqlen = gt_querysubstringmatchiterator_query_seqlen(qsmi);
      queryunitnum = gt_querysubstringmatchiterator_queryunitnum(qsmi);
      if (query_files == NULL || gt_str_array_size(query_files) == 0)
      {
        GT_SEQORENCSEQ_INIT_ENCSEQ(&query_seqorencseq,query_encseq);
        query_seqstart = gt_encseq_seqstartpos(query_encseq,queryunitnum);
      } else
      {
        GT_SEQORENCSEQ_INIT_SEQ(&query_seqorencseq,
                                gt_querysubstringmatchiterator_query(qsmi),
                                gt_querysubstringmatchiterator_desc(qsmi),
                                query_seqlen,
                                NULL,
                                0,
                                true);
        query_seqstart = 0;
      }
      querystart = gt_querysubstringmatchiterator_querystart(qsmi);
      if (eqmf != NULL)
      {
        gt_querymatch_init(exactseed,
                           matchlength,
                           dbseqnum,
                           dbstart - db_seqstart,
                           db_seqstart,
                           dbseqlen,
                           0, /* score */
                           0, /* edist */
                           0, /* mismatches */
                           selfmatch,
                           queryunitnum,
                           matchlength,
                           querystart,
                           query_seqstart,
                           query_seqlen,
                           NULL,
                           NULL);
        eqmf(eqmf_data,suffixarray.encseq,exactseed,&query_seqorencseq,
             same_encseq);
      } else
      {
        GtSeqorEncseq dbes;

        GT_SEQORENCSEQ_INIT_ENCSEQ(&dbes,suffixarray.encseq);
        if (gt_querymatch_complete(exactseed,
                                   out_display_flag,
                                   matchlength,
                                   dbseqnum,
                                   dbstart - db_seqstart,
                                   db_seqstart,
                                   dbseqlen,
                                   0, /* score */
                                   0, /* edist */
                                   0, /* mismatches */
                                   selfmatch,
                                   queryunitnum,
                                   matchlength,
                                   querystart,
                                   &dbes,
                                   &query_seqorencseq,
                                   query_seqstart,
                                   query_seqlen,
                                   dbstart - db_seqstart,
                                   querystart - query_seqstart,
                                   matchlength,
                                   false))
        {
          /* for exact matches we do not output evalues and bitscores */
          gt_querymatch_prettyprint(DBL_MAX,DBL_MAX,out_display_flag,exactseed);
        }
      }
    }
    if (retval == -1)
    {
      haserr = true;
    }
    gt_querymatch_delete(exactseed);
  }
  gt_querysubstringmatchiterator_delete(qsmi);
  gt_freesuffixarray(&suffixarray);
  if (query_encseq_own)
  {
    gt_encseq_delete(query_encseq);
  }
  return haserr ? -1 : 0;
}

static int gt_repfind_runner(int argc,const char **argv, int parsed_args,
                             void *tool_arguments, GtError *err)
{
  bool haserr = false;
  GtMaxpairsoptions *arguments = (GtMaxpairsoptions *) tool_arguments;
  GtLogger *logger = NULL;
  GtXdropmatchinfo *xdropmatchinfo = NULL;
  GtGreedyextendmatchinfo *greedyextendmatchinfo = NULL;
  GtTimer *repfindtimer = NULL;
  GtExtendCharAccess cam_a = GT_EXTEND_CHAR_ACCESS_ANY,
                     cam_b = GT_EXTEND_CHAR_ACCESS_ANY;
  GtFtPolishing_info *pol_info = NULL;
  GtQuerymatchoutoptions *querymatchoutoptions;
  GtKarlinAltschulStat *karlin_altschul_stat = NULL;
  Gt_extend_querymatch_func eqmf = NULL;
  void *eqmf_data = NULL;
  int mode;
  const int modes[] = {GT_READMODE_FORWARD,
                       GT_READMODE_REVERSE,
                       GT_READMODE_REVCOMPL};
  const bool flags[] = {arguments->forward,
                        arguments->reverse,
                        arguments->reverse_complement};
  GtSeedExtendDisplayFlag *out_display_flag = NULL;
  GtFtTrimstat *trimstat = NULL;

  gt_error_check(err);
  logger = gt_logger_new(arguments->beverbose, GT_LOGGER_DEFLT_PREFIX, stdout);
  if (gt_showtime_enabled())
  {
    repfindtimer = gt_timer_new();
    gt_timer_start(repfindtimer);
  }
  if (parsed_args < argc)
  {
    gt_error_set(err,"superfluous arguments: \"%s\"",argv[argc-1]);
    haserr = true;
  }
  if (!haserr)
  {
    GtSeedExtendDisplaySetMode setmode;

    if (gt_option_is_set(arguments->refextendxdropoption) ||
        gt_option_is_set(arguments->refextendgreedyoption))
    {
      setmode = GT_SEED_EXTEND_DISPLAY_SET_STANDARD;
    } else
    {
      setmode = GT_SEED_EXTEND_DISPLAY_SET_EXACT;
    }
    out_display_flag = gt_querymatch_display_flag_new(arguments->display_args,
                                                      setmode,err);
    if (out_display_flag == NULL)
    {
      haserr = true;
    } else
    {
      gt_querymatch_Options_output(stdout,argc,argv,true,arguments->minidentity,
                                   arguments->historysize);
      gt_querymatch_Fields_output(stdout,out_display_flag);
    }
  }
  if (!haserr)
  {
    if (gt_querymatch_evalue_display(out_display_flag) ||
        gt_querymatch_bitscore_display(out_display_flag) ||
        arguments->evalue_threshold != DBL_MAX)
    {
      GtEncseqMetadata *emd
        = gt_encseq_metadata_new(gt_str_get(arguments->indexname),err);
      if (emd == NULL)
      {
        haserr = true;
      } else
      {
        karlin_altschul_stat = gt_karlin_altschul_stat_new_gapped(
                                     gt_encseq_metadata_total_length(emd),
                                     gt_encseq_metadata_num_of_sequences(emd),
                                     NULL);
        gt_encseq_metadata_delete(emd);
      }
    }
  }
  if (!haserr && gt_option_is_set(arguments->refextendxdropoption))
  {
    xdropmatchinfo
      = gt_xdrop_matchinfo_new(arguments->userdefinedleastlength,
                               gt_minidentity2errorpercentage(
                                            arguments->minidentity),
                               arguments->evalue_threshold,
                               arguments->xdropbelowscore,
                               arguments->extendxdrop);
    gt_assert(xdropmatchinfo != NULL);
  }
  if (!haserr)
  {
    if (gt_option_is_set(arguments->refextendgreedyoption) ||
        gt_querymatch_alignment_display(out_display_flag) ||
        gt_option_is_set(arguments->refextendxdropoption))
    {
      if (gt_greedy_extend_char_access(&cam_a,
                                       &cam_b,
                                       gt_str_get(arguments->cam_string),err)
         != 0)
      {
        haserr = true;
      }
    }
  }
  if (!haserr && gt_option_is_set(arguments->refextendgreedyoption))
  {
    GtUword errorpercentage = gt_minidentity2errorpercentage(
                                           arguments->minidentity);
    pol_info = polishing_info_new_with_bias(errorpercentage,
                                            GT_DEFAULT_MATCHSCORE_BIAS,
                                            arguments->historysize);
    greedyextendmatchinfo
      = gt_greedy_extend_matchinfo_new(arguments->maxalignedlendifference,
                                       arguments->historysize,
                                       arguments->perc_mat_history,
                                       arguments->userdefinedleastlength,
                                       errorpercentage,
                                       arguments->evalue_threshold,
                                       cam_a,
                                       cam_b,
                                       false,
                                       arguments->extendgreedy,
                                       pol_info);
    if (arguments->check_extend_symmetry)
    {
      gt_greedy_extend_matchinfo_check_extend_symmetry_set(
                                                greedyextendmatchinfo);
    }
    if (arguments->trimstat_on)
    {
      trimstat = gt_ft_trimstat_new();
      gt_greedy_extend_matchinfo_trimstat_set(greedyextendmatchinfo,trimstat);
    }
  }
  if (!haserr)
  {
    GtEncseq *encseq_for_desc = NULL;
    GtProcessinfo_and_querymatchspaceptr info_querymatch
      = Initializer_GtProcessinfo_and_querymatchspaceptr;
    info_querymatch.karlin_altschul_stat = karlin_altschul_stat;
    info_querymatch.out_display_flag = out_display_flag;
    if (gt_querymatch_alignment_display(out_display_flag) ||
        gt_querymatch_trace_display(out_display_flag) ||
        gt_querymatch_dtrace_display(out_display_flag) ||
        gt_querymatch_cigar_display(out_display_flag) ||
        gt_querymatch_cigarX_display(out_display_flag) ||
        (gt_option_is_set(arguments->refextendxdropoption) &&
         !arguments->noxpolish))
    {
      querymatchoutoptions
        = gt_querymatchoutoptions_new(out_display_flag,
                                      gt_str_get(arguments->indexname),err);
      if (querymatchoutoptions == NULL)
      {
        haserr = true;
      }
      if (!haserr && (gt_option_is_set(arguments->refextendxdropoption) ||
                      gt_option_is_set(arguments->refextendgreedyoption)))
      {
        const bool cam_generic = false;
        const bool weakends = false;
        const GtUword sensitivity
          = gt_option_is_set(arguments->refextendgreedyoption)
              ? arguments->extendgreedy
              : 100;
        gt_querymatchoutoptions_extend(querymatchoutoptions,
                                       gt_minidentity2errorpercentage(
                                               arguments->minidentity),
                                      arguments->evalue_threshold,
                                      arguments->maxalignedlendifference,
                                      arguments->historysize,
                                      arguments->perc_mat_history,
                                      cam_a,
                                      cam_b,
                                      cam_generic,
                                      weakends,
                                      sensitivity,
                                      GT_DEFAULT_MATCHSCORE_BIAS,
                                      true,
                                      out_display_flag);
      }
    } else
    {
      querymatchoutoptions = NULL;
    }
    if (!haserr)
    {
      info_querymatch.querymatchspaceptr = gt_querymatch_new();
      if (querymatchoutoptions != NULL)
      {
        gt_querymatch_outoptions_set(info_querymatch.querymatchspaceptr,
                                     querymatchoutoptions);
      }
      if (arguments->verify_alignment)
      {
        gt_querymatch_verify_alignment_set(info_querymatch.querymatchspaceptr);
      }
      if (gt_option_is_set(arguments->refextendxdropoption))
      {
        eqmf = gt_rf_xdrop_extend_querymatch_with_output;
        info_querymatch.processinfo = xdropmatchinfo;
        eqmf_data = (void *) &info_querymatch;
      } else
      {
        if (gt_option_is_set(arguments->refextendgreedyoption))
        {
          eqmf = gt_rf_greedy_extend_querymatch_with_output;
          info_querymatch.processinfo = greedyextendmatchinfo;
          eqmf_data = (void *) &info_querymatch;
        }
      }
      if (gt_querymatch_subjectid_display(out_display_flag))
      {
        GtEncseqLoader *encseq_loader = gt_encseq_loader_new();
        gt_encseq_loader_require_des_tab(encseq_loader);
        gt_encseq_loader_require_sds_tab(encseq_loader);
        encseq_for_desc
          = gt_encseq_loader_load(encseq_loader,
                                  gt_str_get(arguments->indexname),
                                  err);
        gt_encseq_loader_delete(encseq_loader);
        if (encseq_for_desc == NULL)
        {
          haserr = true;
        }
      }
    }
    if (!haserr && gt_str_array_size(arguments->query_files) == 0 &&
        gt_str_length(arguments->query_indexname) == 0)
    {
      if (arguments->samples > 0)
      {
        if (gt_testmaxpairs(gt_str_get(arguments->indexname),
                            arguments->samples,
                            arguments->seedlength,
                            (GtUword)
                            (100 * arguments->seedlength),
                            logger,
                            err) != 0)
        {
          haserr = true;
        }
      } else
      {
        if (!haserr && arguments->forward)
        {
          GtProcessmaxpairs processmaxpairs;
          void *processmaxpairsdata;

          if (arguments->searchspm)
          {
            processmaxpairs = gt_suffix_prefix_match_with_output;
            processmaxpairsdata = NULL;
          } else
          {
            if (gt_option_is_set(arguments->refextendxdropoption))
            {
              processmaxpairs = gt_generic_extend_selfmatch_xdrop_with_output;
              info_querymatch.processinfo = (void *) xdropmatchinfo;
            } else
            {
              if (gt_option_is_set(arguments->refextendgreedyoption))
              {
                processmaxpairs = gt_generic_simplegreedyselfmatchoutput;
                info_querymatch.processinfo = (void *) greedyextendmatchinfo;
              } else
              {
                processmaxpairs = gt_exact_selfmatch_with_output;
              }
            }
            processmaxpairsdata = (void *) &info_querymatch;
          }
          if (gt_callenummaxpairs(gt_str_get(arguments->indexname),
                                  arguments->seedlength,
                                  arguments->maxfreq,
                                  arguments->scanfile,
                                  processmaxpairs,
                                  processmaxpairsdata,
                                  logger,
                                  err) != 0)
          {
            haserr = true;
          }
        }
        if (!haserr)
        {
          /* only handle reverse and reverse_complement, as forward
             is handled as selfmatch-mode */
          for (mode = 1; !haserr && mode < 3; mode++)
          {
            if (flags[mode])
            {
              gt_querymatch_query_readmode_set(
                   info_querymatch.querymatchspaceptr,
                   modes[mode]);
              if (gt_callenumquerymatches(true,
                                          gt_str_get(arguments->indexname),
                                          NULL,
                                          NULL,
                                          modes[mode],
                                          arguments->seedlength,
                                          querymatchoutoptions,
                                          eqmf,
                                          eqmf_data,
                                          out_display_flag,
                                          logger,
                                          err) != 0)
              {
                haserr = true;
              }
            }
          }
        }
      }
    } else
    {
      for (mode = 0; !haserr && mode < 3; mode++)
      {
        if (flags[mode])
        {
          gt_querymatch_query_readmode_set(info_querymatch.querymatchspaceptr,
                                           modes[mode]);
          if (gt_callenumquerymatches(false,
                                  gt_str_get(arguments->indexname),
                                  arguments->query_files,
                                  arguments->query_indexname,
                                  modes[mode],
                                  arguments->seedlength,
                                  querymatchoutoptions,
                                  eqmf,
                                  eqmf_data,
                                  out_display_flag,
                                  logger,
                                  err) != 0)
          {
            haserr = true;
          }
        }
      }
    }
    gt_encseq_delete(encseq_for_desc);
    gt_querymatchoutoptions_delete(querymatchoutoptions);
    gt_querymatch_delete(info_querymatch.querymatchspaceptr);
  }
  gt_karlin_altschul_stat_delete(karlin_altschul_stat);
  gt_xdrop_matchinfo_delete(xdropmatchinfo);
  gt_greedy_extend_matchinfo_delete(greedyextendmatchinfo);
  gt_ft_trimstat_delete(trimstat);
  polishing_info_delete(pol_info);
  gt_logger_delete(logger);
  gt_querymatch_display_flag_delete(out_display_flag);
  if (repfindtimer != NULL)
  {
    char *keystring = gt_seed_extend_params_keystring(
                             gt_option_is_set(arguments->refextendgreedyoption),
                             gt_option_is_set(arguments->refextendxdropoption),
                             arguments->seedlength,
                             arguments->userdefinedleastlength,
                             arguments->minidentity,
                             arguments->maxalignedlendifference,
                             arguments->perc_mat_history,
                             arguments->extendgreedy,
                             arguments->extendxdrop,
                             arguments->xdropbelowscore);
    printf("# TIME repfind-%s",keystring);
    gt_free(keystring);
    gt_timer_show_formatted(repfindtimer," overall " GT_WD ".%02ld\n",stdout);
    gt_timer_delete(repfindtimer);
  }
  return haserr ? -1 : 0;
}

GtTool* gt_repfind(void)
{
  return gt_tool_new(gt_repfind_arguments_new,
                     gt_repfind_arguments_delete,
                     gt_repfind_option_parser_new,
                     gt_repfind_arguments_check,
                     gt_repfind_runner);
}
