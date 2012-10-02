/*
  Copyright (c) 2007      Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
  Copyright (c)      2010 Sascha Steinbiss <steinbiss@zbh.uni-hamburg.de>
  Copyright (c)      2012 Dirk Willrodt <willrodt@zbh.uni-hamburg.de>
  Copyright (c) 2007-2012 Center for Bioinformatics, University of Hamburg

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

#include <inttypes.h>
#include <stdbool.h>
#include <string.h>

#include "core/basename_api.h"
#include "core/encseq_options.h"
#include "core/error.h"
#include "core/grep_api.h"
#include "core/logger.h"
#include "core/ma_api.h"
#include "core/option_api.h"
#include "core/spacecalc.h"
#include "core/str.h"
#include "core/str_array.h"
#include "core/unused_api.h"
#include "core/versionfunc.h"
#ifndef S_SPLINT_S /* splint reports too many errors for the following */
#include "match/eis-bwtseq-param.h"
#endif
#include "match/index_options.h"

typedef enum {
  GT_INDEX_OPTIONS_UNDEFINED,
  GT_INDEX_OPTIONS_ESA,
  GT_INDEX_OPTIONS_PACKED
} GtIndexOptionsIndexType;

struct GtIndexOptions
{
  unsigned int numofparts,
               prefixlength;
  unsigned long maximumspace;
  GtStrArray *algbounds;
  GtReadmode readmode;
  bool outsuftab,
       outlcptab,
       outbwttab,
       outbcktab,
       outkystab,
       outkyssort,
       lcpdist,
       swallow_tail;
  GtStr *kysargumentstring,
        *indexname,
        *dir,
        *memlimit;
  Sfxstrategy sfxstrategy;
  GtEncseqOptions *encopts;
  GtIndexOptionsIndexType type;
  GtOption *option,
           *optionoutsuftab,
           *optionoutlcptab,
           *optionoutbwttab,
           *optionoutbcktab,
           *optionprefixlength,
           *optioncmpcharbychar,
           *optionnoshortreadsort,
           *optionspmopt,
           *optionstorespecialcodes,
           *optionmaxwidthrealmedian,
           *optionalgbounds,
           *optionparts,
           *optionmemlimit,
           *optiondifferencecover,
           *optionuserdefinedsortmaxdepth,
           *optionkys;
#ifndef S_SPLINT_S /* splint reports too many errors for the following and so
                      we exclude it */
  struct bwtOptions bwtIdxParams;
#endif
};

static GtIndexOptions* gt_index_options_new(void)
{
  GtIndexOptions *oi = gt_malloc(sizeof *oi);
  oi->algbounds = gt_str_array_new();
  oi->dir = gt_str_new_cstr("fwd");
  oi->indexname = NULL;
  oi->kysargumentstring = gt_str_new();
  oi->lcpdist = false;
  oi->maximumspace = 0UL; /* in bytes */
  oi->memlimit = gt_str_new();
  oi->numofparts = 1U;
  oi->option = NULL;
  oi->optionalgbounds = NULL;
  oi->optioncmpcharbychar = NULL;
  oi->optiondifferencecover = NULL;
  oi->optionmaxwidthrealmedian = NULL;
  oi->optionmemlimit = NULL;
  oi->optionoutbcktab = NULL;
  oi->optionoutbwttab = NULL;
  oi->optionoutlcptab = NULL;
  oi->optionoutsuftab = NULL;
  oi->optionparts = NULL;
  oi->optionprefixlength = NULL;
  oi->optionspmopt = NULL;
  oi->optionstorespecialcodes = NULL;
  oi->outbcktab = false;
  oi->outbwttab = false;
  oi->outkyssort = false;
  oi->outkystab = false;
  oi->outlcptab = false;
  oi->outsuftab = false; /* only defined for GT_INDEX_OPTIONS_ESA */
  oi->prefixlength = GT_PREFIXLENGTH_AUTOMATIC;
  oi->swallow_tail = false;
  oi->type = GT_INDEX_OPTIONS_UNDEFINED;
  return oi;
}

#define GT_IDXOPTS_READMAXBOUND(COMP,IDX)\
        if (!haserr)\
        {\
          arg = gt_str_array_get(algbounds, IDX);\
          if (sscanf(arg,"%ld", &readint) != 1 || readint <= 0)\
          {\
            gt_error_set(err,"option -algbds: all arguments must be positive " \
                             "numbers");\
            haserr = true;\
          }\
          sfxstrategy->COMP = (unsigned long) readint;\
        }

int gt_parse_algbounds(Sfxstrategy *sfxstrategy,
                       const GtStrArray *algbounds,
                       GtError *err)
{
  bool haserr = false;
  const char *arg;
  long readint;

  if (gt_str_array_size(algbounds) != 3UL)
  {
    gt_error_set(err,"option -algbds must have exactly 3 arguments");
    haserr = true;
  }
  GT_IDXOPTS_READMAXBOUND(maxinsertionsort, 0);
  GT_IDXOPTS_READMAXBOUND(maxbltriesort, 1UL);
  if (sfxstrategy->maxinsertionsort > sfxstrategy->maxbltriesort)
  {
    gt_error_set(err,"first argument of option -algbds must not be larger "
                     "than second argument");
    haserr = true;
  }
  GT_IDXOPTS_READMAXBOUND(maxcountingsort, 2UL);
  if (sfxstrategy->maxbltriesort > sfxstrategy->maxcountingsort)
  {
    gt_error_set(err,"second argument of option -algbds must not be larger "
                     "than third argument");
    haserr = true;
  }
  return haserr ? -1 : 0;
}

static int gt_index_options_check_set_create_opts(void *oip, GtError *err)
{
  int had_err = 0;
  GtIndexOptions *oi = (GtIndexOptions*) oip;
  gt_assert(oi != NULL && oi->type != GT_INDEX_OPTIONS_UNDEFINED);
  gt_error_check(err);
  /* readmode needs to be initialized to fwd */
  oi->readmode = (GtReadmode) 0;
  if (!had_err)
  {
    if (gt_option_is_set(oi->optionalgbounds))
    {
      if (gt_parse_algbounds(&oi->sfxstrategy,oi->algbounds,err) != 0)
      {
        had_err = -1;
      }
    } else
    {
      oi->sfxstrategy.maxinsertionsort = MAXINSERTIONSORTDEFAULT;
      oi->sfxstrategy.maxbltriesort = MAXBLTRIESORTDEFAULT;
      oi->sfxstrategy.maxcountingsort = MAXCOUNTINGSORTDEFAULT;
    }
  }
  if (!had_err && oi->sfxstrategy.kmerswithencseqreader &&
      oi->sfxstrategy.spmopt_minlength > 0)
  {
    gt_error_set(err,"options -spmopt  and -kmerswithencseqreader are "
                     "not compatible");
    had_err = -1;
  }
  if (!had_err
        && oi->optionmemlimit != NULL
        && gt_option_is_set(oi->optionmemlimit))
  {
    had_err = gt_option_parse_spacespec(&oi->maximumspace,"memlimit",
                                        oi->memlimit,err);
  }
  if (!had_err)
  {
    if (oi->sfxstrategy.maxinsertionsort > oi->sfxstrategy.maxbltriesort)
    {
      gt_error_set(err,"first argument of option -algbds must not be larger "
                       "than second argument");
      had_err = -1;
    } else
    {
      if (oi->sfxstrategy.maxbltriesort > oi->sfxstrategy.maxcountingsort)
      {
        gt_error_set(err,"second argument of option -algbds must not be larger "
                         "than third argument");
        had_err = -1;
      }
    }
  }
  return had_err;
}

static int gt_index_options_check_set_out_opts(void *oip, GtError *err)
{
  int had_err = 0;
  GtIndexOptions *oi = (GtIndexOptions*) oip;
  gt_assert(oi != NULL && oi->type != GT_INDEX_OPTIONS_UNDEFINED);
  gt_error_check(err);

  if (!had_err) {
    int retval;
    retval = gt_readmode_parse(gt_str_get(oi->dir), err);
    if (retval < 0) {
      had_err = -1;
    } else {
      oi->readmode = (GtReadmode) retval;
      if (oi->type == GT_INDEX_OPTIONS_PACKED
            && (oi->readmode == GT_READMODE_COMPL
                  || oi->readmode == GT_READMODE_REVCOMPL)) {
        gt_error_set(err,"construction of packed index not possible for "
                         "complemented and for reverse complemented sequences");
        had_err = -1;
      }
    }
  }
  if (!had_err && oi->type == GT_INDEX_OPTIONS_PACKED) {
#ifndef S_SPLINT_S
    gt_computePackedIndexDefaults(&oi->bwtIdxParams, BWTBaseFeatures);
#endif
  }
  if (!had_err && gt_option_is_set(oi->optionkys)) {
    oi->outkystab = true;
    if (strcmp(gt_str_get(oi->kysargumentstring), "sort") == 0) {
      oi->outkyssort = true;
    } else {
      if (strcmp(gt_str_get(oi->kysargumentstring),"nosort") != 0) {
        gt_error_set(err,"illegal argument to option -kys: either use no "
                         "argument or argument \"sort\"");
        had_err = -1;
      }
    }
  }

  return had_err;
}

static GtIndexOptions*
gt_index_options_register_generic_output(GtOptionParser *op,
                                         GtIndexOptions *idxo,
                                         GtStr *indexname,
                                         GtEncseqOptions *encopts)
{
  gt_assert(idxo != NULL);
  gt_assert(op != NULL && idxo->type != GT_INDEX_OPTIONS_UNDEFINED &&
            encopts != NULL);
  idxo->encopts = encopts;
  idxo->indexname = indexname != NULL ? gt_str_ref(indexname) : NULL;
  idxo->optionkys = gt_option_new_string("kys",
                                   "output/sort according to keys of the form "
                                   "|key| in fasta header",
                                   idxo->kysargumentstring,
                                   "nosort");
  gt_option_argument_is_optional(idxo->optionkys);
  gt_option_imply(idxo->optionkys, gt_encseq_options_sds_option(idxo->encopts));
  gt_option_parser_add_option(op, idxo->optionkys);

  gt_encseq_options_add_readmode_option(op, idxo->dir);

  if (idxo->type == GT_INDEX_OPTIONS_ESA)
  {
    idxo->optionoutsuftab = gt_option_new_bool("suf",
                                   "output suffix array (suftab) to file",
                                   &idxo->outsuftab,
                                   false);
    gt_option_parser_add_option(op, idxo->optionoutsuftab);

    idxo->optionoutlcptab = gt_option_new_bool("lcp",
                                   "output lcp table (lcptab) to file",
                                   &idxo->outlcptab,
                                   false);
    gt_option_parser_add_option(op, idxo->optionoutlcptab);

    idxo->option = gt_option_new_bool("lcpdist",
                              "output distributions of values in lcptab",
                              &idxo->lcpdist,
                              false);
    gt_option_is_extended_option(idxo->option);
    gt_option_imply(idxo->option, idxo->optionoutlcptab);
    gt_option_parser_add_option(op, idxo->option);

    idxo->option = gt_option_new_bool("swallow-tail",
                              "swallow the tail of the suffix array and lcptab",
                              &idxo->swallow_tail,
                              false);
    gt_option_is_development_option(idxo->option);
    gt_option_parser_add_option(op, idxo->option);

    idxo->optionoutbwttab = gt_option_new_bool("bwt",
                                   "output Burrows-Wheeler Transformation "
                                   "(bwttab) to file",
                                   &idxo->outbwttab,
                                   false);
    gt_option_exclude(idxo->optionspmopt, idxo->optionoutbwttab);
    gt_option_parser_add_option(op, idxo->optionoutbwttab);

    idxo->optionoutbcktab = gt_option_new_bool("bck",
                                "output bucket table to file",
                                &idxo->outbcktab,
                                false);
    gt_option_parser_add_option(op, idxo->optionoutbcktab);
  } else {
    idxo->optionoutsuftab
      = idxo->optionoutlcptab = idxo->optionoutbwttab = NULL;
    idxo->sfxstrategy.spmopt_minlength = 0;
#ifndef S_SPLINT_S
    gt_registerPackedIndexOptions(op,
                                  &idxo->bwtIdxParams,
                                  BWTDEFOPT_CONSTRUCTION,
                                  idxo->indexname);
#endif
  }

  gt_option_parser_register_hook(op, gt_index_options_check_set_out_opts, idxo);

  return idxo;
}

static GtIndexOptions* gt_index_options_register_generic_create(
                                                      GtOptionParser *op,
                                                      GtIndexOptionsIndexType t)
{
  GtIndexOptions *idxo;
  gt_assert(op != NULL && t != GT_INDEX_OPTIONS_UNDEFINED);
  idxo = gt_index_options_new();
  idxo->type = t;
  idxo->optionprefixlength = gt_option_new_uint_min("pl",
                                    "specify prefix length for bucket sort\n"
                                    "recommendation: use without argument;\n"
                                    "then a reasonable prefix length is "
                                    "automatically determined.",
                                    &idxo->prefixlength,
                                    GT_PREFIXLENGTH_AUTOMATIC,
                                    1U);
  gt_option_argument_is_optional(idxo->optionprefixlength);
  gt_option_parser_add_option(op, idxo->optionprefixlength);

  idxo->optionuserdefinedsortmaxdepth
    = gt_option_new_uint_min("sortmaxdepth","sort only up to the given depth.",
                             &idxo->sfxstrategy.userdefinedsortmaxdepth,
                             0,
                             1U);
  gt_option_parser_add_option(op, idxo->optionuserdefinedsortmaxdepth);
  gt_option_is_development_option(idxo->optionuserdefinedsortmaxdepth);

  idxo->optiondifferencecover = gt_option_new_uint_min("dc",
                                    "specify difference cover value",
                                    &idxo->sfxstrategy.differencecover,
                                    0,
                                    4U);
  gt_option_parser_add_option(op, idxo->optiondifferencecover);
  gt_option_exclude(idxo->optionuserdefinedsortmaxdepth,
                    idxo->optiondifferencecover);

  idxo->optioncmpcharbychar = gt_option_new_bool("cmpcharbychar",
                                           "compare suffixes character "
                                           "by character",
                                           &idxo->sfxstrategy.cmpcharbychar,
                                           false);
  gt_option_is_development_option(idxo->optioncmpcharbychar);
  gt_option_parser_add_option(op, idxo->optioncmpcharbychar);

  idxo->optionnoshortreadsort = gt_option_new_bool("noshortreadsort",
                                           "do not use short read sort",
                                           &idxo->sfxstrategy.noshortreadsort,
                                           false);
  gt_option_is_development_option(idxo->optionnoshortreadsort);
  gt_option_parser_add_option(op, idxo->optionnoshortreadsort);

  idxo->optionmaxwidthrealmedian = gt_option_new_ulong("maxwidthrealmedian",
                                                 "compute real median for "
                                                 "intervals of at most the "
                                                 "given widthprefixes",
                                                 &idxo->sfxstrategy.
                                                      maxwidthrealmedian,
                                                 1UL);
  gt_option_is_development_option(idxo->optionmaxwidthrealmedian);
  gt_option_parser_add_option(op, idxo->optionmaxwidthrealmedian);

  idxo->optionalgbounds
    = gt_option_new_string_array("algbds",
                                 "length boundaries for the different "
                                 "algorithms to sort buckets of suffixes\n"
                                 "first number: maxbound for insertion sort\n"
                                 "second number: maxbound for blindtrie sort\n"
                                 "third number: maxbound for counting sort",
                                 idxo->algbounds);
  gt_option_is_development_option(idxo->optionalgbounds);
  gt_option_parser_add_option(op, idxo->optionalgbounds);

  idxo->optionstorespecialcodes
    = gt_option_new_bool("storespecialcodes",
                         "store special codes (this may speed up the program)",
                         &idxo->sfxstrategy.storespecialcodes,false);

  gt_option_is_development_option(idxo->optionstorespecialcodes);
  gt_option_parser_add_option(op, idxo->optionstorespecialcodes);

  idxo->optionparts
    = gt_option_new_uint_max("parts",
                             "specify number of parts in which the index "
                             "construction is performed",
                             &idxo->numofparts, 1U, (unsigned) ((1 << 22) - 1));
  gt_option_is_development_option(idxo->optionparts);
  gt_option_parser_add_option(op, idxo->optionparts);

  if (idxo->type == GT_INDEX_OPTIONS_ESA)
  {
    idxo->optionspmopt = gt_option_new_uint_min("spmopt",
                                            "optimize esa-construction for "
                                            "suffix-prefix matching",
                                            &idxo->sfxstrategy.spmopt_minlength,
                                            0,1U);
    gt_option_parser_add_option(op, idxo->optionspmopt);
    gt_option_exclude(idxo->optionspmopt, idxo->optiondifferencecover);
    idxo->optionmemlimit = gt_option_new_string("memlimit",
                           "specify maximal amount of memory to be used during "
                           "index construction (in bytes, the keywords 'MB' "
                           "and 'GB' are allowed)",
                           idxo->memlimit, NULL);
    gt_option_parser_add_option(op, idxo->optionmemlimit);
    gt_option_exclude(idxo->optionmemlimit, idxo->optionparts);
  }

  idxo->option = gt_option_new_bool("iterscan",
                                   "use iteratorbased-kmer scanning",
                                   &idxo->sfxstrategy.iteratorbasedkmerscanning,
                                   false);
  gt_option_is_development_option(idxo->option);
  gt_option_parser_add_option(op, idxo->option);

  idxo->option = gt_option_new_bool("samplewithprefixlengthnull",
                                  "sort sample with prefixlength=0",
                                  &idxo->sfxstrategy.samplewithprefixlengthnull,
                                  false);
  gt_option_is_development_option(idxo->option);
  gt_option_parser_add_option(op, idxo->option);

  idxo->option = gt_option_new_bool("suftabuint",
                                    "use uint32_t for suftab",
                                    &idxo->sfxstrategy.suftabuint,
                                    false);
  gt_option_is_development_option(idxo->option);
  gt_option_parser_add_option(op, idxo->option);

  idxo->option = gt_option_new_bool("onlybucketinsertion",
                                    "perform only bucket insertion",
                                    &idxo->sfxstrategy.onlybucketinsertion,
                                    false);
  gt_option_is_development_option(idxo->option);
  gt_option_parser_add_option(op, idxo->option);

  idxo->option = gt_option_new_bool("kmerswithencseqreader",
                               "always perform kmerscanning with encseq-reader",
                               &idxo->sfxstrategy.kmerswithencseqreader,
                               false);
  gt_option_is_development_option(idxo->option);
  gt_option_parser_add_option(op, idxo->option);

  idxo->option = gt_option_new_bool("dccheck",
                               "check intermediate results in difference cover",
                               &idxo->sfxstrategy.dccheck,
                               false);
  gt_option_is_development_option(idxo->option);
  gt_option_parser_add_option(op, idxo->option);

  idxo->option = gt_option_new_bool("withradixsort",
                                    "use radixsort to sort the buckets",
                                    &idxo->sfxstrategy.withradixsort,
                                    false);
  gt_option_is_development_option(idxo->option);
  gt_option_parser_add_option(op, idxo->option);

  gt_option_parser_register_hook(op, gt_index_options_check_set_create_opts,
                                 idxo);
  return idxo;
}

GtIndexOptions* gt_index_options_register_esa(GtOptionParser *op,
                                              GtEncseqOptions *encopts)
{
  GtIndexOptions *idxo;
  gt_assert(op != NULL);

  idxo = gt_index_options_register_generic_create(op,
                                                  GT_INDEX_OPTIONS_ESA);
  return gt_index_options_register_generic_output(op, idxo, NULL, encopts);
}

GtIndexOptions* gt_index_options_register_esa_noout(GtOptionParser *op)
{
  gt_assert(op != NULL);

  return gt_index_options_register_generic_create(op,
                                                  GT_INDEX_OPTIONS_ESA);
}

GtIndexOptions* gt_index_options_register_packedidx(GtOptionParser *op,
                                                    GtStr *indexname,
                                                    GtEncseqOptions *encopts)
{
  GtIndexOptions *idxo;
  gt_assert(op != NULL);
  idxo = gt_index_options_register_generic_create(op,
                                                  GT_INDEX_OPTIONS_PACKED);
  return gt_index_options_register_generic_output(op, idxo, indexname, encopts);
}

void gt_index_options_delete(GtIndexOptions *oi)
{
  if (oi == NULL) return;
  gt_str_delete(oi->kysargumentstring);
  gt_str_delete(oi->indexname);
  gt_str_delete(oi->dir);
  gt_str_delete(oi->memlimit);
  gt_str_array_delete(oi->algbounds);
  gt_free(oi);
}

/* XXX: clean this up */
#ifndef GT_INDEX_OPTS_GETTER_DEFS_DEFINED
#define GT_INDEX_OPTS_GETTER_DEF_OPT(VARNAME) \
GtOption* gt_index_options_##VARNAME##_option(GtIndexOptions *i) \
{ \
  gt_assert(i != NULL); \
  return i->option##VARNAME; \
}
#define GT_INDEX_OPTS_GETTER_DEF_VAL(VARNAME, TYPE) \
TYPE gt_index_options_##VARNAME##_value(GtIndexOptions *i) \
{ \
  gt_assert(i != NULL); \
  return i->VARNAME; \
}
#define GT_INDEX_OPTS_GETTER_DEF(VARNAME,TYPE) \
GT_INDEX_OPTS_GETTER_DEF_OPT(VARNAME) \
GT_INDEX_OPTS_GETTER_DEF_VAL(VARNAME, TYPE)
#define GT_INDEX_OPTS_GETTER_DEFS_DEFINED
#endif

/* these are available as options and values */
GT_INDEX_OPTS_GETTER_DEF(algbounds, GtStrArray*);
GT_INDEX_OPTS_GETTER_DEF(outbcktab, bool);
GT_INDEX_OPTS_GETTER_DEF(outbwttab, bool);
GT_INDEX_OPTS_GETTER_DEF(outlcptab, bool);
GT_INDEX_OPTS_GETTER_DEF(outsuftab, bool);
GT_INDEX_OPTS_GETTER_DEF(prefixlength, unsigned int);
GT_INDEX_OPTS_GETTER_DEF_OPT(spmopt);
/* these are available as values only, set _after_ option processing */
GT_INDEX_OPTS_GETTER_DEF_VAL(lcpdist, bool);
GT_INDEX_OPTS_GETTER_DEF_VAL(maximumspace, unsigned long);
GT_INDEX_OPTS_GETTER_DEF_VAL(numofparts, unsigned int);
GT_INDEX_OPTS_GETTER_DEF_VAL(outkyssort, bool);
GT_INDEX_OPTS_GETTER_DEF_VAL(outkystab, bool);
GT_INDEX_OPTS_GETTER_DEF_VAL(readmode, GtReadmode);
GT_INDEX_OPTS_GETTER_DEF_VAL(sfxstrategy, Sfxstrategy);
GT_INDEX_OPTS_GETTER_DEF_VAL(swallow_tail, bool);
#ifndef S_SPLINT_S
GT_INDEX_OPTS_GETTER_DEF_VAL(bwtIdxParams, struct bwtOptions);
#endif
