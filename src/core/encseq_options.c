/*
  Copyright (c) 2010 Sascha Steinbiss <steinbiss@zbh.uni-hamburg.de>
  Copyright (c) 2010 Center for Bioinformatics, University of Hamburg

  Permission to use, copy, modify, and distribute this software for any
  purpose with or without fee is hereby granted, provided that the above
  copyright notice and this permission notice appear in all copies.

  THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES
  WITH REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF
  MERCHANTABILITY AND FITNESS. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR
  ANY SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
  WHAciOEVER RESULTING FROM LOSS OF USE, DATA OR PROFIci, WHETHER IN AN
  ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF
  OR IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.
*/

#include "core/assert_api.h"
#include "core/basename_api.h"
#include "core/encseq.h"
#include "core/encseq_options.h"
#include "core/ma.h"
#include "core/option_api.h"
#include "core/str_api.h"
#include "core/str_array.h"

typedef enum {
  GT_ENCSEQ_OPTS_ENCODE,
  GT_ENCSEQ_OPTS_LOAD,
  GT_ENCSEQ_OPTS_UNDEFINED
} GtEncseqOptionsCreationStrategy;

struct GtEncseqOptions {
  GtStr *indexname,
        *sat,
        *smap,
        *dir;
  GtOption *optiondb,
           *optionindexname,
           *optionsat,
           *optionssp,
           *optiondes,
           *optionsds,
           *optionlossless,
           *optiontis,
           *optionmd5,
           *optiondir,
           *optiondna,
           *optionplain,
           *optionprotein,
           *optionsmap,
           *optionmirrored;
  GtStrArray *db;
  bool des,
       ssp,
       sds,
       lossless,
       dna,
       tis,
       md5,
       protein,
       plain,
       mirrored,
       withdb,
       withindexname;
};

static GtEncseqOptions* gt_encseq_options_new(void)
{
  GtEncseqOptions *oi = gt_malloc(sizeof *oi);
  oi->indexname = gt_str_new();
  oi->sat = gt_str_new();
  oi->smap = gt_str_new();
  oi->dir = gt_str_new();
  oi->db = gt_str_array_new();
  oi->des = false;
  oi->ssp = false;
  oi->sds = false;
  oi->md5 = false;
  oi->lossless = false;
  oi->dna = false;
  oi->protein = false;
  oi->plain = false;
  oi->mirrored = false;
  oi->optiondb = NULL;
  oi->optionindexname = NULL;
  oi->optionsat = NULL;
  oi->optionssp = NULL;
  oi->optiondes = NULL;
  oi->optionlossless = NULL;
  oi->optionsds = NULL;
  oi->optiontis = NULL;
  oi->optionmd5 = NULL;
  oi->optiondna = NULL;
  oi->optionplain = NULL;
  oi->optionprotein = NULL;
  oi->optionsmap = NULL;
  oi->optionmirrored = NULL;
  oi->withdb = false;
  oi->withindexname = false;
  return oi;
}

static int gt_encseq_options_check(void *oip, GtError *err)
{
  int had_err = 0;
  GtEncseqOptions *oi = (GtEncseqOptions*) oip;
  gt_assert(oi != NULL);
  gt_error_check(err);

  if (oi->withdb) {
    gt_assert(oi->db != NULL);
    /* set default indexname  */
    if (oi->optiondb != NULL
          && gt_str_array_size(oi->db) == 0UL) {
      if (gt_option_is_set(oi->optiondb)) {
        gt_error_set(err, "missing argument to option -db");
        had_err = -1;
      }
    } else {
      if (oi->optionindexname != NULL
            && !gt_option_is_set(oi->optionindexname)) {
        if (gt_str_array_size(oi->db) > 1UL) {
          gt_error_set(err,"if more than one input file is given, then "
                           "option -indexname is mandatory");
          had_err = -1;
        } else {
          char *basenameptr;
          basenameptr = gt_basename(gt_str_array_get(oi->db, 0));
          gt_str_set(oi->indexname, basenameptr);
          gt_free(basenameptr);
        }
      }
    }
  }
  if (!had_err) {
    if (!oi->des && oi->sds) {
      gt_error_set(err, "option \"-sds yes\" requires \"-des yes\"");
      had_err = -1;
    }
  }
  if (!had_err) {
    if (oi->optionplain != NULL && gt_option_is_set(oi->optionplain)) {
      if (oi->optiondna != NULL && !gt_option_is_set(oi->optiondna) &&
          oi->optionprotein != NULL && !gt_option_is_set(oi->optionprotein) &&
          oi->optionsmap != NULL && !gt_option_is_set(oi->optionsmap)) {
        gt_error_set(err,"if option -plain is used, then any of the options "
                         "-dna, -protein, or -smap is mandatory");
        had_err = -1;
      }
    }
  }
  return had_err;
}

static GtEncseqOptions*
gt_encseq_options_register_generic(GtOptionParser *op,
                                   GtEncseqOptionsCreationStrategy strategy,
                                   GtStr *indexname,
                                   GtStrArray *inputdb)
{
  gt_assert(op != NULL && strategy != GT_ENCSEQ_OPTS_UNDEFINED);
  GtEncseqOptions *oi = gt_encseq_options_new();
  oi->withdb = (inputdb != NULL);
  oi->withindexname = (indexname != NULL);

  if (oi->withindexname) {
    gt_str_delete(oi->indexname);
    oi->indexname = gt_str_ref(indexname);
  }
  if (oi->withdb) {
    gt_str_array_delete(oi->db);
    oi->db = gt_str_array_ref(inputdb);
  }

  /* encoding options */
  if (strategy == GT_ENCSEQ_OPTS_ENCODE) {

    oi->optiontis = gt_option_new_bool("tis",
                                       "output transformed and encoded input "
                                       "sequence to file (deprecated, kept for "
                                       "compatibility reasons)",
                                       &oi->tis,
                                       true);
    gt_option_parser_add_option(op, oi->optiontis);
    gt_option_is_development_option(oi->optiontis);

    oi->optionssp = gt_option_new_bool("ssp",
                                       "output sequence separator positions "
                                       "to file",
                                       &oi->ssp,
                                       true);
    gt_option_parser_add_option(op, oi->optionssp);

    oi->optiondes = gt_option_new_bool("des",
                                       "output sequence descriptions to file",
                                       &oi->des,
                                       true);
    gt_option_parser_add_option(op, oi->optiondes);

    oi->optionsds = gt_option_new_bool("sds",
                                       "output sequence description separator "
                                       "positions to file",
                                       &oi->sds,
                                       true);
    gt_option_parser_add_option(op, oi->optionsds);
    gt_option_imply(oi->optionsds, oi->optiondes);

    oi->optionmd5 = gt_option_new_bool("md5",
                                       "output MD5 sums to file",
                                       &oi->md5,
                                       true);
    gt_option_parser_add_option(op, oi->optionmd5);

    oi->optionsat = gt_option_new_string("sat",
                                         "specify kind of sequence "
                                         "representation\n"
                                         "by one of the keywords direct, "
                                         "bytecompress, eqlen, bit, uchar, "
                                         "ushort, uint32",
                                         oi->sat, NULL);
    gt_option_parser_add_option(op, oi->optionsat);

    oi->optiondna = gt_option_new_bool("dna","input is DNA sequence",
                                       &oi->dna, false);
    gt_option_parser_add_option(op, oi->optiondna);

    oi->optionprotein = gt_option_new_bool("protein","input is "
                                           "protein sequence",
                                           &oi->protein, false);
    gt_option_parser_add_option(op, oi->optionprotein);

    oi->optionplain = gt_option_new_bool("plain","process as plain text",
                                         &oi->plain, false);
    gt_option_parser_add_option(op, oi->optionplain);
    gt_option_is_extended_option(oi->optionplain);

    oi->optiondb = gt_option_new_filename_array("db","specify database files",
                                                oi->db);

    oi->optionsmap = gt_option_new_string("smap",
                                          "specify file containing a "
                                          "symbol mapping",
                                          oi->smap, NULL);

    /* only add -indexname if not present already (e.g. suffixerator) */
    if (gt_option_parser_get_option(op, "indexname") == NULL) {
      oi->optionindexname = gt_option_new_string("indexname",
                                                 "specify name for index to "
                                                 "be generated",
                                                 oi->indexname, NULL);
      if (oi->withindexname) {
        gt_option_parser_add_option(op, oi->optionindexname);
      }
    }

    if (oi->withdb) {
      gt_option_parser_add_option(op, oi->optiondb);
    }

    gt_option_parser_add_option(op, oi->optionsmap);
    gt_option_exclude(oi->optionsmap, oi->optiondna);
    gt_option_exclude(oi->optionsmap, oi->optionprotein);
    gt_option_exclude(oi->optiondna, oi->optionprotein);
  }

  /* decoding options */
  if (strategy == GT_ENCSEQ_OPTS_LOAD) {
    oi->optionmirrored = gt_option_new_bool("mirrored", "virtually append the "
                                            "reverse complement of each "
                                            "sequence",
                                            &oi->mirrored,
                                            false);
    gt_option_parser_add_option(op, oi->optionmirrored);
  }

  /* only add -lossless if not present already (e.g. suffixerator) */
  if (gt_option_parser_get_option(op, "lossless") == NULL) {
    oi->optionlossless = gt_option_new_bool("lossless",
                                            "allow lossless original sequence "
                                            "retrieval",
                                            &oi->lossless,
                                            false);
    gt_option_parser_add_option(op, oi->optionlossless);
  }

  gt_option_parser_register_hook(op, gt_encseq_options_check, oi);

  return oi;
}

GtEncseqOptions* gt_encseq_options_register_encoding(GtOptionParser *op,
                                                     GtStr *indexname,
                                                     GtStrArray *indb)
{
  gt_assert(op != NULL);
  return gt_encseq_options_register_generic(op,
                                            GT_ENCSEQ_OPTS_ENCODE,
                                            indexname,
                                            indb);
}

GtEncseqOptions* gt_encseq_options_register_loading(GtOptionParser *op,
                                                    GtStr *indexname)
{
  gt_assert(op != NULL);
  return gt_encseq_options_register_generic(op,
                                            GT_ENCSEQ_OPTS_LOAD,
                                            indexname,
                                            NULL);
}

void gt_encseq_options_add_readmode_option(GtOptionParser *op, GtStr *readmode)
{
  GtOption *optiondir;
  gt_assert(op && readmode);
  optiondir = gt_option_new_string("dir", "specify reading direction "
                                          "(fwd, cpl, rev, rcl)",
                                   readmode, "fwd");
  gt_option_parser_add_option(op, optiondir);
}

#ifndef GT_ENCSEQ_OPTS_GETTER_DEFS_DEFINED

#define GT_ENCSEQ_OPTS_GETTER_DEF_OPT(VARNAME) \
GtOption* gt_encseq_options_##VARNAME##_option(GtEncseqOptions *i) \
{ \
  gt_assert(i != NULL); \
  return i->option##VARNAME; \
}

#define GT_ENCSEQ_OPTS_GETTER_DEF_VAL(VARNAME, TYPE) \
TYPE gt_encseq_options_##VARNAME##_value(GtEncseqOptions *i) \
{ \
  gt_assert(i != NULL); \
  return i->VARNAME; \
}

#define GT_ENCSEQ_OPTS_GETTER_DEF(VARNAME,TYPE) \
GT_ENCSEQ_OPTS_GETTER_DEF_OPT(VARNAME) \
GT_ENCSEQ_OPTS_GETTER_DEF_VAL(VARNAME, TYPE)

#define GT_ENCSEQ_OPTS_GETTER_DEFS_DEFINED
#endif

GT_ENCSEQ_OPTS_GETTER_DEF(db, GtStrArray*);
GT_ENCSEQ_OPTS_GETTER_DEF(des, bool);
GT_ENCSEQ_OPTS_GETTER_DEF(dna, bool);
GT_ENCSEQ_OPTS_GETTER_DEF(indexname, GtStr*);
GT_ENCSEQ_OPTS_GETTER_DEF(lossless, bool);
GT_ENCSEQ_OPTS_GETTER_DEF(md5, bool);
GT_ENCSEQ_OPTS_GETTER_DEF(mirrored, bool);
GT_ENCSEQ_OPTS_GETTER_DEF(plain, bool);
GT_ENCSEQ_OPTS_GETTER_DEF(protein, bool);
GT_ENCSEQ_OPTS_GETTER_DEF(sat, GtStr*);
GT_ENCSEQ_OPTS_GETTER_DEF(sds, bool);
GT_ENCSEQ_OPTS_GETTER_DEF(smap, GtStr*);
GT_ENCSEQ_OPTS_GETTER_DEF(ssp, bool);
GT_ENCSEQ_OPTS_GETTER_DEF(tis, bool);
GT_ENCSEQ_OPTS_GETTER_DEF_OPT(dir);

void gt_encseq_options_delete(GtEncseqOptions *oi)
{
  if (oi == NULL) return;
  gt_str_delete(oi->indexname);
  gt_str_delete(oi->smap);
  gt_str_delete(oi->sat);
  gt_str_delete(oi->dir);
  gt_str_array_delete(oi->db);
  if (!oi->withdb) {
    gt_option_delete(oi->optiondb);
  }
  if (!oi->withindexname) {
    gt_option_delete(oi->optionindexname);
  }
  gt_free(oi);
}
