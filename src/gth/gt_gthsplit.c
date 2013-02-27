/*
  Copyright (c) 2005-2011 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2005-2008 Center for Bioinformatics, University of Hamburg

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

#include <math.h>
#include "core/cstr_api.h"
#include "core/output_file_api.h"
#include "core/unused_api.h"
#include "gth/gthxml.h"
#include "gth/gthverbosefunc.h"
#include "gth/intermediate.h"
#include "gth/plugins.h"
#include "gth/xml_inter_sa_visitor.h"
#include "gth/gt_gthsplit.h"

#define RANGE_OPT_CSTR  "range"
#define DEFAULT_RANGE   5

typedef enum {
  ALIGNMENTSCORE_SPLIT,
  COVERAGE_SPLIT,
  UNDEF_SPLIT
} Splitmode;

typedef struct {
  Splitmode splitmode;
  unsigned int range;
  GtFileMode file_mode;
  GtStrArray *consensusfiles;
  bool force;
  GthSAFilter *sa_filter; /* owned */
  GthShowVerbose showverbose;
} Gthsplitinfo;

typedef struct {
  Gthsplitinfo *gthsplitinfo;
  char *current_outputfilename;
  GtFile **subset_files;
  GtStr **subset_filenames;
  unsigned long *subset_file_sa_counter,
                num_of_subset_files;
  GthSAFilter *sa_filter; /* reference only */
  GthSAVisitor *sa_visitor;
} Store_in_subset_file_data;

static void close_output_files(Store_in_subset_file_data
                               *store_in_subset_file_data)
{
  unsigned long i;
  GtStr *buf;

  buf = gt_str_new();
  for (i = 0; i < store_in_subset_file_data->num_of_subset_files; i++) {
    if (store_in_subset_file_data->subset_files[i]) {
      if (store_in_subset_file_data->gthsplitinfo->showverbose) {
        gt_str_reset(buf);
        gt_str_append_cstr(buf, "split file created: ");
        gt_str_append_str(buf, store_in_subset_file_data->subset_filenames[i]);
        gt_str_append_cstr(buf, " (size=");
        gt_str_append_ulong(buf,
                          store_in_subset_file_data->subset_file_sa_counter[i]);
        gt_str_append_cstr(buf, ")");
        store_in_subset_file_data->gthsplitinfo->showverbose(gt_str_get(buf));
      }
      gt_assert(store_in_subset_file_data->subset_filenames[i]);
      /* put XML trailer in file before closing it */
      gth_xml_show_trailer(true, store_in_subset_file_data->subset_files[i]);
      gt_file_delete(store_in_subset_file_data->subset_files[i]);
      gt_str_delete(store_in_subset_file_data->subset_filenames[i]);
      store_in_subset_file_data->subset_files[i]           = NULL;
      store_in_subset_file_data->subset_file_sa_counter[i] = 0;
    }
  }
  gt_str_delete(buf);
}

static int store_in_subset_file(void *data, GthSA *sa,
                                const char *outputfilename, GtError *err)
{
  Store_in_subset_file_data *store_in_subset_file_data =
    (Store_in_subset_file_data*) data;
  double split_determing_percentage = 0.0;
  unsigned long filenum;
  char filenamesuffix[4];
  int had_err = 0;

  gt_error_check(err);

  /* filter before we do any further processing */
  if (gth_sa_filter_filter_sa(store_in_subset_file_data->sa_filter, sa)) {
    /* and free it afterwards */
    gth_sa_delete(sa);
    /* discard */
    return 0;
  }

  /* check whether we got a new output file to process */
  if (!store_in_subset_file_data->current_outputfilename) {
    store_in_subset_file_data->current_outputfilename =
      gt_cstr_dup(outputfilename);
  }
  else if (strcmp(store_in_subset_file_data->current_outputfilename,
                  outputfilename)) {
    /* close current output files */
    close_output_files(store_in_subset_file_data);
    gt_free(store_in_subset_file_data->current_outputfilename);
 }

  /* determine in which file the current sa needs to be put */
  switch (store_in_subset_file_data->gthsplitinfo->splitmode) {
    case ALIGNMENTSCORE_SPLIT:
      split_determing_percentage = gth_sa_score(sa);
      strcpy(filenamesuffix, "scr");
      break;
    case COVERAGE_SPLIT:
      split_determing_percentage = gth_sa_coverage(sa);
      strcpy(filenamesuffix, "cov");
      break;
    default: gt_assert(0);
  }
  gt_assert(split_determing_percentage >= 0.0);
  /* XXX: change into an assertion when coverage problem is fixed */
  if (split_determing_percentage > 1.0)
    split_determing_percentage = 1.0;

  if (split_determing_percentage == 1.0)
    filenum = store_in_subset_file_data->num_of_subset_files - 1;
  else {
    filenum =  floor(split_determing_percentage * 100.0 /
                           store_in_subset_file_data->gthsplitinfo->range);
  }
  gt_assert(filenum < store_in_subset_file_data->num_of_subset_files);

  /* make sure the file exists and is open */
  if (!store_in_subset_file_data->subset_files[filenum]) {
    gt_assert(store_in_subset_file_data->subset_filenames[filenum] == NULL);
    store_in_subset_file_data->subset_filenames[filenum] = gt_str_new();
    gt_str_append_cstr_nt(store_in_subset_file_data->subset_filenames[filenum],
                          outputfilename,
                          gt_file_basename_length(outputfilename));
    gt_str_append_char(store_in_subset_file_data->subset_filenames[filenum],
                       '.');
    gt_str_append_cstr(store_in_subset_file_data->subset_filenames[filenum],
                       filenamesuffix);
    gt_str_append_ulong(store_in_subset_file_data->subset_filenames[filenum],
                        filenum *
                        store_in_subset_file_data->gthsplitinfo->range);
    gt_str_append_char(store_in_subset_file_data->subset_filenames[filenum],
                       '-');
    gt_str_append_ulong(store_in_subset_file_data->subset_filenames[filenum],
                     (filenum + 1) *
                     store_in_subset_file_data->gthsplitinfo->range);
    gt_str_append_cstr(store_in_subset_file_data->subset_filenames[filenum],
                       gt_file_mode_suffix(store_in_subset_file_data
                                           ->gthsplitinfo->file_mode));

    /* if not disabled by -force, check if file already exists */
    if (!store_in_subset_file_data->gthsplitinfo->force) {
      store_in_subset_file_data->subset_files[filenum] =
        gt_file_open(store_in_subset_file_data->gthsplitinfo->file_mode,
                     gt_str_get(store_in_subset_file_data
                                ->subset_filenames[filenum]), "r", NULL);
      if (store_in_subset_file_data->subset_files[filenum]) {
        gt_error_set(err, "file \"%s\" exists already. use option -%s to "
                     "overwrite", gt_str_get(store_in_subset_file_data
                                             ->subset_filenames[filenum]),
                     GT_FORCE_OPT_CSTR);
        had_err = -1;
      }
    }
    if (!had_err) {
      /* open split file for writing */
      store_in_subset_file_data->subset_files[filenum] =
          gt_file_xopen_file_mode(store_in_subset_file_data->gthsplitinfo
                                  ->file_mode,
                                  gt_str_get(store_in_subset_file_data
                                             ->subset_filenames[filenum]), "w");
      /* store XML header in file */
      gth_xml_show_leader(true,
                          store_in_subset_file_data->subset_files[filenum]);
    }
  }

  /* put it there */
  if (!had_err) {
    gth_xml_inter_sa_visitor_set_outfp(store_in_subset_file_data->sa_visitor,
                                       store_in_subset_file_data
                                       ->subset_files[filenum]);
    gth_sa_visitor_visit_sa(store_in_subset_file_data->sa_visitor, sa);
  }

  /* adjust counter */
  if (!had_err)
    store_in_subset_file_data->subset_file_sa_counter[filenum]++;

  /* and free it afterwards */
  gth_sa_delete(sa);

  return had_err;
}

static void initGthsplitinfo(Gthsplitinfo *gthsplitinfo)
{
  gthsplitinfo->splitmode      = UNDEF_SPLIT;
  gthsplitinfo->range          = DEFAULT_RANGE;
  gthsplitinfo->file_mode      = GT_FILE_MODE_UNCOMPRESSED;
  gthsplitinfo->showverbose    = NULL;
  gthsplitinfo->force          = false;
  gthsplitinfo->sa_filter      = gth_sa_filter_new();
  gthsplitinfo->consensusfiles = gt_str_array_new();
}

static void freeGthsplitinfo(Gthsplitinfo *gthsplitinfo)
{
  if (!gthsplitinfo) return;
  gt_str_array_delete(gthsplitinfo->consensusfiles);
  gth_sa_filter_delete(gthsplitinfo->sa_filter);
}

static GtOPrval gthsplit_parse_options(int *parsed_args,
                                       Gthsplitinfo *gthsplitinfo,
                                       int argc, const char **argv,
                                       const GthPlugins *plugins, GtError *err)
{
  GtOptionParser *op;
  GtOption *optalignmentscore, *optcoverage, *optrange, *optverbose, *optgzip,
           *optbzip2, *optforce;
  bool alignmentscore, coverage, verbose, gzip, bzip2;
  GtOPrval oprval;

  gt_error_check(err);

  op = gt_option_parser_new("-alignmentscore | -coverage [option ...] "
                            "[file ...]", "Split GenomeThreader output files "
                            "containing intermediate results.");

  /* specify all options with a corresponding help-text */
  optalignmentscore = gt_option_new_bool("alignmentscore", "split according to "
                                      "the overall alignment score (scr)",
                                      &alignmentscore, false);
  gt_option_parser_add_option(op, optalignmentscore);

  optcoverage = gt_option_new_bool("coverage", "split according to coverage "
                                   "(cov)", &coverage, false);
  gt_option_parser_add_option(op, optcoverage);

  optrange = gt_option_new_uint_max(RANGE_OPT_CSTR, "set the percentage range "
                                 "used to create the sets",
                                 &gthsplitinfo->range, DEFAULT_RANGE, 100);
  gt_option_parser_add_option(op, optrange);

  /* add sa_filter options */
  gth_sa_filter_register_options(op, gthsplitinfo->sa_filter, false);

  /* -v */
  optverbose = gt_option_new_verbose(&verbose);
  gt_option_parser_add_option(op, optverbose);

  optgzip = gt_option_new_bool("gzip", "write gzip compressed output file(s)",
                               &gzip, false);
  gt_option_parser_add_option(op, optgzip);

  optbzip2 = gt_option_new_bool("bzip2", "write bzip2 compressed output "
                                "file(s)", &bzip2, false);
  gt_option_parser_add_option(op, optbzip2);

  optforce = gt_option_new_bool(GT_FORCE_OPT_CSTR,"force writing to split "
                                "files", &gthsplitinfo->force, false);
  gt_option_parser_add_option(op, optforce);

  gt_option_exclude(optalignmentscore, optcoverage);
  gt_option_exclude(optgzip, optbzip2);
  gt_option_is_mandatory_either(optalignmentscore, optcoverage);

  gt_option_parser_set_mail_address(op, "<gremme@zbh.uni-hamburg.de>");
  oprval = gt_option_parser_parse(op, parsed_args, argc, argv,
                                  plugins->gth_version_func, err);

  if (oprval == GT_OPTION_PARSER_OK && alignmentscore)
    gthsplitinfo->splitmode = ALIGNMENTSCORE_SPLIT;
  if (oprval == GT_OPTION_PARSER_OK && coverage)
    gthsplitinfo->splitmode = COVERAGE_SPLIT;
  if (oprval == GT_OPTION_PARSER_OK && 100 % gthsplitinfo->range) {
    gt_error_set(err, "argument to option %s must divide 100 without rest",
              RANGE_OPT_CSTR);
    oprval = GT_OPTION_PARSER_ERROR;
  }
  if (oprval == GT_OPTION_PARSER_OK && verbose)
    gthsplitinfo->showverbose = gth_show_on_stdout;
  if (oprval == GT_OPTION_PARSER_OK && gzip)
    gthsplitinfo->file_mode = GT_FILE_MODE_GZIP;
  if (oprval == GT_OPTION_PARSER_OK && bzip2)
    gthsplitinfo->file_mode = GT_FILE_MODE_BZIP2;

  /* save consensus files */
  if (oprval == GT_OPTION_PARSER_OK) {
    while (*parsed_args < argc) {
      gt_str_array_add_cstr(gthsplitinfo->consensusfiles, argv[*parsed_args]);
      (*parsed_args)++;
    }
  }

  if (oprval == GT_OPTION_PARSER_OK &&
      !gt_str_array_size(gthsplitinfo->consensusfiles) &&
      (gt_option_is_set(optgzip) || gt_option_is_set(optbzip2))) {
    gt_error_set(err, "to use compression, at least on input file has to be "
                      "supplied");
    oprval = GT_OPTION_PARSER_ERROR;
  }

  gt_option_parser_delete(op);

  return oprval;
}

static int gthsplit_process_files(Gthsplitinfo *gthsplitinfo,
                                  const GthPlugins *plugins, GtError *err)
{
  Store_in_subset_file_data store_in_subset_file_data;
  GthInput *inputinfo;
  unsigned long i;
  int had_err;

  gt_error_check(err);

  /* initialization */
  inputinfo = gth_input_new(plugins->file_preprocessor, plugins->seq_con_new);
  store_in_subset_file_data.gthsplitinfo           = gthsplitinfo;
  store_in_subset_file_data.num_of_subset_files    = 100 / gthsplitinfo->range;
  store_in_subset_file_data.sa_filter              = gthsplitinfo->sa_filter;
  store_in_subset_file_data.current_outputfilename = NULL;
  store_in_subset_file_data.subset_files =
    gt_malloc(sizeof (GtFile*) *
              store_in_subset_file_data.num_of_subset_files);
  store_in_subset_file_data.subset_filenames =
    gt_malloc(sizeof (GtStr*) *
              store_in_subset_file_data.num_of_subset_files);
  store_in_subset_file_data.subset_file_sa_counter =
    gt_malloc(sizeof (unsigned long) *
              store_in_subset_file_data.num_of_subset_files);
  for (i = 0; i < store_in_subset_file_data.num_of_subset_files; i++) {
    store_in_subset_file_data.subset_files[i]           = NULL;
    store_in_subset_file_data.subset_filenames[i]       = NULL;
    store_in_subset_file_data.subset_file_sa_counter[i] = 0;
  }
  store_in_subset_file_data.sa_visitor = gth_xml_inter_sa_visitor_new(inputinfo,
                                                                      0, NULL);

  if (gthsplitinfo->showverbose)
    gthsplitinfo->showverbose("process all intermediate output files");

  /* split up intermediate files */
  had_err = gth_process_intermediate_files(inputinfo,
                                           gthsplitinfo->consensusfiles,
                                           store_in_subset_file,
                                           &store_in_subset_file_data,
                                           gthsplitinfo->showverbose, err);

  /* close the split files */
  close_output_files(&store_in_subset_file_data);

  /* free */
  gth_sa_visitor_delete(store_in_subset_file_data.sa_visitor);
  gth_input_delete_complete(inputinfo);
  gt_free(store_in_subset_file_data.current_outputfilename);
  gt_free(store_in_subset_file_data.subset_files);
  gt_free(store_in_subset_file_data.subset_filenames);
  gt_free(store_in_subset_file_data.subset_file_sa_counter);

  return had_err;
}

int gt_gthsplit(int argc, const char **argv, const GthPlugins *plugins,
                GtError *err)
{
  Gthsplitinfo gthsplitinfo;
  int had_err, parsed_args;
  gt_error_check(err);

  /* init data structures */
  initGthsplitinfo(&gthsplitinfo);

  switch (gthsplit_parse_options(&parsed_args, &gthsplitinfo, argc, argv,
                                 plugins, err)) {
    case GT_OPTION_PARSER_OK: break;
    case GT_OPTION_PARSER_ERROR:
      freeGthsplitinfo(&gthsplitinfo);
      return -1;
    case GT_OPTION_PARSER_REQUESTS_EXIT:
      freeGthsplitinfo(&gthsplitinfo);
      return 0;
  }
  gt_assert(parsed_args == argc);

  /* process files */
  had_err = gthsplit_process_files(&gthsplitinfo, plugins, err);

  /* free */
  freeGthsplitinfo(&gthsplitinfo);

  return had_err;
}
