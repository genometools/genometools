/*
  Copyright (c) 2003-2010 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2003-2008 Center for Bioinformatics, University of Hamburg

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

#include "gth/gthxml.h"
#include "gth/call_info.h"
#include "gth/gt_gth.h"
#include "gth/gthverbosefunc.h"
#include "gth/gthverbosefuncvm.h"
#include "gth/input.h"
#include "gth/parse_options.h"
#include "gth/similarity_filter.h"
#include "gth/stat.h"
#include "gth/run_header.h"

int gt_gth(int argc, const char **argv, const GthPlugins *plugins, GtError *err)
{
  GthCallInfo *callinfo;
  GthInput *input;
  GthStat *stat;
  int parsed_args, had_err = 0;

  gt_error_check(err);
  gt_assert(plugins && plugins->file_preprocessor && plugins->seq_col_new);
  gt_assert(plugins->gth_version_func);

  /* init data structures */
  callinfo = gth_call_info_new(argv[0]);
  input = gth_input_new(plugins->file_preprocessor, plugins->seq_col_new);
  stat = gth_stat_new();

  /* parse the options */
  switch (gth_parse_options(callinfo, input, &parsed_args, argc,
                            (const char**) argv, false, NULL, stat,
                            gth_show_on_stdout, gth_show_on_stdout_vmatch,
                            plugins->gth_version_func, err)) {
    case GT_OPTION_PARSER_OK: break;
    case GT_OPTION_PARSER_ERROR:
      gth_stat_delete(stat);
      gth_input_delete_complete(input);
      gth_call_info_delete(callinfo);
      return -1;
    case GT_OPTION_PARSER_REQUESTS_EXIT:
      gth_stat_delete(stat);
      gth_input_delete_complete(input);
      gth_call_info_delete(callinfo);
      return 0;
  }
  gt_assert(parsed_args == argc);

  /* show XML leader */
  if (callinfo->out->xmlout)
    gth_xml_show_leader(callinfo->intermediate, callinfo->out->outfp);

  /* show header for this run */
  gth_run_header_show(callinfo, input, plugins->gth_version,
                      INITIAL_XML_INDENTLEVEL + 1, argv + 1);

  /* make sure the necessary indices of all input files are created */
  if (!had_err) {
    had_err = gth_input_preprocess(input, false,
                                   callinfo->simfilterparam.noautoindex,
                                   callinfo->simfilterparam.skipindexcheck,
                                   callinfo->simfilterparam.maskpolyAtails,
                                   callinfo->simfilterparam.online,
                                   callinfo->simfilterparam.inverse,
                                   callinfo->progname,
                                   gt_str_get(callinfo->scorematrixfile),
                                   callinfo->translationtable, callinfo->out,
                                   err);
  }

  if (!had_err) {
    if (callinfo->simfilterparam.createindicesonly) {
      if (callinfo->out->showverbose)
        callinfo->out->showverbose("the indices have been created => stopping");
    }
    else {
      /* set and check position values */
      had_err = gth_input_set_and_check_substring_spec(input, err);

      if (!had_err && callinfo->out->showverbose)
        callinfo->out->showverbose("invoking similarity filter");

      if (!had_err) {
        had_err = gth_similarity_filter(callinfo, input, stat,
                                        INITIAL_XML_INDENTLEVEL, plugins, err);
      }
    }
  }

  /* free space */
  gth_stat_delete(stat);
  gth_input_delete_complete(input);
  gth_call_info_delete(callinfo);

  return had_err;
}
