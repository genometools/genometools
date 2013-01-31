/*
  Copyright (c) 2005-2012 Gordon Gremme <gremme@zbh.uni-hamburg.de>
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

#include "gth/gthxml.h"
#include "gth/gthverbosefunc.h"
#include "gth/gthverbosefuncvm.h"
#include "gth/intermediate.h"
#include "gth/parse_options.h"
#include "gth/proc_sa_collection.h"
#include "gth/gt_gthconsensus.h"

static int process_consensus_files(GtStrArray *consensusfiles,
                                   GthCallInfo *call_info, GthInput *input,
                                   GthStat *stat, unsigned long indentlevel,
                                   GtError *err)
{
  GthSACollection *sa_collection;
  int had_err;

  gt_error_check(err);

  /* initialization */
  sa_collection = gth_sa_collection_new(call_info->duplicate_check);

  if (call_info->out->showverbose)
    call_info->out->showverbose("process all intermediate output files");

  /* build tree of alignments from intermediate files */
  had_err = gth_build_sa_collection(sa_collection, input, consensusfiles,
                                    call_info->sa_filter, stat,
                                    call_info->out->showverbose, err);

  /* make sure the necessary indices of all input files are created */
  if (!had_err) {
    had_err = gth_input_preprocess(input, true,
                                   call_info->simfilterparam.noautoindex,
                                   call_info->simfilterparam.createindicesonly,
                                   call_info->simfilterparam.skipindexcheck,
                                   call_info->simfilterparam.maskpolyAtails,
                                   call_info->simfilterparam.online,
                                   call_info->simfilterparam.inverse,
                                   call_info->progname,
                                   gt_str_get(call_info->scorematrixfile),
                                   call_info->translationtable,
                                   call_info->duplicate_check, call_info->out,
                                   err);
  }

  /* set reference MD5s, if necessary */
  gth_sa_collection_set_md5s(sa_collection, input);

  /* process the alignments */
  if (!had_err)
    gth_proc_sa_collection(sa_collection, call_info, input, stat, indentlevel);

  /* show XML trailer */
  if (!had_err && call_info->out->xmlout)
    gth_xml_show_trailer(call_info->intermediate, call_info->out->outfp);

  /* free collection of spliced alignments */
  gth_sa_collection_delete(sa_collection);

  return had_err;
}

int gt_gthconsensus(int argc, const char **argv, const GthPlugins *plugins,
                    GtError *err)
{
  GthCallInfo *call_info;
  GthInput *input;
  GthStat *stat;
  GtStrArray *consensusfiles;
  int parsed_args, had_err = 0;

  gt_error_check(err);

  /* init data structures */
  call_info = gth_call_info_new(argv[0]);
  input = gth_input_new(plugins->file_preprocessor, plugins->seq_con_new);
  stat = gth_stat_new();
  consensusfiles = gt_str_array_new();

  /* parse the options */
  switch (gth_parse_options(call_info, input, &parsed_args, argc,
                            (const char**) argv, true, consensusfiles, stat,
                            gth_show_on_stdout, gth_show_on_stdout_vmatch,
                            plugins->gth_version_func, plugins->jump_table_new,
                            err)) {
    case GT_OPTION_PARSER_OK: break;
    case GT_OPTION_PARSER_ERROR:
      gt_str_array_delete(consensusfiles);
      gth_stat_delete(stat);
      gth_input_delete_complete(input);
      gth_call_info_delete(call_info);
      return -1;
    case GT_OPTION_PARSER_REQUESTS_EXIT:
      gt_str_array_delete(consensusfiles);
      gth_stat_delete(stat);
      gth_input_delete_complete(input);
      gth_call_info_delete(call_info);
      return 0;
  }
  gt_assert(parsed_args == argc);

  /* show XML leader */
  if (call_info->out->xmlout)
    gth_xml_show_leader(call_info->intermediate, call_info->out->outfp);

  /* process consensus files */
  if (!had_err) {
    had_err = process_consensus_files(consensusfiles, call_info, input, stat,
                                      INITIAL_XML_INDENTLEVEL, err);
  }

  /* output statistics */
  if (!had_err && !call_info->out->gff3out)
    gth_stat_show(stat, false, call_info->out->xmlout, call_info->out->outfp);

  /* free space */
  gt_str_array_delete(consensusfiles);
  gth_stat_delete(stat);
  gth_input_delete_complete(input);
  gth_call_info_delete(call_info);

  return had_err;
}
