/*
  Copyright (c) 2007-2010 Gordon Gremme <gordon@gremme.org>
  Copyright (c) 2007-2008 Center for Bioinformatics, University of Hamburg

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

#include "core/cstr_array.h"
#include "core/option_api.h"
#include "core/tool.h"
#include "core/toolbox.h"
#include "core/versionfunc.h"
#include "gth/gt_gthbssmbuild.h"
#include "gth/gt_gthbssmfileinfo.h"
#include "gth/gt_gthbssmprint.h"
#include "gth/gt_gthbssmrmsd.h"
#include "gth/gt_gthbssmtrain.h"
#include "gth/gt_gthmkbssmfiles.h"
#include "tools/gt_compressedbits.h"
#include "tools/gt_consensus_sa.h"
#include "tools/gt_dev.h"
#include "tools/gt_extracttarget.h"
#include "tools/gt_gdiffcalc.h"
#include "tools/gt_guessprot.h"
#include "tools/gt_idxlocali.h"
#include "tools/gt_magicmatch.h"
#include "tools/gt_mergeesa.h"
#include "tools/gt_paircmp.h"
#include "tools/gt_parsexrf.h"
#include "tools/gt_patternmatch.h"
#include "tools/gt_readreads.h"
#include "tools/gt_regioncov.h"
#include "tools/gt_sain.h"
#include "tools/gt_sam_interface.h"
#include "tools/gt_seqcorrect.h"
#include "tools/gt_seqlensort.h"
#include "tools/gt_sfxmap.h"
#include "tools/gt_skproto.h"
#include "tools/gt_sortbench.h"
#include "tools/gt_trieins.h"

static void* gt_dev_arguments_new(void)
{
  GtToolbox *dev_toolbox = gt_toolbox_new();
  /* add development tools here with a function call like this:
     gt_toolbox_add(dev_toolbox, "devtool", gt_devtool); */
  gt_toolbox_add(dev_toolbox, "gthbssmbuild", gt_gthbssmbuild);
  gt_toolbox_add(dev_toolbox, "gthbssmfileinfo", gt_gthbssmfileinfo);
  gt_toolbox_add(dev_toolbox, "gthbssmprint", gt_gthbssmprint);
  gt_toolbox_add(dev_toolbox, "gthmkbssmfiles", gt_gthmkbssmfiles);
  gt_toolbox_add(dev_toolbox, "guessprot", gt_guessprot);
  gt_toolbox_add(dev_toolbox, "mergeesa", gt_mergeesa);
  gt_toolbox_add(dev_toolbox, "paircmp", gt_paircmp);
  gt_toolbox_add(dev_toolbox, "patternmatch", gt_patternmatch);
  gt_toolbox_add(dev_toolbox, "regioncov", gt_regioncov);
  gt_toolbox_add(dev_toolbox, "trieins", gt_trieins);
  gt_toolbox_add_tool(dev_toolbox, "compbits", gt_compressedbits());
  gt_toolbox_add_tool(dev_toolbox, "consensus_sa", gt_consensus_sa_tool());
  gt_toolbox_add_tool(dev_toolbox, "extracttarget", gt_extracttarget());
  gt_toolbox_add_tool(dev_toolbox, "gdiffcalc", gt_gdiffcalc());
  gt_toolbox_add_tool(dev_toolbox, "gthbssmrmsd", gt_gthbssmrmsd());
  gt_toolbox_add_tool(dev_toolbox, "gthbssmtrain", gt_gthbssmtrain());
  gt_toolbox_add_tool(dev_toolbox, "idxlocali", gt_idxlocali());
  gt_toolbox_add_tool(dev_toolbox, "magicmatch", gt_magicmatch());
  gt_toolbox_add_tool(dev_toolbox, "parsexrf", gt_parsexrf());
  gt_toolbox_add_tool(dev_toolbox, "readreads", gt_readreads());
  gt_toolbox_add_tool(dev_toolbox, "sain", gt_sain());
  gt_toolbox_add_tool(dev_toolbox, "sambam", gt_sam_interface());
  gt_toolbox_add_tool(dev_toolbox, "seqcorrect", gt_seqcorrect());
  gt_toolbox_add_tool(dev_toolbox, "seqlensort", gt_seqlensort());
  gt_toolbox_add_tool(dev_toolbox, "sfxmap", gt_sfxmap());
  gt_toolbox_add_tool(dev_toolbox, "skproto", gt_skproto());
  gt_toolbox_add_tool(dev_toolbox, "sortbench", gt_sortbench());
  return dev_toolbox;
}

static void gt_dev_arguments_delete(void *tool_arguments)
{
  GtToolbox *dev_toolbox = tool_arguments;
  if (!dev_toolbox) return;
  gt_toolbox_delete(dev_toolbox);
}

static GtOptionParser* gt_dev_option_parser_new(void *tool_arguments)
{
  GtToolbox *dev_toolbox = tool_arguments;
  GtOptionParser *op;
  gt_assert(dev_toolbox);
  op = gt_option_parser_new("[option ...] dev_tool_name [argument ...]",
                            "Call development tool with name dev_tool_name and "
                            "pass argument(s) to it.");
  gt_option_parser_set_comment_func(op, gt_toolbox_show, dev_toolbox);
  gt_option_parser_set_min_args(op, 1);
  return op;
}

static int gt_dev_runner(int argc, const char **argv, int parsed_args,
                         void *tool_arguments, GtError *err)
{
  GtToolbox *dev_toolbox = tool_arguments;
  GtToolfunc toolfunc;
  GtTool *tool = NULL;
  int had_err = 0;
  char **nargv = NULL;

  gt_error_check(err);
  gt_assert(dev_toolbox);

  /* get development tools */
  if (!gt_toolbox_has_tool(dev_toolbox, argv[parsed_args])) {
    gt_error_set(err, "development tool '%s' not found; option -help lists "
                      "possible tools", argv[parsed_args]);
    had_err = -1;
  }

  /* call development tool */
  if (!had_err) {
    if (!(toolfunc = gt_toolbox_get(dev_toolbox, argv[parsed_args]))) {
      tool = gt_toolbox_get_tool(dev_toolbox, argv[parsed_args]);
      gt_assert(tool);
    }
    nargv = gt_cstr_array_prefix_first(argv + parsed_args,
                                       gt_error_get_progname(err));
    gt_error_set_progname(err, nargv[0]);
    if (toolfunc)
      had_err = toolfunc(argc - parsed_args, (const char**) nargv, err);
    else
      had_err = gt_tool_run(tool, argc - parsed_args, (const char**) nargv,
                            err);
  }

  /* free */
  gt_cstr_array_delete(nargv);

  return had_err;
}

GtTool* gt_dev(void)
{
  GtTool *tool = gt_tool_new(gt_dev_arguments_new,
                             gt_dev_arguments_delete,
                             gt_dev_option_parser_new,
                             NULL,
                             gt_dev_runner);
  gt_tool_set_toolbox_new(tool, (GtToolToolboxNew) gt_dev_arguments_new);
  return tool;
}
