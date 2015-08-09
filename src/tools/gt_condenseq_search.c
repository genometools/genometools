/*
  Copyright (c) 2014 Florian Markowsky <moltenboron@web.de>
  Copyright (c) 2014 Dirk Willrodt <willrodt@zbh.uni-hamburg.de>
  Copyright (c) 2014 Center for Bioinformatics, University of Hamburg

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

#include <errno.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <sys/types.h>
#ifndef _WIN32
#include <sys/wait.h>
#endif

#include "core/cstr_array.h"
#include "core/error.h"
#include "core/ma_api.h"
#include "core/toolbox_api.h"
#include "core/unused_api.h"

#include "tools/gt_condenseq_blast.h"
#include "tools/gt_condenseq_hmmsearch.h"

#include "tools/gt_condenseq_search.h"

static void* gt_condenseq_search_toolbox_new(void)
{
  GtToolbox *condenseq_search_toolbox = gt_toolbox_new();
  gt_toolbox_add_tool(condenseq_search_toolbox,
                      "blast", gt_condenseq_blast());
  gt_toolbox_add_tool(condenseq_search_toolbox,
                      "hmmsearch", gt_condenseq_hmmsearch());
  return condenseq_search_toolbox;
}

static void gt_condenseq_search_toolbox_delete(void *toolbox)
{
  GtToolbox *search_toolbox = toolbox;
  gt_toolbox_delete(search_toolbox);
}

static GtOptionParser*
  gt_condenseq_search_option_parser_new(void *toolbox)
{
  GtToolbox *search_toolbox = toolbox;
  GtOptionParser *option_parser;
  gt_assert(search_toolbox);

  /* init */
  option_parser = gt_option_parser_new("tool [option ...]",
                            "Call one of the CONDENSER search tools to query "
                            "redundancy compressed genomic data.");
  gt_option_parser_set_comment_func(option_parser, gt_toolbox_show,
                                    search_toolbox);
  gt_option_parser_set_min_args(option_parser, 1U);
  return option_parser;
}

static int gt_condenseq_search_runner(GT_UNUSED int argc,
                                      GT_UNUSED const char **argv,
                                      GT_UNUSED int parsed_args,
                                      void *toolbox,
                                      GtError *err)
{
  GtToolbox *search_toolbox = toolbox;
  GtTool *tool;
  int had_err = 0;
  char **nargv = NULL;

  gt_error_check(err);
  gt_assert(search_toolbox);

  /* get search tools */
  tool = gt_toolbox_get_tool(search_toolbox, argv[parsed_args]);
  if (!tool) {
    gt_error_set(err,
                 "condenseq search tool '%s' not found; option -help lists "
                 "possible tools", argv[parsed_args]);
    had_err = -1;
  }
  /* call search tool */
  if (!had_err) {
    gt_assert(tool);
    nargv = gt_cstr_array_prefix_first(argv + parsed_args,
                                       gt_error_get_progname(err));
    gt_error_set_progname(err, nargv[0]);
    had_err = gt_tool_run(tool, argc - parsed_args, (const char**) nargv, err);
  }

  /* free */
  gt_cstr_array_delete(nargv);

  return had_err;
}

GtTool* gt_condenseq_search(void)
{
  return gt_tool_new(gt_condenseq_search_toolbox_new,
                     gt_condenseq_search_toolbox_delete,
                     gt_condenseq_search_option_parser_new,
                     NULL,
                     gt_condenseq_search_runner);
}
