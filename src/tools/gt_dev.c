/*
  Copyright (c) 2007-2008 Gordon Gremme <gremme@zbh.uni-hamburg.de>
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

#include "libgtcore/cstr.h"
#include "libgtcore/option.h"
#include "libgtcore/versionfunc.h"
#include "libgtext/toolbox.h"
#include "tools/gt_dev.h"
#include "tools/gt_guessprot.h"
#include "tools/gt_regioncov.h"
#include "tools/gt_seqiterator.h"
#include "tools/gt_sfxmap.h"
#include "tools/gt_trieins.h"
#include "tools/gt_mergeesa.h"
#include "tools/gt_maxpairs.h"
#include "tools/gt_paircmp.h"
#include "tools/gt_patternmatch.h"
#include "tools/gt_skproto.h"

static OPrval parse_options(int *parsed_args, int argc, const char **argv,
                            Toolbox *dev_toolbox, Error *err)
{
  OptionParser *op;
  OPrval oprval;
  error_check(err);
  op = option_parser_new("[option ...] dev_tool_name [argument ...]",
                         "Call development tool with name dev_tool_name and "
                         "pass argument(s) to it.");
  option_parser_set_comment_func(op, toolbox_show, dev_toolbox);
  option_parser_set_min_args(op, 1);
  oprval = option_parser_parse(op, parsed_args, argc, argv, versionfunc, err);
  option_parser_delete(op);
  return oprval;
}

void register_devtools(Toolbox *dev_toolbox)
{
  assert(dev_toolbox);
  /* add development tools here with a function call like this:
     toolbox_add(dev_toolbox, "devtool", gt_devtool); */
  toolbox_add(dev_toolbox, "guessprot", gt_guessprot);
  toolbox_add(dev_toolbox, "regioncov", gt_regioncov);
  toolbox_add(dev_toolbox, "sfxmap", gt_sfxmap);
  toolbox_add(dev_toolbox, "seqiterator", gt_seqiterator);
  toolbox_add(dev_toolbox, "trieins", gt_trieins);
  toolbox_add(dev_toolbox, "mergeesa", gt_mergeesa);
  toolbox_add(dev_toolbox, "skproto", gt_skproto);
  toolbox_add(dev_toolbox, "maxpairs", gt_maxpairs);
  toolbox_add(dev_toolbox, "patternmatch", gt_patternmatch);
  toolbox_add(dev_toolbox, "paircmp", gt_paircmp);
}

int gt_dev(int argc, const char **argv, Error *err)
{
  Toolbox *dev_toolbox;
  Toolfunc devtoolfunc;
  int parsed_args, had_err = 0;
  char **nargv = NULL;
  error_check(err);

  /* option parsing */
  dev_toolbox = toolbox_new();
  register_devtools(dev_toolbox);
  switch (parse_options(&parsed_args, argc, argv, dev_toolbox, err)) {
    case OPTIONPARSER_OK: break;
    case OPTIONPARSER_ERROR:
      toolbox_delete(dev_toolbox);
      return -1;
    case OPTIONPARSER_REQUESTS_EXIT:
      toolbox_delete(dev_toolbox);
      return 0;
  }
  assert(parsed_args < argc);

  /* get development tools */
  if (!(devtoolfunc = toolbox_get(dev_toolbox, argv[1]))) {
    error_set(err, "development tool '%s' not found; option -help lists "
                   "possible tools", argv[1]);
    had_err = -1;
  }

  /* call development tool */
  if (!had_err) {
    nargv = cstr_array_prefix_first(argv+parsed_args, argv[0]);
    error_set_progname(err, nargv[0]);
    had_err = devtoolfunc(argc-parsed_args, (const char**) nargv, err);
  }

  /* free */
  cstr_array_delete(nargv);
  toolbox_delete(dev_toolbox);

  return had_err;
}
