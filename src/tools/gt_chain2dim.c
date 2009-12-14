/*
  Copyright (c) 2009 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
  Copyright (c) 2009 Center for Bioinformatics, University of Hamburg

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

#include "core/option.h"
#include "core/ma.h"
#include "core/unused_api.h"
#include "core/tool.h"
#include "match/chaindef.h"

static void *gt_chain2dim_arguments_new (void)
{
  return gt_malloc (sizeof (Chaincalloptions));
}

static void gt_chain2dim_arguments_delete (void *tool_arguments)
{
  Chaincalloptions *arguments = tool_arguments;

  if (!arguments)
  {
    return;
  }
  gt_str_delete (arguments->matchfile);
  gt_str_array_delete (arguments->globalargs);
  gt_free (arguments);
}

static GtOptionParser *gt_chain2dim_option_parser_new (void *tool_arguments)
{
  Chaincalloptions *arguments = tool_arguments;
  GtOptionParser *op;
  GtOption *option;

  gt_assert (arguments != NULL);
  arguments->matchfile = gt_str_new ();
  arguments->globalargs = gt_str_array_new();

  op = gt_option_parser_new
    ("[options] -m matchfile","chain pairwise matches");

  gt_option_parser_set_mailaddress (op, "<kurtz@zbh.uni-hamburg.de>");
  option = gt_option_new_filename("m","Specify file containing the matches",
                                  arguments->matchfile);
  gt_option_parser_add_option(op, option);
  gt_option_is_mandatory (option);

  optionglobal = gt_option_new_stringarray("global","perform global chaining",
                                           arguments->globalargs);
  gt_option_parser_add_option(op, optionglobal);

  optionlocal = gt_option_new_stringarray("local","perform local chaining",
                                          arguments->localargs);
  gt_option_parser_add_option(op, optionlocal);

  return op;
}
