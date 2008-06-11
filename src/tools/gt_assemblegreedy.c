/*
  Copyright (c) 2008 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2008 Center for Bioinformatics, University of Hamburg

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

#include "libgtcore/bioseq.h"
#include "libgtcore/ma.h"
#include "libgtcore/unused.h"
#include "libgtexercise/fragment_overlaps.h"
#include "libgtexercise/greedy_assembly.h"
#include "tools/gt_assemblegreedy.h"

typedef struct {
  unsigned long minlength;
  bool showoverlaps,
       showpath;
} AssemblegreedyArguments;

static void* gt_assemblegreedy_arguments_new(void)
{
  return ma_malloc(sizeof (AssemblegreedyArguments));
}

static void gt_assemblegreedy_arguments_delete(void *tool_arguments)
{
  AssemblegreedyArguments *arguments = tool_arguments;
  if (!arguments) return;
  ma_free(arguments);
}

static OptionParser* gt_assemblegreedy_option_parser_new(void *tool_arguments)
{
  OptionParser *op;
  Option *option, *showoverlaps_option, *showpath_option;
  AssemblegreedyArguments *arguments = tool_arguments;
  assert(arguments);
  op = option_parser_new("[option ...] fragment_file",
                         "Assemble fragments given in fragment_file in greedy "
                         "fashion.");
  option = option_new_ulong("minlength", "set the minimum length an overlap "
                            "must have to be considered", &arguments->minlength,
                            5);
  option_parser_add_option(op, option);
  showoverlaps_option = option_new_bool("showoverlaps", "show only the "
                                        "overlaps between the fragments",
                                        &arguments->showoverlaps, false);
  option_parser_add_option(op, showoverlaps_option);
  showpath_option = option_new_bool("showpath", "show the assembled fragment "
                                    "path instead of the assembled sequence",
                                    &arguments->showpath, false);
  option_parser_add_option(op, showpath_option);
  option_exclude(showoverlaps_option, showpath_option);
  option_parser_set_min_max_args(op, 1, 1);
  return op;
}

static int gt_assemblegreedy_runner(UNUSED int argc, const char **argv,
                                 int parsed_args, void *tool_arguments,
                                 Error *err)
{
  AssemblegreedyArguments *arguments = tool_arguments;
  Bioseq *fragments;
  int had_err = 0;

  error_check(err);
  assert(arguments);

  /* init */
  fragments = bioseq_new(argv[parsed_args], err);
  if (!fragments)
     had_err = -1;

  if (!had_err) {
    FragmentOverlaps *fragment_overlaps;
    fragment_overlaps = fragment_overlaps_new(fragments, arguments->minlength);

    /* greedy assembly */
    if (arguments->showoverlaps)
      fragment_overlaps_show(fragment_overlaps);
    else {
      GreedyAssembly *greedy_assembly;
      fragment_overlaps_sort(fragment_overlaps);
      greedy_assembly = greedy_assembly_new(fragments, fragment_overlaps);
      if (arguments->showpath)
        greedy_assembly_show_path(greedy_assembly);
      else
        greedy_assembly_show(greedy_assembly, fragments);
      greedy_assembly_delete(greedy_assembly);
    }
    fragment_overlaps_delete(fragment_overlaps);
  }

  /* free */
  bioseq_delete(fragments);

  return had_err;
}

Tool* gt_assemblegreedy(void)
{
  return tool_new(gt_assemblegreedy_arguments_new,
                  gt_assemblegreedy_arguments_delete,
                  gt_assemblegreedy_option_parser_new,
                  NULL,
                  gt_assemblegreedy_runner);
}
