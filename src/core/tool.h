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

#ifndef TOOL_H
#define TOOL_H

#include "core/option.h"

/* the tool class */
typedef struct GtTool GtTool;

/* tool functions */
typedef void*           (*GtToolArgumentsNew)(void);
typedef void            (*GtToolArgumentsDelete)(void *tool_arguments);
typedef GtOptionParser* (*GtToolOptionParserNew)(void *tool_arguments);
typedef int             (*GtToolArgumentsCheck)(int rest_argc,
                                                void *tool_arguments, GtError*);
typedef int             (*GtToolRunner)(int argc, const char **argv,
                                        int parsed_args, void *tool_arguments,
                                        GtError*);

/* the type of a tool constructor */
typedef GtTool*         (*GtToolConstructor)(void);

/*
   Create a new tool object, with
   - a tool argument constructor <gt_tool_arguments_new> (optional),
   - a tool argument destructor <gt_tool_arguments_delete> (optional).
   - a tool option parser constructor <gt_tool_option_parser_new) (required),
   - a tool argument checker <gt_tool_arguments_check> (optional),
   - a tool runner <gt_tool_runner> (required), and
   <tool_arguments_new> and <tool_arguments_delete> imply each other.
   Returns a new GtTool object.
*/
GtTool* gt_tool_new(GtToolArgumentsNew tool_arguments_new,
                    GtToolArgumentsDelete tool_arguments_delete,
                    GtToolOptionParserNew tool_option_parser_new,
                    GtToolArgumentsCheck tool_arguments_check,
                    GtToolRunner tool_runner);

/*
  Run the given <tool> as follows:
  - Create a tool arguments object, if necessary.
  - Create a new option parser and pass the tool arguments along.
  - Parse the options (<argc> and <argv>) with the created option parser.
  - Return upon error, continue otherwise.
  - Check the tool arguments, if necessary.
  - Return upon error, continue otherwise.
  - Run the actual tool with the given arguments (the tool arguments object
    is passed along).
  - Delete the tool arguments object, if one was created.
  Returns -1 and sets <err> on error, returns 0 otherwise.
*/
int   gt_tool_run(GtTool*, int argc, const char **argv, GtError *err);

/* Delete the given <tool>. */
void  gt_tool_delete(GtTool*);

#endif
