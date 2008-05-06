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

#include "libgtcore/option.h"

/* the tool class */
typedef struct Tool Tool;

/* tool functions */
typedef void*         (*ToolArgumentsNew)(void);
typedef void          (*ToolArgumentsDelete)(void *tool_arguments);
typedef OptionParser* (*ToolOptionParserNew)(void *tool_arguments);
typedef int           (*ToolArgumentsCheck)(int rest_argc,
                                            void *tool_arguments, Error*);
typedef int           (*ToolRunner)(int argc, const char **argv,
                                    int parsed_args, void *tool_arguments,
                                    Error*);

/* the type of a tool constructor */
typedef Tool*         (*ToolConstructor)(void);

/*
   Create a new tool object, with
   - a tool argument constructor <tool_arguments_new> (optional),
   - a tool argument destructor <tool_arguments_delete> (optional).
   - a tool option parser constructor <tool_option_parser_new) (required),
   - a tool argument checker <tool_arguments_check> (optional),
   - a tool runner <tool_runner> (required), and
   <tool_arguments_new> and <tool_arguments_check> imply each other.
   Returns a new Tool object.
*/
Tool* tool_new(ToolArgumentsNew tool_arguments_new,
               ToolArgumentsDelete tool_arguments_delete,
               ToolOptionParserNew tool_option_parser_new,
               ToolArgumentsCheck tool_arguments_check,
               ToolRunner tool_runner);

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
int   tool_run(Tool*, int argc, const char **argv, Error *err);

/* Delete the given <tool>. */
void  tool_delete(Tool*);

#endif
