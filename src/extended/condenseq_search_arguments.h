/*
  Copyright (c) 2015 Dirk Willrodt <willrodt@zbh.uni-hamburg.de>
  Copyright (c) 2015 Center for Bioinformatics, University of Hamburg

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

#ifndef CONDENSEQ_SEARCH_ARGUMENTS_H
#define CONDENSEQ_SEARCH_ARGUMENTS_H

#include "core/logger_api.h"
#include "core/option_api.h"
#include "core/str_api.h"
#include "extended/condenseq.h"

typedef struct GtCondenseqSearchArguments GtCondenseqSearchArguments;

/* Returns a new <GtCondenseqSearchArguments> object */
GtCondenseqSearchArguments* gt_condenseq_search_arguments_new(void);

/* Deletes <condenseq_search_arguments> and frees all associated memory. */
void                        gt_condenseq_search_arguments_delete(
                        GtCondenseqSearchArguments *condenseq_search_arguments);

/* register the options -db for the mandatory input archive and -verbose for
   verbose output */
void                        gt_condenseq_search_register_options(
                         GtCondenseqSearchArguments *condenseq_search_arguments,
                         GtOptionParser *option_parser);
/* Returns the <GtCondenseq> object read from file given by -db option */
GtCondenseq*                gt_condenseq_search_arguments_read_condenseq(
                   const GtCondenseqSearchArguments *condenseq_search_arguments,
                   GtLogger *logger,
                   GtError *err);
/* Returns true if -verbose was set to true/yes */
bool                        gt_condenseq_search_arguments_verbose(
                  const GtCondenseqSearchArguments *condenseq_search_arguments);

/* Returns the path to a db-file with suffix <suffix>. No checks are done if
   file exists, so this can be used for additional new file names. Caller owns
   the returned <GtStr>! */
GtStr*                      gt_condenseq_search_arguments_db_filename(
                         GtCondenseqSearchArguments *condenseq_search_arguments,
                         const char *suffix);
#endif
