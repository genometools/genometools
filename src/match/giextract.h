/*
  Copyright (c) 2007-2008 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
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

#ifndef GIEXTRACT_H
#define GIEXTRACT_H
#include <stdbool.h>
#include "core/str_api.h"
#include "core/fileutils_api.h"
#include "core/error.h"
#include "core/logger.h"

int gt_extractkeysfromfastafile(bool verbose,
                                GtFile *outfp,
                                GtUword width,
                                const GtStr *fileofkeystoextract,
                                GtStrArray *referencefiletab,
                                GtError *err);

int gt_extractkeysfromfastaindex(const char *indexname,
                                 const GtStr *fileofkeystoextract,
                                 GtUword linewidth,GtError *err);

int gt_extractkeysfromdesfile(const char *indexname,
                              bool sortkeys,
                              GtLogger *logger,
                              GtError *err);

bool gt_deskeysfileexists(const char *indexname);

#endif
