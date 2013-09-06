/*
  Copyright (c) 2009 Gordon Gremme <gordon@gremme.org>

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

#ifndef MATCHER_H
#define MATCHER_H

#include "gth/input.h"

typedef void* (*GthMatcherArgumentsNew)(bool checksubstrspec,
                                        GthInput *input,
                                        const char *queryfilename,
                                        const char *indexfilename,
                                        bool directmatches,
                                        bool refseqisdna,
                                        const char *progname,
                                        char *proteinsmap,
                                        bool exact,
                                        bool edist,
                                        bool hamming,
                                        GtUword hammingdistance,
                                        GtUword minmatchlength,
                                        GtUword seedlength,
                                        GtUword exdrop,
                                        GtUword prminmatchlen,
                                        GtUword prseedlength,
                                        GtUword prhdist,
                                        GtUword translationtable,
                                        bool online,
                                        bool noautoindex,
                                        bool usepolyasuffix,
                                        bool dbmaskmatch);
typedef void  (*GthMatcherArgumentsDelete)(void *matcher_arguments);
typedef void  (*GthMatcherRunner)(void *matcher_arguments, GthShowVerbose,
                                  GthShowVerboseVM, void *match_processor_data);

#endif
