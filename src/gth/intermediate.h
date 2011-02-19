/*
  Copyright (c) 2004-2011 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2004-2008 Center for Bioinformatics, University of Hamburg

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

#ifndef INTERMEDIATE_H
#define INTERMEDIATE_H

#include "gth/sa_filter.h"
#include "gth/sa_collection.h"

/* XXX: impelement hashing method */
#define GTH_UNDEFINED_HASH  "undefined"

typedef int (*GthSAProcessFunc)(void *data, GthSA*,
                                const char *outputfilename, GtError*);

bool gth_intermediate_output_is_correct(char *outputfilename,
                                        GthSACollection *orig_sacollection,
                                        GthInput*, GtFile **outfp, GtError*);
int  gth_process_intermediate_files(GthInput*, GtStrArray *consensusfiles,
                                    GthSAProcessFunc, void *data,
                                    GthShowVerbose, GtError*);

/* The following function builds a tree of alignments from a set of consensus
  files. If no consensus file is given, stdin is used as input. */
int  gth_build_sa_collection(GthSACollection*, GthInput*,
                             GtStrArray *consensusfiles, GthSAFilter*, GthStat*,
                             GthShowVerbose, GtError*);

#endif
