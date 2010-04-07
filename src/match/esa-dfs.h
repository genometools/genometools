/*
  Copyright (c) 2007 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
  Copyright (c) 2007 Center for Bioinformatics, University of Hamburg

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

#ifndef ESA_DFS_H
#define ESA_DFS_H
#include <stdbool.h>
#include "core/error_api.h"
#include "core/unused_api.h"

#include "esa-seqread.h"

typedef struct Dfsinfo Dfsinfo;
typedef struct Dfsstate Dfsstate;

int gt_depthfirstesa(Sequentialsuffixarrayreader *ssar,
                  Dfsinfo *(*allocateDfsinfo)(Dfsstate *),
                  void(*freeDfsinfo)(Dfsinfo *,Dfsstate *),
                  int(*processleafedge)(bool,unsigned long,Dfsinfo *,
                                        unsigned long,Dfsstate *,
                                        GtError *),
                  int(*processbranchedge)(bool,
                                          unsigned long,
                                          Dfsinfo *,
                                          Dfsinfo *,
                                          Dfsstate *,
                                          GtError *),
                  int(*processcompletenode)(unsigned long,
                                            Dfsinfo *,
                                            unsigned long,
                                            Dfsstate *,
                                            GtError *),
                  void(*assignleftmostleaf)(Dfsinfo *,
                                            unsigned long,
                                            Dfsstate *),
                  void(*assignrightmostleaf)(Dfsinfo *,
                                             unsigned long,
                                             unsigned long,
                                             unsigned long,
                                             Dfsstate *),
                  Dfsstate *state,
                  GT_UNUSED GtLogger *logger,
                  GtError *err);

#endif
