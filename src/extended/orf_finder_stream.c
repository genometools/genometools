/*
  Copyright (c) 2010      Sascha Kastens <sascha.kastens@studium.uni-hamburg.de>
  Copyright (c)      2011 Sascha Steinbiss <steinbiss@zbh.uni-hamburg.de>
  Copyright (c) 2010-2011 Center for Bioinformatics, University of Hamburg

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

#include <string.h>
#include "core/class_alloc_lock.h"
#include "core/log.h"
#include "core/ma.h"
#include "core/mathsupport.h"
#include "core/range.h"
#include "core/str.h"
#include "core/unused_api.h"
#include "extended/node_stream_api.h"
#include "extended/feature_node.h"
#include "extended/feature_node_iterator_api.h"
#include "extended/reverse_api.h"
#include "extended/orf_finder_stream.h"
#include "extended/orf_finder_visitor.h"

struct GtORFFinderStream {
  const GtNodeStream parent_instance;
  GtNodeStream *in_stream;
  GtORFFinderVisitor *lv;
};

#define gt_orf_finder_stream_cast(GS)\
        gt_node_stream_cast(gt_orf_finder_stream_class(), GS)

static int gt_orf_finder_stream_next(GtNodeStream *gs, GtGenomeNode **gn,
                                    GtError *err)
{
  GtORFFinderStream *ls;
  int had_err;

  gt_error_check(err);
  ls = gt_orf_finder_stream_cast(gs);

  had_err = gt_node_stream_next(ls->in_stream, gn, err);
  if (!had_err && *gn) {
    (void) gt_genome_node_accept(*gn, (GtNodeVisitor*) ls->lv, err);
  }
  if (had_err) {
    gt_genome_node_delete(*gn);
    *gn = NULL;
  }
  return had_err;
}

static void gt_orf_finder_stream_free(GtNodeStream *gs)
{
  GtORFFinderStream *ls = gt_orf_finder_stream_cast(gs);
  gt_node_visitor_delete((GtNodeVisitor*) ls->lv);
  gt_node_stream_delete(ls->in_stream);
}

const GtNodeStreamClass* gt_orf_finder_stream_class(void)
{
  static const GtNodeStreamClass *gsc;
  gt_class_alloc_lock_enter();
  if (!gsc) {
    gsc = gt_node_stream_class_new(sizeof (GtORFFinderStream),
                                   gt_orf_finder_stream_free,
                                   gt_orf_finder_stream_next);
  }
  gt_class_alloc_lock_leave();
  return gsc;
}

GtNodeStream* gt_orf_finder_stream_new(GtNodeStream *in_stream,
                                       GtEncseq *encseq,
                                       GtHashmap *types,
                                       unsigned int min,
                                       unsigned int max,
                                       bool all,
                                       GtError *err)
{
  GtNodeStream *gs;
  GtORFFinderStream *ls;
  gs = gt_node_stream_create(gt_orf_finder_stream_class(), false);
  ls = gt_orf_finder_stream_cast(gs);
  ls->in_stream = gt_node_stream_ref(in_stream);
  ls->lv = (GtORFFinderVisitor*) gt_orf_finder_visitor_new(encseq,types,
                                                           min, max, all, err);
  return gs;
}
