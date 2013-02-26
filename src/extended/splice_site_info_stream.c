/*
  Copyright (c) 2007-2011 Gordon Gremme <gremme@zbh.uni-hamburg.de>
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

#include "core/assert_api.h"
#include "core/class_alloc_lock.h"
#include "extended/genome_node.h"
#include "extended/node_stream_api.h"
#include "extended/splice_site_info_stream.h"
#include "extended/splice_site_info_visitor.h"

struct GtSpliceSiteInfoStream {
  const GtNodeStream parent_instance;
  GtNodeStream *in_stream;
  GtNodeVisitor *splice_site_info_visitor;
};

#define gt_splice_site_info_stream_cast(GS)\
        gt_node_stream_cast(gt_splice_site_info_stream_class(), GS)

static int gt_splice_site_info_stream_next(GtNodeStream *ns, GtGenomeNode **gn,
                                           GtError *err)
{
  GtSpliceSiteInfoStream *ssis;
  int had_err;
  gt_error_check(err);
  ssis = gt_splice_site_info_stream_cast(ns);
  had_err = gt_node_stream_next(ssis->in_stream, gn, err);
  if (!had_err) {
    gt_assert(ssis->splice_site_info_visitor);
    if (*gn) {
      had_err = gt_genome_node_accept(*gn, ssis->splice_site_info_visitor, err);
      if (had_err) {
        /* we own the node -> delete it */
        gt_genome_node_delete(*gn);
        *gn = NULL;
      }
    }
  }
  return had_err;
}

static void gt_splice_site_info_stream_free(GtNodeStream *ns)
{
  GtSpliceSiteInfoStream *ssis = gt_splice_site_info_stream_cast(ns);
  gt_node_visitor_delete(ssis->splice_site_info_visitor);
  gt_node_stream_delete(ssis->in_stream);
}

const GtNodeStreamClass* gt_splice_site_info_stream_class(void)
{
  static const GtNodeStreamClass *nsc = NULL;
  gt_class_alloc_lock_enter();
  if (!nsc) {
    nsc = gt_node_stream_class_new(sizeof (GtSpliceSiteInfoStream),
                                   gt_splice_site_info_stream_free,
                                   gt_splice_site_info_stream_next);
  }
  gt_class_alloc_lock_leave();
  return nsc;
}

GtNodeStream* gt_splice_site_info_stream_new(GtNodeStream *in_stream,
                                             GtRegionMapping *rm)
{
  GtNodeStream *ns = gt_node_stream_create(gt_splice_site_info_stream_class(),
                                           false);
  GtSpliceSiteInfoStream *ssis = gt_splice_site_info_stream_cast(ns);
  ssis->in_stream = gt_node_stream_ref(in_stream);
  ssis->splice_site_info_visitor = gt_splice_site_info_visitor_new(rm);
  return ns;
}

bool gt_splice_site_info_stream_show(GtNodeStream *ns, GtFile *outfp)
{
  GtSpliceSiteInfoStream *ssis;
  gt_assert(ns);
  ssis = gt_splice_site_info_stream_cast(ns);
  return gt_splice_site_info_visitor_show(ssis->splice_site_info_visitor,
                                          outfp);
}

bool gt_splice_site_info_stream_intron_processed(GtNodeStream *ns)
{
  GtSpliceSiteInfoStream *ssis;
  gt_assert(ns);
  ssis = gt_splice_site_info_stream_cast(ns);
  return gt_splice_site_info_visitor_intron_processed(ssis
                                                    ->splice_site_info_visitor);
}

bool gt_splice_site_info_stream_show_canonical(GtNodeStream *ns, bool show_gc)
{
  GtSpliceSiteInfoStream *ssis;
  gt_assert(ns);
  ssis = gt_splice_site_info_stream_cast(ns);
  return gt_splice_site_info_visitor_show_canonical(ssis
                                                    ->splice_site_info_visitor,
                                                    show_gc);
}
