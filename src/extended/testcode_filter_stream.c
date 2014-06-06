/*
  Copyright (c) 2014 Sascha Steinbiss <ss34@sanger.ac.uk>
  Copyright (c) 2014 Genome Research Ltd.

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

#include <ctype.h>
#include "core/array_api.h"
#include "core/assert_api.h"
#include "core/class_alloc_lock.h"
#include "core/log.h"
#include "core/ma.h"
#include "core/parseutils.h"
#include "core/queue.h"
#include "core/symbol_api.h"
#include "core/undef_api.h"
#include "core/unused_api.h"
#include "extended/feature_node.h"
#include "extended/feature_type.h"
#include "extended/merge_stream.h"
#include "extended/testcode_filter_stream.h"

struct GtTestcodeFilterStream {
  const GtNodeStream parent_instance;
  GtNodeStream *in_stream;
  double threshold;
  GtQueue *outqueue;
};

#define gt_testcode_filter_stream_cast(GS)\
        gt_node_stream_cast(gt_testcode_filter_stream_class(), GS);

static int testcode_filter_stream_next(GtNodeStream *ns, GtGenomeNode **gn,
                                        GtError *err)
{
  GtTestcodeFilterStream *tfs;
  int had_err = 0;
  GtGenomeNode *mygn = NULL;
  GtFeatureNode *fn = NULL;
  const char *cds = gt_symbol(gt_ft_CDS);
  gt_error_check(err);
  tfs = gt_testcode_filter_stream_cast(ns);

  while (!had_err) {
    had_err = gt_node_stream_next(tfs->in_stream, &mygn, err);
    if (!had_err && mygn) {
      if ((fn = gt_feature_node_try_cast(mygn))) {
        const char *type = gt_feature_node_get_type(fn);
        if (type == cds) {
          const char *tcode;
          if ((tcode = gt_feature_node_get_attribute(fn, "tcode_mean"))) {
            double tcode_val;
            if (gt_parse_double(&tcode_val, tcode) == 0) {
              if (tcode_val >= tfs->threshold) {
                *gn = mygn;
                return 0;
              }
            }
          }
        }
        gt_genome_node_delete(mygn);
      } else {
        *gn = mygn;
        return 0;
      }
    } else break;
  }

  return had_err;
}

static void testcode_filter_stream_free(GtNodeStream *ns)
{
  GtTestcodeFilterStream *tfs;
  if (!ns) return;
  tfs = gt_testcode_filter_stream_cast(ns);
  while (gt_queue_size(tfs->outqueue) > 0) {
    gt_genome_node_delete((GtGenomeNode*) gt_queue_get(tfs->outqueue));
  }
  gt_node_stream_delete(tfs->in_stream);
  gt_queue_delete(tfs->outqueue);

}

const GtNodeStreamClass* gt_testcode_filter_stream_class(void)
{
  static const GtNodeStreamClass *nsc = NULL;
  gt_class_alloc_lock_enter();
  if (!nsc) {
    nsc = gt_node_stream_class_new(sizeof (GtTestcodeFilterStream),
                                   testcode_filter_stream_free,
                                   testcode_filter_stream_next);
  }
  gt_class_alloc_lock_leave();
  return nsc;
}

GtNodeStream* gt_testcode_filter_stream_new(GtNodeStream *in_stream,
                                            double threshold)
{
  GtTestcodeFilterStream *tfs;
  GtNodeStream *ns;
  gt_assert(in_stream);
  ns = gt_node_stream_create(gt_testcode_filter_stream_class(), true);
  tfs = gt_testcode_filter_stream_cast(ns);
  tfs->in_stream = gt_node_stream_ref(in_stream);
  tfs->threshold = threshold;
  tfs->outqueue = gt_queue_new();
  return ns;
}
