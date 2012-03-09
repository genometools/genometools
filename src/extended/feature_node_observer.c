/*
  Copyright (c) 2009 Sascha Steinbiss <steinbiss@zbh.uni-hamburg.de>
  Copyright (c) 2009 Center for Bioinformatics, University of Hamburg

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

#include "core/ma.h"
#include "extended/feature_node_observer.h"

GtFeatureNodeObserver* gt_feature_node_observer_new()
{
  GtFeatureNodeObserver* fno = gt_calloc(1, sizeof (GtFeatureNodeObserver));
  return fno;
}

GtFeatureNodeObserver* gt_feature_node_observer_ref(GtFeatureNodeObserver *o)
{
  gt_assert(o);
  o->reference_count++;
  return o;
}

void gt_feature_node_observer_delete(GtFeatureNodeObserver *o)
{
  if (!o) return;
  if (o->reference_count) {
    o->reference_count--;
    return;
  }
  if (o->data)
    gt_free(o->data);
  gt_free(o);
}
