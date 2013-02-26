/*
  Copyright (c) 2006-2012 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2006-2008 Center for Bioinformatics, University of Hamburg

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

#include "core/class_alloc_lock.h"
#include "extended/add_ids_stream.h"
#include "extended/cds_check_stream.h"
#include "extended/gff3_in_stream.h"
#include "extended/gff3_in_stream_plain.h"
#include "extended/multi_sanitizer_visitor.h"
#include "extended/node_stream_api.h"
#include "extended/tidy_region_node_stream.h"
#include "extended/visitor_stream_api.h"

struct GtGFF3InStream {
  const GtNodeStream parent_instance;
  GtNodeStream *gff3_in_stream_plain,
               *add_ids_stream,
               *cds_check_stream,
               *fix_region_stream,
               *last_stream,
               *multi_sanitize_stream;
};

#define gff3_in_stream_cast(NS)\
        gt_node_stream_cast(gt_gff3_in_stream_class(), NS)

static int gff3_in_stream_next(GtNodeStream *ns, GtGenomeNode **gn,
                               GtError *err)
{
  GtGFF3InStream *is;
  gt_error_check(err);
  is = gff3_in_stream_cast(ns);
  return gt_node_stream_next(is->last_stream, gn, err);
}

static void gff3_in_stream_free(GtNodeStream *ns)
{
  GtGFF3InStream *gff3_in_stream = gff3_in_stream_cast(ns);
  gt_node_stream_delete(gff3_in_stream->cds_check_stream);
  gt_node_stream_delete(gff3_in_stream->add_ids_stream);
  gt_node_stream_delete(gff3_in_stream->gff3_in_stream_plain);
  gt_node_stream_delete(gff3_in_stream->fix_region_stream);
  gt_node_stream_delete(gff3_in_stream->multi_sanitize_stream);
}

const GtNodeStreamClass* gt_gff3_in_stream_class(void)
{
  static const GtNodeStreamClass *nsc = NULL;
  gt_class_alloc_lock_enter();
  if (!nsc) {
    nsc = gt_node_stream_class_new(sizeof (GtGFF3InStream),
                                   gff3_in_stream_free,
                                   gff3_in_stream_next);
  }
  gt_class_alloc_lock_leave();
  return nsc;
}

void gt_gff3_in_stream_check_id_attributes(GtGFF3InStream *is)
{
  gt_assert(is);
  gt_gff3_in_stream_plain_check_id_attributes((GtGFF3InStreamPlain*)
                                              is->gff3_in_stream_plain);
}

void gt_gff3_in_stream_show_progress_bar(GtGFF3InStream *is)
{
  gt_assert(is);
  gt_gff3_in_stream_plain_show_progress_bar((GtGFF3InStreamPlain*)
                                            is->gff3_in_stream_plain);
}

void gt_gff3_in_stream_set_type_checker(GtNodeStream *ns,
                                        GtTypeChecker *type_checker)
{
  GtGFF3InStream *is = gff3_in_stream_cast(ns);
  gt_assert(is);
  gt_gff3_in_stream_plain_set_type_checker(is->gff3_in_stream_plain,
                                           type_checker);
}

GtStrArray* gt_gff3_in_stream_get_used_types(GtNodeStream *ns)
{
  GtGFF3InStream *is = gff3_in_stream_cast(ns);
  gt_assert(is);
  return gt_gff3_in_stream_plain_get_used_types(is->gff3_in_stream_plain);
}

void gt_gff3_in_stream_set_offset(GtNodeStream *ns, long offset)
{
  GtGFF3InStream *is = gff3_in_stream_cast(ns);
  gt_assert(is);
  gt_gff3_in_stream_plain_set_offset(is->gff3_in_stream_plain, offset);
}

int gt_gff3_in_stream_set_offsetfile(GtNodeStream *ns, GtStr *offsetfile,
                                     GtError *err)
{
  GtGFF3InStream *is = gff3_in_stream_cast(ns);
  gt_assert(is);
  return gt_gff3_in_stream_plain_set_offsetfile(is->gff3_in_stream_plain,
                                                offsetfile, err);
}

void gt_gff3_in_stream_disable_add_ids(GtNodeStream *ns)
{
  GtGFF3InStream *is = gff3_in_stream_cast(ns);
  gt_assert(is && is->add_ids_stream);
  gt_add_ids_stream_disable(is->add_ids_stream);
}

void gt_gff3_in_stream_enable_strict_mode(GtGFF3InStream *is)
{
  gt_assert(is);
  gt_gff3_in_stream_plain_enable_strict_mode(is->gff3_in_stream_plain);
}

void gt_gff3_in_stream_enable_tidy_mode(GtGFF3InStream *is)
{
  gt_assert(is);
  gt_gff3_in_stream_plain_enable_tidy_mode(is->gff3_in_stream_plain);
  gt_cds_check_stream_enable_tidy_mode((GtCDSCheckStream*)
                                       is->cds_check_stream);
}

void gt_gff3_in_stream_fix_region_boundaries(GtGFF3InStream *is)
{
  gt_assert(is);
  gt_gff3_in_stream_plain_do_not_check_region_boundaries(
                               (GtGFF3InStreamPlain*) is->gff3_in_stream_plain);
  is->last_stream = is->fix_region_stream =
             gt_tidy_region_node_stream_new(is->last_stream);
}

GtNodeStream* gt_gff3_in_stream_new_unsorted(int num_of_files,
                                             const char **filenames)
{
  GtNodeStream *ns = gt_node_stream_create(gt_gff3_in_stream_class(), false);
  GtGFF3InStream *is = gff3_in_stream_cast(ns);
  is->fix_region_stream = NULL;
  is->last_stream = is->gff3_in_stream_plain =
                  gt_gff3_in_stream_plain_new_unsorted(num_of_files, filenames);
  gt_gff3_in_stream_plain_check_region_boundaries(
                               (GtGFF3InStreamPlain*) is->gff3_in_stream_plain);
  is->last_stream = is->add_ids_stream = gt_add_ids_stream_new(is->last_stream);
  is->last_stream = is->multi_sanitize_stream =
       gt_visitor_stream_new(is->last_stream, gt_multi_sanitizer_visitor_new());
  is->last_stream = is->cds_check_stream =
                                       gt_cds_check_stream_new(is->last_stream);
  return ns;
}

GtNodeStream* gt_gff3_in_stream_new_sorted(const char *filename)
{
  GtNodeStream *ns = gt_node_stream_create(gt_gff3_in_stream_class(), true);
  GtGFF3InStream *is = gff3_in_stream_cast(ns);
  is->fix_region_stream = NULL;
  is->last_stream = is->gff3_in_stream_plain =
                                   gt_gff3_in_stream_plain_new_sorted(filename);
  gt_gff3_in_stream_plain_check_region_boundaries(
                               (GtGFF3InStreamPlain*) is->gff3_in_stream_plain);
  is->last_stream = is->add_ids_stream = gt_add_ids_stream_new(is->last_stream);
  is->last_stream = is->multi_sanitize_stream =
       gt_visitor_stream_new(is->last_stream, gt_multi_sanitizer_visitor_new());
  is->last_stream = is->cds_check_stream =
                                       gt_cds_check_stream_new(is->last_stream);
  return ns;
}
