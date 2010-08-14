/*
  Copyright (c) 2006-2010 Gordon Gremme <gremme@zbh.uni-hamburg.de>
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

#include "extended/cds_check_stream.h"
#include "extended/gff3_in_stream.h"
#include "extended/gff3_in_stream_plain.h"
#include "extended/node_stream_api.h"

struct GtGFF3InStream {
  const GtNodeStream parent_instance;
  GtNodeStream *gff3_in_stream_plain,
               *cds_check_stream;
};

#define gff3_in_stream_cast(NS)\
        gt_node_stream_cast(gt_gff3_in_stream_class(), NS)

static int gff3_in_stream_next(GtNodeStream *ns, GtGenomeNode **gn,
                               GtError *err)
{
  GtGFF3InStream *is = gff3_in_stream_cast(ns);
  return gt_node_stream_next(is->cds_check_stream, gn, err);
}

static void gff3_in_stream_free(GtNodeStream *ns)
{
  GtGFF3InStream *gff3_in_stream = gff3_in_stream_cast(ns);
  gt_node_stream_delete(gff3_in_stream->cds_check_stream);
  gt_node_stream_delete(gff3_in_stream->gff3_in_stream_plain);
}

const GtNodeStreamClass* gt_gff3_in_stream_class(void)
{
  static const GtNodeStreamClass *nsc = NULL;
  if (!nsc) {
    nsc = gt_node_stream_class_new(sizeof (GtGFF3InStream),
                                   gff3_in_stream_free,
                                   gff3_in_stream_next);
  }
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

void gt_gff3_in_stream_enable_tidy_mode(GtGFF3InStream *is)
{
  gt_assert(is);
  gt_gff3_in_stream_plain_enable_tidy_mode(is->gff3_in_stream_plain);
  gt_cds_check_stream_enable_tidy_mode((GtCDSCheckStream*)
                                       is->cds_check_stream);
}

GtNodeStream* gt_gff3_in_stream_new_unsorted(int num_of_files,
                                             const char **filenames)
{
  GtNodeStream *ns = gt_node_stream_create(gt_gff3_in_stream_class(), false);
  GtGFF3InStream *is = gff3_in_stream_cast(ns);
  is->gff3_in_stream_plain = gt_gff3_in_stream_plain_new_unsorted(num_of_files,
                                                                  filenames);
  is->cds_check_stream = gt_cds_check_stream_new(is->gff3_in_stream_plain);
  return ns;
}

GtNodeStream* gt_gff3_in_stream_new_sorted(const char *filename)
{
  GtNodeStream *ns = gt_node_stream_create(gt_gff3_in_stream_class(), true);
  GtGFF3InStream *is = gff3_in_stream_cast(ns);
  is->gff3_in_stream_plain = gt_gff3_in_stream_plain_new_sorted(filename);
  is->cds_check_stream = gt_cds_check_stream_new(is->gff3_in_stream_plain);
  return ns;
}
