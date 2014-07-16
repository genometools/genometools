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
#include <string.h>
#include "core/array_api.h"
#include "core/assert_api.h"
#include "core/class_alloc_lock.h"
#include "core/file_api.h"
#include "core/log.h"
#include "core/ma.h"
#include "core/parseutils_api.h"
#include "core/queue.h"
#include "core/qsort_r_api.h"
#include "core/splitter_api.h"
#include "core/symbol_api.h"
#include "core/undef_api.h"
#include "core/unused_api.h"
#include "extended/feature_node.h"
#include "extended/feature_type.h"
#include "extended/gff3_defines.h"
#include "extended/gff3_linesorted_out_stream.h"
#include "extended/gff3_visitor.h"

struct GtGFF3LinesortedOutStream {
  const GtNodeStream parent_instance;
  GtNodeStream *in_stream;
  GtArray *cur_node_set;
  GtRange cur_node_range;
  GtQueue *outqueue;
  GtFile *outfp;
  GtStr *buf;
  GtSplitter *splitter;
  GtNodeVisitor *gff3vis;
  char **outstrings;
  GtUword outstrings_length;
};

#define gt_gff3_linesorted_out_stream_cast(GS)\
        gt_node_stream_cast(gt_gff3_linesorted_out_stream_class(), GS);

static int gt_linesorted_gff3_cmp(const void *val1, const void *val2,
                                  GT_UNUSED void *data)
{
  GtUword p1s = 0, p1e = 0,
          p2s = 0, p2e = 0;
  GT_UNUSED int p1scanned = 0,
                p2scanned = 0;
  const char *s1 = *(const char**) val1;
  const char *s2 = *(const char**) val2;

  if (s1[0] == '#' || s2[0] == '\0')
    return 1;
  if (s2[0] == '#' || s1[0] == '\0')
    return -1;

  p1scanned = sscanf(s1, "%*s\t%*s\t%*s\t"GT_WU"\t"GT_WU, &p1s, &p1e);
  gt_assert(p1s != 0);
  gt_assert(p1scanned == 2);
  p2scanned = sscanf(s2, "%*s\t%*s\t%*s\t"GT_WU"\t"GT_WU, &p2s, &p2e);
  gt_assert(p2s != 0);
  gt_assert(p2scanned == 2);

  if (p1s == p2s)
    return 0;
  if (p1s > p2s)
    return 1;
  else
    return -1;
}

static int gff3_linesorted_out_stream_process_current_cluster(
                                               GtGFF3LinesortedOutStream *lsos,
                                               GtError *err)
{
  int had_err = 0;
  GtUword i;
  bool shown_sep = false;
  GtUword nof_nodes = gt_array_size(lsos->cur_node_set), nof_lines;
  gt_error_check(err);

  /* do not waste time on empty clusters */
  if (nof_nodes == 0) return had_err;

  /* collect output */
  gt_str_reset(lsos->buf);
  for (i = 0; !had_err && i < nof_nodes; i++) {
    GtGenomeNode *n = *(GtGenomeNode**) gt_array_get(lsos->cur_node_set, i);
    had_err = gt_genome_node_accept(n, lsos->gff3vis, err);
    gt_queue_add(lsos->outqueue, n);
  }

  /* split buffered lines */
  gt_splitter_split(lsos->splitter, gt_str_get(lsos->buf),
                    gt_str_length(lsos->buf), '\n');
  nof_lines = gt_splitter_size(lsos->splitter);

  /* sort lines by start positions */
  lsos->outstrings = gt_splitter_get_tokens(lsos->splitter);
  gt_qsort_r(lsos->outstrings, nof_lines, sizeof (char*),
             NULL, gt_linesorted_gff3_cmp);

  /* output */
  for (i = 0; i < nof_lines-1; i++) {
    if (strlen(lsos->outstrings[i]) > 0) {
      if (strcmp(GT_GFF_TERMINATOR, lsos->outstrings[i]) == 0) {
        if (shown_sep)
          continue;
        else
          shown_sep = true;
      }
      gt_file_xprintf(lsos->outfp, "%s\n", lsos->outstrings[i]);
    }
  }

  /* cleanup */
  gt_splitter_reset(lsos->splitter);
  gt_array_reset(lsos->cur_node_set);
  return had_err;
}

static int gff3_linesorted_out_stream_next(GtNodeStream *ns, GtGenomeNode **gn,
                                           GtError *err)
{
  GtGFF3LinesortedOutStream *lsos;
  int had_err = 0;
  bool complete_cluster = false;
  GtGenomeNode *mygn = NULL;
  GtFeatureNode *fn = NULL;
  gt_error_check(err);
  lsos = gt_gff3_linesorted_out_stream_cast(ns);

  /* if there are still nodes left in the buffer, output them */
  if (gt_queue_size(lsos->outqueue) > 0) {
    *gn = (GtGenomeNode*) gt_queue_get(lsos->outqueue);
    return had_err;
  } else complete_cluster = false;

  while (!had_err && !complete_cluster) {
    had_err = gt_node_stream_next(lsos->in_stream, &mygn, err);

    /* stop if stream is at the end */
    if (had_err || !mygn) {
      /* do not forget to finish last cluster */
      if (!had_err) {
        had_err = gff3_linesorted_out_stream_process_current_cluster(lsos, err);
      }
      break;
    }

    /* process all feature nodes */
    if ((fn = gt_feature_node_try_cast(mygn))) {
      GtGenomeNode *addgn;
      GtRange new_rng = gt_genome_node_get_range(mygn);
      if (gt_array_size(lsos->cur_node_set) == 0UL) {
        /* new overlapping node cluster */
        addgn = gt_genome_node_ref(mygn);
        gt_array_add(lsos->cur_node_set, addgn);
        lsos->cur_node_range = gt_genome_node_get_range(mygn);
      } else {
        if (gt_range_overlap(&new_rng, &lsos->cur_node_range)) {
          /* node overlaps with current one, add to cluster */
          addgn = gt_genome_node_ref(mygn);
          gt_array_add(lsos->cur_node_set, addgn);
          lsos->cur_node_range = gt_range_join(&lsos->cur_node_range, &new_rng);
        } else {
          /* finish current cluster and start a new one */
          had_err = gff3_linesorted_out_stream_process_current_cluster(lsos,
                                                                       err);
          if (!had_err) {
            gt_assert(gt_array_size(lsos->cur_node_set) == 0);
            addgn = gt_genome_node_ref(mygn);
            gt_array_add(lsos->cur_node_set, addgn);
            lsos->cur_node_range = gt_genome_node_get_range(mygn);
          }
          if (gt_queue_size(lsos->outqueue) > 0) {
            *gn = (GtGenomeNode*) gt_queue_get(lsos->outqueue);
            complete_cluster = true;
          }
        }
      }
      /* from now on, nodes are kept in clusters only */
      gt_genome_node_delete(mygn);
    } else {
      /* other nodes */
      had_err = gff3_linesorted_out_stream_process_current_cluster(lsos, err);
      if (!had_err) {
        gt_str_reset(lsos->buf);
        had_err = gt_genome_node_accept(mygn, lsos->gff3vis, err);
      }
      if (!had_err) {
        gt_file_xprintf(lsos->outfp, "%s", gt_str_get(lsos->buf));
        gt_queue_add(lsos->outqueue, mygn);
      }
      if (gt_queue_size(lsos->outqueue) > 0) {
        *gn = (GtGenomeNode*) gt_queue_get(lsos->outqueue);
        complete_cluster = true;
      }
    }
  }
  return had_err;
}

static void gff3_linesorted_out_stream_free(GtNodeStream *ns)
{
  GtUword i;
  GtGFF3LinesortedOutStream *lsos;
  if (!ns) return;
  lsos = gt_gff3_linesorted_out_stream_cast(ns);
  while (gt_queue_size(lsos->outqueue) > 0) {
    gt_genome_node_delete((GtGenomeNode*) gt_queue_get(lsos->outqueue));
  }
  gt_queue_delete(lsos->outqueue);
  for (i = 0; i < gt_array_size(lsos->cur_node_set); i++) {
    gt_genome_node_delete(*(GtGenomeNode**)
                                           gt_array_get(lsos->cur_node_set, i));
  }
  gt_node_stream_delete(lsos->in_stream);
  gt_str_delete(lsos->buf);
  gt_node_visitor_delete(lsos->gff3vis);
  gt_splitter_delete(lsos->splitter);
  gt_array_delete(lsos->cur_node_set);
}

const GtNodeStreamClass* gt_gff3_linesorted_out_stream_class(void)
{
  static const GtNodeStreamClass *nsc = NULL;
  gt_class_alloc_lock_enter();
  if (!nsc) {
    nsc = gt_node_stream_class_new(sizeof (GtGFF3LinesortedOutStream),
                                   gff3_linesorted_out_stream_free,
                                   gff3_linesorted_out_stream_next);
  }
  gt_class_alloc_lock_leave();
  return nsc;
}

GtNodeStream* gt_gff3_linesorted_out_stream_new(GtNodeStream *in_stream,
                                                GtFile *outfp)
{
  GtGFF3LinesortedOutStream *lsos;
  GtNodeStream *ns;
  gt_assert(in_stream);
  ns = gt_node_stream_create(gt_gff3_linesorted_out_stream_class(), true);
  lsos = gt_gff3_linesorted_out_stream_cast(ns);
  lsos->cur_node_set = gt_array_new(sizeof (GtFeatureNode*));
  lsos->in_stream = gt_node_stream_ref(in_stream);
  lsos->cur_node_range.start = lsos->cur_node_range.end = GT_UNDEF_UWORD;
  lsos->outqueue = gt_queue_new();
  lsos->outfp = outfp;
  lsos->outstrings = NULL;
  lsos->outstrings_length = 0;
  lsos->splitter = gt_splitter_new();
  lsos->buf = gt_str_new();
  lsos->gff3vis = gt_gff3_visitor_new_to_str(lsos->buf);
  return ns;
}

void gt_gff3_linesorted_out_stream_set_fasta_width(GtGFF3LinesortedOutStream
                                                               *gff3_out_stream,
                                                   GtUword fasta_width)
{
  gt_assert(gff3_out_stream);
  gt_gff3_visitor_set_fasta_width((GtGFF3Visitor*)
                                  gff3_out_stream->gff3vis, fasta_width);
}

void gt_gff3_linesorted_out_stream_retain_id_attributes(
                                     GtGFF3LinesortedOutStream *gff3_out_stream)
{
  gt_assert(gff3_out_stream);
  gt_gff3_visitor_retain_id_attributes((GtGFF3Visitor*)
                                       gff3_out_stream->gff3vis);
}