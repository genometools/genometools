/*
  Copyright (c) 2014-2015 Genome Research Ltd.
  Copyright (c) 2022 Sascha Steinbiss <sascha@steinbiss.name>

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
#include "core/ma_api.h"
#include "core/minmax_api.h"
#include "core/parseutils_api.h"
#include "core/qsort_r_api.h"
#include "core/splitter_api.h"
#include "core/symbol_api.h"
#include "core/undef_api.h"
#include "core/unused_api.h"
#include "extended/feature_node.h"
#include "extended/feature_node_iterator_api.h"
#include "extended/feature_type_api.h"
#include "extended/gff3_defines.h"
#include "extended/gff3_linesorted_out_stream.h"
#include "extended/gff3_visitor.h"

#define GT_LINESORTED_SEP '\t'

struct GtGFF3LinesortedOutStream {
  const GtNodeStream parent_instance;
  GtNodeStream *in_stream;
  GtArray *cur_node_set;
  GtFile *outfp;
  GtStr *buf;
  GtSplitter *splitter;
  GtNodeVisitor *gff3vis;
  GtStr *last_seqid;
  char **outstrings;
  GtUword outstrings_length;
};

#define gt_gff3_linesorted_out_stream_cast(GS)\
        gt_node_stream_cast(gt_gff3_linesorted_out_stream_class(), GS);

/* this is parsing self-generated content, so we reasonably expect enough
   fields -- fail assertion if there is a parsing problem */
static inline void gt_linesorted_gff3_next_token_pair(const char **s1,
                                                      const char **s2,
                                                      char *buf1, size_t len1,
                                                      char *buf2, size_t len2,
                                                      char sep)
{
  const char *tokend1, *tokend2;
  size_t idlength1, idlength2;

  tokend1 = strchr(*s1, sep);
  gt_assert(tokend1 && tokend1 > *s1);
  tokend2 = strchr(*s2, sep);
  gt_assert(tokend2 && tokend2 > *s2);
  idlength1 = GT_MIN((tokend1 - *s1), len1 - 1);
  gt_assert(idlength1 > 0);
  idlength2 = GT_MIN((tokend2 - *s2), len2 - 1);
  gt_assert(idlength2 > 0);
  (void) strncpy(buf1, *s1, idlength1);
  (void) strncpy(buf2, *s2, idlength2);
  buf1[idlength1] = '\0';
  buf2[idlength2] = '\0';
  *s1 = tokend1 + 1;
  *s2 = tokend2 + 1;
}

/* this is comparing self-generated content, so we reasonably expect enough
   fields -- fail assertion if there is a parsing problem */
static int gt_linesorted_gff3_cmp(const void *val1, const void *val2,
                                  GT_UNUSED void *data)
{
  GtUword p1 = 0, p2 = 0;
  const char *s1 = *(const char**) val1,
             *s2 = *(const char**) val2,
             *reststart1 = s1, *reststart2 = s2;
  int str_cmp_result,
      GT_UNUSED rval = 0;
  char buf1[BUFSIZ], buf2[BUFSIZ];

  if (s1[0] == '#' || s2[0] == '\0')
    return 1;
  if (s2[0] == '#' || s1[0] == '\0')
    return -1;

  /* parse seqid */
  gt_linesorted_gff3_next_token_pair(&reststart1, &reststart2, buf1, BUFSIZ,
                                     buf2, BUFSIZ, GT_LINESORTED_SEP);
  str_cmp_result = strcmp(buf1, buf2);
  if (str_cmp_result != 0)
    return str_cmp_result;

  /* parse source */
  gt_linesorted_gff3_next_token_pair(&reststart1, &reststart2, buf1, BUFSIZ,
                                     buf2, BUFSIZ, GT_LINESORTED_SEP);
  /* parse type */
  gt_linesorted_gff3_next_token_pair(&reststart1, &reststart2, buf1, BUFSIZ,
                                     buf2, BUFSIZ, GT_LINESORTED_SEP);

  /* parse start */
  gt_linesorted_gff3_next_token_pair(&reststart1, &reststart2, buf1, BUFSIZ,
                                     buf2, BUFSIZ, GT_LINESORTED_SEP);
  rval = gt_parse_uword(&p1, buf1);
  gt_assert(rval == 0);
  rval = gt_parse_uword(&p2, buf2);
  gt_assert(rval == 0);

  if (p1 == p2) {
    /* parse end */
    gt_linesorted_gff3_next_token_pair(&reststart1, &reststart2, buf1, BUFSIZ,
                                       buf2, BUFSIZ, GT_LINESORTED_SEP);
    rval = gt_parse_uword(&p1, buf1);
    gt_assert(rval == 0);
    rval = gt_parse_uword(&p2, buf2);
    gt_assert(rval == 0);
    if (p1 == p2)
      return 0;
    if (p1 > p2)
      return 1;
    else
      return -1;
  }
  if (p1 > p2)
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

  /* collect and free output */
  gt_str_reset(lsos->buf);
  for (i = 0; !had_err && i < nof_nodes; i++) {
    GtGenomeNode *n = *(GtGenomeNode**) gt_array_get(lsos->cur_node_set, i);
    had_err = gt_genome_node_accept(n, lsos->gff3vis, err);
    gt_genome_node_delete(n);
  }
  gt_array_reset(lsos->cur_node_set);

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
  return had_err;
}

static int gff3_linesorted_out_stream_next(GtNodeStream *ns, GtGenomeNode **gn,
                                           GtError *err)
{
  GtGFF3LinesortedOutStream *lsos;
  int had_err = 0;
  GtGenomeNode *node;
  GtFeatureNode *fn = NULL;
  gt_error_check(err);
  lsos = gt_gff3_linesorted_out_stream_cast(ns);

  /* we do not need to pass on any nodes, this is an output stream */
  *gn = NULL;

  while (!(had_err = gt_node_stream_next(lsos->in_stream, &node,
                                          err))) {
     if (node == NULL) {
       /* do not forget to finish last cluster */
       if (!had_err) {
         had_err = gff3_linesorted_out_stream_process_current_cluster(lsos, err);
       }
       break;
     }
     if ((fn = gt_feature_node_try_cast(node))) {
       if (lsos->last_seqid == NULL) {
         lsos->last_seqid = gt_str_clone(gt_genome_node_get_seqid(node));
       }
       if (gt_str_cmp(lsos->last_seqid, gt_genome_node_get_seqid(node)) != 0) {
         /* new sequence reached, finish old cluster and update current seqid;
            we can do this as GenomeTools enforces all features within a CC to
            share the same seqid */
         had_err = gff3_linesorted_out_stream_process_current_cluster(lsos, err);
         if (had_err) {
           return had_err;
         }
         gt_str_reset(lsos->last_seqid);
         gt_str_append_str(lsos->last_seqid, gt_genome_node_get_seqid(node));
       }
       gt_array_add(lsos->cur_node_set, fn);
     } else {
       if (!had_err) {
        gt_str_reset(lsos->buf);
        had_err = gt_genome_node_accept(node, lsos->gff3vis, err);
      }
      if (!had_err) {
        gt_file_xprintf(lsos->outfp, "%s", gt_str_get(lsos->buf));
        gt_genome_node_delete(node);
      } else {
        return had_err;
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
  for (i = 0; i < gt_array_size(lsos->cur_node_set); i++) {
    gt_genome_node_delete(*(GtGenomeNode**)
                                           gt_array_get(lsos->cur_node_set, i));
  }
  gt_node_stream_delete(lsos->in_stream);
  gt_str_delete(lsos->buf);
  gt_str_delete(lsos->last_seqid);
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
  lsos->outfp = outfp;
  lsos->outstrings = NULL;
  lsos->outstrings_length = 0;
  lsos->last_seqid = NULL;
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
