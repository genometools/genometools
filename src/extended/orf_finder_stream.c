/*
  Copyright (c) 2010      Sascha Kastens <mail@skastens.de>
  Copyright (c) 2011-2014 Sascha Steinbiss <ss34@sanger.ac.uk>
  Copyright (c) 2010-2011 Center for Bioinformatics, University of Hamburg
  Copyright (c)      2014 Genome Research Ltd.

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
#include "core/codon_iterator_api.h"
#include "core/codon_iterator_simple_api.h"
#include "core/cstr.h"
#include "core/log.h"
#include "core/ma.h"
#include "core/mathsupport.h"
#include "core/range.h"
#include "core/seq_iterator_sequence_buffer_api.h"
#include "core/str.h"
#include "core/symbol_api.h"
#include "core/trans_table.h"
#include "core/translator.h"
#include "core/unused_api.h"
#include "extended/node_stream_api.h"
#include "extended/feature_node.h"
#include "extended/feature_node_iterator_api.h"
#include "extended/reverse_api.h"
#include "extended/orf_finder_stream.h"
#include "extended/orf_finder_visitor.h"
#include "extended/orf_iterator_api.h"

struct GtORFFinderStream {
  const GtNodeStream parent_instance;
  GtNodeStream *in_stream;
  GtCodonIterator* ci;
  GtTranslator* translator;
  GtORFIterator* orfi;
  GtSeqIterator *seqit;
  GtStr *desc_str;
  const GtUchar *seq;
  char *seqr, *desc;
  const char *type;
  GtUword seqlen;
  unsigned int min, max;
  bool reverse, all;
  GtORFFinderVisitor *lv;
};

#define gt_orf_finder_stream_cast(GS)\
        gt_node_stream_cast(gt_orf_finder_stream_class(), GS)

static int gt_orf_finder_stream_next_v(GtNodeStream *gs, GtGenomeNode **gn,
                                       GtError *err)
{
  GtORFFinderStream *ls;
  int had_err;

  gt_error_check(err);
  ls = gt_orf_finder_stream_cast(gs);

  had_err = gt_node_stream_next(ls->in_stream, gn, err);
  if (!had_err && *gn) {
    had_err = gt_genome_node_accept(*gn, (GtNodeVisitor*) ls->lv, err);
  }
  if (had_err) {
    gt_genome_node_delete(*gn);
    *gn = NULL;
  }
  return had_err;
}

static int gt_orf_finder_stream_next_s(GtNodeStream *gs, GtGenomeNode **gn,
                                       GtError *err)
{
  GtORFFinderStream *ls;
  int had_err = 0;
  bool found = false;
  GtRange orf_rng;
  GtORFIteratorStatus state;
  unsigned int orf_frame;

  gt_error_check(err);
  ls = gt_orf_finder_stream_cast(gs);
  gt_assert(ls->seqit);

  while (!had_err && !found) {
    while ((state = gt_orf_iterator_next(ls->orfi, &orf_rng, &orf_frame,
                                                 err)) == GT_ORF_ITERATOR_END) {
      /* no more ORFs, try next sequence  */
      if (!ls->reverse) {
        gt_free(ls->seqr);
        ls->seqr = gt_cstr_dup((const char*) ls->seq);
        had_err = gt_reverse_complement(ls->seqr, ls->seqlen, err);
        ls->seq = (const GtUchar*) ls->seqr;
        ls->reverse = true;
      } else {
        int rval;
        rval = gt_seq_iterator_next(ls->seqit, &ls->seq, &ls->seqlen,
                                    &ls->desc, err);
        if (rval < 0) {
          had_err = -1;
          return had_err;
        }
        if (rval == 0) {
          *gn = NULL;
          return 0;
        }
        gt_assert(rval == 1);
        ls->reverse = false;
      }
      if (!had_err) {
        gt_assert(ls->seq && ls->seqlen > 0);
        if (ls->desc) {
          gt_str_reset(ls->desc_str);
          gt_str_append_cstr(ls->desc_str, ls->desc);
        }
        gt_codon_iterator_delete(ls->ci);
        ls->ci = gt_codon_iterator_simple_new((char*) ls->seq, ls->seqlen, err);
        gt_assert(ls->ci);
        gt_translator_delete(ls->translator);
        ls->translator = gt_translator_new(ls->ci);
        gt_assert(ls->translator);
        gt_orf_iterator_delete(ls->orfi);
        ls->orfi = gt_orf_iterator_new(ls->ci, ls->translator);
        if (ls->all)
          gt_orf_iterator_find_all(ls->orfi);
        gt_assert(ls->orfi);
      }
    }
    if (!had_err && state == GT_ORF_ITERATOR_OK) {
        GtUword rlen = gt_range_length(&orf_rng);
        GtRange orig_rng = orf_rng;
        GtStrand new_strand;
        new_strand = ls->reverse ? GT_STRAND_REVERSE : GT_STRAND_FORWARD;
        if (ls->reverse) {
          orf_rng.end = ls->seqlen - 1 - orig_rng.start;
          orf_rng.start = ls->seqlen - 1 - orig_rng.end;
        }
        if (rlen <= ls->max && rlen >= ls->min) {
          *gn = gt_feature_node_new(ls->desc_str, gt_symbol(ls->type),
                                    orf_rng.start + 1,
                                    orf_rng.end + 1,
                                    new_strand);
          found = true;
        }
    }

    if (!had_err && state == GT_ORF_ITERATOR_ERROR)
      had_err = -1;
  }
  return had_err;
}

static int gt_orf_finder_stream_next(GtNodeStream *gs, GtGenomeNode **gn,
                                     GtError *err)
{
  GtORFFinderStream *ls;
  gt_assert(gs);
  ls = gt_orf_finder_stream_cast(gs);

  if (ls->seqit)
    return gt_orf_finder_stream_next_s(gs, gn, err);
  else
    return gt_orf_finder_stream_next_v(gs, gn, err);
}

static void gt_orf_finder_stream_free(GtNodeStream *gs)
{
  GtORFFinderStream *ls = gt_orf_finder_stream_cast(gs);
  gt_node_visitor_delete((GtNodeVisitor*) ls->lv);
  gt_node_stream_delete(ls->in_stream);
  gt_translator_delete(ls->translator);
  gt_orf_iterator_delete(ls->orfi);
  gt_codon_iterator_delete(ls->ci);
  gt_str_delete(ls->desc_str);
  gt_free(ls->seqr);
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
                                       GtRegionMapping *rmap,
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
  ls->lv = (GtORFFinderVisitor*) gt_orf_finder_visitor_new(rmap,types,
                                                           min, max, all, err);
  ls->seqit = NULL;
  return gs;
}

void gt_orf_finder_stream_set_type(GtORFFinderStream *orfs, const char *type)
{
  gt_assert(orfs && type);
  orfs->type = type;
}

GtNodeStream* gt_orf_finder_stream_new_from_seq(GtSeqIterator *seqit,
                                                unsigned int min,
                                                unsigned int max,
                                                bool all,
                                                GtError *err)
{
  GtNodeStream *gs;
  GtORFFinderStream *ls;
  int rval;
  gs = gt_node_stream_create(gt_orf_finder_stream_class(), false);
  ls = gt_orf_finder_stream_cast(gs);
  ls->in_stream = NULL;
  ls->lv = NULL;
  ls->seqit = seqit;
  ls->ci = NULL;
  ls->orfi = NULL;
  ls->min = min;
  ls->max = max;
  ls->seqr = NULL;
  ls->reverse = false;
  ls->all = all;
  ls->type = "open_reading_frame";
  ls->desc_str = gt_str_new();
  while (true) {
    rval = gt_seq_iterator_next(ls->seqit, &ls->seq, &ls->seqlen,
                                &ls->desc, err);
    if (rval < 0) {
      gt_node_stream_delete(gs);
      return NULL;
    }
    if (rval == 0)
      break;
    if (ls->desc) {
      gt_str_reset(ls->desc_str);
      gt_str_append_cstr(ls->desc_str, ls->desc);
    }
    if (rval == 1) {
      if (ls->seqlen < 3) {
        continue;
      } else {
        break;
      }
    }
  }
  if (rval == 0) {
    gt_error_set(err, "no translatable sequence given");
    gt_node_stream_delete(gs);
    return NULL;
  }
  ls->ci = gt_codon_iterator_simple_new((char*) ls->seq, ls->seqlen, err);
  gt_assert(ls->ci);
  ls->translator = gt_translator_new(ls->ci);
  ls->orfi = gt_orf_iterator_new(ls->ci, ls->translator);
  if (ls->all)
    gt_orf_iterator_find_all(ls->orfi);
  gt_assert(ls->orfi);

  return gs;
}
