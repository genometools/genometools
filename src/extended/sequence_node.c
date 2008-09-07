/*
  Copyright (c) 2008 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2008 Center for Bioinformatics, University of Hamburg

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

#include "core/unused.h"
#include "extended/sequence_node.h"
#include "extended/genome_node_rep.h"

struct SequenceNode
{
  const GT_GenomeNode parent_instance;
  Str *description,
      *sequence;
};

#define sequence_node_cast(GN)\
        gt_genome_node_cast(sequence_node_class(), GN)

static void sequence_node_free(GT_GenomeNode *gn)
{
  SequenceNode *sn = sequence_node_cast(gn);
  str_delete(sn->sequence);
  str_delete(sn->description);
}

static Str* sequence_node_get_seqid(GT_GenomeNode *gn)
{
  SequenceNode *sn = sequence_node_cast(gn);
  return sn->description;
}

static GT_Range sequence_node_get_range(UNUSED GT_GenomeNode *gn)
{
  GT_Range range;
  range.start = 0;
  range.end = 0;
  return range;
}

static void sequence_node_change_seqid(GT_GenomeNode *gn, Str *seqid)
{
  SequenceNode *sn = sequence_node_cast(gn);
  assert(sn && seqid);
  str_delete(sn->description);
  sn->description = str_ref(seqid);
}

static int sequence_node_accept(GT_GenomeNode *gn, GenomeVisitor *gv, GT_Error *err)
{
  SequenceNode *sn;
  gt_error_check(err);
  sn = sequence_node_cast(gn);
  return genome_visitor_visit_sequence_node(gv, sn, err);
}

const GT_GenomeNodeClass* sequence_node_class()
{
  static const GT_GenomeNodeClass gnc = { sizeof (SequenceNode),
                                       sequence_node_free,
                                       sequence_node_get_seqid,
                                       sequence_node_get_seqid,
                                       sequence_node_get_range,
                                       NULL,
                                       sequence_node_change_seqid,
                                       sequence_node_accept };
  return &gnc;
}

GT_GenomeNode* sequence_node_new(const char *description, Str *sequence)
{
  GT_GenomeNode *gn = gt_genome_node_create(sequence_node_class());
  SequenceNode *sn = sequence_node_cast(gn);
  assert(description && sequence);
  sn->description = str_new_cstr(description);
  sn->sequence = sequence;
  return gn;
}

const char* sequence_node_get_description(const SequenceNode *sn)
{
  assert(sn);
  return str_get(sn->description);
}

const char* sequence_node_get_sequence(const SequenceNode *sn)
{
  assert(sn);
  return str_get(sn->sequence);
}

unsigned long sequence_node_get_sequence_length(const SequenceNode *sn)
{
  assert(sn);
  return str_length(sn->sequence);
}
