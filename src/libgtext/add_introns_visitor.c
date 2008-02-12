/*
  Copyright (c) 2007-2008 Gordon Gremme <gremme@zbh.uni-hamburg.de>
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

#include <assert.h>
#include "libgtcore/undef.h"
#include "libgtcore/unused.h"
#include "libgtext/add_introns_visitor.h"
#include "libgtext/genome_visitor_rep.h"

struct AddIntronsVisitor {
  const GenomeVisitor parent_instance;
  GenomeFeature *parent_feature,
                *previous_exon_feature;
};

#define add_introns_visitor_cast(GV)\
        genome_visitor_cast(add_introns_visitor_class(), GV)

static int add_introns_in_children(GenomeNode *gn, void *data,
                                   UNUSED Error *err)
{
  AddIntronsVisitor *v = (AddIntronsVisitor*) data;
  GenomeFeature *current_feature;
  GenomeNode *intron_node;
  Range previous_range, current_range, intron_range;
  Strand previous_strand, current_strand, intron_strand;
  Str *parent_seqid;
  error_check(err);
  current_feature = genome_node_cast(genome_feature_class(), gn);
  assert(current_feature);
  if (genome_feature_get_type(current_feature) == gft_exon) {
    if (v->previous_exon_feature) {
      /* determine intron range */
      previous_range = genome_node_get_range((GenomeNode*)
                                             v->previous_exon_feature);
      current_range = genome_node_get_range(gn);
      assert(previous_range.end < current_range.start);
      intron_range.start = previous_range.end + 1;
      intron_range.end = current_range.start - 1;

      /* determine intron strand */
      previous_strand = genome_feature_get_strand(v->previous_exon_feature);
      current_strand = genome_feature_get_strand(current_feature);
      assert(previous_strand == current_strand);
      intron_strand = previous_strand;

      /* determine sequence id */
      parent_seqid = genome_node_get_seqid((GenomeNode*) v->parent_feature);
      assert(!str_cmp(parent_seqid,
             genome_node_get_seqid((GenomeNode*) v->previous_exon_feature)));
      assert(!str_cmp(parent_seqid,
             genome_node_get_seqid((GenomeNode*) current_feature)));

      /* create intron */
      intron_node = genome_feature_new(gft_intron, intron_range, intron_strand,
                                       NULL, UNDEF_ULONG);
      genome_node_set_seqid(intron_node, parent_seqid);
      genome_node_is_part_of_genome_node((GenomeNode*) v->parent_feature,
                                         intron_node);
    }
    v->previous_exon_feature = current_feature;
  }
  return 0;
}

static int add_introns_if_necessary(GenomeNode *gn, void *data, Error *e)
{
  AddIntronsVisitor *v = (AddIntronsVisitor*) data;
  GenomeFeature *gf;
  error_check(e);
  gf = genome_node_cast(genome_feature_class(), gn);
  assert(gf);
  v->parent_feature = gf;
  v->previous_exon_feature = NULL;
  return genome_node_traverse_direct_children(gn, v, add_introns_in_children,
                                              e);
}

static int add_introns_visitor_genome_feature(GenomeVisitor *gv,
                                             GenomeFeature *gf, Error *e)
{
  AddIntronsVisitor *v;
  error_check(e);
  v = add_introns_visitor_cast(gv);
  return genome_node_traverse_children((GenomeNode*) gf, v,
                                       add_introns_if_necessary, false, e);
}

const GenomeVisitorClass* add_introns_visitor_class()
{
  static const GenomeVisitorClass gvc = { sizeof (AddIntronsVisitor),
                                          NULL,
                                          NULL,
                                          add_introns_visitor_genome_feature,
                                          NULL };
  return &gvc;
}

GenomeVisitor* add_introns_visitor_new(void)
{
  return genome_visitor_create(add_introns_visitor_class());
}
