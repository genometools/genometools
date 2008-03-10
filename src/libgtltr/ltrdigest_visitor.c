/*
  Copyright (c) 2008 Sascha Steinbiss <ssteinbiss@zbh.uni-hamburg.de>
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

#include <assert.h>
#include "libgtcore/unused.h"
#include "libgtcore/range.h"
#include "libgtext/genome_feature_type.h"
#include "libgtext/genome_visitor_rep.h"
#include "libgtext/sequence_region.h"
#include "libgtltr/ltrdigest_visitor.h"

struct LTRdigestVisitor {
  const GenomeVisitor parent_instance;
  LTRElement *element;
};

#define ltrdigest_visitor_cast(GV)\
        genome_visitor_cast(ltrdigest_visitor_class(), GV)

static int ltrdigest_visitor_genome_feature(GenomeVisitor *gv, GenomeFeature *gf,
                                            UNUSED Error *err)
{
  LTRdigestVisitor *lv = ltrdigest_visitor_cast(gv);
  GenomeFeatureType gft;
  Range node_range;
  assert(lv);
  error_check(err);

  gft = genome_feature_get_type(gf);

  switch(gft)
  {
    case gft_LTR_retrotransposon:
      lv->element->mainnode = gf;
      break;
    case gft_long_terminal_repeat:
      /* XXX: check order if unsorted! */
      if(lv->element->leftLTR == NULL)
      {
        node_range = genome_node_get_range((GenomeNode*) gf);
        lv->element->leftLTR = gf;
        lv->element->leftLTR_5 = node_range.start;
        lv->element->leftLTR_3 = node_range.end;
      }
      else
      {
        node_range = genome_node_get_range((GenomeNode*) gf);
        lv->element->rightLTR = gf;
        lv->element->rightLTR_5 = node_range.start;
        lv->element->rightLTR_3 = node_range.end;
      }
      break;
    case gft_target_site_duplication:
      /* XXX: check order if unsorted! */
      if(lv->element->leftTSD == NULL)
      {
        node_range = genome_node_get_range((GenomeNode*) gf);
        lv->element->leftTSD = gf;
      }
      else
      {
        node_range = genome_node_get_range((GenomeNode*) gf);
        lv->element->rightTSD = gf;
      }
      break;
    default:
      break;
  }
  return 0;
}

int ltrdigest_visitor_sequence_region(GenomeVisitor *gv,
                                             SequenceRegion *sr,
                                             UNUSED Error *err)
{
  LTRdigestVisitor *v = ltrdigest_visitor_cast(gv);
  assert(v && sr);
  error_check(err);

  return 0;
}

const GenomeVisitorClass* ltrdigest_visitor_class()
{
  static const GenomeVisitorClass gvc = { sizeof (LTRdigestVisitor),
                                          NULL,
                                          NULL,
                                          ltrdigest_visitor_genome_feature,
                                          NULL };
  return &gvc;
}

GenomeVisitor* ltrdigest_visitor_new(LTRElement *element)
{
  GenomeVisitor *gv;
  LTRdigestVisitor *lv;
  assert(element);
  gv = genome_visitor_create(ltrdigest_visitor_class());
  lv = ltrdigest_visitor_cast(gv);
  lv->element = element;
  assert(lv);
  return gv;
}
