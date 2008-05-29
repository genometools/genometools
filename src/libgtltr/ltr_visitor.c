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
#include "libgtcore/cstr.h"
#include "libgtcore/hashtable.h"
#include "libgtcore/ma.h"
#include "libgtcore/range.h"
#include "libgtcore/unused.h"
#include "libgtext/genome_feature_type.h"
#include "libgtext/genome_visitor_rep.h"
#include "libgtext/sequence_region.h"
#include "libgtltr/ltr_visitor.h"

struct LTRVisitor {
  const GenomeVisitor parent_instance;
  LTRElement *element;
};

#define ltr_visitor_cast(GV)\
        genome_visitor_cast(ltr_visitor_class(), GV)

static int ltr_visitor_genome_feature(GenomeVisitor *gv, GenomeFeature *gf,
                                      UNUSED Error *err)
{
  LTRVisitor *lv;
  Range node_range;
  Array *pdomarr = NULL;
  const char* pfamname;
  lv = ltr_visitor_cast(gv);
  assert(lv);
  error_check(err);

  switch (genome_feature_get_type(gf))
  {
    case gft_LTR_retrotransposon:
      lv->element->mainnode = gf;
      break;
    case gft_long_terminal_repeat:
      /* XXX: check order if unsorted! */
      if (lv->element->leftLTR == NULL)
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
      if (lv->element->leftTSD == NULL)
      {
        lv->element->leftTSD = gf;
      }
      else
      {
        lv->element->rightTSD = gf;
      }
      break;
    case gft_RR_tract:
      if (lv->element->ppt == NULL)
      {
        lv->element->ppt = gf;
      }
      break;
    case gft_primer_binding_site:
      if (lv->element->pbs == NULL)
      {
        lv->element->pbs = gf;
      }
      break;
    case gft_protein_match:
      if (!lv->element->pdoms)
      {
        lv->element->pdoms = hashtable_new(HASH_STRING,
                                           ma_free_func,
                                           (FreeFunc) array_delete);
      }
      pfamname = genome_feature_get_attribute((GenomeNode*) gf,
                                              "pfamname");
      if (!(pdomarr = (Array*) hashtable_get(lv->element->pdoms, pfamname)))
      {
        char *pfamcpy = cstr_dup(pfamname);
        pdomarr = array_new(sizeof (GenomeFeature*));
        hashtable_add(lv->element->pdoms, pfamcpy, pdomarr);
        array_add(lv->element->pdomorder, pfamcpy);
      }
      array_add(pdomarr, gf);
      break;
    default:
      break;
  }
  return 0;
}

int ltr_visitor_sequence_region(GenomeVisitor *gv,
                                SequenceRegion *sr,
                                UNUSED Error *err)
{
  LTRVisitor *v = ltr_visitor_cast(gv);
  assert(v && sr);
  error_check(err);

  return 0;
}

const GenomeVisitorClass* ltr_visitor_class()
{
  static const GenomeVisitorClass gvc = { sizeof (LTRVisitor),
                                          NULL,
                                          NULL,
                                          ltr_visitor_genome_feature,
                                          NULL };
  return &gvc;
}

GenomeVisitor* ltr_visitor_new(LTRElement *element)
{
  GenomeVisitor *gv;
  LTRVisitor *lv;
  assert(element);
  gv = genome_visitor_create(ltr_visitor_class());
  lv = ltr_visitor_cast(gv);
  lv->element = element;
  assert(lv);
  return gv;
}
