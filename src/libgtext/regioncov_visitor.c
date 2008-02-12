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
#include "libgtcore/cstr.h"
#include "libgtcore/hashtable.h"
#include "libgtcore/ma.h"
#include "libgtcore/minmax.h"
#include "libgtcore/unused.h"
#include "libgtext/genome_visitor_rep.h"
#include "libgtext/regioncov_visitor.h"

struct RegionCovVisitor {
  const GenomeVisitor parent_instance;
  unsigned long max_feature_dist;
  Hashtable *region2rangelist;
};

#define regioncov_visitor_cast(GV)\
        genome_visitor_cast(regioncov_visitor_class(), GV)

static void regioncov_visitor_free(GenomeVisitor *gv)
{
  RegionCovVisitor *regioncov_visitor = regioncov_visitor_cast(gv);
  hashtable_delete(regioncov_visitor->region2rangelist);
}

static int regioncov_visitor_genome_feature(GenomeVisitor *gv,
                                            GenomeFeature *gf,
                                            UNUSED Error *err)
{
  Range *old_range_ptr, old_range, new_range;
  Array *ranges;
  RegionCovVisitor *regioncov_visitor;
  error_check(err);
  regioncov_visitor = regioncov_visitor_cast(gv);
  ranges = hashtable_get(regioncov_visitor->region2rangelist,
                         str_get(genome_node_get_seqid((GenomeNode*) gf)));
  assert(ranges);
  new_range = genome_node_get_range((GenomeNode*) gf);
  if (!array_size(ranges))
    array_add(ranges, new_range);
  else {
    old_range_ptr = array_get_last(ranges);
    old_range = *old_range_ptr;
    old_range.end += regioncov_visitor->max_feature_dist;
    if (range_overlap(old_range, new_range)) {
      old_range_ptr->end = MAX(old_range_ptr->end, new_range.end);
    }
    else
      array_add(ranges, new_range);
  }
  return 0;
}

static int regioncov_visitor_sequence_region(GenomeVisitor *gv,
                                             SequenceRegion *sr,
                                             UNUSED Error *err)
{
  RegionCovVisitor *regioncov_visitor;
  Array *rangelist;
  error_check(err);
  regioncov_visitor = regioncov_visitor_cast(gv);
  rangelist = array_new(sizeof (Range));
  hashtable_add(regioncov_visitor->region2rangelist,
                cstr_dup(str_get(genome_node_get_seqid((GenomeNode*) sr))),
                rangelist);
  return 0;
}

const GenomeVisitorClass* regioncov_visitor_class()
{
  static const GenomeVisitorClass gvc = { sizeof (RegionCovVisitor),
                                          regioncov_visitor_free,
                                          NULL,
                                          regioncov_visitor_genome_feature,
                                          regioncov_visitor_sequence_region };
  return &gvc;
}

GenomeVisitor* regioncov_visitor_new(unsigned long max_feature_dist)
{
  GenomeVisitor *gv = genome_visitor_create(regioncov_visitor_class());
  RegionCovVisitor *regioncov_visitor = regioncov_visitor_cast(gv);
  regioncov_visitor->max_feature_dist = max_feature_dist;
  regioncov_visitor->region2rangelist = hashtable_new(HASH_STRING,
                                                      ma_free_func,
                                                      (FreeFunc) array_delete);
  return gv;
}

static int show_rangelist(void *key, void *value, UNUSED void *data,
                          UNUSED Error *err)
{
  unsigned long i;
  Array *rangelist;
  Range *rangeptr;
  error_check(err);
  assert(key && value);
  rangelist = (Array*) value;
  if (array_size(rangelist)) {
    assert(ranges_are_sorted_and_do_not_overlap(rangelist));
    printf("%s:\n", (char*) key);
    for (i = 0; i < array_size(rangelist); i++) {
      rangeptr = array_get(rangelist, i);
      printf("%lu, %lu\n", rangeptr->start, rangeptr->end);
    }
  }
  return 0;
}

void regioncov_visitor_show_coverage(GenomeVisitor *gv)
{
  RegionCovVisitor *regioncov_visitor = regioncov_visitor_cast(gv);
  int had_err;
  had_err = hashtable_foreach_ao(regioncov_visitor->region2rangelist,
                                 show_rangelist, NULL, NULL);
  assert(!had_err); /* show_rangelist() is sane */
}
