/*
  Copyright (c) 2010-2011 Giorgio Gonnella <gonnella@zbh.uni-hamburg.de>
  Copyright (c) 2010-2011 Center for Bioinformatics, University of Hamburg

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

#ifndef RDJ_STRGRAPH_COUNTS_CHAR_DEF_H
#define RDJ_STRGRAPH_COUNTS_CHAR_DEF_H

#include <limits.h>
#include "core/hashmap-generic.h"

/* --- for exclusive use in rdj-strgraph.c (__ = private) --- */

#define GT_STRGRAPH_COUNTS_REPRESENTATION "char array + hash table"

typedef unsigned char GtStrgraphCount__Small;
typedef unsigned long GtStrgraphCount__Large;
typedef GtStrgraphCount__Large GtStrgraphCount;
#define FormatGtStrgraphCount       "%lu"
#define PRINTGtStrgraphCountcast(X) (X)
#define SCANGtStrgraphCountcast(X)  (X)
#define GT_STRGRAPH_COUNT_MAX       ULONG_MAX

/* hash map for large counts: vertex number -> count */
DECLARE_HASHMAP(GtStrgraphVnum, v, GtStrgraphCount__Large, c_,
                static, inline);
DEFINE_HASHMAP(GtStrgraphVnum, v, GtStrgraphCount__Large, c_,
               gt_ht_ul_elem_hash, gt_ht_ul_elem_cmp, NULL_DESTRUCTOR,
               NULL_DESTRUCTOR, static, inline);

#define GT_STRGRAPH_DECLARE_COUNTS\
  GtStrgraphCount__Small *__small_counts;\
  GtHashtable            *__large_counts;\
  GtStrgraphVnum       __n_counts

#define GT_STRGRAPH__COUNT_IS_LARGE ((GtStrgraphCount__Small)255)

#define GT_STRGRAPH_ALLOC_COUNTS(STRGRAPH, NOFVERTICES)\
  (STRGRAPH)->__small_counts = gt_calloc((size_t)(NOFVERTICES),\
      sizeof (*(STRGRAPH)->__small_counts));\
  (STRGRAPH)->__large_counts = v_c__gt_hashmap_new();\
  (STRGRAPH)->__n_counts = (NOFVERTICES)

#define GT_STRGRAPH_FREE_COUNTS(STRGRAPH)\
  gt_free((STRGRAPH)->__small_counts);\
  (STRGRAPH)->__small_counts = NULL;\
  if ((STRGRAPH)->__large_counts != NULL)\
    v_c__gt_hashmap_delete((STRGRAPH)->__large_counts);\
  (STRGRAPH)->__large_counts = NULL

#define GT_STRGRAPH_GET_COUNT(STRGRAPH, POSITION) \
  (((STRGRAPH)->__small_counts[(POSITION)] < GT_STRGRAPH__COUNT_IS_LARGE)\
    ? (GtStrgraphCount)((STRGRAPH)->__small_counts[(POSITION)]) \
    : (GtStrgraphCount)*v_c__gt_hashmap_get((STRGRAPH)->__large_counts,\
      (POSITION)))

#define GT_STRGRAPH__SET_COUNT(STRGRAPH, POSITION, VALUE) \
  if ((VALUE) < (GtStrgraphCount)GT_STRGRAPH__COUNT_IS_LARGE)\
  {\
    ((STRGRAPH)->__small_counts[(POSITION)] = (GtStrgraphCount__Small)(VALUE));\
  }\
  else\
  {\
    ((STRGRAPH)->__small_counts[(POSITION)] = GT_STRGRAPH__COUNT_IS_LARGE);\
    v_c__gt_hashmap_add((STRGRAPH)->__large_counts, (POSITION), (VALUE));\
  }\

#define GT_STRGRAPH_INC_COUNT(STRGRAPH, POSITION) \
  if ((STRGRAPH)->__small_counts[(POSITION)] < GT_STRGRAPH__COUNT_IS_LARGE -\
      (GtStrgraphCount__Small)1)\
    ((STRGRAPH)->__small_counts[(POSITION)])++; \
  else if ((STRGRAPH)->__small_counts[(POSITION)] == \
      GT_STRGRAPH__COUNT_IS_LARGE - (GtStrgraphCount__Small)1) \
  { \
    (STRGRAPH)->__small_counts[(POSITION)] = GT_STRGRAPH__COUNT_IS_LARGE;\
    gt_assert(v_c__gt_hashmap_get((STRGRAPH)->__large_counts,\
          (POSITION)) == NULL);\
    v_c__gt_hashmap_add((STRGRAPH)->__large_counts, (POSITION),\
        (GtStrgraphCount__Large)GT_STRGRAPH__COUNT_IS_LARGE);\
  } \
  else\
  {\
    gt_assert((STRGRAPH)->__small_counts[(POSITION)] == \
        GT_STRGRAPH__COUNT_IS_LARGE);\
    gt_assert(v_c__gt_hashmap_get((STRGRAPH)->__large_counts,\
          (POSITION)) != NULL);\
    (*v_c__gt_hashmap_get((STRGRAPH)->__large_counts, (POSITION)))++;\
  }

enum iterator_op gt_strgraph__save_large_count(GtStrgraphVnum vnum,
   GtStrgraphCount__Large count, GtFile *outfp, GT_UNUSED GtError *err)
{
  gt_file_xwrite(outfp, &vnum, sizeof (GtStrgraphVnum));
  gt_file_xwrite(outfp, &count, sizeof (GtStrgraphCount__Large));
  return CONTINUE_ITERATION;
}

#define GT_STRGRAPH_SERIALIZE_COUNTS(STRGRAPH, FP) \
  GT_STRGRAPH_SERIALIZE_DATA((FP), (STRGRAPH)->__n_counts,\
      (STRGRAPH)->__small_counts);\
  (void)v_c__gt_hashmap_foreach((STRGRAPH)->__large_counts,\
      (v_c__gt_hashmap_iteratorfunc)gt_strgraph__save_large_count,\
      (FP), NULL)

#define GT_STRGRAPH__DESERIALIZE_LARGE_COUNTS(STRGRAPH, FP) \
  {\
    GtStrgraphVnum          __large_count_vnum;\
    GtStrgraphCount__Large  __large_count = 0;\
    GT_UNUSED int           __large_count_read_retval;\
    \
    while (gt_file_xread((FP), &__large_count_vnum,\
        sizeof (__large_count_vnum)) == (int)sizeof (__large_count_vnum))\
    {\
      __large_count_read_retval = gt_file_xread((FP), &__large_count, \
          sizeof (__large_count));\
      gt_assert(__large_count_read_retval == (int)sizeof (__large_count));\
      GT_STRGRAPH__SET_COUNT((STRGRAPH), __large_count_vnum, __large_count);\
    }\
  }

#define GT_STRGRAPH_DESERIALIZE_COUNTS(STRGRAPH, FP) \
  GT_STRGRAPH_DESERIALIZE_DATA((FP), \
      (STRGRAPH)->__n_counts, (STRGRAPH)->__small_counts);\
  GT_STRGRAPH__DESERIALIZE_LARGE_COUNTS(STRGRAPH, FP)

#endif
