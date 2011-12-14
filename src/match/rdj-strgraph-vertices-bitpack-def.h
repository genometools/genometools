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

#ifndef RDJ_STRGRAPH_VERTICES_BITPACK_DEF_H
#define RDJ_STRGRAPH_VERTICES_BITPACK_DEF_H

/* --- for exclusive use in rdj-strgraph.c (__ = private) --- */

#define GT_STRGRAPH_VERTICES_REPRESENTATION "bitpack"

typedef uint64_t GtStrgraphVEdgenum;
#define FormatGtStrgraphVEdgenum       Formatuint64_t
#define PRINTGtStrgraphVEdgenumcast(X) PRINTuint64_tcast(X)
#define SCANGtStrgraphVEdgenumcast(X)  SCANuint64_tcast(X)
#define GT_STRGRAPH_V_EDGENUM_MAX      (UINT64_MAX - 1)
/* - 1 is not necessary, but avoids a compiler warning */

typedef uint64_t GtStrgraphEdgenum;
#define FormatGtStrgraphEdgenum       Formatuint64_t
#define PRINTGtStrgraphEdgenumcast(X) PRINTuint64_tcast(X)
#define SCANGtStrgraphEdgenumcast(X)  SCANuint64_tcast(X)
#define GT_STRGRAPH_EDGENUM_MAX       (UINT64_MAX - 1)
/* - 1 is not necessary, but avoids a compiler warning */

#define GT_STRGRAPH_DECLARE_VERTICES\
  BitPackArray       *__v_mark;\
  BitPackArray       *__v_offset;\
  GtStrgraphEdgenum  __offset_max;\
  BitPackArray       *__v_outdeg;\
  GtStrgraphVEdgenum __outdeg_max;\
  GtStrgraphVnum     __n_vertices

/* vertex accessor macros */

#define GT_STRGRAPH_SET_NOFVERTICES(STRGRAPH, N) \
  ((STRGRAPH)->__n_vertices = (N))

#define GT_STRGRAPH_NOFVERTICES(STRGRAPH) \
  ((STRGRAPH)->__n_vertices)

#define GT_STRGRAPH__ALLOC_VMARKS(STRGRAPH)\
  (STRGRAPH)->__v_mark = bitpackarray_new((unsigned int)GT_STRGRAPH_VMARK_BITS,\
      (BitOffset)(GT_STRGRAPH_NOFVERTICES(STRGRAPH) + (GtStrgraphVnum)1), true)

#define GT_STRGRAPH__OFFSET_BITS(STRGRAPH)\
  (gt_requiredUInt64Bits((STRGRAPH)->__offset_max))

#define GT_STRGRAPH__ALLOC_OFFSETS(STRGRAPH)\
  (STRGRAPH)->__offset_max = gt_strgraph_counts_sum(STRGRAPH);\
  (STRGRAPH)->__v_offset = bitpackarray_new(GT_STRGRAPH__OFFSET_BITS(STRGRAPH),\
      (BitOffset)(GT_STRGRAPH_NOFVERTICES(STRGRAPH) + (GtStrgraphVnum)1), true)

#define GT_STRGRAPH__OUTDEG_BITS(STRGRAPH)\
  (gt_requiredUInt64Bits((STRGRAPH)->__outdeg_max))

#define GT_STRGRAPH__ALLOC_OUTDEGS(STRGRAPH)\
  gt_assert(sizeof (GtStrgraphVEdgenum) >= sizeof (GtStrgraphCount));\
  (STRGRAPH)->__outdeg_max = (GtStrgraphVEdgenum)\
      gt_strgraph_largest_count(STRGRAPH);\
  (STRGRAPH)->__v_outdeg = bitpackarray_new(GT_STRGRAPH__OUTDEG_BITS(STRGRAPH),\
      (BitOffset)(GT_STRGRAPH_NOFVERTICES(STRGRAPH) + (GtStrgraphVnum)1), true)

#define GT_STRGRAPH_ALLOC_VERTICES(STRGRAPH)\
  GT_STRGRAPH__ALLOC_VMARKS(STRGRAPH);\
  GT_STRGRAPH__ALLOC_OFFSETS(STRGRAPH);\
  GT_STRGRAPH__ALLOC_OUTDEGS(STRGRAPH)

#define GT_STRGRAPH__SIZEOF_OUTDEGS(STRGRAPH)\
  (sizeofbitarray(GT_STRGRAPH__OUTDEG_BITS(STRGRAPH),\
                  (BitOffset)(\
                    GT_STRGRAPH_NOFVERTICES(STRGRAPH) + (GtStrgraphVnum)1)))

#define GT_STRGRAPH__SIZEOF_OFFSETS(STRGRAPH)\
  (sizeofbitarray(GT_STRGRAPH__OFFSET_BITS(STRGRAPH),\
                  (BitOffset)(\
                  GT_STRGRAPH_NOFVERTICES(STRGRAPH) + (GtStrgraphVnum)1)))

#define GT_STRGRAPH__SIZEOF_VMARKS(STRGRAPH)\
  (sizeofbitarray(GT_STRGRAPH_VMARK_BITS,\
                  (BitOffset)(\
                  GT_STRGRAPH_NOFVERTICES(STRGRAPH) + (GtStrgraphVnum)1)))

#define GT_STRGRAPH_SIZEOF_VERTICES(STRGRAPH) \
  (sizeof ((STRGRAPH)->__v_mark)+\
  GT_STRGRAPH__SIZEOF_VMARKS(STRGRAPH) +\
  sizeof ((STRGRAPH)->__v_offset)+\
  sizeof ((STRGRAPH)->__offset_max)+\
  GT_STRGRAPH__SIZEOF_OFFSETS(STRGRAPH) +\
  sizeof ((STRGRAPH)->__v_outdeg)+\
  sizeof ((STRGRAPH)->__outdeg_max)+\
  GT_STRGRAPH__SIZEOF_OUTDEGS(STRGRAPH))

#define GT_STRGRAPH_SERIALIZE_VERTICES(STRGRAPH, FP)\
  GT_STRGRAPH_SERIALIZE_DATA((FP),\
      bitElemsAllocSize(GT_STRGRAPH_VMARK_BITS * \
        GT_STRGRAPH_NOFVERTICES(STRGRAPH) + 1), (STRGRAPH)->__v_mark);\
  GT_STRGRAPH_SERIALIZE_DATA((FP), \
      bitElemsAllocSize(GT_STRGRAPH__OUTDEG_BITS(STRGRAPH) * \
        GT_STRGRAPH_NOFVERTICES(STRGRAPH) + 1), (STRGRAPH)->__v_outdeg);\
  GT_STRGRAPH_SERIALIZE_DATA((FP), \
      bitElemsAllocSize(GT_STRGRAPH__OFFSET_BITS(STRGRAPH) * \
        GT_STRGRAPH_NOFVERTICES(STRGRAPH) + 1), (STRGRAPH)->__v_offset)

#define GT_STRGRAPH_DESERIALIZE_VERTICES(STRGRAPH, FP)\
  GT_STRGRAPH_DESERIALIZE_DATA((FP),\
      bitElemsAllocSize(GT_STRGRAPH_VMARK_BITS * \
        GT_STRGRAPH_NOFVERTICES(STRGRAPH) + 1), (STRGRAPH)->__v_mark);\
  GT_STRGRAPH_DESERIALIZE_DATA((FP), \
      bitElemsAllocSize(GT_STRGRAPH__OUTDEG_BITS(STRGRAPH) * \
        GT_STRGRAPH_NOFVERTICES(STRGRAPH) + 1), (STRGRAPH)->__v_outdeg);\
  GT_STRGRAPH_DESERIALIZE_DATA((FP), \
      bitElemsAllocSize(GT_STRGRAPH__OFFSET_BITS(STRGRAPH) * \
        GT_STRGRAPH_NOFVERTICES(STRGRAPH) + 1), (STRGRAPH)->__v_offset)

#define GT_STRGRAPH_FREE_VERTICES(STRGRAPH)\
  bitpackarray_delete((STRGRAPH)->__v_mark);\
  bitpackarray_delete((STRGRAPH)->__v_outdeg);\
  bitpackarray_delete((STRGRAPH)->__v_offset);\
  (STRGRAPH)->__v_mark = NULL;\
  (STRGRAPH)->__v_outdeg = NULL;\
  (STRGRAPH)->__v_offset = NULL

#define GT_STRGRAPH_V_SET_OFFSET(STRGRAPH, V, VALUE) \
  bitpackarray_store_uint64((STRGRAPH)->__v_offset, (BitOffset)(V),\
     (uint64_t)(VALUE))

#define GT_STRGRAPH_V_OFFSET(STRGRAPH, V) \
  ((GtStrgraphEdgenum)\
   bitpackarray_get_uint64((STRGRAPH)->__v_offset, (BitOffset)(V)))

#define GT_STRGRAPH_V_OUTDEG(STRGRAPH, V) \
  ((GtStrgraphVEdgenum)\
   bitpackarray_get_uint64((STRGRAPH)->__v_outdeg, (BitOffset)(V)))

#define GT_STRGRAPH_V_INC_OUTDEG(STRGRAPH, V) \
  bitpackarray_store_uint64((STRGRAPH)->__v_outdeg, (BitOffset)(V), \
    (uint64_t)(GT_STRGRAPH_V_OUTDEG(STRGRAPH, V) + (GtStrgraphVEdgenum)1))

#define GT_STRGRAPH_V_DEC_OUTDEG(STRGRAPH, V) \
  bitpackarray_store_uint64((STRGRAPH)->__v_outdeg, (BitOffset)(V), \
    (uint64_t)(GT_STRGRAPH_V_OUTDEG(STRGRAPH, V) - (GtStrgraphVEdgenum)1))

#define GT_STRGRAPH_V_SET_MARK(STRGRAPH, V, VALUE) \
  bitpackarray_store_uint64((STRGRAPH)->__v_mark, (BitOffset)(V),\
     (uint64_t)(VALUE))

#define GT_STRGRAPH_V_MARK(STRGRAPH, V) \
  ((GtStrgraphVmark)\
   bitpackarray_get_uint64((STRGRAPH)->__v_mark, (BitOffset)(V)))

#endif
