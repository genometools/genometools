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

#ifndef RDJ_STRGRAPH_EDGES_BITPACK_DEF_H
#define RDJ_STRGRAPH_EDGES_BITPACK_DEF_H

#include <limits.h>
#include "core/bitpackarray.h"
#include "core/intbits.h"

#define GT_STRGRAPH_EDGES_REPRESENTATION "bitpack"

typedef uint64_t GtStrgraphLength;
#define FormatGtStrgraphLength        Formatuint64_t
#define PRINTGtStrgraphLengthcast(X)  PRINTuint64_tcast(X)
#define SCANGtStrgraphLengthcast(X)   SCANuint64_tcast(X)

/* LENGTH_MAX is limited by __e_len */
#define GT_STRGRAPH_LENGTH_MAX   (GtStrgraphLength)UINT64_MAX

#if UINT64_MAX >= ULONG_MAX
/* N_READS_MAX is limited by read numbers in the encseq */
#define GT_STRGRAPH_N_READS_MAX ((ULONG_MAX >> 1) - 1)
#else
/* N_READS_MAX is limited by __e_dest */
#define GT_STRGRAPH_N_READS_MAX (unsigned long)((UINT64_MAX >> 1) - 1)
#endif

#define GT_STRGRAPH_DECLARE_EDGES\
  BitPackArray       *__e_len;\
  BitPackArray       *__e_dest;\
  GtBitsequence      *__e_mark;\
  GtStrgraphLength   __len_max;\
  GtStrgraphEdgenum   __n_edges

#define GT_STRGRAPH_SET_NOFEDGES(STRGRAPH, VALUE) \
  (STRGRAPH)->__n_edges = (VALUE)

#define GT_STRGRAPH_NOFEDGES(STRGRAPH) \
  ((STRGRAPH)->__n_edges)

#define GT_STRGRAPH__LEN_BITS(STRGRAPH)\
  (gt_requiredUInt64Bits((STRGRAPH)->__len_max))

#define GT_STRGRAPH__UNDEF_EDGE_LEN(STRGRAPH) \
  ((STRGRAPH)->__len_max)

#define GT_STRGRAPH__EDGE_REDUCED(STRGRAPH) \
  GT_STRGRAPH__UNDEF_EDGE_LEN(STRGRAPH)

#define GT_STRGRAPH__ALLOC_LENGTHSLIST(STRGRAPH)\
  (STRGRAPH)->__len_max = gt_strgraph_longest_read(STRGRAPH) - \
    (STRGRAPH)->minmatchlen + 1;\
  (STRGRAPH)->__e_len = bitpackarray_new(GT_STRGRAPH__LEN_BITS(STRGRAPH),\
      (BitOffset)GT_STRGRAPH_NOFEDGES(STRGRAPH), true)

#define GT_STRGRAPH__DEST_BITS(STRGRAPH)\
  (gt_requiredUInt64Bits(GT_STRGRAPH_NOFVERTICES(STRGRAPH)))

#define GT_STRGRAPH__ALLOC_NEIGHBOURSLIST(STRGRAPH)\
  (STRGRAPH)->__e_dest = bitpackarray_new(GT_STRGRAPH__DEST_BITS(STRGRAPH),\
      (BitOffset)GT_STRGRAPH_NOFEDGES(STRGRAPH), true)

#define GT_STRGRAPH_ALLOC_EDGES(STRGRAPH)\
  GT_STRGRAPH__ALLOC_NEIGHBOURSLIST(STRGRAPH);\
  GT_STRGRAPH__ALLOC_LENGTHSLIST(STRGRAPH);\
  GT_INITBITTAB((STRGRAPH)->__e_mark, GT_STRGRAPH_NOFEDGES(STRGRAPH))

#define GT_STRGRAPH__SIZEOF_LENGTHSLIST(STRGRAPH)\
  (sizeofbitarray(GT_STRGRAPH__LEN_BITS(STRGRAPH),\
                  (BitOffset)GT_STRGRAPH_NOFEDGES(STRGRAPH)))

#define GT_STRGRAPH__SIZEOF_NEIGHBOURSLIST(STRGRAPH)\
  (sizeofbitarray(GT_STRGRAPH__DEST_BITS(STRGRAPH),\
                  (BitOffset)GT_STRGRAPH_NOFEDGES(STRGRAPH)))

#define GT_STRGRAPH__SIZEOF_TOREDUCELIST(STRGRAPH)\
  (sizeof (GtBitsequence) * \
   GT_NUMOFINTSFORBITS(GT_STRGRAPH_NOFEDGES(STRGRAPH)))

#define GT_STRGRAPH_SIZEOF_EDGES(STRGRAPH) \
  (sizeof ((STRGRAPH)->__e_len)+\
      GT_STRGRAPH__SIZEOF_LENGTHSLIST(STRGRAPH) +\
      sizeof ((STRGRAPH)->__len_max)+\
      sizeof ((STRGRAPH)->__e_dest)+\
      GT_STRGRAPH__SIZEOF_NEIGHBOURSLIST(STRGRAPH) +\
      sizeof ((STRGRAPH)->__e_mark)+\
      GT_STRGRAPH__SIZEOF_TOREDUCELIST(STRGRAPH)+\
      sizeof ((STRGRAPH)->__n_edges))

#define GT_STRGRAPH_SERIALIZE_EDGES(STRGRAPH, FP)\
  do {\
    GT_STRGRAPH_SERIALIZE_DATA((FP),\
        bitElemsAllocSize(GT_STRGRAPH__DEST_BITS(STRGRAPH) * \
          GT_STRGRAPH_NOFEDGES(STRGRAPH)), (STRGRAPH)->__e_dest);\
    GT_STRGRAPH_SERIALIZE_DATA((FP), \
        bitElemsAllocSize(GT_STRGRAPH__LEN_BITS(STRGRAPH) * \
          GT_STRGRAPH_NOFEDGES(STRGRAPH)), (STRGRAPH)->__e_len);\
    GT_STRGRAPH_SERIALIZE_DATA((FP),\
        GT_NUMOFINTSFORBITS(GT_STRGRAPH_NOFEDGES(STRGRAPH)),\
        (STRGRAPH)->__e_mark);\
  } while (false)

#define GT_STRGRAPH_DESERIALIZE_EDGES(STRGRAPH, FP)\
  do {\
    GT_STRGRAPH_DESERIALIZE_DATA((FP),\
        bitElemsAllocSize(GT_STRGRAPH__DEST_BITS(STRGRAPH) * \
          GT_STRGRAPH_NOFEDGES(STRGRAPH)), (STRGRAPH)->__e_dest);\
    GT_STRGRAPH_DESERIALIZE_DATA((FP), \
        bitElemsAllocSize(GT_STRGRAPH__LEN_BITS(STRGRAPH) * \
          GT_STRGRAPH_NOFEDGES(STRGRAPH)), (STRGRAPH)->__e_len);\
    GT_STRGRAPH_DESERIALIZE_DATA((FP),\
        GT_NUMOFINTSFORBITS(GT_STRGRAPH_NOFEDGES(STRGRAPH)),\
        (STRGRAPH)->__e_mark);\
  } while (false)

#define GT_STRGRAPH_SHRINK_EDGES(STRGRAPH, NEWSIZE)\
  do {\
    gt_assert((NEWSIZE) < GT_STRGRAPH_NOFEDGES(STRGRAPH));\
    GT_STRGRAPH_SET_NOFEDGES(STRGRAPH, NEWSIZE);\
    if ((NEWSIZE) == 0)\
    {\
      GT_STRGRAPH_FREE_EDGES(STRGRAPH);\
    }\
    else\
    {\
      /* XXX: realloc of __e_dest and __e_len not implemented */\
      (STRGRAPH)->__e_mark = gt_realloc((STRGRAPH)->__e_mark,\
          sizeof (*(STRGRAPH)->__e_mark) * \
          GT_NUMOFINTSFORBITS(GT_STRGRAPH_NOFEDGES(STRGRAPH)));\
    }\
  } while (false)

#define GT_STRGRAPH_FREE_EDGES(STRGRAPH)\
  bitpackarray_delete((STRGRAPH)->__e_dest);\
  bitpackarray_delete((STRGRAPH)->__e_len);\
  gt_free((STRGRAPH)->__e_mark);\
  (STRGRAPH)->__e_dest = NULL;\
  (STRGRAPH)->__e_len = NULL;\
  (STRGRAPH)->__e_mark = NULL

#define GT_STRGRAPH_EDGE_SET_DEST(STRGRAPH, V, EDGENUM, DEST) \
  bitpackarray_store_uint64((STRGRAPH)->__e_dest,\
      (BitOffset)GT_STRGRAPH_V_NTH_EDGE_OFFSET(STRGRAPH, V, EDGENUM), (DEST))

#define GT_STRGRAPH_EDGE_DEST(STRGRAPH, V, EDGENUM) \
  ((GtStrgraphVnum)(bitpackarray_get_uint64((STRGRAPH)->__e_dest,\
      (BitOffset)GT_STRGRAPH_V_NTH_EDGE_OFFSET(STRGRAPH, V, EDGENUM))))

#define GT_STRGRAPH_EDGE_SET_LEN(STRGRAPH, V, EDGENUM, LEN) \
  bitpackarray_store_uint64((STRGRAPH)->__e_len,\
      (BitOffset)GT_STRGRAPH_V_NTH_EDGE_OFFSET(STRGRAPH, V, EDGENUM), (LEN))

#define GT_STRGRAPH_EDGE_LEN(STRGRAPH, V, EDGENUM) \
  ((GtStrgraphLength)bitpackarray_get_uint64((STRGRAPH)->__e_len,\
      (BitOffset)GT_STRGRAPH_V_NTH_EDGE_OFFSET(STRGRAPH, V, EDGENUM)))

#define GT_STRGRAPH_EDGE_INIT(STRGRAPH, V, EDGENUM) \
  (GT_UNSETIBIT((STRGRAPH)->__e_mark,\
     GT_STRGRAPH_V_NTH_EDGE_OFFSET(STRGRAPH, V, EDGENUM)))

#define GT_STRGRAPH_EDGE_SET_MARK(STRGRAPH, V, EDGENUM) \
  (GT_SETIBIT((STRGRAPH)->__e_mark,\
     GT_STRGRAPH_V_NTH_EDGE_OFFSET(STRGRAPH, V, EDGENUM)))

#define GT_STRGRAPH_EDGE_HAS_MARK(STRGRAPH, V, EDGENUM) \
  (GT_ISIBITSET((STRGRAPH)->__e_mark,\
     GT_STRGRAPH_V_NTH_EDGE_OFFSET(STRGRAPH, V, EDGENUM))\
   ? true : false)

#define GT_STRGRAPH_EDGE_REDUCE(STRGRAPH, V, EDGENUM) \
  GT_STRGRAPH_EDGE_SET_LEN(STRGRAPH, V, EDGENUM,\
      GT_STRGRAPH__EDGE_REDUCED(STRGRAPH))

#define GT_STRGRAPH_EDGE_IS_REDUCED(STRGRAPH, V, EDGENUM) \
  (GT_STRGRAPH_EDGE_LEN(STRGRAPH, V, EDGENUM) == \
      GT_STRGRAPH__EDGE_REDUCED(STRGRAPH) \
   ? true : false)

#define GT_STRGRAPH_COPY_EDGE(STRGRAPH, OFFSET_SOURCE, OFFSET_DEST)\
  do {\
    bitpackarray_store_uint64((STRGRAPH)->__e_dest,\
        (BitOffset)(OFFSET_DEST), bitpackarray_get_uint64((STRGRAPH->__e_dest),\
          (BitOffset)(OFFSET_SOURCE)));\
    bitpackarray_store_uint64((STRGRAPH)->__e_len,\
        (BitOffset)(OFFSET_DEST), bitpackarray_get_uint64((STRGRAPH->__e_len),\
          (BitOffset)(OFFSET_SOURCE)));\
    if (GT_ISIBITSET((STRGRAPH)->__e_mark, OFFSET_SOURCE))\
    {\
      GT_SETIBIT((STRGRAPH)->__e_mark, OFFSET_DEST);\
    }\
    else\
    {\
      GT_UNSETIBIT((STRGRAPH)->__e_mark, OFFSET_DEST);\
    }\
  } while (false)

/* order by length from < to > */

struct GtStrgraph__EdgeData
{
  GtStrgraphVnum     dest;
  GtStrgraphLength   len;
  bool               to_reduce;
};

int gt_strgraph_edges_compare_by_length(const void *edgea,
    const void *edgeb)
{
  const struct GtStrgraph__EdgeData *a = edgea, *b = edgeb;
  return (int)(a->len > b->len) - (int)(a->len < b->len);
}

#define GT_STRGRAPH_SORT_V_EDGES(STRGRAPH, VNUM)\
  do {\
    GtStrgraphVEdgenum v_nofedges, v_edges_i;\
    struct GtStrgraph__EdgeData *v_edges;\
    v_nofedges = GT_STRGRAPH_V_NOFEDGES(STRGRAPH, VNUM);\
    if (v_nofedges > (GtStrgraphVEdgenum)1)\
    {\
      v_edges = gt_malloc(sizeof (*v_edges) * v_nofedges);\
      for (v_edges_i = 0; v_edges_i < v_nofedges; v_edges_i++)\
      {\
        v_edges[v_edges_i].dest = \
          GT_STRGRAPH_EDGE_DEST(STRGRAPH, VNUM, v_edges_i);\
        v_edges[v_edges_i].len = \
          GT_STRGRAPH_EDGE_LEN(STRGRAPH, VNUM, v_edges_i);\
        v_edges[v_edges_i].to_reduce = \
          GT_STRGRAPH_EDGE_HAS_MARK(STRGRAPH, VNUM, v_edges_i);\
      }\
      qsort(v_edges, v_nofedges, sizeof (*v_edges),\
          gt_strgraph_edges_compare_by_length);\
      for (v_edges_i = 0; v_edges_i < v_nofedges; v_edges_i++)\
      {\
        GT_STRGRAPH_EDGE_SET_DEST(STRGRAPH, VNUM, \
            v_edges_i, v_edges[v_edges_i].dest);\
        GT_STRGRAPH_EDGE_SET_LEN(STRGRAPH, VNUM, v_edges_i, \
            v_edges[v_edges_i].len);\
        GT_STRGRAPH_EDGE_INIT(STRGRAPH, VNUM, v_edges_i);\
        if (v_edges[v_edges_i].to_reduce)\
          GT_STRGRAPH_EDGE_SET_MARK(STRGRAPH, VNUM, v_edges_i);\
      }\
      gt_free(v_edges);\
    }\
  } while (false)

#define GT_STRGRAPH_FIND_LONGEST_EDGE(STRGRAPH, VNUM, LONGEST) \
  do {\
    GtStrgraphVEdgenum ith_edge_before_last, vnofedges;\
    vnofedges = GT_STRGRAPH_V_NOFEDGES(STRGRAPH, VNUM);\
    LONGEST = (STRGRAPH)->__len_max;\
    for (ith_edge_before_last = 0; ith_edge_before_last < vnofedges;\
        ith_edge_before_last++)\
    {\
      if (!GT_STRGRAPH_EDGE_IS_REDUCED(STRGRAPH, VNUM, vnofedges -\
            ith_edge_before_last - 1))\
      {\
        LONGEST = GT_STRGRAPH_EDGE_LEN(STRGRAPH, VNUM, vnofedges -\
            ith_edge_before_last - 1);\
        break;\
      }\
    }\
  } while (false)

#endif
