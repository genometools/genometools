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

#ifndef RDJ_STRGRAPH_EDGES_SHORT_DEF_H
#define RDJ_STRGRAPH_EDGES_SHORT_DEF_H

#include <limits.h>
#include "core/intbits.h"

/*
 * Representation of string graph edges for reads shorter
 * than 255 and less than 2**31 reads (about 2.1E9).
 */

#define GT_STRGRAPH_EDGES_REPRESENTATION "short"

typedef unsigned char GtStrgraphLength;
#define FormatGtStrgraphLength       "%u"
#define PRINTGtStrgraphLengthcast(X) (X)
#define SCANGtStrgraphLengthcast(X)  (X)
#define GT_STRGRAPH__UNDEF_EDGE_LEN  (GtStrgraphLength)UCHAR_MAX
#define GT_STRGRAPH_LENGTH_MAX       (GtStrgraphLength)(UCHAR_MAX-1)

typedef uint32_t GtStrgraphVnum__Short;
#if ULONG_MAX >= UINT32_MAX
#define GT_STRGRAPH_N_READS_MAX (unsigned long)((UINT32_MAX >> 1) - 1)
#else
#define GT_STRGRAPH_N_READS_MAX ((ULONG_MAX >> 1) - 1)
#endif

#define GT_STRGRAPH_DECLARE_EDGES \
  GtStrgraphLength         *__e_len; \
  GtStrgraphVnum__Short    *__e_dest;\
  GtBitsequence            *__e_mark

#define GT_STRGRAPH__EDGE_REDUCED GT_STRGRAPH__UNDEF_EDGE_LEN

#define GT_STRGRAPH_EDGE_SET_DEST(STRGRAPH, V, EDGENUM, DEST) \
  ((STRGRAPH)->__e_dest[GT_STRGRAPH_V_NTH_EDGE_OFFSET(STRGRAPH, V, EDGENUM)] =\
   (DEST))

#define GT_STRGRAPH_EDGE_DEST(STRGRAPH, V, EDGENUM) \
  ((GtStrgraphVnum)\
   ((STRGRAPH)->__e_dest[GT_STRGRAPH_V_NTH_EDGE_OFFSET(STRGRAPH, V, EDGENUM)]))

#define GT_STRGRAPH_EDGE_LEN(STRGRAPH, V, EDGENUM) \
  ((GtStrgraphLength)\
   ((STRGRAPH)->__e_len[GT_STRGRAPH_V_NTH_EDGE_OFFSET(STRGRAPH, V, EDGENUM)]))

#define GT_STRGRAPH_EDGE_SET_LEN(STRGRAPH, V, EDGENUM, LEN) \
  ((STRGRAPH)->__e_len[GT_STRGRAPH_V_NTH_EDGE_OFFSET(STRGRAPH, V, EDGENUM)] =\
   (LEN))

#define GT_STRGRAPH_EDGE_INIT(STRGRAPH, V, EDGENUM) \
  (GT_UNSETIBIT((STRGRAPH)->__e_mark,\
     GT_STRGRAPH_V_NTH_EDGE_OFFSET(STRGRAPH, V, EDGENUM)))

#define GT_STRGRAPH_EDGE_SET_MARK(STRGRAPH, V, EDGENUM) \
  (GT_SETIBIT((STRGRAPH)->__e_mark,\
     GT_STRGRAPH_V_NTH_EDGE_OFFSET(STRGRAPH, V, EDGENUM)))

#define GT_STRGRAPH_EDGE_HAS_MARK(STRGRAPH, V, EDGENUM) \
  (GT_ISIBITSET((STRGRAPH)->__e_mark,\
     GT_STRGRAPH_V_NTH_EDGE_OFFSET(STRGRAPH, V, EDGENUM)))

#define GT_STRGRAPH_EDGE_REDUCE(STRGRAPH, V, EDGENUM) \
  GT_STRGRAPH_EDGE_SET_LEN(STRGRAPH, V, EDGENUM, GT_STRGRAPH__EDGE_REDUCED)

#define GT_STRGRAPH_EDGE_IS_REDUCED(STRGRAPH, V, EDGENUM) \
  (GT_STRGRAPH_EDGE_LEN(STRGRAPH, V, EDGENUM) == GT_STRGRAPH__EDGE_REDUCED)

#define GT_STRGRAPH_NOFEDGES(STRGRAPH) \
  GT_STRGRAPH_V_OFFSET((STRGRAPH), GT_STRGRAPH_NOFVERTICES(STRGRAPH))

#define GT_STRGRAPH_SET_NOFEDGES(STRGRAPH, VALUE) \
  GT_STRGRAPH_V_SET_OFFSET((STRGRAPH),\
      GT_STRGRAPH_NOFVERTICES(STRGRAPH), (VALUE))

#define GT_STRGRAPH_ALLOC_EDGES(STRGRAPH)\
  (STRGRAPH)->__e_dest = gt_calloc(GT_STRGRAPH_NOFEDGES(STRGRAPH),\
      sizeof (*(STRGRAPH)->__e_dest));\
  (STRGRAPH)->__e_len = gt_calloc(GT_STRGRAPH_NOFEDGES(STRGRAPH),\
      sizeof (*(STRGRAPH)->__e_len));\
  GT_INITBITTAB((STRGRAPH)->__e_mark, GT_STRGRAPH_NOFEDGES(STRGRAPH))

#define GT_STRGRAPH_SIZEOF_EDGES(STRGRAPH) \
  (((sizeof (*(STRGRAPH)->__e_len) + sizeof (*(STRGRAPH)->__e_dest)) * \
    GT_STRGRAPH_NOFEDGES(STRGRAPH)) + (sizeof (GtBitsequence) * \
    GT_NUMOFINTSFORBITS(GT_STRGRAPH_NOFEDGES(STRGRAPH))))

#define GT_STRGRAPH_SHRINK_EDGES(STRGRAPH, NEWSIZE)\
  gt_assert((NEWSIZE) < GT_STRGRAPH_NOFEDGES(STRGRAPH));\
  GT_STRGRAPH_SET_NOFEDGES(STRGRAPH, NEWSIZE);\
  if ((NEWSIZE) == 0)\
  {\
    GT_STRGRAPH_FREE_EDGES(STRGRAPH);\
  }\
  else\
  {\
    (STRGRAPH)->__e_dest = gt_realloc((STRGRAPH)->__e_dest,\
        sizeof (*(STRGRAPH)->__e_dest) * GT_STRGRAPH_NOFEDGES(STRGRAPH));\
    (STRGRAPH)->__e_len = gt_realloc((STRGRAPH)->__e_len,\
        sizeof (*(STRGRAPH)->__e_len) * GT_STRGRAPH_NOFEDGES(STRGRAPH));\
    (STRGRAPH)->__e_mark = gt_realloc((STRGRAPH)->__e_mark,\
        sizeof (*(STRGRAPH)->__e_mark) * \
        GT_NUMOFINTSFORBITS(GT_STRGRAPH_NOFEDGES(STRGRAPH)));\
  }

#define GT_STRGRAPH_FREE_EDGES(STRGRAPH)\
  gt_free((STRGRAPH)->__e_dest);\
  gt_free((STRGRAPH)->__e_len);\
  gt_free((STRGRAPH)->__e_mark);\
  (STRGRAPH)->__e_dest = NULL;\
  (STRGRAPH)->__e_len = NULL;\
  (STRGRAPH)->__e_mark = NULL

#define GT_STRGRAPH_COPY_EDGE(STRGRAPH, OFFSET_SOURCE, OFFSET_DEST)\
  (STRGRAPH)->__e_dest[OFFSET_DEST] = (STRGRAPH)->__e_dest[OFFSET_SOURCE];\
  (STRGRAPH)->__e_len[OFFSET_DEST] = (STRGRAPH)->__e_len[OFFSET_SOURCE];\
  if (GT_ISIBITSET((STRGRAPH)->__e_mark, OFFSET_SOURCE))\
  {\
    GT_SETIBIT((STRGRAPH)->__e_mark, OFFSET_DEST);\
  }\
  else\
  {\
    GT_UNSETIBIT((STRGRAPH)->__e_mark, OFFSET_DEST);\
  }

/* order by length from < to > */

struct GtStrgraphEdgeData
{
  GtStrgraphVnum__Short  n;
  GtStrgraphLength         len;
  bool                     to_reduce;
};

int gt_strgraph_edges_compare_by_length(const void *edgea,
    const void *edgeb)
{
  const struct GtStrgraphEdgeData *a = edgea, *b = edgeb;
  return (int)(a->len > b->len) - (int)(a->len < b->len);
}

#define GT_STRGRAPH_SORT_V_EDGES(STRGRAPH, VNUM)\
{\
  uint64_t v_nofedges, v_edges_i;\
  struct GtStrgraphEdgeData *v_edges;\
  v_nofedges = GT_STRGRAPH_V_NOFEDGES(STRGRAPH, VNUM);\
  if (v_nofedges > (uint64_t)1)\
  {\
    v_edges = gt_malloc(sizeof (*v_edges) * v_nofedges);\
    for (v_edges_i = 0; v_edges_i < v_nofedges; v_edges_i++)\
    {\
      v_edges[v_edges_i].n = GT_STRGRAPH_EDGE_DEST(STRGRAPH, VNUM, v_edges_i);\
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
          v_edges_i, v_edges[v_edges_i].n);\
      GT_STRGRAPH_EDGE_SET_LEN(STRGRAPH, VNUM, v_edges_i, \
          v_edges[v_edges_i].len);\
      GT_STRGRAPH_EDGE_INIT(STRGRAPH, VNUM, v_edges_i);\
      if (v_edges[v_edges_i].to_reduce)\
        GT_STRGRAPH_EDGE_SET_MARK(STRGRAPH, VNUM, v_edges_i);\
    }\
    gt_free(v_edges);\
  }\
}

#define GT_STRGRAPH_SERIALIZE_EDGES(STRGRAPH, FP)\
  GT_STRGRAPH_SERIALIZE_DATA((FP), \
      GT_STRGRAPH_NOFEDGES(STRGRAPH), (STRGRAPH)->__e_dest);\
  GT_STRGRAPH_SERIALIZE_DATA((FP),\
      GT_STRGRAPH_NOFEDGES(STRGRAPH), (STRGRAPH)->__e_len);\
  GT_STRGRAPH_SERIALIZE_DATA((FP),\
      GT_NUMOFINTSFORBITS(GT_STRGRAPH_NOFEDGES(STRGRAPH)),\
      (STRGRAPH)->__e_mark)

#define GT_STRGRAPH_DESERIALIZE_EDGES(STRGRAPH, FP)\
  GT_STRGRAPH_DESERIALIZE_DATA((FP),\
      GT_STRGRAPH_NOFEDGES(STRGRAPH), (STRGRAPH)->__e_dest);\
  GT_STRGRAPH_DESERIALIZE_DATA((FP),\
      GT_STRGRAPH_NOFEDGES(STRGRAPH), (STRGRAPH)->__e_len);\
  GT_STRGRAPH_DESERIALIZE_DATA((FP),\
      GT_NUMOFINTSFORBITS(GT_STRGRAPH_NOFEDGES(STRGRAPH)),\
      (STRGRAPH)->__e_mark)

#define GT_STRGRAPH_FIND_LONGEST_EDGE(STRGRAPH, VNUM, LONGEST) \
        {\
          unsigned long ith_edge_before_last, vnofedges;\
          vnofedges = GT_STRGRAPH_V_NOFEDGES(STRGRAPH, VNUM);\
          LONGEST = GT_STRGRAPH__UNDEF_EDGE_LEN;\
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
        }

#endif
