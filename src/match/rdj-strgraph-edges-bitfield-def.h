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

#ifndef RDJ_STRGRAPH_EDGES_BITFIELD_DEF_H
#define RDJ_STRGRAPH_EDGES_BITFIELD_DEF_H

#define GT_STRGRAPH_EDGES_REPRESENTATION "bitfield"

typedef unsigned long GtStrgraphLength;
#define FormatGtStrgraphLength       "%lu"
#define PRINTGtStrgraphLengthcast(X) (X)
#define SCANGtStrgraphLengthcast(X)  (X)
#define GT_STRGRAPH__LEN_BITS        18
#define GT_STRGRAPH_LENGTH_MAX \
  ((GtStrgraphLength)(1UL << GT_STRGRAPH__LEN_BITS) - 1)

#define GT_STRGRAPH__DEST_BITS 44
#ifndef S_SPLINT_S
#define GT_STRGRAPH__DEST_MAX ((1ULL << GT_STRGRAPH__DEST_BITS) - 1)
#else
#define GT_STRGRAPH__DEST_MAX ULONG_MAX
#endif

#if GT_STRGRAPH__DEST_MAX >= ULONG_MAX
#define GT_STRGRAPH_N_READS_MAX (unsigned long)(GT_STRGRAPH__DEST_MAX >> 1)
#else
#define GT_STRGRAPH_N_READS_MAX (ULONG_MAX >> 1)
#endif

typedef struct {
  GtStrgraphVnum    __dest      :GT_STRGRAPH__DEST_BITS;
  GtStrgraphLength  __len       :GT_STRGRAPH__LEN_BITS;
  bool              __reduced   :1;
  bool              __mark      :1;
} GtStrgraph__Edge;

#define GT_STRGRAPH_DECLARE_EDGES \
  GtStrgraphEdgenum  __n_edges;\
  GtStrgraph__Edge  *__e

#define GT_STRGRAPH__EDGE_NOT_REDUCED   false
#define GT_STRGRAPH__EDGE_REDUCED       true
#define GT_STRGRAPH__EDGE_NOT_MARKED    false
#define GT_STRGRAPH__EDGE_MARKED        true

#define GT_STRGRAPH__V_NTH_EDGE(STRGRAPH, V, N) \
  ((STRGRAPH)->__e + GT_STRGRAPH_V_NTH_EDGE_OFFSET(STRGRAPH, V, N))

/* edge accessor macros */

#define GT_STRGRAPH_EDGE_SET_DEST(STRGRAPH, V, EDGENUM, DEST) \
  GT_STRGRAPH__V_NTH_EDGE(STRGRAPH, V, EDGENUM)->__dest = (DEST)

#define GT_STRGRAPH_EDGE_DEST(STRGRAPH, V, EDGENUM) \
  ((GtStrgraphVnum)(GT_STRGRAPH__V_NTH_EDGE(STRGRAPH, V, EDGENUM)->__dest))

#define GT_STRGRAPH_EDGE_SET_LEN(STRGRAPH, V, EDGENUM, LEN) \
  GT_STRGRAPH__V_NTH_EDGE(STRGRAPH, V, EDGENUM)->__len = (LEN)

#define GT_STRGRAPH_EDGE_LEN(STRGRAPH, V, EDGENUM) \
  ((GtStrgraphLength)(GT_STRGRAPH__V_NTH_EDGE(STRGRAPH, V, EDGENUM)->__len))

#define GT_STRGRAPH_EDGE_INIT(STRGRAPH, V, EDGENUM) \
  GT_STRGRAPH__V_NTH_EDGE(STRGRAPH, V, EDGENUM)->__reduced = \
      GT_STRGRAPH__EDGE_NOT_REDUCED; \
  GT_STRGRAPH__V_NTH_EDGE(STRGRAPH, V, EDGENUM)->__mark = \
      GT_STRGRAPH__EDGE_NOT_MARKED

#define GT_STRGRAPH_EDGE_SET_MARK(STRGRAPH, V, EDGENUM) \
  GT_STRGRAPH__V_NTH_EDGE(STRGRAPH, V, EDGENUM)->__mark = \
      GT_STRGRAPH__EDGE_MARKED

#define GT_STRGRAPH_EDGE_HAS_MARK(STRGRAPH, V, EDGENUM) \
  ((GT_STRGRAPH__V_NTH_EDGE(STRGRAPH, V, EDGENUM)->__mark == \
   GT_STRGRAPH__EDGE_MARKED) ? true : false)

#define GT_STRGRAPH_EDGE_REDUCE(STRGRAPH, V, EDGENUM) \
  (GT_STRGRAPH__V_NTH_EDGE(STRGRAPH, V, EDGENUM)->__reduced = \
   GT_STRGRAPH__EDGE_REDUCED)

#define GT_STRGRAPH_EDGE_IS_REDUCED(STRGRAPH, V, EDGENUM) \
  (GT_STRGRAPH__V_NTH_EDGE(STRGRAPH, V, EDGENUM)->__reduced == \
   GT_STRGRAPH__EDGE_REDUCED)

#define GT_STRGRAPH_NOFEDGES(STRGRAPH) \
  ((STRGRAPH)->__n_edges)

#define GT_STRGRAPH_SET_NOFEDGES(STRGRAPH, VALUE) \
  (STRGRAPH)->__n_edges = (VALUE)

#define GT_STRGRAPH__SIZEOF_EDGE \
  (sizeof (GtStrgraph__Edge))

#define GT_STRGRAPH_SIZEOF_EDGES(STRGRAPH) \
  (GT_STRGRAPH__SIZEOF_EDGE * GT_STRGRAPH_NOFEDGES(STRGRAPH) +\
   sizeof ((STRGRAPH)->__n_edges))

#define GT_STRGRAPH_ALLOC_EDGES(STRGRAPH)\
  (STRGRAPH)->__e = gt_calloc((size_t)GT_STRGRAPH_NOFEDGES(STRGRAPH),\
      sizeof (*(STRGRAPH)->__e))

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
      (STRGRAPH)->__e = \
          gt_realloc((STRGRAPH)->__e, GT_STRGRAPH_NOFEDGES(STRGRAPH) * \
          sizeof (*(STRGRAPH)->__e));\
    }\
  } while (false)

#define GT_STRGRAPH_FREE_EDGES(STRGRAPH)\
  gt_free((STRGRAPH)->__e);\
  (STRGRAPH)->__e = NULL

#define GT_STRGRAPH_COPY_EDGE(STRGRAPH, OFFSET_SOURCE, OFFSET_DEST)\
  (STRGRAPH)->__e[OFFSET_DEST] = (STRGRAPH)->__e[OFFSET_SOURCE]

/* order by length from < to > */
int gt_strgraph_edges_compare_by_length(const void *edgea,
    const void *edgeb)
{
  const GtStrgraph__Edge *a = edgea, *b = edgeb;
  return (int)(a->__len > b->__len) - (int)(a->__len < b->__len);
}

#define GT_STRGRAPH_SORT_V_EDGES(STRGRAPH, VNUM)\
  do {\
    GtStrgraphVEdgenum v_nofedges;\
    v_nofedges = GT_STRGRAPH_V_NOFEDGES(STRGRAPH, VNUM);\
    if (v_nofedges > (GtStrgraphVEdgenum)1)\
      qsort((STRGRAPH)->__e + GT_STRGRAPH_V_OFFSET(STRGRAPH, VNUM),\
          (size_t)v_nofedges, sizeof (*(STRGRAPH)->__e),\
          gt_strgraph_edges_compare_by_length);\
  } while (false)

#define GT_STRGRAPH_SERIALIZE_EDGES(STRGRAPH, FP)\
  GT_STRGRAPH_SERIALIZE_DATA((FP), \
      GT_STRGRAPH_NOFEDGES(STRGRAPH), (STRGRAPH)->__e)

#define GT_STRGRAPH_DESERIALIZE_EDGES(STRGRAPH, FP)\
  GT_STRGRAPH_DESERIALIZE_DATA((FP), \
      GT_STRGRAPH_NOFEDGES(STRGRAPH), (STRGRAPH)->__e)

#define GT_STRGRAPH_FIND_LONGEST_EDGE(STRGRAPH, VNUM, LONGEST) \
        LONGEST = GT_STRGRAPH_EDGE_LEN(STRGRAPH, VNUM, \
          GT_STRGRAPH_V_NOFEDGES(STRGRAPH, VNUM) - 1);

#endif
