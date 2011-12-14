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

#ifndef RDJ_STRGRAPH_VERTICES_BITFIELD_DEF_H
#define RDJ_STRGRAPH_VERTICES_BITFIELD_DEF_H

/* --- for exclusive use in rdj-strgraph.c (__ = private) --- */

#define GT_STRGRAPH_VERTICES_REPRESENTATION "bitfield"

typedef unsigned long GtStrgraphVEdgenum;
#define FormatGtStrgraphVEdgenum       "%lu"
#define PRINTGtStrgraphVEdgenumcast(X) (X)
#define SCANGtStrgraphVEdgenumcast(X)  (X)
#define GT_STRGRAPH__OUTDEG_BITS       24
#define GT_STRGRAPH_V_EDGENUM_MAX \
  ((GtStrgraphVEdgenum)((1UL << GT_STRGRAPH__OUTDEG_BITS) - 1))

typedef uint64_t GtStrgraphEdgenum;
#define FormatGtStrgraphEdgenum       Formatuint64_t
#define PRINTGtStrgraphEdgenumcast(X) PRINTuint64_tcast(X)
#define SCANGtStrgraphEdgenumcast(X)  SCANuint64_tcast(X)
#define GT_STRGRAPH__OFFSET_BITS      38
#define GT_STRGRAPH_EDGENUM_MAX \
  ((GtStrgraphEdgenum)((1ULL << GT_STRGRAPH__OFFSET_BITS) - 1))

typedef struct {
  GtStrgraphEdgenum  offset :GT_STRGRAPH__OFFSET_BITS;
  GtStrgraphVEdgenum outdeg :GT_STRGRAPH__OUTDEG_BITS;
  GtStrgraphVmark    mark   :GT_STRGRAPH_VMARK_BITS;
} GtStrgraph__Vertex;

#define GT_STRGRAPH_DECLARE_VERTICES\
  GtStrgraphVnum       __n_vertices;\
  GtStrgraph__Vertex   *__v

#define GT_STRGRAPH_SET_NOFVERTICES(STRGRAPH, N) \
  ((STRGRAPH)->__n_vertices = (N))

#define GT_STRGRAPH_NOFVERTICES(STRGRAPH) \
  ((STRGRAPH)->__n_vertices)

#define GT_STRGRAPH_ALLOC_VERTICES(STRGRAPH)\
  (STRGRAPH)->__v = \
      gt_calloc((size_t)((GT_STRGRAPH_NOFVERTICES(STRGRAPH)) + 1),\
          sizeof (GtStrgraph__Vertex))

#define GT_STRGRAPH_SIZEOF_VERTICES(STRGRAPH) \
  ((sizeof (GtStrgraph__Vertex) * (GT_STRGRAPH_NOFVERTICES(STRGRAPH) + 1)) +\
   sizeof ((STRGRAPH)->__n_vertices))

#define GT_STRGRAPH_FREE_VERTICES(STRGRAPH)\
  gt_free((STRGRAPH)->__v);\
  (STRGRAPH)->__v = NULL

#define GT_STRGRAPH_SERIALIZE_VERTICES(STRGRAPH, FP)\
  GT_STRGRAPH_SERIALIZE_DATA((FP),\
      GT_STRGRAPH_NOFVERTICES(STRGRAPH) + 1, (STRGRAPH)->__v)

#define GT_STRGRAPH_DESERIALIZE_VERTICES(STRGRAPH, FP)\
  GT_STRGRAPH_DESERIALIZE_DATA((FP),\
      GT_STRGRAPH_NOFVERTICES(STRGRAPH) + 1, (STRGRAPH)->__v)

#define GT_STRGRAPH_V_SET_OFFSET(STRGRAPH, V, VALUE) \
  ((STRGRAPH)->__v[(V)].offset = (VALUE))

#define GT_STRGRAPH_V_OFFSET(STRGRAPH, V) \
  ((GtStrgraphEdgenum)((STRGRAPH)->__v[(V)].offset))

#define GT_STRGRAPH_V_INC_OUTDEG(STRGRAPH, V) \
  (((STRGRAPH)->__v[(V)].outdeg)++)

#define GT_STRGRAPH_V_DEC_OUTDEG(STRGRAPH, V) \
  (((STRGRAPH)->__v[(V)].outdeg)--)

#define GT_STRGRAPH_V_OUTDEG(STRGRAPH, V) \
  ((GtStrgraphVEdgenum)((STRGRAPH)->__v[(V)].outdeg))

#define GT_STRGRAPH_V_SET_MARK(STRGRAPH, V, VALUE) \
  ((STRGRAPH)->__v[(V)].mark = (VALUE))

#define GT_STRGRAPH_V_MARK(STRGRAPH, V) \
  ((STRGRAPH)->__v[(V)].mark)

#endif
