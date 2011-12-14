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

#ifndef RDJ_STRGRAPH_VERTICES_COMMON_DEF_H
#define RDJ_STRGRAPH_VERTICES_COMMON_DEF_H

/* --- for exclusive use in rdj-strgraph.c --- */

#define GT_STRGRAPH_NOFSPM_MAX \
  (GT_STRGRAPH_EDGENUM_MAX >> 1)

#define GT_STRGRAPH_NOFREADS(STRGRAPH) \
  (GT_STRGRAPH_NOFVERTICES(STRGRAPH) >> 1)

#define GT_STRGRAPH_V_NTH_EDGE_OFFSET(STRGRAPH, V, N)\
  (GT_STRGRAPH_V_OFFSET(STRGRAPH, V) + (N))

#define GT_STRGRAPH_V_INDEG(STRGRAPH, V) \
  GT_STRGRAPH_V_OUTDEG(STRGRAPH, GT_STRGRAPH_V_OTHER(V))

#define GT_STRGRAPH_V_IS_INTERNAL(STRGRAPH, V) \
  ((GT_STRGRAPH_V_OUTDEG(STRGRAPH, V) == (GtStrgraphVEdgenum)1) && \
   (GT_STRGRAPH_V_INDEG(STRGRAPH, V) == (GtStrgraphVEdgenum)1))

#define GT_STRGRAPH_V_IS_JUNCTION(STRGRAPH, I) \
  ((GT_STRGRAPH_V_OUTDEG(STRGRAPH, I) > (GtStrgraphVEdgenum)1 && \
    GT_STRGRAPH_V_INDEG(STRGRAPH, I) > 0) || \
   (GT_STRGRAPH_V_OUTDEG(STRGRAPH, I) == (GtStrgraphVEdgenum)1 && \
    GT_STRGRAPH_V_INDEG(STRGRAPH, I) > (GtStrgraphVEdgenum)1))

#define GT_STRGRAPH_V_NOFEDGES(STRGRAPH, V) \
   (GtStrgraphVEdgenum)(GT_STRGRAPH_V_OFFSET(STRGRAPH, (V)+1) \
   - GT_STRGRAPH_V_OFFSET(STRGRAPH, V))

/* begin/end vertices: */

#define GT_STRGRAPH_V_B(READNUM) \
  ((GtStrgraphVnum)(READNUM) << 1)

#define GT_STRGRAPH_V_E(READNUM) \
  (GT_STRGRAPH_V_B(READNUM) + (GtStrgraphVnum)1)

#define GT_STRGRAPH_V_READNUM(V) \
  (unsigned long)((V) >> 1)

#define GT_STRGRAPH_V_IS_E(V) \
  (((V) & (GtStrgraphVnum)1) == (GtStrgraphVnum)1)

#define GT_STRGRAPH_V_IS_B(V) \
  (!(GT_STRGRAPH_V_IS_E(V)))

#define GT_STRGRAPH_V_OTHER(V) \
  (GT_STRGRAPH_V_IS_E(V) ? (V) - (GtStrgraphVnum)1 : (V) + (GtStrgraphVnum)1)

#define GT_STRGRAPH_V_CHAR(V) \
  (GT_STRGRAPH_V_IS_E(V) ? 'E' : 'B')

#endif
