/*
  Copyright (c) 2004-2011 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2004-2008 Center for Bioinformatics, University of Hamburg

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

#ifndef PGL_H
#define PGL_H

#include "core/array.h"
#include "core/error.h"
#include "extended/node_visitor_api.h"
#include "gth/sa.h"

typedef struct {
  GthSA *representative; /* the representative of this cluster */
  GtArray *members;      /* additional members of this cluster */
} GthSACluster;

typedef struct GthPGLObject GthPGLObject;

typedef struct {
  GtRange maxrange;      /* contains the leftmost and rightmost sequence
                            position of all AGSs in this PGL */
  GtArray *assemblies;   /* contains all assemblies (AGSs) */

  /* this goes into the assembly phase */
  GtArray *alignments;   /* contains the underlying alignments */

  /* this is computed in the assembly phase */
  GtArray *saclusters;

  GthPGLObject *pglo;
} GthPGL;

GthPGL*        gth_pgl_new(bool forward);
void           gth_pgl_delete(GthPGL*);
void           gth_pgl_add_sa(GthPGL*, GthSA*);
/* Returns alternative gene structure <i> from <pgl>. */
struct GthAGS* gth_pgl_get_ags(const GthPGL*, unsigned long i);
/* Returns the number of alternative gene structures in <pgl>. */
unsigned long  gth_pgl_num_of_ags(const GthPGL*);
/* Set the maximum number of AGSs (<maxagsnum>) which are allowed for <pgl>. */
void           gth_pgl_set_max_ags(GthPGL *pgl, unsigned int maxagsnum);
bool           gth_pgl_is_forward(const GthPGL*);
unsigned long  gth_pgl_filenum(const GthPGL*);
unsigned long  gth_pgl_seqnum(const GthPGL*);
unsigned long  gth_pgl_total_length(const GthPGL*);
unsigned long  gth_pgl_genomic_offset(const GthPGL*);
GtRange        gth_pgl_genomic_range(const GthPGL*);
GtStrand       gth_pgl_genomic_strand(const GthPGL*);
const char*    gth_pgl_gen_id(const GthPGL*);

#endif
