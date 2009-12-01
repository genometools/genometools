/*
  Copyright (c) 2003-2009 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2003-2008 Center for Bioinformatics, University of Hamburg

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

#ifndef AGS_H
#define AGS_H

#include "core/array.h"
#include "core/error.h"
#include "gth/gthoutput.h"
#include "gth/bssm_param.h"
#include "gth/pgl.h"
#include "gth/sa.h"

#define SHOWGENPOSAGS(P)\
        SHOWGENPOS(gth_ags_is_forward(ags), gth_ags_total_length(ags),\
                   gth_ags_genomic_offset(ags), P)

typedef struct {
  GtRange range; /* the borders of the exon in the genomic sequence */
  GthDbl score;
} GthExonAGS;

typedef struct {
  GthFlt donorsiteprob,
         acceptorsiteprob;
} GthSpliceSiteProb;

typedef struct GthAGSObject GthAGSObject;

/* the alternative gene structure (AGS) class */
typedef struct GthAGS {
  GtStr *gen_id;                 /* the id of the genomic sequence this two
                                    values could also be included in the PGL
                                    structure, but including them here is better
                                    for the stand alone version of the assemble
                                    program */
  GtArray *exons,                /* contains the actual structure of the AGS and
                                    the exonscores */
          *splicesiteprobs,      /* contains the splice site probabilities */
          *alignments;           /* pointer to the generating spliced alignments
                                    I.e., the alignments which constitute this
                                    alternative gene structure. */
  unsigned long numofstoredsaclusters; /* number of stored SA clusters */
                                 /* (needed in assembly phase) */
  GthDbl overallscore; /* overall score used for sorting of AGSs */

  GthAGSObject *agso;
} GthAGS;

GthAGS*       gth_ags_new(const GthPGL*);
void          gth_ags_delete(GthAGS*);
bool          gth_ags_is_forward(const GthAGS*);
unsigned long gth_ags_filenum(const GthAGS*);
unsigned long gth_ags_total_length(const GthAGS*);
unsigned long gth_ags_genomic_offset(const GthAGS*);
GtStr*        gth_ags_get_gen_id(const GthAGS*);
GthExonAGS*   gth_ags_get_exon(const GthAGS *ags, unsigned long exon);
unsigned long gth_ags_num_of_exons(const GthAGS *ags);
GtStrand      gth_ags_genomic_strand(const GthAGS*);
GtRange       gth_ags_donor_site_range(const GthAGS*, unsigned long intron);
GtRange       gth_ags_acceptor_site_range(const GthAGS*, unsigned long intron);
double        gth_ags_donor_site_prob(const GthAGS*, unsigned long intron);
double        gth_ags_acceptor_site_prob(const GthAGS*, unsigned long intron);

#endif
