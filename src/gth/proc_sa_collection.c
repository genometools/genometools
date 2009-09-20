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

#include "extended/gff3_visitor.h"
#include "gth/indent.h"
#include "gth/proc_sa_collection.h"
#include "gth/gthsadistri.h"
#include "gth/pgl_collection.h"
#include "gth/gff3_pgl_visitor.h"
#include "gth/gff3_sa_visitor.h"
#include "gth/txt_pgl_visitor.h"
#include "gth/txt_sa_visitor.h"
#include "gth/xml_final_sa_visitor.h"
#include "gth/xml_inter_sa_visitor.h"
#include "gth/xml_pgl_visitor.h"

static void calc_sa_distributions(GthStat *stat, GthSACollection *sa_collection)
{
  gt_assert(stat);
  /* compute exon and intron length distributions */
  if (gth_stat_get_exondistri(stat) || gth_stat_get_introndistri(stat)) {
    gthcalcSAdistributions(gth_stat_get_exondistri(stat),
                           gth_stat_get_introndistri(stat),
                           gth_stat_get_exondistribution(stat),
                           gth_stat_get_introndistribution(stat),
                           sa_collection);
  }
}

void proc_sa_collection(GthSACollection *sa_collection, GthCallInfo *callinfo,
                        GthInput *input, GthStat *stat,
                        unsigned int indentlevel)
{
  if (callinfo->out->showverbose)
    callinfo->out->showverbose("output spliced alignments");

  /* output alignments */
  if (!callinfo->out->skipalignmentout) {
    GthSAVisitor *sa_visitor;
    if (callinfo->out->xmlout) {
      if (callinfo->intermediate) {
        sa_visitor = gth_xml_inter_sa_visitor_new(input, indentlevel + 1,
                                                  callinfo->out->outfp);
      }
      else {
        sa_visitor = gth_xml_final_sa_visitor_new(input,
                                                  callinfo->dp_options_core
                                                  ->dpminintronlength,
                                                  callinfo->translationtable,
                                                  indentlevel,
                                                  callinfo->out->outfp);
      }
    }
    else if (callinfo->out->gff3out) {
      sa_visitor = gth_gff3_sa_visitor_new(input, callinfo->out->outfp);
    }
    else {
      sa_visitor = gth_txt_sa_visitor_new(input, callinfo->out->gs2out,
                                      callinfo->dp_options_core
                                      ->dpminintronlength,
                                      callinfo->out->widthforgenpos,
                                      callinfo->out->showintronmaxlen,
                                      callinfo->translationtable,
                                      callinfo->out->showseqnums,
                                      callinfo->out->outfp);
    }
    gth_sa_collection_traverse(sa_collection, sa_visitor);
    gth_sa_visitor_delete(sa_visitor);
  }

  if (gth_sa_collection_contains_sa(sa_collection) &&
      !callinfo->intermediate) { /* do not stop after SA computation */
    GthPGLCollection *pgl_collection;
    GthPGLVisitor *pgl_visitor;
    if (callinfo->out->showverbose)
      callinfo->out->showverbose("compute predicted gene locations");

    /* compute PGLs */
    pgl_collection = gth_pgl_collection_new(sa_collection,
                                            callinfo->disableclustersas);
    if (callinfo->out->sortags)
      gth_pgl_collection_sortAGSs(pgl_collection, callinfo->out->sortagswf);

    if (callinfo->out->showverbose)
      callinfo->out->showverbose("output predicted gene locations");

    /* save memory statistics for PGLs */
    gth_stat_increase_numofPGLs_stored(stat,
                                       gth_pgl_collection_size(pgl_collection));

    /* output PGLs */
    if (callinfo->out->xmlout) {
      pgl_visitor = gth_xml_pgl_visitor_new(input, callinfo->translationtable,
                                            indentlevel + 1, callinfo->out);
    }
    else if (callinfo->out->gff3out) {
      pgl_visitor = gth_gff3_pgl_visitor_new(input, callinfo->out->outfp);
    }
    else {
      pgl_visitor = gth_txt_pgl_visitor_new(input, callinfo->translationtable,
                                            indentlevel + 1, callinfo->out);
    }
    gth_pgl_collection_traverse(pgl_collection, pgl_visitor);
    gth_pgl_visitor_delete(pgl_visitor);

    /* free*/
    gth_pgl_collection_delete(pgl_collection);
  }

  /* compute exon and intron length distributions */
  calc_sa_distributions(stat, sa_collection);
}
