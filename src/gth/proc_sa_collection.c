/*
  Copyright (c) 2003-2010 Gordon Gremme <gremme@zbh.uni-hamburg.de>
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

#include "core/undef_api.h"
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

void gth_proc_sa_collection(GthSACollection *sa_collection,
                            GthCallInfo *call_info,
                            GthInput *input, GthStat *stat,
                            unsigned int indentlevel)
{
  if (call_info->out->showverbose)
    call_info->out->showverbose("output spliced alignments");

  gth_input_delete_current(input);

  /* output alignments */
  if (!call_info->out->skipalignmentout) {
    GthSAVisitor *sa_visitor;
    if (call_info->out->xmlout) {
      if (call_info->intermediate) {
        sa_visitor = gth_xml_inter_sa_visitor_new(input, indentlevel + 1,
                                                  call_info->out->outfp);
      }
      else {
        sa_visitor = gth_xml_final_sa_visitor_new(input,
                                                  call_info->dp_options_core
                                                  ->dpminintronlength,
                                                  call_info->translationtable,
                                                  indentlevel,
                                                  call_info->out->outfp);
      }
    }
    else if (call_info->out->gff3out) {
      sa_visitor = gth_gff3_sa_visitor_new(input,
                                           call_info->out->gff3descranges,
                                           call_info->out->outfp);
    }
    else {
      sa_visitor = gth_txt_sa_visitor_new(input, call_info->out->gs2out,
                                      call_info->dp_options_core
                                      ->dpminintronlength,
                                      call_info->out->widthforgenpos,
                                      call_info->out->showintronmaxlen,
                                      call_info->translationtable,
                                      call_info->out->showseqnums,
                                      call_info->out->outfp);
    }
    gth_sa_collection_traverse(sa_collection, sa_visitor, input);
    gth_sa_visitor_delete(sa_visitor);
  }

  if (gth_sa_collection_contains_sa(sa_collection) &&
      !call_info->intermediate) { /* do not stop after SA computation */
    GthPGLCollection *pgl_collection;
    GthPGLVisitor *pgl_visitor;
    if (call_info->out->showverbose)
      call_info->out->showverbose("compute predicted gene locations");

    /* compute PGLs */
    pgl_collection = gth_pgl_collection_new(sa_collection,
                                            call_info->disableclustersas);
    if (call_info->out->sortags)
      gth_pgl_collection_sortAGSs(pgl_collection, call_info->out->sortagswf);

    if (call_info->out->maxagsnum != GT_UNDEF_UINT)
      gth_pgl_collection_set_max_ags(pgl_collection, call_info->out->maxagsnum);

    if (call_info->out->showverbose)
      call_info->out->showverbose("output predicted gene locations");

    /* save memory statistics for PGLs */
    gth_stat_increase_numofPGLs_stored(stat,
                                       gth_pgl_collection_size(pgl_collection));

    /* output PGLs */
    if (call_info->out->xmlout) {
      pgl_visitor = gth_xml_pgl_visitor_new(input, call_info->translationtable,
                                            indentlevel + 1, call_info->out);
    }
    else if (call_info->out->gff3out) {
      pgl_visitor = gth_gff3_pgl_visitor_new(input,
                                             call_info->out->gff3descranges,
                                             call_info->out->minORFlength,
                                             call_info->out->start_codon,
                                             call_info->out->final_stop_codon,
                                             call_info->out->outfp);
    }
    else {
      pgl_visitor = gth_txt_pgl_visitor_new(input, call_info->translationtable,
                                            indentlevel + 1, call_info->out);
    }
    gth_pgl_collection_traverse(pgl_collection, pgl_visitor, input,
                                call_info->out->gff3descranges);
    gth_pgl_visitor_delete(pgl_visitor);

    /* free*/
    gth_pgl_collection_delete(pgl_collection);
  }

  /* compute exon and intron length distributions */
  calc_sa_distributions(stat, sa_collection);
}
