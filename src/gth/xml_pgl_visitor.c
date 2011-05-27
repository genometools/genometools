/*
  Copyright (c) 2004-2009 Gordon Gremme <gremme@zbh.uni-hamburg.de>
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

#include "core/undef_api.h"
#include "core/unused_api.h"
#include "gth/indent.h"
#include "gth/gthtrans.h"
#include "gth/pgl_visitor_rep.h"
#include "gth/xml_pgl_visitor.h"

struct GthXMLPGLVisitor {
  const GthPGLVisitor parent_instance;
  GthInput *input;
  unsigned long translationtable;
  unsigned int indentlevel;
  GthOutput *out;
};

#define xml_pgl_visitor_cast(GV)\
        gth_pgl_visitor_cast(gth_xml_pgl_visitor_class(), GV)

static void xml_outputAGSline(const GthAGS *ags, unsigned long agsnum,
                              unsigned int indentlevel, GtFile *outfp)
{
  GthExonAGS *exon;
  unsigned long i;

  gth_indent(outfp, indentlevel);
  gt_file_xprintf(outfp, "<AGS_line AGS_serial=\"%lu\">\n",
                  agsnum + OUTPUTOFFSET);
  indentlevel++;
  gth_indent(outfp, indentlevel);
  gt_file_xprintf(outfp, "<exon_coordinates>\n");
  indentlevel++;

  for (i = 0; i < gth_ags_num_of_exons(ags); i++) {
    exon = gth_ags_get_exon(ags, i);
    gth_indent(outfp, indentlevel);
    gt_file_xprintf(outfp, "<exon e_start=\"%lu\" e_stop=\"%lu\"/>\n",
                    SHOWGENPOSAGS(exon->range.start),
                    SHOWGENPOSAGS(exon->range.end));
  }

  indentlevel--;
  gth_indent(outfp, indentlevel);
  gt_file_xprintf(outfp, "</exon_coordinates>\n");
  indentlevel--;
  gth_indent(outfp, indentlevel);
  gt_file_xprintf(outfp, "</AGS_line>\n");
}

static void xml_outputSCRline(const GthAGS *ags, unsigned int indentlevel,
                              GtFile *outfp)
{
  GthSpliceSiteProb *splicesiteprob;
  unsigned long i;

  gth_indent(outfp, indentlevel);
  gt_file_xprintf(outfp, "<SCR_line>\n");
  indentlevel++;

  for (i = 0; i < gt_array_size(ags->exons) - 1; i++) {
    splicesiteprob = (GthSpliceSiteProb*) gt_array_get(ags->splicesiteprobs, i);
    gth_indent(outfp, indentlevel);
    gt_file_xprintf(outfp, "<exon-intron don_prob=\"%.3f\" "
                       "acc_prob=\"%.3f\" e_score=\"%.3f\"/>\n",
                       splicesiteprob->donorsiteprob,
                       splicesiteprob->acceptorsiteprob,
                       ((GthExonAGS*) gt_array_get(ags->exons, i))->score);
  }

  gth_indent(outfp, indentlevel);
  gt_file_xprintf(outfp, "<exon-only e_score=\"%.3f\"/>\n",
                  ((GthExonAGS*) gt_array_get(ags->exons, i))->score);
  indentlevel--;
  gth_indent(outfp, indentlevel);
  gt_file_xprintf(outfp, "</SCR_line>\n");
}

static void xml_output_exon_intron_lines(const GthAGS *ags,
                                         unsigned int indentlevel,
                                         GtFile *outfp)
{
  GthSpliceSiteProb *splicesiteprob;
  GthExonAGS *exon;
  unsigned long i, leftexonborder, rightexonborder, exonlength,
                leftintronborder = GT_UNDEF_ULONG, rightintronborder,
                intronlength;
  GthDbl exonscore;
  GthFlt donorsiteprob, acceptorsiteprob;

  gth_indent(outfp, indentlevel);
  gt_file_xprintf(outfp,
                  "<exon-intron_info xmlns=\"http://www.genomethreader.org/"
                  "GTH_output/PGL_module/predicted_gene_location/"
                  "AGS_information/exon-intron_info/\">\n");
  indentlevel++;

  for (i = 0; i < gt_array_size(ags->exons); i++) {
    exon            = (GthExonAGS*) gt_array_get(ags->exons, i);
    leftexonborder  = exon->range.start;
    rightexonborder = exon->range.end;
    exonlength      = rightexonborder - leftexonborder + 1;
    exonscore       = exon->score;

    if (i > 0) {
      rightintronborder = leftexonborder - 1;
      intronlength      = rightintronborder - leftintronborder + 1;
      splicesiteprob    = (GthSpliceSiteProb*)
                          gt_array_get(ags->splicesiteprobs, i-1);
      donorsiteprob     = splicesiteprob->donorsiteprob;
      acceptorsiteprob  = splicesiteprob->acceptorsiteprob;

      /* output intron */
      gth_indent(outfp, indentlevel);
      gt_file_xprintf(outfp, "<intron i_serial=\"%lu\" don_prob=\"%.3f\" "
                      "acc_prob=\"%.3f\">\n",  i - 1 + OUTPUTOFFSET,
                      donorsiteprob, acceptorsiteprob);
      indentlevel++;
      gth_indent(outfp, indentlevel);
      gt_file_xprintf(outfp,
                      "<gDNA_intron_boundary i_start=\"%lu\" i_stop=\"%lu\" "
                      "i_length=\"%lu\"/>\n",
                      SHOWGENPOSAGS(leftintronborder),
                      SHOWGENPOSAGS(rightintronborder),  intronlength);
      indentlevel--;
      gth_indent(outfp, indentlevel);
      gt_file_xprintf(outfp, "</intron>\n");
    }
    leftintronborder = rightexonborder + 1;

    /* output exon */
    gth_indent(outfp, indentlevel);
    gt_file_xprintf(outfp, "<exon e_serial=\"%lu\" e_score=\"%.3f\">\n",
                    i + OUTPUTOFFSET, exonscore);
    indentlevel++;
    gth_indent(outfp, indentlevel);
    gt_file_xprintf(outfp,
                    "<gDNA_exon_boundary e_start=\"%lu\" e_stop=\"%lu\" "
                    "e_length=\"%lu\"/>\n", SHOWGENPOSAGS(leftexonborder),
                    SHOWGENPOSAGS(rightexonborder), exonlength);
    indentlevel--;
    gth_indent(outfp, indentlevel);
    gt_file_xprintf(outfp, "</exon>\n");
  }

  indentlevel--;
  gth_indent(outfp, indentlevel);
  gt_file_xprintf(outfp, "</exon-intron_info>\n");
}

static void xml_outputPGSlines(GtArray *alignments, unsigned int indentlevel,
                               GtFile *outfp)
{
  unsigned long i, j;
  GthSA *sa;

  gth_indent(outfp, indentlevel);
  gt_file_xprintf(outfp, "<supporting_evidence xmlns=\""
                  "http://www.genomethreader.org/"
                  "GTH_output/PGL_module/predicted_gene_location/"
                  "AGS_information/supporting_evidence/\">\n");
  indentlevel++;

  for (i = 0; i < gt_array_size(alignments); i++) {
    sa = *(GthSA**) gt_array_get(alignments, i);

    gth_indent(outfp, indentlevel);
    gt_file_xprintf(outfp, "<PGS_line>\n");
    indentlevel++;
    gth_indent(outfp, indentlevel);
    gt_file_xprintf(outfp, "<gDNA_exon_coordinates>\n");
    indentlevel++;

    for (j = 0; j < gth_sa_num_of_exons(sa); j++) {
      gth_indent(outfp, indentlevel);
      gt_file_xprintf(outfp, "<exon start=\"%lu\" stop=\"%lu\"/>\n",
                      gth_sa_left_genomic_exon_border(sa, j),
                      gth_sa_right_genomic_exon_border(sa, j));
    }
    indentlevel--;
    gth_indent(outfp, indentlevel);
    gt_file_xprintf(outfp, "</gDNA_exon_coordinates>\n");
    gth_indent(outfp, indentlevel);
    if (gth_sa_alphatype(sa) == DNA_ALPHA) {
      gt_file_xprintf(outfp, "<referenceDNA id=\"%s\" strand=\"%c\"/>\n",
                         gth_sa_ref_id(sa),
                         gth_sa_ref_strand_char(sa));
    }
    else {
      gt_file_xprintf(outfp, "<referenceProtein id=\"%s\"/>\n",
                         gth_sa_ref_id(sa));
    }
    indentlevel--;
    gth_indent(outfp, indentlevel);
    gt_file_xprintf(outfp, "</PGS_line>\n");
  }

  indentlevel--;
  gth_indent(outfp, indentlevel);
  gt_file_xprintf(outfp, "</supporting_evidence>\n");
}

static void xml_show_ags(const GthAGS *ags, unsigned long pglnum,
                         unsigned long agsnum, unsigned long translationtable,
                         GthInput *input, unsigned int indentlevel,
                         GthOutput *out)
{
  gth_indent(out->outfp, indentlevel);
  gt_file_xprintf(out->outfp, "<AGS_information>\n");
  indentlevel++;

  /* output AGS line */
  xml_outputAGSline(ags, agsnum, indentlevel, out->outfp);

  /* output SCR line */
  xml_outputSCRline(ags, indentlevel, out->outfp);

  /* output exon/intron lines */
  xml_output_exon_intron_lines(ags, indentlevel, out->outfp);

  /* output PGS lines */
  xml_outputPGSlines(ags->alignments, indentlevel, out->outfp);

  /* output 3-phase translation */
  gt_outputtranslationandorf(pglnum, ags, agsnum, translationtable, input,
                          indentlevel, out);

  indentlevel--;
  gth_indent(out->outfp, indentlevel);
  gt_file_xprintf(out->outfp, "</AGS_information>\n");
}

static void xml_show_pgl(GthPGL *pgl, unsigned long pglnum,
                         unsigned long translationtable, GthInput *input,
                         unsigned int indentlevel, GthOutput *out)
{
  unsigned long i;

  gth_indent(out->outfp, indentlevel);
  gt_file_xprintf(out->outfp, "<predicted_gene_location>\n");
  indentlevel++;
  gth_indent(out->outfp, indentlevel);
  gt_file_xprintf(out->outfp,
                     "<PGL_line PGL_serial=\"%lu\" PGL_strand=\"%c\" "
                     "PGL_start=\"%lu\" PGL_stop=\"%lu\"/>\n",
                     pglnum + OUTPUTOFFSET,
                     SHOWSTRAND(gth_pgl_is_forward(pgl)),
                     SHOWGENPOS(gth_pgl_is_forward(pgl),
                                gth_pgl_total_length(pgl),
                                gth_pgl_genomic_offset(pgl),
                                pgl->maxrange.start),
                     SHOWGENPOS(gth_pgl_is_forward(pgl),
                                gth_pgl_total_length(pgl),
                                gth_pgl_genomic_offset(pgl),
                                pgl->maxrange.end));

  for (i = 0; i < gth_pgl_num_of_ags(pgl); i++) {
    xml_show_ags(gth_pgl_get_ags(pgl, i), pglnum, i, translationtable, input,
                 indentlevel, out);
  }

  indentlevel--;
  gth_indent(out->outfp, indentlevel);
  gt_file_xprintf(out->outfp, "</predicted_gene_location>\n");
}

static void xml_pgl_visitor_preface(GthPGLVisitor *pgl_visitor,
                                    GT_UNUSED unsigned long num_of_pgls)
{
  GthXMLPGLVisitor *visitor = xml_pgl_visitor_cast(pgl_visitor);
  gth_indent(visitor->out->outfp, visitor->indentlevel);
  gt_file_xprintf(visitor->out->outfp,
                     "<PGL_module xmlns=\"http://www.genomethreader."
                     "org/GTH_output/PGL_module/\">\n");
}

static void xml_pgl_visitor_visit_pgl(GthPGLVisitor *pgl_visitor,
                                      GthPGL *pgl, unsigned long pglnum)
{
  GthXMLPGLVisitor *visitor = xml_pgl_visitor_cast(pgl_visitor);
  gt_assert(pgl);
  xml_show_pgl(pgl, pglnum, visitor->translationtable, visitor->input,
           visitor->indentlevel, visitor->out);
}

static void xml_pgl_visitor_trailer(GthPGLVisitor *pgl_visitor)
{
  GthXMLPGLVisitor *visitor = xml_pgl_visitor_cast(pgl_visitor);
  gth_indent(visitor->out->outfp, visitor->indentlevel);
  gt_file_xprintf(visitor->out->outfp, "</PGL_module>\n");
}

const GthPGLVisitorClass* gth_xml_pgl_visitor_class()
{
  static const GthPGLVisitorClass pglvc = { sizeof (GthXMLPGLVisitor),
                                            NULL,
                                            xml_pgl_visitor_preface,
                                            NULL,
                                            xml_pgl_visitor_visit_pgl,
                                            xml_pgl_visitor_trailer };
  return &pglvc;
}

GthPGLVisitor* gth_xml_pgl_visitor_new(GthInput *input,
                                       unsigned long translationtable,
                                       unsigned int indentlevel, GthOutput *out)
{
  GthPGLVisitor *pgl_visitor =
    gth_pgl_visitor_create(gth_xml_pgl_visitor_class());
  GthXMLPGLVisitor *visitor = xml_pgl_visitor_cast(pgl_visitor);
  visitor->input = input;
  visitor->translationtable = translationtable;
  visitor->indentlevel = indentlevel;
  visitor->out = out;
  return pgl_visitor;
}
