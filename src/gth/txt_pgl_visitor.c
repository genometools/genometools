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

#include "core/undef_api.h"
#include "gth/indent.h"
#include "gth/gthtrans.h"
#include "gth/pgl_visitor_rep.h"
#include "gth/txt_pgl_visitor.h"

#define PGLS_DELIMITERCHAR     '-'
#define PGL_DELIMITERLINECHAR  '*'

struct GthTxtPGLVisitor {
  const GthPGLVisitor parent_instance;
  GthInput *input;
  unsigned long translationtable;
  unsigned int indentlevel;
  GthOutput *out;
};

#define txt_pgl_visitor_cast(GV)\
        gth_pgl_visitor_cast(gth_txt_pgl_visitor_class(), GV)

static void outputAGSline(const GthAGS *ags, unsigned long agsnum,
                          GtFile *outfp)
{
  GthExonAGS *exon;
  unsigned long i;

  gt_file_xprintf(outfp, "AGS-%lu (",  agsnum + OUTPUTOFFSET);
  for (i = 0; i < gth_ags_num_of_exons(ags); i++) {
    exon = gth_ags_get_exon(ags, i);
    if (i > 0)
      gt_file_xfputc(',', outfp);
    gt_file_xprintf(outfp, "%lu  %lu", SHOWGENPOSAGS(exon->range.start),
                    SHOWGENPOSAGS(exon->range.end));
  }
  gt_file_xprintf(outfp, ")\n");
}

static void outputSCRline(const GthAGS *ags, GtFile *outfp)
{
  GthSpliceSiteProb *splicesiteprob;
  unsigned long i;

  gt_file_xprintf(outfp, "SCR   (");
  for (i = 0; i < gt_array_size(ags->exons) - 1; i++) {
    splicesiteprob = (GthSpliceSiteProb*) gt_array_get(ags->splicesiteprobs, i);
    gt_file_xprintf(outfp, "e %5.3f  d %5.3f a %5.3f,",
                    ((GthExonAGS*) gt_array_get(ags->exons, i))->score,
                    splicesiteprob->donorsiteprob,
                    splicesiteprob->acceptorsiteprob);
  }
  gt_file_xprintf(outfp, "e %5.3f)\n",
                  ((GthExonAGS*) gt_array_get(ags->exons, i))->score);
  gt_file_xfputc('\n', outfp);
}

static void output_exon_intron_lines(const GthAGS *ags, int widthforgenpos,
                                     GtFile *outfp)
{
  GthSpliceSiteProb *splicesiteprob;
  GthExonAGS *exon;
  unsigned long i, leftexonborder, rightexonborder, exonlength,
                leftintronborder = GT_UNDEF_ULONG, rightintronborder,
                intronlength;
  GthDbl exonscore;
  GthFlt donorsiteprob, acceptorsiteprob;

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
      gt_file_xprintf(outfp,"    Intron %2lu %*lu %*lu (%4lu n);           "
                      "Pd: %5.3f  Pa: %5.3f\n",  i - 1 + OUTPUTOFFSET,
                      widthforgenpos, SHOWGENPOSAGS(leftintronborder),
                      widthforgenpos, SHOWGENPOSAGS(rightintronborder),
                      intronlength, donorsiteprob, acceptorsiteprob);
    }
    leftintronborder = rightexonborder + 1;

    /* output exon */
    gt_file_xprintf(outfp, "  Exon %2lu %*lu %*lu (%4lu n); score: %5.3f\n",
                    i + OUTPUTOFFSET, widthforgenpos,
                    SHOWGENPOSAGS(leftexonborder),
                    widthforgenpos,
                    SHOWGENPOSAGS(rightexonborder), exonlength, exonscore);
  }
  gt_file_xfputc('\n', outfp);
}

static void outputPGSlines(GtArray *alignments, GtFile *outfp)
{
  unsigned long i, j;
  GthSA *sa;

  for (i = 0; i < gt_array_size(alignments); i++) {
    sa = *(GthSA**) gt_array_get(alignments, i);

    gt_file_xprintf(outfp, "  PGS (");
    for (j = 0; j < gth_sa_num_of_exons(sa); j++) {
      if (j > 0)
        gt_file_xfputc(',', outfp);
      gt_file_xprintf(outfp, "%lu  %lu",
                      gth_sa_left_genomic_exon_border(sa, j),
                      gth_sa_right_genomic_exon_border(sa, j));
    }
    gt_file_xprintf(outfp, ")\t%s%c\n", gth_sa_ref_id(sa),
                       gth_sa_ref_strand_char(sa));
  }

  gt_file_xfputc('\n', outfp);
}

static void show_ags(const GthAGS *ags, unsigned long pglnum,
                     unsigned long agsnum, unsigned long translationtable,
                     GthInput *input, unsigned int indentlevel, GthOutput *out)
{
  GtFile *outfp = out->outfp;

  /* output AGS line */
  outputAGSline(ags, agsnum, out->outfp);

  /* output SCR line */
  outputSCRline(ags, out->outfp);

  /* output exon/intron lines */
  output_exon_intron_lines(ags, out->widthforgenpos, out->outfp);

  /* output PGS lines */
  outputPGSlines(ags->alignments, out->outfp);

  /* output 3-phase translation */
  gt_outputtranslationandorf(pglnum, ags, agsnum, translationtable, input,
                             indentlevel, out);

  /* output three final newlines */
  gt_file_xprintf(outfp, "\n\n\n");
}

static void show_pgl(GthPGL *pgl, unsigned long pglnum,
                     unsigned long translationtable, GthInput *input,
                     unsigned int indentlevel, GthOutput *out)
{
  unsigned long i;
  GtFile *outfp = out->outfp;

  gt_assert(!out->gff3out);

  if (out->xmlout) {
    gth_indent(outfp, indentlevel);
    gt_file_xprintf(outfp, "<predicted_gene_location>\n");
    indentlevel++;
    gth_indent(outfp, indentlevel);
    gt_file_xprintf(outfp, "<PGL_line PGL_serial=\"%lu\" PGL_strand=\"%c\" "
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
  }
  else {
    gt_file_xprintf(outfp, "PGL %3lu (%c strand):      %lu     %lu",
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
    if (out->pglgentemplate)
      gt_file_xprintf(outfp, " (genomic template '%s')", gth_pgl_gen_id(pgl));
    gt_file_xfputc('\n', outfp);
  }

  for (i = 0; i < gt_array_size(pgl->assemblies); i++) {
    show_ags(gth_pgl_get_ags(pgl, i), pglnum, i, translationtable, input,
             indentlevel, out);
  }

  if (out->xmlout) {
    indentlevel--;
    gth_indent(outfp, indentlevel);
    gt_file_xprintf(outfp, "</predicted_gene_location>\n");
  }
}

static void txt_pgl_visitor_preface(GthPGLVisitor *pgl_visitor,
                                    unsigned long num_of_pgls)
{
  unsigned long i;
  GthTxtPGLVisitor *visitor = txt_pgl_visitor_cast(pgl_visitor);
  for (i = 0; i < DELIMITERLINELENGTH; i++)
    gt_file_xfputc(PGLS_DELIMITERCHAR, visitor->out->outfp);
  gt_file_xprintf(visitor->out->outfp, "\n\n");
  gt_file_xprintf(visitor->out->outfp,
                     "Predicted gene locations (%lu):\n\n\n", num_of_pgls);
}

static void txt_pgl_visitor_visit_pgl(GthPGLVisitor *pgl_visitor,
                                      GthPGL *pgl, unsigned long pglnum)
{
  GthTxtPGLVisitor *visitor = txt_pgl_visitor_cast(pgl_visitor);
  gt_assert(pgl);
  show_pgl(pgl, pglnum, visitor->translationtable, visitor->input,
           visitor->indentlevel, visitor->out);
}

const GthPGLVisitorClass* gth_txt_pgl_visitor_class()
{
  static const GthPGLVisitorClass pglvc = { sizeof (GthTxtPGLVisitor),
                                            NULL,
                                            txt_pgl_visitor_preface,
                                            NULL,
                                            txt_pgl_visitor_visit_pgl,
                                            NULL };
  return &pglvc;
}

GthPGLVisitor* gth_txt_pgl_visitor_new(GthInput *input,
                                       unsigned long translationtable,
                                       unsigned int indentlevel,
                                       GthOutput *out)
{
  GthPGLVisitor *pgl_visitor =
    gth_pgl_visitor_create(gth_txt_pgl_visitor_class());
  GthTxtPGLVisitor *visitor = txt_pgl_visitor_cast(pgl_visitor);
  visitor->input = input;
  visitor->translationtable = translationtable;
  visitor->indentlevel = indentlevel;
  visitor->out = out;
  return pgl_visitor;
}
