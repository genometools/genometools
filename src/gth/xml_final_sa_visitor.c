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

#include "core/chardef.h"
#include "core/undef_api.h"
#include "core/unused_api.h"
#include "gth/indent.h"
#include "gth/sa_visitor_rep.h"
#include "gth/xml_final_sa_visitor.h"

struct GthXMLFinalSAVisitor {
  const GthSAVisitor parent_instance;
  unsigned long minintronlength,
                translationtable;
  GthInput *input;
  unsigned int indentlevel;
  GtFile *outfp;
};

#define xml_final_sa_visitor_cast(GV)\
        gth_sa_visitor_cast(gth_xml_final_sa_visitor_class(), GV)

static void xml_showgthreferenceinformation(GthSA *sa,
                                            GthInput *input,
                                            unsigned int indentlevel,
                                            GtFile *outfp)
{
  gt_assert(gth_sa_ref_file_num(sa) != GT_UNDEF_ULONG);

  gth_indent(outfp, indentlevel);

  switch (gth_sa_alphatype(sa)) {
    case DNA_ALPHA:
      gt_file_xprintf(outfp, "<reference ref_file=\"%s\" ref_id=\"%s\" "
                                "ref_strand=\"%c\" ref_description=\"",
                         gth_input_get_reference_filename(input,
                                                  gth_sa_ref_file_num(sa)),
                         gth_sa_ref_id(sa),
                         gth_sa_ref_strand_char(sa));
      break;
    case PROTEIN_ALPHA:
      gt_file_xprintf(outfp, "<reference ref_file=\"%s\" ref_id=\"%s\" "
                                "ref_description=\"",
                         gth_input_get_reference_filename(input,
                                                  gth_sa_ref_file_num(sa)),
                         gth_sa_ref_id(sa));
      break;
    default: gt_assert(0);
  }

  gth_input_echo_reference_description(input, gth_sa_ref_file_num(sa),
                                       gth_sa_ref_seq_num(sa), outfp);

  gt_file_xprintf(outfp, "\">\n");
}

static void xml_showgthgenomicinformation(GthSA *sa,
                                          GthInput *input,
                                          unsigned int indentlevel,
                                          GtFile *outfp)
{
  gt_assert(gth_sa_gen_file_num(sa) != GT_UNDEF_ULONG);

  gth_indent(outfp, indentlevel);
  gt_file_xprintf(outfp, "<gDNA_segment>\n");
  indentlevel++;
  gth_indent(outfp, indentlevel);
  gt_file_xprintf(outfp, "<template temp_file=\"%s\" temp_id=\"%s\" "
                            "temp_strand=\"%c\" temp_description=\"",
                     gth_input_get_genomic_filename(input,
                                                    gth_sa_gen_file_num(sa)),
                     gth_sa_gen_id(sa),
                     gth_sa_gen_strand_char(sa));

  gth_input_echo_genomic_description(input, gth_sa_gen_file_num(sa),
                                     gth_sa_gen_seq_num(sa), outfp);

  gt_file_xprintf(outfp, "\">\n");
  indentlevel++;
  gth_indent(outfp, indentlevel);
  gt_file_xprintf(outfp, "<position start=\"%lu\" stop=\"%lu\"/>\n",
                     gth_sa_gen_dp_start_show(sa),
                     gth_sa_gen_dp_end_show(sa));
  indentlevel--;
  gth_indent(outfp, indentlevel);
  gt_file_xprintf(outfp, "</template>\n");
  indentlevel--;
  gth_indent(outfp, indentlevel);
  gt_file_xprintf(outfp, "</gDNA_segment>\n");
}

/*
  The following function prints a PPA line, which shows the start and end
  position of the poly-A tail in the cDNA (iff a poly-A tail could be
  determined).
*/
static void xml_showppaline(GthSA *sa, unsigned int indentlevel,
                           GtFile *outfp)
{
  if (gth_sa_polyAtail_start(sa) ||
      gth_sa_polyAtail_stop(sa)) {
    gth_indent(outfp, indentlevel);
    gt_file_xprintf(outfp,
                       "<PPA_line polyA_start=\"%lu\" polyA_stop=\"%lu\"/>\n",
                       gth_sa_polyAtail_start(sa) + OUTPUTOFFSET,
                       gth_sa_polyAtail_stop(sa) + OUTPUTOFFSET);
  }
}

/* The following function prints the "classic" GeneSeqer2 MATCH line */
static void xml_showmatchline(GthSA *sa, unsigned int indentlevel,
                              GtFile *outfp)
{
  gth_indent(outfp, indentlevel);
  gt_file_xprintf(outfp, "<MATCH_line gen_id=\"%s\" gen_strand=\"%c\" ",
                     gth_sa_gen_id(sa),
                     gth_sa_gen_strand_char(sa));
  if (gth_sa_alphatype(sa) == DNA_ALPHA) {
    gt_file_xprintf(outfp, "ref_id=\"%s\" ref_strand=\"%c\">\n",
                       gth_sa_ref_id(sa),
                       gth_sa_ref_strand_char(sa));
  }
  else
    gt_file_xprintf(outfp, "ref_id=\"%s\">\n", gth_sa_ref_id(sa));

  indentlevel++;
  gth_indent(outfp, indentlevel);
  gt_file_xprintf(outfp,
                     "<total_alignment_score>%.3f</total_alignment_score>\n",
                     gth_sa_score(sa));
  gth_indent(outfp, indentlevel);
  gt_file_xprintf(outfp, "<cumulative_length_of_scored_exons>%lu"
                     "</cumulative_length_of_scored_exons>\n",
                     gth_sa_cumlen_scored_exons(sa));
  gth_indent(outfp, indentlevel);
  gt_file_xprintf(outfp, "<coverage percentage=\"%.3f\" high_type=\"",
                     gth_sa_coverage(sa));
  gt_file_xfputc(gth_sa_coverage_char(sa), outfp);

  gt_file_xprintf(outfp, "\"/>\n");
  indentlevel--;
  gth_indent(outfp, indentlevel);
  gt_file_xprintf(outfp, "</MATCH_line>\n");
}

/*
  The following function prints the "classic" GeneSeqer2 PGS line
*/
static void xml_showpgsline(GthSA *sa, unsigned int indentlevel,
                            GtFile *outfp)
{
  unsigned long i, numofexons;
  gt_assert(sa);
  numofexons = gth_sa_num_of_exons(sa);
  gth_indent(outfp, indentlevel);
  gt_file_xprintf(outfp, "<PGS_line>\n");
  indentlevel++;
  gth_indent(outfp, indentlevel);
  gt_file_xprintf(outfp, "<gDNA gen_id=\"%s\" gen_strand=\"%c\"/>\n",
                     gth_sa_gen_id(sa), gth_sa_gen_strand_char(sa));
  gth_indent(outfp, indentlevel);
  if (gth_sa_alphatype(sa) == DNA_ALPHA) {
    gt_file_xprintf(outfp, "<rDNA rDNA_id=\"%s\" rDNA_strand=\"%c\"/>\n",
                       gth_sa_ref_id(sa), gth_sa_ref_strand_char(sa));
  }
  else {
    gt_file_xprintf(outfp, "<rProt rProt_id=\"%s\"/>\n",
                       gth_sa_ref_id(sa));
  }
  gth_indent(outfp, indentlevel);
  gt_file_xprintf(outfp, "<gDNA_exon_coordinates>\n");
  indentlevel++;

  for (i = 0; i < numofexons; i++) {
    gth_indent(outfp, indentlevel);
    gt_file_xprintf(outfp, "<exon e_start=\"%lu\" e_stop=\"%lu\"/>\n",
                    gth_sa_left_genomic_exon_border(sa, i),
                    gth_sa_right_genomic_exon_border(sa, i));
  }

  indentlevel--;
  gth_indent(outfp, indentlevel);
  gt_file_xprintf(outfp, "</gDNA_exon_coordinates>\n");
  indentlevel--;
  gth_indent(outfp, indentlevel);
  gt_file_xprintf(outfp, "</PGS_line>\n");
}

static void xml_showalignmentheader(GthSA *sa,
                                    unsigned long minintronlength,
                                    unsigned int indentlevel, GtFile *outfp)
{
  unsigned long i, leftreferenceexonborder, rightreferenceexonborder,
                referenceexonlength;
  GthDbl exonscore, donorsitescore, acceptorsitescore;
  GthFlt donorsiteprobability, acceptorsiteprobability;
  Exoninfo *exoninfo;
  Introninfo *introninfo;

  gth_indent(outfp, indentlevel);
  gt_file_xprintf(outfp, "<predicted_gene_structure>\n");
  indentlevel++;
  gth_indent(outfp, indentlevel);
  gt_file_xprintf(outfp, "<exon-intron_info>\n");
  indentlevel++;

  for (i = 0; i < gth_sa_num_of_exons(sa); i++) {
    exoninfo = gth_sa_get_exon(sa, i);
    leftreferenceexonborder  = exoninfo->leftreferenceexonborder;
    rightreferenceexonborder = exoninfo->rightreferenceexonborder;
    referenceexonlength      = rightreferenceexonborder
                               - leftreferenceexonborder + 1;
    exonscore                = exoninfo->exonscore;

    if (i > 0) {
      introninfo = gth_sa_get_intron(sa, i-1);
      donorsiteprobability    = introninfo->donorsiteprobability;
      donorsitescore          = introninfo->donorsitescore;
      acceptorsiteprobability = introninfo->acceptorsiteprobability;
      acceptorsitescore       = introninfo->acceptorsitescore;

      gth_indent(outfp, indentlevel);
      gt_file_xprintf(outfp, "<intron i_serial=\"%lu\">\n",
                      i - 1 + OUTPUTOFFSET);
      indentlevel++;
      gth_indent(outfp, indentlevel);
      gt_file_xprintf(outfp,
                      "<gDNA_intron_boundary i_start=\"%lu\" i_stop=\"%lu\" "
                      "i_length=\"%lu\">\n",
                      gth_sa_left_intron_border(sa, i-1),
                      gth_sa_right_intron_border(sa, i-1),
                      gth_sa_intron_length(sa, i-1));
      indentlevel++;
      gth_indent(outfp, indentlevel);
      gt_file_xprintf(outfp, "<donor d_prob=\"%.3f\"", donorsiteprobability);
      if (gth_sa_alphatype(sa) == DNA_ALPHA)
        gt_file_xprintf(outfp, " d_score=\"%.2f\"", donorsitescore);
      gt_file_xprintf(outfp, "/>\n");
      gth_indent(outfp, indentlevel);
      gt_file_xprintf(outfp, "<acceptor a_prob=\"%.3f\"",
                      acceptorsiteprobability);
      if (gth_sa_alphatype(sa) == DNA_ALPHA)
        gt_file_xprintf(outfp, " a_score=\"%.2f\"", acceptorsitescore);
      gt_file_xprintf(outfp, "/>\n");
      indentlevel--;
      gth_indent(outfp, indentlevel);
      gt_file_xprintf(outfp, "</gDNA_intron_boundary>\n");

      /* if the intron is shorter or equal than the minimal intron length an
         additional tag is shown */
      if (gth_sa_intron_length(sa, i-1) <= minintronlength) {
        gth_indent(outfp, indentlevel);
        gt_file_xprintf(outfp, "<shorter_than_min_intron_len/>\n");
      }

      indentlevel--;
      gth_indent(outfp, indentlevel);
      gt_file_xprintf(outfp, "</intron>\n");
    }

    gth_indent(outfp, indentlevel);
    gt_file_xprintf(outfp, "<exon e_serial=\"%lu\">\n", i + OUTPUTOFFSET);
    indentlevel++;
    gth_indent(outfp, indentlevel);
    gt_file_xprintf(outfp, "<gDNA_exon_boundary g_start=\"%lu\" g_stop="
                    "\"%lu\" g_length=\"%lu\"/>\n",
                    gth_sa_left_genomic_exon_border(sa, i),
                    gth_sa_right_genomic_exon_border(sa, i),
                    gth_sa_genomic_exon_length(sa, i));
    gth_indent(outfp, indentlevel);
    gt_file_xprintf(outfp,
                    "<reference_exon_boundary r_type=\"%s\" r_start=\"%lu\" "
                    "r_stop=\"%lu\" r_length=\"%lu\" r_score=\"%5.3f\"/>\n",
                    gth_sa_alphastring(sa),
                    leftreferenceexonborder  + OUTPUTOFFSET ,
                    rightreferenceexonborder + OUTPUTOFFSET ,
                    referenceexonlength, exonscore);
    indentlevel--;
    gth_indent(outfp, indentlevel);
    gt_file_xprintf(outfp, "</exon>\n");
  }

  indentlevel--;
  gth_indent(outfp, indentlevel);
  gt_file_xprintf(outfp, "</exon-intron_info>\n");

  /* showing PPA line (if an poly-A tail was determined) */
  if (gth_sa_alphatype(sa) == DNA_ALPHA)
    xml_showppaline(sa, indentlevel, outfp);

  /* showing MATCH line */
  xml_showmatchline(sa, indentlevel, outfp);

  /* showing PGS line */
  xml_showpgsline(sa, indentlevel, outfp);
}

static void showconcreteline(const unsigned char *alignmentline,
                             unsigned long cols, GtFile *outfp)
{
  unsigned long i;

  for (i = 0; i < cols; i++) {
    switch (alignmentline[i]) {
      case ABSTRACTGAPSYMBOL:
        gt_file_xfputc(CONCRETEGAPSYMBOL, outfp);
        break;
      case ABSTRACTINTRONSYMBOL:
        gt_file_xfputc(CONCRETEINTRONSYMBOL, outfp);
        break;
      default:
        gt_file_xfputc(alignmentline[i], outfp);
    }
  }
}

static void xml_final_show_spliced_alignment(GthSA *sa, GthInput *input,
                                             unsigned long minintronlength,
                                             unsigned long translationtable,
                                             unsigned int indentlevel,
                                             GtFile *outfp)
{
  unsigned char *first_line, *second_line, *third_line;
  GT_UNUSED bool reverse_subject_pos = false;
  unsigned long cols;

  gt_assert(sa && input);

  gth_indent(outfp, indentlevel);
  gt_file_xprintf(outfp,
                  "<spliced_alignment xmlns=\"http://www.genomethreader.org/"
                  "GTH_output/alignment_module/spliced_alignment/\">\n");
  indentlevel++;

  /* If the reverse complement of the genomic DNA is considered, this
     opition is needed for correct output of the genomic sequence positions
     by the function showalignmentgeneric() */
  if (!gth_sa_gen_strand_forward(sa))
    reverse_subject_pos = true;

  xml_showgthreferenceinformation(sa, input, indentlevel, outfp);

  indentlevel++;
  gth_indent(outfp, indentlevel);
  gt_file_xprintf(outfp, "<seq>");
  gth_sa_echo_reference_sequence(sa, input, false, outfp);
  gt_file_xprintf(outfp, "</seq>\n");
  indentlevel--;

  gth_indent(outfp, indentlevel);
  gt_file_xprintf(outfp, "</reference>\n");
  xml_showgthgenomicinformation(sa, input, indentlevel, outfp);

  xml_showalignmentheader(sa, minintronlength, indentlevel, outfp);

  indentlevel++;
  gth_indent(outfp, indentlevel);
  gt_file_xprintf(outfp, "<alignment>\n");

  /* compute the alignment lines */
  cols = gth_sa_get_alignment_lines(sa, &first_line, &second_line, &third_line,
                                    translationtable, input);

  indentlevel++;
  gth_indent(outfp, indentlevel);
  gt_file_xprintf(outfp, "<genome_strand>");
  showconcreteline(first_line, cols, outfp);
  gt_file_xprintf(outfp, "</genome_strand>\n");
  gth_indent(outfp, indentlevel);
  switch (gth_sa_alphatype(sa)) {
    case DNA_ALPHA:
      gt_file_xprintf(outfp, "<mrna_strand>");
      showconcreteline(second_line, cols, outfp);
      gt_file_xprintf(outfp, "</mrna_strand>\n");
      break;
    case PROTEIN_ALPHA:
      gt_file_xprintf(outfp, "<genomeProt>");
      showconcreteline(second_line, cols, outfp);
      gt_file_xprintf(outfp, "</genomeProt>\n");
      gth_indent(outfp, indentlevel);
      gt_file_xprintf(outfp, "<queryProt>");
      showconcreteline(third_line, cols, outfp);
      gt_file_xprintf(outfp, "</queryProt>\n");

      break;
    default: gt_assert(0);
  }

  indentlevel--;
  gth_indent(outfp, indentlevel);
  gt_file_xprintf(outfp, "</alignment>\n");
  indentlevel--;
  gth_indent(outfp, indentlevel);
  gt_file_xprintf(outfp, "</predicted_gene_structure>\n");
  indentlevel--;
  gth_indent(outfp, indentlevel);
  gt_file_xprintf(outfp, "</spliced_alignment>\n");

  /* free */
  gt_free(first_line);
  gt_free(second_line);
  gt_free(third_line);
}

static void xml_final_sa_visitor_preface(GthSAVisitor *sa_visitor)
{
  GthXMLFinalSAVisitor *visitor = xml_final_sa_visitor_cast(sa_visitor);
  gth_indent(visitor->outfp, visitor->indentlevel);
  gt_file_xprintf(visitor->outfp, "<alignment_module>\n");
}

static void xml_final_sa_visitor_visit_sa(GthSAVisitor *sa_visitor,
                                          GthSA *sa)
{
  GthXMLFinalSAVisitor *visitor = xml_final_sa_visitor_cast(sa_visitor);
  gt_assert(sa);
  xml_final_show_spliced_alignment(sa, visitor->input, visitor->minintronlength,
                                   visitor->translationtable,
                                   visitor->indentlevel, visitor->outfp);
}

static void  xml_final_sa_visitor_trailer(GthSAVisitor *sa_visitor,
                                          unsigned long num_of_sas)
{
  GthXMLFinalSAVisitor *visitor = xml_final_sa_visitor_cast(sa_visitor);
  visitor->indentlevel++;
  gth_indent(visitor->outfp, visitor->indentlevel);
  gt_file_xprintf(visitor->outfp, "<total_number_ESTs_reported>%lu"
                  "</total_number_ESTs_reported>\n", num_of_sas);
  gth_indent(visitor->outfp, visitor->indentlevel);
  gt_file_xprintf(visitor->outfp, "</alignment_module>\n");
  visitor->indentlevel--;
}

const GthSAVisitorClass* gth_xml_final_sa_visitor_class()
{
  static const GthSAVisitorClass savc = { sizeof (GthXMLFinalSAVisitor),
                                          NULL,
                                          xml_final_sa_visitor_preface,
                                          xml_final_sa_visitor_visit_sa,
                                          xml_final_sa_visitor_trailer
                                        };
  return &savc;
}

GthSAVisitor* gth_xml_final_sa_visitor_new(GthInput *input,
                                           unsigned long minintronlength,
                                           unsigned long translationtable,
                                           unsigned int indentlevel,
                                           GtFile *outfp)
{
  GthSAVisitor *sa_visitor =
    gth_sa_visitor_create(gth_xml_final_sa_visitor_class());
  GthXMLFinalSAVisitor *visitor = xml_final_sa_visitor_cast(sa_visitor);
  visitor->input = input;
  visitor->minintronlength = minintronlength;
  visitor->translationtable = translationtable;
  visitor->indentlevel = indentlevel;
  visitor->outfp = outfp;
  return sa_visitor;
}
