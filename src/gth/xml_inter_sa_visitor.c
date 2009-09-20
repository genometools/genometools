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

#include "gth/indent.h"
#include "gth/showbool.h"
#include "gth/intermediate.h"
#include "gth/sa_visitor_rep.h"
#include "gth/xml_inter_sa_visitor.h"

#define PRECISION  16

struct GthXMLInterSAVisitor {
  const GthSAVisitor parent_instance;
  GthInput *input;
  unsigned int indentlevel;
  GtFile *outfp;
};

#define xml_inter_sa_visitor_cast(GV)\
        gth_sa_visitor_cast(gth_xml_inter_sa_visitor_class(), GV)

static void showgenomicfilename(GthSA *sa, GthInput *input,
                                unsigned int indentlevel, GtFile *outfp)
{
  gth_indent(outfp, indentlevel);
  gt_file_xprintf(outfp, "<genomicfile>\n");
  indentlevel++;

  gth_indent(outfp, indentlevel);
  gt_file_xprintf(outfp, "<genomicfilename>%s</genomicfilename>\n",
                     gth_input_get_genomic_filename(input,
                                                   gth_sa_gen_file_num(sa)));
  gth_indent(outfp, indentlevel);
  gt_file_xprintf(outfp, "<genomicfilehash>%s</genomicfilehash>\n",
                     GTH_UNDEFINED_HASH);

  indentlevel--;
  gth_indent(outfp, indentlevel);
  gt_file_xprintf(outfp, "</genomicfile>\n");
}

static void showreferencefilename(GthSA *sa,
                                  GthInput *input,
                                  unsigned int indentlevel, GtFile *outfp)
{
  gth_indent(outfp, indentlevel);
  gt_file_xprintf(outfp, "<referencefile>\n");
  indentlevel++;

  gth_indent(outfp, indentlevel);
  gt_file_xprintf(outfp, "<referencefilename>%s</referencefilename>\n",
                     gth_input_get_reference_filename(input,
                                                 gth_sa_ref_file_num(sa)));
  gth_indent(outfp, indentlevel);
  gt_file_xprintf(outfp, "<referencefilehash>%s</referencefilehash>\n",
                     GTH_UNDEFINED_HASH);

  indentlevel--;
  gth_indent(outfp, indentlevel);
  gt_file_xprintf(outfp, "</referencefile>\n");
}

static void showalignmentcutoffs(GthSA *sa, unsigned int indentlevel,
                                 GtFile *outfp)
{
  gth_indent(outfp, indentlevel);
  gt_file_xprintf(outfp, "<cutoffs>\n");
  indentlevel++;
  gth_indent(outfp, indentlevel);
  gt_file_xprintf(outfp, "<cutoffsstart>\n");
  indentlevel++;
  gth_indent(outfp, indentlevel);
  gt_file_xprintf(outfp, "<genomiccutoff>%lu</genomiccutoff>\n",
                     gth_sa_genomiccutoff_start(sa));
  gth_indent(outfp, indentlevel);
  gt_file_xprintf(outfp, "<referencecutoff>%lu</referencecutoff>\n",
                     gth_sa_referencecutoff_start(sa));
  gth_indent(outfp, indentlevel);
  gt_file_xprintf(outfp, "<eopcutoff>%lu</eopcutoff>\n",
                     gth_sa_eopcutoff_start(sa));
  indentlevel--;
  gth_indent(outfp, indentlevel);
  gt_file_xprintf(outfp, "</cutoffsstart>\n");
  gth_indent(outfp, indentlevel);
  gt_file_xprintf(outfp, "<cutoffsend>\n");
  indentlevel++;
  gth_indent(outfp, indentlevel);
  gt_file_xprintf(outfp, "<genomiccutoff>%lu</genomiccutoff>\n",
                     gth_sa_genomiccutoff_end(sa));
  gth_indent(outfp, indentlevel);
  gt_file_xprintf(outfp, "<referencecutoff>%lu</referencecutoff>\n",
                     gth_sa_referencecutoff_end(sa));
  gth_indent(outfp, indentlevel);
  gt_file_xprintf(outfp, "<eopcutoff>%lu</eopcutoff>\n",
                     gth_sa_eopcutoff_end(sa));
  indentlevel--;
  gth_indent(outfp, indentlevel);
  gt_file_xprintf(outfp, "</cutoffsend>\n");
  indentlevel--;
  gth_indent(outfp, indentlevel);
  gt_file_xprintf(outfp, "</cutoffs>\n");
}

static void showexons(GthSA *sa, unsigned int indentlevel,
                      GtFile *outfp)
{
  Exoninfo *exoninfo;
  unsigned long i;

  for (i = 0; i < gth_sa_num_of_exons(sa); i++) {
    exoninfo = gth_sa_get_exon(sa, i);

    gth_indent(outfp, indentlevel);
    gt_file_xprintf(outfp, "<exoninfo>\n");
    indentlevel++;

    gth_indent(outfp, indentlevel);
    gt_file_xprintf(outfp,
                    "<leftgenomicexonborder>%lu</leftgenomicexonborder>\n",
                    exoninfo->leftgenomicexonborder);
    gth_indent(outfp, indentlevel);
    gt_file_xprintf(outfp,
                    "<rightgenomicexonborder>%lu</rightgenomicexonborder>\n",
                    exoninfo->rightgenomicexonborder);
    gth_indent(outfp, indentlevel);
    gt_file_xprintf(outfp,
                    "<leftreferenceexonborder>%lu</leftreferenceexonborder>\n",
                    exoninfo->leftreferenceexonborder);
    gth_indent(outfp, indentlevel);
    gt_file_xprintf(outfp, "<rightreferenceexonborder>%lu</"
                    "rightreferenceexonborder>\n",
                    exoninfo->rightreferenceexonborder);
    gth_indent(outfp, indentlevel);
    gt_file_xprintf(outfp, "<exonscore>%.*f</exonscore>\n", PRECISION,
                    exoninfo->exonscore);

    indentlevel--;
    gth_indent(outfp, indentlevel);
    gt_file_xprintf(outfp, "</exoninfo>\n");
  }
}

static void showintrons(GthSA *sa, bool dnaalpha,
                        unsigned int indentlevel, GtFile *outfp)
{
  Introninfo *introninfo;
  unsigned long i;

  for (i = 0; i < gth_sa_num_of_introns(sa); i++) {
    introninfo = gth_sa_get_intron(sa, i);

    gth_indent(outfp, indentlevel);
    gt_file_xprintf(outfp, "<introninfo>\n");
    indentlevel++;

    gth_indent(outfp, indentlevel);
    gt_file_xprintf(outfp,
                    "<donorsiteprobability>%.*f</donorsiteprobability>\n",
                    PRECISION, introninfo->donorsiteprobability);

    gth_indent(outfp, indentlevel);
    gt_file_xprintf(outfp,
                    "<acceptorsiteprobability>%.*f</acceptorsiteprobability>\n",
                    PRECISION, introninfo->acceptorsiteprobability);

    gth_indent(outfp, indentlevel);
    gt_file_xprintf(outfp, "<donorsitescore>%.*f</donorsitescore>\n",
                       PRECISION,
                       dnaalpha ? introninfo->donorsitescore
                                : UNDEFINED_SPLICE_SITE_SCORE);

    gth_indent(outfp, indentlevel);
    gt_file_xprintf(outfp, "<acceptorsitescore>%.*f</acceptorsitescore>\n",
                       PRECISION,
                       dnaalpha ? introninfo->acceptorsitescore
                                : UNDEFINED_SPLICE_SITE_SCORE);

    indentlevel--;
    gth_indent(outfp, indentlevel);
    gt_file_xprintf(outfp, "</introninfo>\n");
  }
}

static void showpolyAtailpos(GthSA *sa, unsigned int indentlevel,
                             GtFile *outfp)
{
  gth_indent(outfp, indentlevel);
  gt_file_xprintf(outfp, "<polyAtailpos>\n");
  indentlevel++;

  gth_indent(outfp, indentlevel);
  gt_file_xprintf(outfp, "<polyAstart>%lu</polyAstart>\n",
                     gth_sa_polyAtail_start(sa));
  gth_indent(outfp, indentlevel);
  gt_file_xprintf(outfp, "<polyAstop>%lu</polyAstop>\n",
                     gth_sa_polyAtail_stop(sa));

  indentlevel--;
  gth_indent(outfp, indentlevel);
  gt_file_xprintf(outfp, "</polyAtailpos>\n");
}

static void xml_inter_show_spliced_alignment(GthSA *sa, GthInput *input,
                                             unsigned int indentlevel,
                                             GtFile *outfp)
{
  bool dnaalpha = true;

  gth_indent(outfp, indentlevel);
  gt_file_xprintf(outfp,
                  "<spliced_alignment xmlns=\"http://www.GenomeThreader.org/"
                  "SplicedAlignment/spliced_alignment/\">\n");
  indentlevel++;

  gth_indent(outfp, indentlevel);
  gt_file_xprintf(outfp, "<referencealphatype>");
  switch (gth_sa_alphatype(sa)) {
    case DNA_ALPHA:
      gt_file_xprintf(outfp, "DNA_ALPHA");
      break;
    case PROTEIN_ALPHA:
      gt_file_xprintf(outfp, "PROTEIN_ALPHA");
      dnaalpha = false;
      break;
    default: gt_assert(0);
  }
  gt_file_xprintf(outfp, "</referencealphatype>\n");

  gth_indent(outfp, indentlevel);
  gt_file_xprintf(outfp, "<editoperations>\n");
  indentlevel++;
  gth_backtrace_path_show_complete(gth_sa_backtrace_path(sa), true, indentlevel,
                                   outfp);
  indentlevel--;
  gth_indent(outfp, indentlevel);
  gt_file_xprintf(outfp, "</editoperations>\n");

  gth_indent(outfp, indentlevel);
  gt_file_xprintf(outfp, "<indelcount>%lu</indelcount>\n",
                     gth_sa_indelcount(sa));

  gth_indent(outfp, indentlevel);
  gt_file_xprintf(outfp, "<genomiclengthDP>%lu</genomiclengthDP>\n",
                     gth_sa_gen_dp_length(sa));

  gth_indent(outfp, indentlevel);
  gt_file_xprintf(outfp, "<genomiclengthtotal>%lu</genomiclengthtotal>\n",
                     gth_sa_gen_total_length(sa));

  gth_indent(outfp, indentlevel);
  gt_file_xprintf(outfp, "<genomicoffset>%lu</genomicoffset>\n",
                     gth_sa_gen_offset(sa));

  gth_indent(outfp, indentlevel);
  gt_file_xprintf(outfp, "<referencelength>%lu</referencelength>\n",
                     gth_sa_ref_total_length(sa));

  gth_indent(outfp, indentlevel);
  gt_file_xprintf(outfp, "<dpstartpos>%lu</dpstartpos>\n",
                     gth_sa_gen_dp_start(sa));

  gth_indent(outfp, indentlevel);
  gt_file_xprintf(outfp, "<dpendpos>%lu</dpendpos>\n",
                     gth_sa_gen_dp_end(sa));

  showgenomicfilename(sa, input, indentlevel, outfp);

  gth_indent(outfp, indentlevel);
  gt_file_xprintf(outfp, "<genomicseqnum>%lu</genomicseqnum>\n",
                     gth_sa_gen_seq_num(sa));

  showreferencefilename(sa, input, indentlevel, outfp);

  gth_indent(outfp, indentlevel);
  gt_file_xprintf(outfp, "<referenceseqnum>%lu</referenceseqnum>\n",
                     gth_sa_ref_seq_num(sa));

  gth_indent(outfp, indentlevel);
  gt_file_xprintf(outfp, "<genomicid>%s</genomicid>\n", gth_sa_gen_id(sa));

  gth_indent(outfp, indentlevel);
  gt_file_xprintf(outfp, "<referenceid>%s</referenceid>\n",
                  gth_sa_ref_id(sa));

  gth_indent(outfp, indentlevel);
  gt_file_xprintf(outfp,
                  "<genomicstrandisforward>%s</genomicstrandisforward>\n",
                  GTH_SHOWBOOL(gth_sa_gen_strand_forward(sa)));

  gth_indent(outfp, indentlevel);
  gt_file_xprintf(outfp,
                    "<referencestrandisforward>%s</referencestrandisforward>\n",
                    GTH_SHOWBOOL(gth_sa_ref_strand_forward(sa)));

  showalignmentcutoffs(sa, indentlevel, outfp);

  showexons(sa, indentlevel, outfp);

  showintrons(sa, dnaalpha, indentlevel, outfp);

  showpolyAtailpos(sa, indentlevel, outfp);

  gth_indent(outfp, indentlevel);
  gt_file_xprintf(outfp, "<alignmentscore>%.*f</alignmentscore>\n",
                  PRECISION, gth_sa_score(sa));

  gth_indent(outfp, indentlevel);
  gt_file_xprintf(outfp, "<coverage>%.*f</coverage>\n", PRECISION,
                     gth_sa_coverage(sa));

  gth_indent(outfp, indentlevel);
  gt_file_xprintf(outfp, "<coverageofgenomicsegmentishighest>%s"
                  "</coverageofgenomicsegmentishighest>\n",
                  GTH_SHOWBOOL(gth_sa_genomic_cov_is_highest(sa)));

  gth_indent(outfp, indentlevel);
  gt_file_xprintf(outfp, "<cumulativelengthofscoredexons>%lu"
                     "</cumulativelengthofscoredexons>\n",
                     gth_sa_cumlen_scored_exons(sa));

  indentlevel--;
  gth_indent(outfp, indentlevel);
  gt_file_xprintf(outfp, "</spliced_alignment>\n");
}

static void xml_inter_sa_visitor_visit_sa(GthSAVisitor *sa_visitor,
                                          GthSA *sa)
{
  GthXMLInterSAVisitor *visitor = xml_inter_sa_visitor_cast(sa_visitor);
  gt_assert(sa);
  xml_inter_show_spliced_alignment(sa, visitor->input, visitor->indentlevel,
                                   visitor->outfp);
}

const GthSAVisitorClass* gth_xml_inter_sa_visitor_class()
{
  static const GthSAVisitorClass savc = { sizeof (GthXMLInterSAVisitor),
                                          NULL,
                                          NULL,
                                          xml_inter_sa_visitor_visit_sa,
                                          NULL };
  return &savc;
}

GthSAVisitor* gth_xml_inter_sa_visitor_new(GthInput *input,
                                           unsigned int indentlevel,
                                           GtFile *outfp)
{
  GthSAVisitor *sa_visitor =
    gth_sa_visitor_create(gth_xml_inter_sa_visitor_class());
  GthXMLInterSAVisitor *visitor = xml_inter_sa_visitor_cast(sa_visitor);
  visitor->input = input;
  visitor->indentlevel = indentlevel;
  visitor->outfp = outfp;
  return sa_visitor;
}

void gth_xml_inter_sa_visitor_set_outfp(GthSAVisitor *sa_visitor,
                                        GtFile *outfp)
{
  GthXMLInterSAVisitor *visitor = xml_inter_sa_visitor_cast(sa_visitor);
  visitor->outfp = outfp;
}
