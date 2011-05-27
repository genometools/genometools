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
#include "gth/sa_visitor_rep.h"
#include "gth/txt_sa_visitor.h"

#define SA_DELIMITERLINECHAR  '*'

struct GthTxtSAVisitor {
  const GthSAVisitor parent_instance;
  GthInput *input;
  bool gs2out,
       showseqnums;
  unsigned long minintronlength,
                widthforgenpos,
                showintronmaxlen,
                translationtable;
  GtFile *outfp;
};

#define txt_sa_visitor_cast(GV)\
        gth_sa_visitor_cast(gth_txt_sa_visitor_class(), GV)

static void showgthreferenceinformation(GthSA *sa, GthInput *input,
                                        bool showseqnums,
                                        GtFile *outfp)
{
  gt_assert(gth_sa_ref_file_num(sa) != GT_UNDEF_ULONG);

  switch (gth_sa_alphatype(sa)) {
    case DNA_ALPHA:
      gt_file_xprintf(outfp,
                         "EST Sequence: file=%s, strand=%c, description=",
                         gth_input_get_reference_filename(input,
                                                  gth_sa_ref_file_num(sa)),
                         gth_sa_ref_strand_char(sa));
      break;
    case PROTEIN_ALPHA:
      gt_file_xprintf(outfp, "Protein Sequence: file=%s, description=",
                         gth_input_get_reference_filename(input,
                                                 gth_sa_ref_file_num(sa)));
      break;
    default: gt_assert(0);
  }

  gth_sa_echo_reference_description(sa, input, outfp);

  if (showseqnums)
    gt_file_xprintf(outfp, ", seqnum=%lu",  gth_sa_ref_seq_num(sa));

  gt_file_xfputc('\n', outfp);
  gt_file_xfputc('\n', outfp);
}

static void showgs2referenceinformation(GthSA *sa, GtFile *outfp)
{
  gt_file_xprintf(outfp, "EST sequence %4lu %cstrand (File: %s%c)\n\n",
                     gth_sa_call_number(sa),
                     gth_sa_ref_strand_char(sa),
                     gth_sa_ref_id(sa),
                     gth_sa_ref_strand_char(sa));
}

static void showgthgenomicinformation(GthSA *sa, GthInput *input,
                                      bool showseqnums, GtFile *outfp)
{
  gt_assert(gth_sa_gen_file_num(sa) != GT_UNDEF_ULONG);

  gt_file_xprintf(outfp, "Genomic Template: file=%s, strand=%c, from=%lu, "
                            "to=%lu, description=",
                     gth_input_get_genomic_filename(input,
                                                    gth_sa_gen_file_num(sa)),
                     gth_sa_gen_strand_char(sa),
                     gth_sa_gen_dp_start_show(sa),
                     gth_sa_gen_dp_end_show(sa));

  gth_sa_echo_genomic_description(sa, input, outfp);

  if (showseqnums)
    gt_file_xprintf(outfp, ", seqnum=%lu",  gth_sa_gen_seq_num(sa));

  gt_file_xfputc('\n', outfp);
  gt_file_xfputc('\n', outfp);
}

/*
  The following function prints a PPA line, which shows the start and end
  position of the poly-A tail in the cDNA (iff a poly-A tail could be
  determined).
*/
static void showppaline(GthSA *sa, GtFile *outfp)
{
  if (gth_sa_polyAtail_start(sa) ||
      gth_sa_polyAtail_stop(sa)) {
    gt_file_xprintf(outfp,
                       " PPA                                  cDNA %6lu %6lu\n",
                       gth_sa_polyAtail_start(sa) + OUTPUTOFFSET,
                       gth_sa_polyAtail_stop(sa) + OUTPUTOFFSET);
  }
}

/* The following function prints the "classic" GeneSeqer2 MATCH line */
static void showmatchline(GthSA *sa, GtFile *outfp)
{
  gt_file_xprintf(outfp, "MATCH\t%s%c\t%s%c\t%5.3f\t%lu\t%5.3f\t%c\n",
                     gth_sa_gen_id(sa),
                     gth_sa_gen_strand_char(sa),
                     gth_sa_ref_id(sa),
                     gth_sa_ref_strand_char(sa),
                     gth_sa_score(sa),
                     gth_sa_cumlen_scored_exons(sa),
                     gth_sa_coverage(sa),
                     gth_sa_coverage_char(sa));
}

/*
  The following function prints the "classic" GeneSeqer2 PGS line
*/
static void showpgsline(GthSA *sa, GtFile *outfp)
{
  unsigned long i, numofexons;
  gt_assert(sa);
  numofexons = gth_sa_num_of_exons(sa);
  gt_file_xprintf(outfp, "PGS_%s%c_%s%c\t(",
                     gth_sa_gen_id(sa),
                     gth_sa_gen_strand_char(sa),
                     gth_sa_ref_id(sa),
                     gth_sa_ref_strand_char(sa));

  for (i = 0; i < numofexons; i++) {
    gt_file_xprintf(outfp, "%lu  %lu",
                    gth_sa_left_genomic_exon_border(sa, i),
                    gth_sa_right_genomic_exon_border(sa, i));
    if (i == numofexons - 1)
      gt_file_xprintf(outfp, ")\n\n");
    else
      gt_file_xfputc(',', outfp);
  }
}

static void showalignmentheader(GthSA *sa, bool gs2out, int widthforgenpos,
                                unsigned long minintronlength, GtFile *outfp)
{
  unsigned long i, leftreferenceexonborder, rightreferenceexonborder,
                referenceexonlength;
  GthDbl exonscore, donorsitescore, acceptorsitescore;
  GthFlt donorsiteprobability, acceptorsiteprobability;
  Exoninfo *exoninfo;
  Introninfo *introninfo;

  gt_file_xprintf(outfp, "Predicted gene structure");
  if (gs2out) {
    gt_file_xprintf(outfp, " (within gDNA segment %lu to %lu):\n",
                       gth_sa_gen_dp_start_show(sa),
                       gth_sa_gen_dp_end_show(sa));
  }
  else
    gt_file_xprintf(outfp, ":\n");
  gt_file_xfputc('\n', outfp);

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

      gt_file_xprintf(outfp, "  Intron %2lu %*lu %*lu (%4lu n); ",
                      i - 1 + OUTPUTOFFSET, widthforgenpos,
                      gth_sa_left_intron_border(sa, i-1),
                      widthforgenpos,
                      gth_sa_right_intron_border(sa, i-1),
                      gth_sa_intron_length(sa, i-1));

      gt_file_xprintf(outfp, "Pd: %5.3f ", donorsiteprobability);
      if (gth_sa_alphatype(sa) == DNA_ALPHA) {
        if (donorsitescore == 0.0)
          gt_file_xprintf(outfp, "(s:    0), ");
        else
          gt_file_xprintf(outfp, "(s: %4.2f), ", donorsitescore);
      }
      else
        gt_file_xprintf(outfp, "  ");
      gt_file_xprintf(outfp, "Pa: %5.3f ", acceptorsiteprobability);
      if (gth_sa_alphatype(sa) == DNA_ALPHA) {
        if (acceptorsitescore == 0.0)
          gt_file_xprintf(outfp, "(s:    0)");
        else
          gt_file_xprintf(outfp, "(s: %4.2f)", acceptorsitescore);
      }
      /* if the intron is shorter or equal than the minimum intron length two
         question marks are shown at the end of the line */
      if (gth_sa_intron_length(sa, i-1) <= minintronlength)
        gt_file_xprintf(outfp, " ??");
      gt_file_xfputc('\n', outfp);
    }

    gt_file_xprintf(outfp,
                    " Exon %2lu %*lu %*lu (%4lu n);  %s %6lu %6lu (%4lu %s); "
                    "score: %5.3f\n", i + OUTPUTOFFSET, widthforgenpos,
                    gth_sa_left_genomic_exon_border(sa, i),
                    widthforgenpos,
                    gth_sa_right_genomic_exon_border(sa, i),
                    gth_sa_genomic_exon_length(sa, i),
                    gth_sa_alphastring(sa),
                    leftreferenceexonborder  + OUTPUTOFFSET,
                    rightreferenceexonborder + OUTPUTOFFSET,
                    referenceexonlength,
                    gth_sa_alphatype(sa) == DNA_ALPHA
                    ? "n" : "aa",
                    exonscore);
  }

  /* showing PPA line (if an poly-A tail was determined) */
  if (gth_sa_alphatype(sa) == DNA_ALPHA)
    showppaline(sa, outfp);
  gt_file_xfputc('\n', outfp);

  /* showing MATCH line */
  showmatchline(sa, outfp);

  /* showing PGS line */
  showpgsline(sa, outfp);
}

static void showdelimiterline(GtFile *outfp)
{
  unsigned long i;
  for (i = 0; i < DELIMITERLINELENGTH; i++)
    gt_file_xfputc(SA_DELIMITERLINECHAR, outfp);
  gt_file_xfputc('\n', outfp);
}

static void show_spliced_alignment(GthSA *sa, GthInput *input, bool gs2out,
                                   unsigned long minintronlength,
                                   unsigned long widthforgenpos,
                                   unsigned long showintronmaxlen,
                                   unsigned long translationtable,
                                   bool showseqnums, GtFile *outfp)
{
  bool wildcardimplosion = false;

  showdelimiterline(outfp);

  if (gs2out) {
    /* all wildcards (N,S,Y,W,R,K,V,B,D,H,M) will be replaced by the wildcard N
       makes only sense for a DNA alphabet */
    wildcardimplosion = true;
    showgs2referenceinformation(sa, outfp);
    gth_sa_echo_reference_sequence(sa, input, true, outfp);
  }
  else {
    showgthreferenceinformation(sa, input, showseqnums, outfp);
    gth_sa_echo_reference_sequence(sa, input, true, outfp);
    showgthgenomicinformation(sa, input, showseqnums, outfp);
  }

  showalignmentheader(sa, gs2out, widthforgenpos, minintronlength, outfp);

  gt_file_xprintf(outfp,
                     "Alignment (genomic DNA sequence = upper lines):\n\n");

  gth_sa_echo_alignment(sa, showintronmaxlen, translationtable,
                        wildcardimplosion, input, outfp);
}

static void txt_sa_visitor_visit_sa(GthSAVisitor *sa_visitor, GthSA *sa)
{
  GthTxtSAVisitor *visitor = txt_sa_visitor_cast(sa_visitor);
  gt_assert(sa);
  show_spliced_alignment(sa, visitor->input, visitor->gs2out,
                         visitor->minintronlength, visitor->widthforgenpos,
                         visitor->showintronmaxlen, visitor->translationtable,
                         visitor->showseqnums, visitor->outfp);
}

const GthSAVisitorClass* gth_txt_sa_visitor_class()
{
  static const GthSAVisitorClass savc = { sizeof (GthTxtSAVisitor),
                                          NULL,
                                          NULL,
                                          txt_sa_visitor_visit_sa,
                                          NULL };
  return &savc;
}

GthSAVisitor* gth_txt_sa_visitor_new(GthInput *input, bool gs2out,
                                     unsigned long minintronlength,
                                     unsigned long widthforgenpos,
                                     unsigned long showintronmaxlen,
                                     unsigned long translationtable,
                                     bool showseqnums, GtFile *outfp)
{
  GthSAVisitor *sa_visitor = gth_sa_visitor_create(gth_txt_sa_visitor_class());
  GthTxtSAVisitor *visitor = txt_sa_visitor_cast(sa_visitor);
  visitor->input = input;
  visitor->gs2out = gs2out;
  visitor->minintronlength = minintronlength;
  visitor->widthforgenpos = widthforgenpos;
  visitor->showintronmaxlen = showintronmaxlen;
  visitor->translationtable = translationtable;
  visitor->showseqnums = showseqnums;
  visitor->outfp = outfp;
  return sa_visitor;
}
