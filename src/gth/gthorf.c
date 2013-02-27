/*
  Copyright (c) 2004-2010 Gordon Gremme <gremme@zbh.uni-hamburg.de>
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

#include "core/codon_api.h"
#include "core/orf.h"
#include "core/trans_table.h"
#include "core/unused_api.h"
#include "gth/gthorf.h"
#include "gth/gthstopcodon.h"
#include "gth/indent.h"

#define CONSOLIDATE_ORFS_MAX_OVERLAP            30
#define CONSOLIDATE_ORFS_WITHIN_LENGTH          30

#define ORFBLOCKWIDTH           10
#define ORFLINEWIDTH            60
#define ORFBPWIDTH              5
#define ORFRESIDUESWIDTH        4

typedef struct {
  GtRange splseqrange;   /* genomic positions refering to spliced seq.
                            (without stopcodon) */
  unsigned long framenum,
                lengthwithstopcodon,
                lengthwithoutstopcodon;
  const char *frame;
  bool stopcodon;
} MaximalORF;

typedef struct {
  GtArray *maximalORFs;
  unsigned long minORFlength;
} SaveORFInfo;

static void saveORF(void *data, GtRange *range, unsigned long framenum,
                    const char *frame, bool ends_with_stop_codon)
{
  MaximalORF orf;
  SaveORFInfo *info = data;
  gt_assert(info && range);
  orf.framenum = framenum;
  orf.lengthwithstopcodon = gt_range_length(range);
  orf.frame = frame + range->start;
  if (ends_with_stop_codon)
    orf.lengthwithoutstopcodon = orf.lengthwithstopcodon - 1;
  else
    orf.lengthwithoutstopcodon = orf.lengthwithstopcodon;
  orf.stopcodon = ends_with_stop_codon;
  orf.splseqrange.start = range->start * GT_CODON_LENGTH + orf.framenum;
  orf.splseqrange.end = orf.splseqrange.start
                        + orf.lengthwithoutstopcodon * GT_CODON_LENGTH - 1;
                        /* -1 to shift from the genomic position _after_ the
                           last codon to the _last_ position of the last codon
                        */
  /* save ORF, if necessary */
  if (orf.lengthwithoutstopcodon >= info->minORFlength)
    gt_array_add(info->maximalORFs, orf);
}

static int compareORFs(const void *dataA, const void *dataB)
{
  MaximalORF *orfA = (MaximalORF*) dataA,
             *orfB = (MaximalORF*) dataB;
  gt_assert(orfA && orfB);
  if (orfA->lengthwithoutstopcodon > orfB->lengthwithoutstopcodon)
    return -1;
  else if (orfA->lengthwithoutstopcodon == orfB->lengthwithoutstopcodon) {
    if (orfA->framenum < orfB->framenum)
      return -1;
    else if (orfA->framenum == orfB->framenum)
      return 0;
    else
      return 1;
  }
  else
    return 1;
}

#ifndef NDEBUG
static bool orfs_are_sorted(GtArray *maximalORFs)
{
  MaximalORF *orfA, *orfB;
  unsigned long i;
  gt_assert(maximalORFs);
  for (i = 1; i < gt_array_size(maximalORFs); i++) {
    orfA = (MaximalORF*) gt_array_get(maximalORFs, i-1);
    orfB = (MaximalORF*) gt_array_get(maximalORFs, i);
    if (orfA->lengthwithoutstopcodon < orfB->lengthwithoutstopcodon)
      return false;
    else if ((orfA->lengthwithoutstopcodon == orfB->lengthwithoutstopcodon) &&
             (orfA->framenum > orfB->framenum)) {
      return false;
    }
  }
  return true;
}
#endif

static void sortmaximalORFs(GtArray *maximalORFs)
{
  gt_assert(maximalORFs);
  qsort(gt_array_get_space(maximalORFs), gt_array_size(maximalORFs),
        sizeof (MaximalORF), compareORFs);
  gt_assert(orfs_are_sorted(maximalORFs));
}

static void showsingleORF(MaximalORF *orf, bool gen_strand_forward,
                          unsigned long gen_total_length,
                          unsigned long gen_offset, const char *gen_id,
                          unsigned long pglnum, unsigned long agsnum,
                          unsigned long ppsnum, GthSplicedSeq *splicedseq,
                          unsigned int indentlevel, GthOutput *out)
{
  unsigned long i;
  GT_UNUSED bool hasborder = false;
  GtFile *outfp = out->outfp;

  if (out->xmlout) {
    gth_indent(outfp, indentlevel);
    gt_file_xprintf(outfp, "<orf_entry>\n");
    indentlevel++;
    gth_indent(outfp, indentlevel);
    gt_file_xprintf(outfp, "<id_line>\n");
    indentlevel++;
    gth_indent(outfp, indentlevel);
    gt_file_xprintf(outfp, "<gDNA id=\"%s\" strand=\"%c\"/>\n", gen_id,
                    SHOWSTRAND(gen_strand_forward));
    gth_indent(outfp, indentlevel);
    gt_file_xprintf(outfp, "<serials PGL_serial=\"%lu\" AGS_serial=\"%lu\" "
                    "PPS_serial=\"%lu\"/>\n", pglnum + OUTPUTOFFSET,
                    agsnum + OUTPUTOFFSET,  ppsnum + OUTPUTOFFSET);
    gth_indent(outfp, indentlevel);
    gt_file_xprintf(outfp, "<orf_info>\n");
    indentlevel++;
    gth_indent(outfp, indentlevel);
    gt_file_xprintf(outfp, "<exon_boundaries>\n");
    indentlevel++;
    gth_indent(outfp, indentlevel);
    gt_file_xprintf(outfp, "<exon start=\"%lu\" ",
                    SHOWGENPOS(gen_strand_forward, gen_total_length,
                               gen_offset,
                               splicedseq
                               ->positionmapping[orf->splseqrange.start]));
  }
  else {
    gt_file_xprintf(outfp, ">%s%c_PGL-%lu_AGS-%lu_PPS_%lu (%lu  ", gen_id,
                    SHOWSTRAND(gen_strand_forward), pglnum + OUTPUTOFFSET,
                    agsnum + OUTPUTOFFSET,  ppsnum + OUTPUTOFFSET,
                    SHOWGENPOS(gen_strand_forward, gen_total_length,
                               gen_offset,
                               splicedseq
                               ->positionmapping[orf->splseqrange.start]));
  }

  for (i = orf->splseqrange.start;
       i < orf->splseqrange.end + orf->stopcodon * GT_CODON_LENGTH; i++) {
    if (gth_spliced_seq_pos_is_border(splicedseq, i)) {
      if (out->xmlout) {
        gt_file_xprintf(outfp, "stop=\"%lu\"/>\n",
                        SHOWGENPOS(gen_strand_forward, gen_total_length,
                                   gen_offset,
                                   splicedseq->positionmapping[i]));
        gth_indent(outfp, indentlevel);
        gt_file_xprintf(outfp, "<exon start=\"%lu\" ",
                        SHOWGENPOS(gen_strand_forward, gen_total_length,
                                   gen_offset,
                                   splicedseq->positionmapping[i+1]));
      }
      else {
        gt_file_xprintf(outfp, "%lu,%lu  ",
                        SHOWGENPOS(gen_strand_forward, gen_total_length,
                                   gen_offset,
                                   splicedseq->positionmapping[i]),
                        SHOWGENPOS(gen_strand_forward, gen_total_length,
                                   gen_offset,
                                   splicedseq->positionmapping[i+1]));
      }

      hasborder = true;
    }
  }

  if (out->xmlout) {
    gt_file_xprintf(outfp, "stop=\"%lu\"/>\n",
                    SHOWGENPOS(gen_strand_forward, gen_total_length,
                               gen_offset, splicedseq
                               ->positionmapping[orf->splseqrange.end +
                                                 orf->stopcodon *
                                                 GT_CODON_LENGTH]));
    indentlevel--;
    gth_indent(outfp, indentlevel);
    gt_file_xprintf(outfp, "</exon_boundaries>\n");
  }
  else {
    gt_file_xprintf(outfp, "%lu)",
                    SHOWGENPOS(gen_strand_forward,
                               gen_total_length, gen_offset, splicedseq
                               ->positionmapping[orf->splseqrange.end +
                                                 orf->stopcodon *
                                                 GT_CODON_LENGTH]));
  }

  if (out->xmlout)
  {

    gth_indent(outfp, indentlevel);
    gt_file_xprintf(outfp, "<frame>%lu</frame>\n",  orf->framenum);
    gth_indent(outfp, indentlevel);
    gt_file_xprintf(outfp, "<number_coding_nucleotides>%lu"
                    "</number_coding_nucleotides>\n",
                    orf->lengthwithoutstopcodon * GT_CODON_LENGTH);
    gth_indent(outfp, indentlevel);
    gt_file_xprintf(outfp, "<number_encoded_amino_acids>%lu"
                    "</number_encoded_amino_acids>\n",
                    orf->lengthwithoutstopcodon);
    indentlevel--;
    gth_indent(outfp, indentlevel);
    gt_file_xprintf(outfp, "</orf_info>\n");
  }
  else {
    gt_file_xprintf(outfp, "\t(frame '%lu'; %*lu bp, %*lu residues)",
                    orf->framenum, ORFBPWIDTH,
                    orf->lengthwithoutstopcodon * GT_CODON_LENGTH,
                    ORFRESIDUESWIDTH,
                    orf->lengthwithoutstopcodon);
  }

  if (out->xmlout) {
    indentlevel--;
    gth_indent(outfp, indentlevel);
    gt_file_xprintf(outfp, "</id_line>\n");
    gth_indent(outfp, indentlevel);
    gt_file_xprintf(outfp, "<predicted_protein_sequence>");
  }

  for (i = 0; i < orf->lengthwithstopcodon; i++) {
    if (!out->xmlout) {
      if (i % ORFBLOCKWIDTH == 0)
        gt_file_xfputc(' ', outfp);
      if (i % ORFLINEWIDTH == 0) {
        gt_file_xprintf(outfp, "\n%*lu  ", out->widthforgenpos,
                        i + OUTPUTOFFSET);
      }
    }

    if (out->gs2out && i + 1 == orf->lengthwithstopcodon &&
        orf->frame[i] == GT_STOP_AMINO) {
      gt_file_xfputc(GS2_STOPCODON, outfp);
    }
    else
      gt_file_xprintf(outfp, "%c", orf->frame[i]);
  }

  if (out->xmlout) {
    gt_file_xprintf(outfp, "</predicted_protein_sequence>\n");
    indentlevel--;
    gth_indent(outfp, indentlevel);
    gt_file_xprintf(outfp, "</orf_entry>\n");
  }
  else
    gt_file_xprintf(outfp, "\n\n");
}

static void outputconsolidatedORFs(GtArray *consolidatedORFs,
                                   bool gen_strand_forward,
                                   unsigned long gen_total_length,
                                   unsigned long gen_offset,
                                   const char *gen_id,
                                   unsigned long pglnum,
                                   unsigned long agsnum,
                                   GthSplicedSeq *splicedseq,
                                   unsigned int indentlevel,
                                   GthOutput *out)
{
  GtFile *outfp = out->outfp;
  unsigned long i;

  if (out->xmlout) {
    gth_indent(outfp, indentlevel);
    gt_file_xprintf(outfp,
               "<probable_ORFs xmlns=\"http://www.genomethreader.org/"
               "GTH_output/PGL_module/predicted_gene_location/"
               "AGS_information/three_phase_translation/probable_ORFs/"
               "\">\n");
    indentlevel++;
  }
  else {
    gt_file_xprintf(outfp,
              "Maximal non-overlapping open reading frames (>= %lu codons)\n",
                out->minORFlength);
  }

  if (gt_array_size(consolidatedORFs) == 0) {
    if (out->xmlout) {
      gth_indent(outfp, indentlevel);
      gt_file_xprintf(outfp, "<none gDNA_id=\"%s\"/>\n", gen_id);
    }
    else if (!out->xmlout)
      gt_file_xprintf(outfp, "none\n");
  }
  else {
    if (!out->xmlout)
      gt_file_xfputc('\n', outfp);

    for (i = 0; i < gt_array_size(consolidatedORFs); i++) {
      showsingleORF(gt_array_get(consolidatedORFs, i), gen_strand_forward,
                    gen_total_length, gen_offset, gen_id, pglnum,
                    agsnum, i, splicedseq, indentlevel, out);
    }
  }

  if (out->xmlout) {
    indentlevel--;
    gth_indent(outfp, indentlevel);
    gt_file_xprintf(outfp, "</probable_ORFs>\n");
  }
}

static bool overlapwithotherranges(GtRange range, GtArray *ranges,
                                   unsigned long delta)
{
  unsigned long i;
  gt_assert(ranges);
  for (i = 0; i < gt_array_size(ranges); i++) {
    if (gt_range_overlap_delta(&range, gt_array_get(ranges, i), delta))
      return true;
  }
  return false;
}

static void consolidateORFs(GtArray *consolidatedORFs, GtArray *maximalORFs)
{
  unsigned long i, maxlength;
  MaximalORF *orf;
  GtArray *ranges;

  gt_assert(consolidatedORFs && maximalORFs);

  ranges = gt_array_new(sizeof (GtRange));

  /* save first ORF */
  orf = gt_array_get_first(maximalORFs);
  gt_array_add(consolidatedORFs, *orf);
  /* and its genomic range */
  gt_array_add(ranges, orf->splseqrange);
  /* and the maximal length */
  maxlength = orf->lengthwithoutstopcodon;

  /* save the other ORFs if they fit the criteria */
  for (i = 1; i < gt_array_size(maximalORFs); i++) {
    orf = gt_array_get(maximalORFs, i);
    if (!overlapwithotherranges(orf->splseqrange, ranges,
                                CONSOLIDATE_ORFS_MAX_OVERLAP)
        || maxlength - orf->lengthwithoutstopcodon
           <= CONSOLIDATE_ORFS_WITHIN_LENGTH) {
      /* save actual ORF and its genomic range */
      gt_array_add(consolidatedORFs, *orf);
      gt_array_add(ranges, orf->splseqrange);
    }
  }

  /* the number of consolidated ORFs equals the number of ranges */
  gt_assert(gt_array_size(consolidatedORFs) == gt_array_size(ranges));

  gt_array_delete(ranges);
}

void gthshowORFs(char *frame0, char *frame1, char *frame2,
                 unsigned long frame0len, unsigned long frame1len,
                 unsigned long frame2len, bool gen_strand_forward,
                 unsigned long gen_total_length, unsigned long gen_offset,
                 const char *gen_id, unsigned long pglnum,
                 unsigned long agsnum, GthSplicedSeq *splicedseq,
                 unsigned int indentlevel, GthOutput *out)
{
  GtArray *maximalORFs, *consolidatedORFs;
  SaveORFInfo info;

  maximalORFs = gt_array_new(sizeof (MaximalORF));
  consolidatedORFs = gt_array_new(sizeof (MaximalORF));

  info.maximalORFs = maximalORFs;
  info.minORFlength = out->minORFlength;

  gt_determine_ORFs(saveORF, &info, 0, frame0, frame0len, out->start_codon,
                    out->final_stop_codon, true, NULL);
  gt_determine_ORFs(saveORF, &info, 1, frame1, frame1len, out->start_codon,
                    out->final_stop_codon, true, NULL);
  gt_determine_ORFs(saveORF, &info, 2, frame2, frame2len, out->start_codon,
                    out->final_stop_codon, true, NULL);

  if (gt_array_size(maximalORFs)) {
    sortmaximalORFs(maximalORFs);
    consolidateORFs(consolidatedORFs, maximalORFs);
  }

  outputconsolidatedORFs(consolidatedORFs, gen_strand_forward,
                         gen_total_length, gen_offset, gen_id, pglnum,
                         agsnum, splicedseq, indentlevel, out);

  gt_array_delete(maximalORFs);
  gt_array_delete(consolidatedORFs);
}
