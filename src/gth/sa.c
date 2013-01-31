/*
  Copyright (c) 2003-2012 Gordon Gremme <gremme@zbh.uni-hamburg.de>
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

#include "core/mathsupport.h"
#include "core/md5_seqid.h"
#include "core/safearith.h"
#include "core/undef_api.h"
#include "core/unused_api.h"
#include "extended/gff3_escaping.h"
#include "gth/default.h"
#include "gth/sa.h"
#include "gth/txt_sa_visitor.h"

#define SA_DELIMITERLINECHAR    '*'
#define CALCPOLYATAILWINDOW     50
#define MINIMUMPOLYATAILLENGTH  10

/* Keep in sync with gth_sas_are_equal()! */
struct GthSA {
  GthBacktracePath *backtrace_path;  /* the edit operations */
  unsigned long gen_total_length,    /* total length of the genomic sequence */
                gen_offset,          /* offset of the genomic sequence where
                                        this SA refers to */
                ref_total_length,    /* total length of reference sequence */
                gen_file_num,        /* genomic file number */
                gen_seq_num,         /* genomic sequence number */
                ref_file_num,        /* reference file number */
                ref_seq_num;         /* reference sequence number */
  GtStr *gen_id,                     /* contains id of genomic sequence */
        *ref_id,                     /* contains id of reference sequence */
        *gen_md5,                    /* contains MD5 of genomic sequence, if
                                        necessary */
        *ref_md5,                    /* contains MD5 of reference sequence, if
                                        necessary */
        *gen_desc,                   /* contains description of genomic
                                        sequence, if necessary */
        *ref_desc;                   /* contains description of reference
                                        sequence, if necessary */
  bool gen_strand_forward,           /* equals true, iff the alignment relates
                                        to the forward strand of the genomic
                                        sequence */
       ref_strand_forward;           /* equals true, iff the alignment relates
                                        to the forward strand of the reference
                                        sequence */
  GtArray *exons,                    /* contains the information about all
                                        exons */
          *introns;                  /* contains the information about all
                                        introns */
  GtRange polyAtailpos;              /* (GS2 = ppa), contains the positions of
                                        the  poly-A tail in the reference
                                        sequence. Iff both positions are equal
                                        to zero, no poly-A tail could be
                                        determined. */
  GthFlt alignmentscore,             /* the score of the alignment
                                        (GS2 = simscore) */
         coverage;                   /* coverage of the genomic DNA or cDNA
                                        segment, whichever is highest
                                        (GS2 = cvrge) */
  bool genomic_cov_is_highest;       /* (GS2 = no equivalent) */
  unsigned long cumlen_scored_exons, /* cumulative length of scored exons
                                        (GS2 = mlgth) */
                call_number;         /* the consecutive number of successful DP
                                        calls. not really usefull, only for
                                        compatibility with GS2 */
  GtStr *gff3_target_attribute;
};

GthSA* gth_sa_new(void)
{
  GthSA *sa;

  /* allocating space */
  sa = gt_calloc(1, sizeof *sa);

  /* initialize strings */
  sa->gen_id = gt_str_new();
  sa->ref_id = gt_str_new();

  /* initialize all arrays */
  sa->backtrace_path = gth_backtrace_path_new(GT_UNDEF_ULONG,
                                              GT_UNDEF_ULONG,
                                              0, /* ref_dp_start */
                                              GT_UNDEF_ULONG);
  sa->exons = gt_array_new(sizeof (Exoninfo));
  sa->introns = gt_array_new(sizeof (Introninfo));

  return sa;
}

GthSA* gth_sa_new_and_set(bool gen_strand_forward,
                          bool ref_strand_forward,
                          GthInput *input,
                          unsigned long gen_file_num,
                          unsigned long gen_seq_num,
                          unsigned long ref_file_num,
                          unsigned long ref_seq_num,
                          unsigned long call_number,
                          unsigned long gen_total_length,
                          unsigned long gen_offset,
                          unsigned long ref_total_length)
{
  GthSA *sa;

  /* alloc and init of arrays */
  sa = gth_sa_new();

  /* setting the strand directions */
  sa->gen_strand_forward = gen_strand_forward;
  sa->ref_strand_forward = ref_strand_forward;

  /* saving sequence ids */
  gth_input_save_gen_id(input, sa->gen_id, gen_file_num, gen_seq_num);
  gth_input_save_ref_id(input, sa->ref_id, ref_file_num, ref_seq_num);

  /* saving MD5s, if necessary */
  gth_input_save_gen_md5(input, &sa->gen_md5, gen_file_num, gen_seq_num);
  gth_input_save_ref_md5(input, &sa->ref_md5, ref_file_num, ref_seq_num);

  /* saving descriptions, if necessary */
  gth_input_save_gen_desc(input, &sa->gen_desc, gen_file_num, gen_seq_num);
  gth_input_save_ref_desc(input, &sa->ref_desc, ref_file_num, ref_seq_num);

  /* save the consecutive call number */
  sa->call_number = call_number;

  /* save total length of genomic sequence */
  sa->gen_total_length = gen_total_length;

  /* save genomic offset */
  sa->gen_offset = gen_offset;

  /* save total length of reference sequence */
  sa->ref_total_length = ref_total_length;
  gth_backtrace_path_set_ref_dp_length(sa->backtrace_path, ref_total_length);

  /* save file and sequence numbers */
  sa->gen_file_num = gen_file_num;
  sa->gen_seq_num = gen_seq_num;
  sa->ref_file_num = ref_file_num;
  sa->ref_seq_num = ref_seq_num;

  return sa;
}

void gth_sa_set(GthSA *sa, GthAlphatype ref_alphatype,
                unsigned long gen_dp_start, unsigned long gen_dp_length)
{
  gth_backtrace_path_set_gen_dp_start(sa->backtrace_path, gen_dp_start);
  gth_backtrace_path_set_gen_dp_length(sa->backtrace_path, gen_dp_length);
  gth_sa_set_score(sa, 0.0);
  gth_sa_set_coverage(sa, 0.0);
  gth_sa_set_highest_cov(sa, true);
  gth_sa_set_cumlen_scored_exons(sa, 0);
  /* reset edit operations */
  gth_backtrace_path_reset(sa->backtrace_path);
  gth_backtrace_path_set_alphatype(sa->backtrace_path, ref_alphatype);
  /* reset arrays */
  gt_array_reset(sa->exons);
  gt_array_reset(sa->introns);
}

void gth_sa_set_gen_dp_length(GthSA *sa, unsigned long gen_dp_length)
{
  gt_assert(sa);
  gth_backtrace_path_set_gen_dp_length(sa->backtrace_path, gen_dp_length);
}

void gth_sa_delete(GthSA *sa)
{
  if (!sa) return;
  gth_backtrace_path_delete(sa->backtrace_path);
  gt_array_delete(sa->exons);
  gt_array_delete(sa->introns);
  gt_str_delete(sa->gen_id);
  gt_str_delete(sa->ref_id);
  gt_str_delete(sa->gen_md5);
  gt_str_delete(sa->ref_md5);
  gt_str_delete(sa->gen_desc);
  gt_str_delete(sa->ref_desc);
  gt_str_delete(sa->gff3_target_attribute);
  gt_free(sa);
}

void gth_sa_show_exons(const GthSA *sa, GtFile *outfp)
{
  Exoninfo *exoninfo;
  unsigned long i;
  gt_assert(sa);
  for (i = 0; i < gt_array_size(sa->exons); i++) {
    exoninfo = (Exoninfo*) gt_array_get(sa->exons, i);
    gt_file_xprintf(outfp, "(%lu,%lu)", exoninfo->leftgenomicexonborder,
                    exoninfo->rightgenomicexonborder);
  }
  gt_file_xfputc('\n', outfp);
}

void gth_sa_get_exons(const GthSA *sa, GtArray *ranges)
{
  Exoninfo *exoninfo;
  unsigned long i;
  GtRange range;
  gt_assert(sa && ranges);
  for (i = 0; i < gt_array_size(sa->exons); i++) {
    exoninfo = gt_array_get(sa->exons, i);
    range.start = exoninfo->leftgenomicexonborder;
    range.end = exoninfo->rightgenomicexonborder;
    gt_array_add(ranges, range);
  }
}

#ifndef NDEBUG
bool gth_sa_exons_are_forward_and_consecutive(const GthSA *sa)
{
  GtArray *ranges;
  gt_assert(sa);
  ranges = gt_array_new(sizeof (GtRange));
  gth_sa_get_exons(sa, ranges);
  if (!gt_ranges_are_consecutive(ranges)) {
    gt_array_delete(ranges);
    return false;
  }
  gt_array_delete(ranges);
  return true;
}
#endif

GtRange gth_sa_range_forward(const GthSA *sa)
{
  GtRange range;
  unsigned long leftgenomicborder, rightgenomicborder;

  gt_assert(sa);

  leftgenomicborder  = ((Exoninfo*) gt_array_get_first(sa->exons))
                       ->leftgenomicexonborder;
  rightgenomicborder = ((Exoninfo*) gt_array_get_last(sa->exons))
                       ->rightgenomicexonborder;

  if (sa->gen_strand_forward) {
    range.start = leftgenomicborder;
    range.end = rightgenomicborder;
  }
  else {
    /* genomic offset is defined */
    gt_assert(sa->gen_offset != GT_UNDEF_ULONG);
    range.start  = sa->gen_total_length - 1
                   - (rightgenomicborder - sa->gen_offset)
                   + sa->gen_offset;
    range.end = sa->gen_total_length - 1
                - (leftgenomicborder - sa->gen_offset)
                + sa->gen_offset;
  }

  return range;
}

GtRange gth_sa_range_actual(const GthSA *sa)
{
  GtRange range;

  gt_assert(sa);

  range.start = ((Exoninfo*) gt_array_get_first(sa->exons))
                ->leftgenomicexonborder;
  range.end = ((Exoninfo*) gt_array_get_last(sa->exons))
              ->rightgenomicexonborder;

  return range;
}

GthFlt gth_sa_average_splice_site_prob(const GthSA *sa)
{
  unsigned long i, numofintrons;
  GthFlt averagepdpa = 0.0;
  Introninfo *introninfo;
  gt_assert(sa);
  numofintrons = gt_array_size(sa->introns);
  if (numofintrons > 0) {
    for (i = 0; i < numofintrons; i++) {
      introninfo = gt_array_get(sa->introns, i);
      averagepdpa += introninfo->donorsiteprobability;
      averagepdpa += introninfo->acceptorsiteprobability;
    }
    averagepdpa /= (2 * numofintrons);
  }
  return averagepdpa;
}

bool gth_sa_is_poor(const GthSA *sa, GthFlt minaveragessp)
{
  unsigned long num_of_introns;
  GthFlt averagepdpa;

  gt_assert(sa);

  num_of_introns = gt_array_size(sa->introns);
  averagepdpa = gth_sa_average_splice_site_prob(sa);

  if (num_of_introns == 0 &&
      sa->polyAtailpos.start > sa->polyAtailpos.end)
    return true;
  if (num_of_introns == 1 &&
      averagepdpa < (GthFlt) (1.5 * minaveragessp)) {
    return true;
  }
  if (num_of_introns > 1 && averagepdpa < minaveragessp)
    return true;

  return false;
}

bool gth_sa_B_is_better_than_A(const GthSA *saA, const GthSA *saB)
{
  GthFlt firstaveragepdpa, secondaveragepdpa;

  gt_assert(saA && saB);

  firstaveragepdpa  = gth_sa_average_splice_site_prob(saA);
  secondaveragepdpa = gth_sa_average_splice_site_prob(saB);

  if (gt_array_size(saA->introns) == 0 ||
      (gt_array_size(saB->introns) > 0 &&
       secondaveragepdpa > firstaveragepdpa)) {
    return true;
  }
  else
    return false;
}

unsigned long gth_sa_left_genomic_exon_border(const GthSA *sa,
                                              unsigned long exon)
{
  Exoninfo *exoninfo;
  gt_assert(sa);
  exoninfo = gth_sa_get_exon(sa, exon);
  return SHOWGENPOS(sa->gen_strand_forward,
                    sa->gen_total_length,
                    sa->gen_offset, exoninfo->leftgenomicexonborder);
}

unsigned long gth_sa_right_genomic_exon_border(const GthSA *sa,
                                               unsigned long exon)
{
  Exoninfo *exoninfo;
  gt_assert(sa);
  exoninfo = gth_sa_get_exon(sa, exon);
  return SHOWGENPOS(sa->gen_strand_forward,
                    sa->gen_total_length,
                    sa->gen_offset, exoninfo->rightgenomicexonborder);
}

double gth_sa_exon_score(const GthSA *sa, unsigned long exon)
{
  Exoninfo *exoninfo;
  gt_assert(sa);
  exoninfo = gth_sa_get_exon(sa, exon);
  return exoninfo->exonscore;
}

GtRange gth_sa_donor_site_range(const GthSA *sa, unsigned long intron)
{
  GtRange range;
  gt_assert(sa);
  range.start = gth_sa_left_intron_border(sa, intron);
  range.end = range.start + 1;
  return range;
}

GtRange gth_sa_acceptor_site_range(const GthSA *sa, unsigned long intron)
{
  GtRange range;
  gt_assert(sa);
  range.end = gth_sa_right_intron_border(sa, intron);
  range.start = range.end - 1;
  return range;
}

float gth_sa_donor_site_prob(const GthSA *sa, unsigned long intron)
{
  float prob;
  gt_assert(sa);
  prob = ((Introninfo*) gt_array_get(sa->introns, intron))
         ->donorsiteprobability;
  if (0.0005 > prob) /* return only three positions after decimal point */
    return 0.0;
  return prob;
}

float gth_sa_acceptor_site_prob(const GthSA *sa, unsigned long intron)
{
  float prob;
  gt_assert(sa);
  prob = ((Introninfo*) gt_array_get(sa->introns, intron))
         ->acceptorsiteprobability;
  if (0.0005 > prob) /* return only three positions after decimal point */
    return 0.0;
  return prob;
}

unsigned long gth_sa_genomic_exon_length(const GthSA *sa, unsigned long exon)
{
  Exoninfo *exoninfo;
  gt_assert(sa);
  exoninfo = gth_sa_get_exon(sa, exon);
  return exoninfo->rightgenomicexonborder - exoninfo->leftgenomicexonborder + 1;
}

unsigned long gth_sa_left_intron_border(const GthSA *sa, unsigned long intron)
{
  Exoninfo *exoninfo;
  gt_assert(sa);
  exoninfo = gth_sa_get_exon(sa, intron);
  return SHOWGENPOS(sa->gen_strand_forward,
                    sa->gen_total_length,
                    sa->gen_offset,
                    exoninfo->rightgenomicexonborder + 1);
}

unsigned long gth_sa_right_intron_border(const GthSA *sa, unsigned long intron)
{
  Exoninfo *exoninfo;
  gt_assert(sa);
  exoninfo = gth_sa_get_exon(sa, intron + 1);
  return SHOWGENPOS(sa->gen_strand_forward,
                    sa->gen_total_length,
                    sa->gen_offset,
                    exoninfo->leftgenomicexonborder - 1);
}

unsigned long gth_sa_intron_length(const GthSA *sa, unsigned long intron)
{
  Exoninfo *left_exon, *right_exon;
  gt_assert(sa);
  left_exon  = gth_sa_get_exon(sa, intron);
  right_exon = gth_sa_get_exon(sa, intron + 1);
  return right_exon->leftgenomicexonborder - left_exon->rightgenomicexonborder
         - 1;
}

GthBacktracePath* gth_sa_backtrace_path(const GthSA *sa)
{
  gt_assert(sa);
  return sa->backtrace_path;
}

Editoperation* gth_sa_get_editoperations(const GthSA *sa)
{
  gt_assert(sa);
  return gth_backtrace_path_get(sa->backtrace_path);
}

unsigned long gth_sa_get_editoperations_length(const GthSA *sa)
{
  gt_assert(sa);
  return gth_backtrace_path_length(sa->backtrace_path);
}

unsigned long gth_sa_indelcount(const GthSA *sa)
{
  gt_assert(sa);
  return gth_backtrace_path_indelcount(sa->backtrace_path);
}

unsigned long gth_sa_gen_dp_length(const GthSA *sa)
{
  gt_assert(sa);
  return gth_backtrace_path_gen_dp_length(sa->backtrace_path);
}

unsigned long gth_sa_gen_total_length(const GthSA *sa)
{
  gt_assert(sa);
  return sa->gen_total_length;
}

void gth_sa_set_gen_total_length(GthSA *sa, unsigned long gen_total_length)
{
  gt_assert(sa);
  sa->gen_total_length = gen_total_length;
}

unsigned long gth_sa_gen_offset(const GthSA *sa)
{
  gt_assert(sa);
  return sa->gen_offset;
}

void gth_sa_set_gen_offset(GthSA *sa, unsigned long gen_offset)
{
  gt_assert(sa);
  sa->gen_offset = gen_offset;
}

unsigned long gth_sa_ref_total_length(const GthSA *sa)
{
  gt_assert(sa);
  return sa->ref_total_length;
}

void gth_sa_set_ref_total_length(GthSA *sa, unsigned long reflen)
{
  gt_assert(sa);
  sa->ref_total_length = reflen;
  /* XXX: ??? */
  gth_backtrace_path_set_ref_dp_length(sa->backtrace_path, reflen);
}

unsigned long gth_sa_gen_dp_start(const GthSA *sa)
{
  gt_assert(sa);
  return gth_backtrace_path_gen_dp_start(sa->backtrace_path);
}

unsigned long gth_sa_gen_dp_start_show(const GthSA *sa)
{
  gt_assert(sa);
  return SHOWGENPOS(sa->gen_strand_forward,
                    sa->gen_total_length,
                    sa->gen_offset,
                    gth_backtrace_path_gen_dp_start(sa->backtrace_path));
}

void gth_sa_set_gen_dp_start(GthSA *sa, unsigned long gen_dp_start)
{
  gt_assert(sa);
  gth_backtrace_path_set_gen_dp_start(sa->backtrace_path, gen_dp_start);
}

unsigned long gth_sa_gen_dp_end(const GthSA *sa)
{
  gt_assert(sa);
  return gth_sa_gen_dp_start(sa) + gth_sa_gen_dp_length(sa) - 1;
}

unsigned long gth_sa_gen_dp_end_show(const GthSA *sa)
{
  gt_assert(sa);
  return SHOWGENPOS(sa->gen_strand_forward,
                    sa->gen_total_length,
                    sa->gen_offset,
                    gth_sa_gen_dp_start(sa) + gth_sa_gen_dp_length(sa) - 1);
}

unsigned long gth_sa_gen_file_num(const GthSA *sa)
{
  gt_assert(sa);
  return sa->gen_file_num;
}

void gth_sa_set_gen_file_num(GthSA *sa, unsigned long filenum)
{
  gt_assert(sa);
  sa->gen_file_num = filenum;
}

unsigned long gth_sa_gen_seq_num(const GthSA *sa)
{
  gt_assert(sa);
  return sa->gen_seq_num;
}

void gth_sa_set_gen_seq_num(GthSA *sa, unsigned long seqnum)
{
  gt_assert(sa);
  sa->gen_seq_num = seqnum;
}

unsigned long gth_sa_ref_file_num(const GthSA *sa)
{
  gt_assert(sa);
  return sa->ref_file_num;
}

void gth_sa_set_ref_file_num(GthSA *sa, unsigned long filenum)
{
  gt_assert(sa);
  sa->ref_file_num = filenum;
}

unsigned long gth_sa_ref_seq_num(const GthSA *sa)
{
  gt_assert(sa);
  return sa->ref_seq_num;
}

void gth_sa_set_ref_seq_num(GthSA *sa, unsigned long seqnum)
{
  gt_assert(sa);
  sa->ref_seq_num = seqnum;
}

const char* gth_sa_gen_id(const GthSA *sa)
{
  gt_assert(sa);
  return gt_str_get(sa->gen_id);
}

GtStr* gth_sa_gen_id_str(const GthSA *sa)
{
  gt_assert(sa);
  return sa->gen_id;
}

void gth_sa_set_gen_id(GthSA *sa, const char *id)
{
  gt_assert(sa);
  gt_str_set(sa->gen_id, id);
}

const char* gth_sa_ref_id(const GthSA *sa)
{
  gt_assert(sa);
  return gt_str_get(sa->ref_id);
}

GtStr* gth_sa_ref_id_str(const GthSA *sa)
{
  gt_assert(sa);
  return sa->ref_id;
}

void gth_sa_set_ref_id(GthSA *sa, const char *id)
{
  gt_assert(sa);
  gt_str_set(sa->ref_id, id);
}

GtStr* gth_sa_gen_md5(const GthSA *sa)
{
  gt_assert(sa && sa->gen_md5);
  return sa->gen_md5;
}

GtStr* gth_sa_ref_md5(const GthSA *sa)
{
  gt_assert(sa && sa->ref_md5);
  return sa->ref_md5;
}

GtStr* gth_sa_gen_desc(const GthSA *sa)
{
  gt_assert(sa && sa->gen_desc);
  return sa->gen_desc;
}

GtStr* gth_sa_ref_desc(const GthSA *sa)
{
  gt_assert(sa && sa->ref_desc);
  return sa->ref_desc;
}

GtStrand gth_sa_gen_strand(const GthSA *sa)
{
  gt_assert(sa);
  return sa->gen_strand_forward ? GT_STRAND_FORWARD : GT_STRAND_REVERSE;
}

bool gth_sa_gen_strand_forward(const GthSA *sa)
{
  gt_assert(sa);
  return sa->gen_strand_forward;
}

char gth_sa_gen_strand_char(const GthSA *sa)
{
  gt_assert(sa);
  return SHOWSTRAND(sa->gen_strand_forward);
}

void gth_sa_set_gen_strand(GthSA *sa, bool forward)
{
  gt_assert(sa);
  sa->gen_strand_forward = forward;
}

bool gth_sa_ref_strand_forward(const GthSA *sa)
{
  gt_assert(sa);
  return sa->ref_strand_forward;
}

char gth_sa_ref_strand_char(const GthSA *sa)
{
  gt_assert(sa);
  return SHOWSTRAND(sa->ref_strand_forward);
}

void gth_sa_set_ref_strand(GthSA *sa, bool forward)
{
  gt_assert(sa);
  sa->ref_strand_forward = forward;
}

unsigned long gth_sa_genomiccutoff_start(const GthSA *sa)
{
  gt_assert(sa);
  return gth_backtrace_path_genomiccutoff_start(sa->backtrace_path);
}

unsigned long gth_sa_referencecutoff_start(const GthSA *sa)
{
  gt_assert(sa);
  return gth_backtrace_path_referencecutoff_start(sa->backtrace_path);
}

unsigned long gth_sa_eopcutoff_start(const GthSA *sa)
{
  gt_assert(sa);
  return gth_backtrace_path_eopcutoff_start(sa->backtrace_path);
}

unsigned long gth_sa_genomiccutoff_end(const GthSA *sa)
{
  gt_assert(sa);
  return gth_backtrace_path_genomiccutoff_end(sa->backtrace_path);
}

unsigned long gth_sa_referencecutoff_end(const GthSA *sa)
{
  gt_assert(sa);
  return gth_backtrace_path_referencecutoff_end(sa->backtrace_path);
}

unsigned long gth_sa_eopcutoff_end(const GthSA *sa)
{
  gt_assert(sa);
  return gth_backtrace_path_eopcutoff_end(sa->backtrace_path);
}

void gth_sa_set_cutoffs_start(GthSA *sa, Cutoffs *cutoffs)
{
 gt_assert(sa && cutoffs);
 gth_backtrace_path_set_cutoffs_start(sa->backtrace_path, cutoffs);
}

void gth_sa_set_cutoffs_end(GthSA *sa, Cutoffs *cutoffs)
{
 gt_assert(sa && cutoffs);
 gth_backtrace_path_set_cutoffs_end(sa->backtrace_path, cutoffs);
}

GthAlphatype gth_sa_alphatype(const GthSA *sa)
{
  gt_assert(sa);
  return gth_backtrace_path_alphatype(sa->backtrace_path);
}

const char* gth_sa_alphastring(const GthSA *sa)
{
  gt_assert(sa);
  return gth_sa_alphatype(sa) == DNA_ALPHA ? "cDNA" : "Protein";
}

void gth_sa_set_alphatype(GthSA *sa, GthAlphatype alphatype)
{
  gt_assert(sa);
  gth_backtrace_path_set_alphatype(sa->backtrace_path, alphatype);
}

Exoninfo* gth_sa_get_exon(const GthSA *sa, unsigned long exon)
{
  gt_assert(sa && sa->exons);
  gt_assert(exon < gt_array_size(sa->exons));
  return gt_array_get(sa->exons, exon);
}

void gth_sa_add_exon(GthSA *sa, Exoninfo *exoninfo)
{
  gt_assert(sa && exoninfo);
  gt_array_add(sa->exons, *exoninfo);
}

unsigned long gth_sa_num_of_exons(const GthSA *sa)
{
  gt_assert(sa && sa->exons);
  return gt_array_size(sa->exons);
}

Introninfo* gth_sa_get_intron(const GthSA *sa, unsigned long intron)
{
  gt_assert(sa && sa->introns);
  gt_assert(intron < gt_array_size(sa->exons));
  return gt_array_get(sa->introns, intron);
}

void gth_sa_add_intron(GthSA *sa, Introninfo *introninfo)
{
  gt_assert(sa && introninfo);
  gt_array_add(sa->introns, *introninfo);
}

unsigned long gth_sa_num_of_introns(const GthSA *sa)
{
  gt_assert(sa && sa->introns);
  return gt_array_size(sa->introns);
}

void gth_sa_calc_polyAtailpos(GthSA *sa, const unsigned char *ref_seq_tran,
                              GtAlphabet *ref_alphabet)
{
  unsigned long ppa, mma, rightreferenceborder, referencelength;
  long i, leftreferenceborder;

  sa->polyAtailpos.start = 0;
  sa->polyAtailpos.end = 0;
  ppa = mma = 0;

  rightreferenceborder = ((Exoninfo*) gt_array_get_last(sa->exons))
                         ->rightreferenceexonborder;
  leftreferenceborder  = ((Exoninfo*) gt_array_get_first(sa->exons))
                         ->leftreferenceexonborder;

  /* setting i */
  referencelength = gth_sa_ref_total_length(sa);
  if ((rightreferenceborder + 1) >=
      (referencelength - 1 - CALCPOLYATAILWINDOW)) {
    i = gt_safe_cast2long(rightreferenceborder + 1);
  }
  else {
    if (referencelength < 1 + CALCPOLYATAILWINDOW)
      i = 0;
    else
      i =  referencelength - 1 - CALCPOLYATAILWINDOW;
  }

  for (/* i already set */; i < gt_safe_cast2long(referencelength); i++) {
    if (ref_seq_tran[i] == gt_alphabet_encode(ref_alphabet, 'A'))
      ppa++;
    else {
      if (ppa > 0 && mma < 1) {
        mma++;
        continue;
      }
      else {
        if (ppa >= MINIMUMPOLYATAILLENGTH)
          break;
        else {
          ppa = mma = 0;
          continue;
        }
      }
    }
  }

  if (ppa >= MINIMUMPOLYATAILLENGTH) {
    sa->polyAtailpos.start = gt_safe_cast2ulong(i - ppa - mma);
    sa->polyAtailpos.end = i - 1;
  }
  else {
    ppa = mma = 0;

    /* setting i */
    if ((leftreferenceborder - 1) <= CALCPOLYATAILWINDOW)
      i = leftreferenceborder - 1;
    else
      i =  CALCPOLYATAILWINDOW - 1;

    for (/* i already set */; i >= 0; i--) {
      if (ref_seq_tran[i] == gt_alphabet_encode(ref_alphabet, 'T'))
        ppa++;
      else {
        if (ppa > 0 && mma < 1) {
          mma++;
          continue;
        }
        else {
          if (ppa >= MINIMUMPOLYATAILLENGTH)
            break;
          else {
            ppa = mma = 0;
            continue;
          }
        }
      }
    }

    if (ppa >= MINIMUMPOLYATAILLENGTH) {
      sa->polyAtailpos.start  = gt_safe_cast2ulong(i + ppa + mma);
      sa->polyAtailpos.end = i + 1;
    }
  }
}

unsigned long gth_sa_polyAtail_start(const GthSA *sa)
{
  gt_assert(sa);
  return sa->polyAtailpos.start;
}

unsigned long gth_sa_polyAtail_stop(const GthSA *sa)
{
  gt_assert(sa);
  return sa->polyAtailpos.end;
}

void gth_sa_set_polyAtail_start(GthSA *sa, unsigned long start)
{
  gt_assert(sa);
  sa->polyAtailpos.start = start;
}

void gth_sa_set_polyAtail_stop(GthSA *sa, unsigned long stop)
{
  gt_assert(sa);
  sa->polyAtailpos.end = stop;
}

GthFlt gth_sa_score(const GthSA *sa)
{
  gt_assert(sa);
  return sa->alignmentscore;
}

void gth_sa_set_score(GthSA *sa, GthFlt score)
{
  gt_assert(sa);
  sa->alignmentscore = score;
}

GthFlt gth_sa_coverage(const GthSA *sa)
{
  gt_assert(sa);
  return sa->coverage;
}

void gth_sa_set_coverage(GthSA *sa, GthFlt coverage)
{
  gt_assert(sa);
  sa->coverage = coverage;
}

bool gth_sa_genomic_cov_is_highest(const GthSA *sa)
{
  gt_assert(sa);
  return sa->genomic_cov_is_highest;
}

char gth_sa_coverage_char(const GthSA *sa)
{
  gt_assert(sa);
  if (gth_sa_genomic_cov_is_highest(sa))
    return 'G';
  else {
    if (gth_sa_alphatype(sa) == DNA_ALPHA)
      return 'C';
    else
      return 'P';
  }
}

void gth_sa_set_highest_cov(GthSA *sa, bool genomic)
{
  gt_assert(sa);
  sa->genomic_cov_is_highest = genomic;
}

unsigned long gth_sa_cumlen_scored_exons(const GthSA *sa)
{
  gt_assert(sa);
  return sa->cumlen_scored_exons;
}

void gth_sa_set_cumlen_scored_exons(GthSA *sa, unsigned long cumlen)
{
  gt_assert(sa);
  sa->cumlen_scored_exons = cumlen;
}

unsigned long gth_sa_call_number(const GthSA *sa)
{
  gt_assert(sa);
  return sa->call_number;
}

static void set_gff3_target_attribute(GthSA *sa, bool md5ids)
{
  gt_assert(sa && !sa->gff3_target_attribute);
  sa->gff3_target_attribute = gt_str_new();
  if (md5ids) {
    gt_assert(sa->ref_md5);
    gt_str_append_cstr(sa->gff3_target_attribute, GT_MD5_SEQID_PREFIX);
    gt_str_append_str(sa->gff3_target_attribute, sa->ref_md5);
    gt_str_append_char(sa->gff3_target_attribute, ':');
  }
  gt_gff3_escape(sa->gff3_target_attribute, gt_str_get(sa->ref_id),
                 gt_str_length(sa->ref_id));
  gt_str_append_char(sa->gff3_target_attribute, ' ');
  gt_str_append_ulong(sa->gff3_target_attribute,
                      gth_sa_referencecutoff_start(sa) + 1); /* XXX: use
                                                                reference
                                                                dpstartpos */
  gt_str_append_char(sa->gff3_target_attribute, ' ');
  gt_str_append_ulong(sa->gff3_target_attribute,
                      gth_sa_ref_total_length(sa) - /* XXX */
                      gth_sa_referencecutoff_end(sa));
  gt_str_append_char(sa->gff3_target_attribute, ' ');
  if (sa->ref_strand_forward) {
    gt_str_append_char(sa->gff3_target_attribute,
                       GT_STRAND_CHARS[GT_STRAND_FORWARD]);
  }
  else {
    gt_str_append_char(sa->gff3_target_attribute,
                       GT_STRAND_CHARS[GT_STRAND_REVERSE]);
  }
}

const char* gth_sa_gff3_target_attribute(GthSA *sa, bool md5ids)
{
  gt_assert(sa);
  if (!sa->gff3_target_attribute && (md5ids || gt_str_length(sa->ref_id)))
    set_gff3_target_attribute(sa, md5ids);
  return gt_str_get(sa->gff3_target_attribute);
}

void gth_sa_determine_cutoffs(GthSA *sa, GthCutoffmode leadcutoffsmode,
                              GthCutoffmode termcutoffsmode,
                              unsigned long cutoffsminexonlen)
{
  gt_assert(sa);
  gth_backtrace_path_determine_cutoffs(sa->backtrace_path, leadcutoffsmode,
                                       termcutoffsmode, cutoffsminexonlen);
}

void gth_sa_cutoff_start(GthSA *sa)
{
  gt_assert(sa);
  gth_backtrace_path_cutoff_start(sa->backtrace_path);
}

void gth_sa_cutoff_end(GthSA *sa)
{
  gt_assert(sa);
  gth_backtrace_path_cutoff_end(sa->backtrace_path);
}

void gth_sa_cutoff_walked_path(GthSA *sa, const GthPathWalker *pw,
                               bool showeops, GtFile *outfp)
{
  gt_assert(sa && pw);
  gth_backtrace_path_cutoff_walked_path(sa->backtrace_path, pw, showeops,
                                        outfp);
}

void gth_sa_prepend(GthSA *sa, const GthBacktracePath *eops)
{
  gt_assert(sa && eops);
  gth_backtrace_path_prepend(sa->backtrace_path, eops);
}

void gth_sa_append(GthSA *sa, const GthBacktracePath *eops)
{
  gt_assert(sa && eops);
  gth_backtrace_path_append(sa->backtrace_path, eops);
}

void gth_sa_remove_zero_base_exons(GthSA *sa, GthStat *stat)
{
  gt_assert(sa);
  gth_backtrace_path_remove_zero_base_exons(sa->backtrace_path, stat);
}

bool gth_sa_contains_no_zero_base_exons(const GthSA *sa)
{
  gt_assert(sa);
  return gth_backtrace_path_contains_no_zero_base_exons(sa->backtrace_path);
}

void gth_sa_echo_genomic_description(const GthSA *sa, GthInput *input,
                                     GtFile *outfp)
{
  gt_assert(sa && input);
  gth_input_echo_genomic_description(input, sa->gen_file_num,
                                     sa->gen_seq_num, outfp);
}

void gth_sa_echo_reference_description(const GthSA *sa, GthInput *input,
                                       GtFile *outfp)
{
  gt_assert(sa && input);
  gth_input_echo_reference_description(input, sa->ref_file_num,
                                       sa->ref_seq_num, outfp);
}

void gth_sa_echo_reference_sequence(const GthSA *sa, GthInput *input,
                                    bool format, GtFile *outfp)
{
  gt_assert(sa && input);
  gth_input_echo_reference_sequence(input, format, sa->ref_file_num,
                                    sa->ref_seq_num, sa->ref_strand_forward,
                                    outfp);
}

void gth_sa_echo_alignment(const GthSA *sa, unsigned long showintronmaxlen,
                           unsigned long translationtable,
                           bool wildcardimplosion, GthInput *input,
                           GtFile *outfp)
{
  unsigned long genomicstartcutoff, genomicendcutoff, genomictotalcutoff,
                referencestartcutoff, referenceendcutoff, referencetotalcutoff;
  bool reverse_subject_pos = false;
  const unsigned char *gen_seq_orig, *ref_seq_orig;
  GthSeqCon *ref_seq_con;
  GtAlphabet *ref_alphabet;

  gt_assert(sa && input);

  /* only for cosmetic reasons */
  genomicstartcutoff   = gth_sa_genomiccutoff_start(sa);
  genomicendcutoff     = gth_sa_genomiccutoff_end(sa);
  genomictotalcutoff   = genomicstartcutoff + genomicendcutoff;
  referencestartcutoff = gth_sa_referencecutoff_start(sa);
  referenceendcutoff   = gth_sa_referencecutoff_end(sa);
  referencetotalcutoff = referencestartcutoff + referenceendcutoff;

  /* make sure that the correct files are loaded */
  gth_input_load_reference_file(input, gth_sa_ref_file_num(sa), false);
  ref_seq_con = gth_input_current_ref_seq_con(input);
  ref_alphabet = gth_input_current_ref_alphabet(input);

  /* If the reverse complement of the genomic DNA is considered, this
     opition is needed for correct output of the genomic sequence positions
     by the function showalignmentgeneric() */
  if (!gth_sa_gen_strand_forward(sa))
    reverse_subject_pos = true;

  /* get genomic sequence */
  gen_seq_orig =
    gth_input_original_genomic_sequence(input, sa->gen_file_num,
                                        sa->gen_strand_forward)
    + gth_sa_gen_dp_start(sa);

  /* get reference sequence */
  if (gth_sa_ref_strand_forward(sa)) {
    ref_seq_orig =
      gth_seq_con_get_orig_seq(ref_seq_con, gth_sa_ref_seq_num(sa));
  }
  else {
    ref_seq_orig =
      gth_seq_con_get_orig_seq_rc(ref_seq_con, gth_sa_ref_seq_num(sa));
  }

  switch (gth_sa_alphatype(sa)) {
    case DNA_ALPHA:
      gthshowalignmentdna(outfp,ALIGNMENTLINEWIDTH,
                          gth_sa_get_editoperations(sa),
                          gth_sa_get_editoperations_length(sa),
                          gth_sa_indelcount(sa),
                          gen_seq_orig + genomicstartcutoff,
                          gth_sa_gen_dp_length(sa) - genomictotalcutoff,
                          ref_seq_orig + referencestartcutoff,
                          gth_sa_ref_total_length(sa) -
                          referencetotalcutoff,
                          gth_sa_gen_dp_start(sa) + genomicstartcutoff -
                          gth_sa_gen_offset(sa), referencestartcutoff,
                          gth_sa_gen_total_length(sa), showintronmaxlen,
                          ref_alphabet, reverse_subject_pos,
                          wildcardimplosion);
      break;
    case PROTEIN_ALPHA:
      gthshowalignmentprotein(outfp, ALIGNMENTLINEWIDTH,
                              gth_sa_get_editoperations(sa),
                              gth_sa_get_editoperations_length(sa),
                              gth_sa_indelcount(sa),
                              gen_seq_orig + genomicstartcutoff,
                              gth_sa_gen_dp_length(sa) - genomictotalcutoff,
                              ref_seq_orig + referencestartcutoff,
                              gth_sa_ref_total_length(sa) -
                              referencetotalcutoff,
                              gth_sa_gen_dp_start(sa) + genomicstartcutoff -
                              gth_sa_gen_offset(sa), referencestartcutoff,
                              gth_sa_gen_total_length(sa), showintronmaxlen,
                              ref_alphabet, translationtable,
                              gth_input_score_matrix(input),
                              gth_input_score_matrix_alpha(input),
                              reverse_subject_pos, wildcardimplosion);
      break;
    default: gt_assert(0);
  }
}

unsigned long gth_sa_get_alignment_lines(const GthSA *sa,
                                         unsigned char **first_line,
                                         unsigned char **second_line,
                                         unsigned char **third_line,
                                         unsigned long translationtable,
                                         GthInput *input)
{
  unsigned long genomicstartcutoff, genomicendcutoff, genomictotalcutoff,
                referencestartcutoff, referenceendcutoff, referencetotalcutoff;
  GT_UNUSED bool reverse_subject_pos = false;

  gt_assert(sa && first_line && second_line && third_line && input);

  /* only for cosmetic reasons */
  genomicstartcutoff   = gth_sa_genomiccutoff_start(sa);
  genomicendcutoff     = gth_sa_genomiccutoff_end(sa);
  genomictotalcutoff   = genomicstartcutoff + genomicendcutoff;
  referencestartcutoff = gth_sa_referencecutoff_start(sa);
  referenceendcutoff   = gth_sa_referencecutoff_end(sa);
  referencetotalcutoff = referencestartcutoff + referenceendcutoff;

  /* sequences */
  unsigned char *gen_seq_orig, *ref_seq_orig;
  unsigned long cols = 0;
  GthSeqCon *ref_seq_con;

  /* make sure that the correct files are loaded */
  gth_input_load_reference_file(input, gth_sa_ref_file_num(sa), false);
  ref_seq_con = gth_input_current_ref_seq_con(input);

  /* If the reverse complement of the genomic DNA is considered, this
     opition is needed for correct output of the genomic sequence positions
     by the function showalignmentgeneric() */
  if (!gth_sa_gen_strand_forward(sa))
    reverse_subject_pos = true;

  /* get genomic sequence */
  gen_seq_orig = (unsigned char*)
    gth_input_original_genomic_sequence(input, gth_sa_gen_file_num(sa),
                                        gth_sa_gen_strand_forward(sa))
    + gth_sa_gen_dp_start(sa);

  /* get reference sequence */
  if (gth_sa_ref_strand_forward(sa)) {
    ref_seq_orig =
      gth_seq_con_get_orig_seq(ref_seq_con, gth_sa_ref_seq_num(sa));
  }
  else {
    ref_seq_orig =
      gth_seq_con_get_orig_seq_rc(ref_seq_con, gth_sa_ref_seq_num(sa));
  }

  switch (gth_sa_alphatype(sa)) {
    case DNA_ALPHA:
      /* compute the two alignment lines */
      cols = gthfillthetwoalignmentlines(first_line,
                                         second_line,
                                         gen_seq_orig +
                                         genomicstartcutoff,
                                         gth_sa_gen_dp_length(sa) -
                                         genomictotalcutoff,
                                         ref_seq_orig +
                                         referencestartcutoff,
                                         gth_sa_ref_total_length(sa) -
                                         referencetotalcutoff,
                                         gth_sa_get_editoperations(sa),
                                         gth_sa_get_editoperations_length(sa),
                                         0,   /* linewidth not important here */
                                         0,   /* no short introns here */
                                         NULL,/* therefore no shortintroninfo */
                                         gth_sa_indelcount(sa));
      *third_line = NULL;
      break;
    case PROTEIN_ALPHA:
      /* compute the three alignment lines */
      cols = gthfillthethreealignmentlines(first_line,
                                           second_line,
                                           third_line,
                                           gth_sa_get_editoperations(sa),
                                           gth_sa_get_editoperations_length(sa),
                                           gth_sa_indelcount(sa),
                                           gen_seq_orig +
                                           genomicstartcutoff,
                                           gth_sa_gen_dp_length(sa) -
                                           genomictotalcutoff,
                                           ref_seq_orig +
                                           referencestartcutoff,
                                           gth_sa_ref_total_length(sa) -
                                           referencetotalcutoff,
                                           translationtable);
      break;
    default: gt_assert(0);
  }

  return cols;
}

bool gth_sa_is_valid(const GthSA *sa)
{
  gt_assert(sa);
  return gth_backtrace_path_is_valid(sa->backtrace_path);
}

void gth_sa_show(GthSA *sa, GthInput *input, GtFile *outfp)
{
  GthSAVisitor *sa_visitor;
  gt_assert(sa && input);
  gth_input_load_genomic_file(input, sa->gen_file_num, false);
  gth_input_load_reference_file(input, sa->ref_file_num, false);
  sa_visitor = gth_txt_sa_visitor_new(input,
                                      GTH_DEFAULT_GS2OUT,
                                      GTH_DEFAULT_DPMININTRONLENGTH,
                                      6, /* XXX */
                                      GTH_DEFAULT_SHOWINTRONMAXLEN,
                                      GTH_DEFAULT_TRANSLATIONTABLE,
                                      GTH_DEFAULT_SHOWSEQNUMS,
                                      outfp);
  gth_sa_visitor_visit_sa(sa_visitor, sa);
  gth_sa_visitor_delete(sa_visitor);
}

void gth_sa_save_ref_md5(GthSA *sa, GthInput *input)
{
  gt_assert(sa && input);
  gth_input_save_ref_md5(input, &sa->ref_md5, sa->ref_file_num,
                         sa->ref_seq_num);
}

bool gth_sas_are_equal(const GthSA *saA, const GthSA *saB)
{
  Exoninfo *exoninfoA, *exoninfoB;
  Introninfo *introninfoA, *introninfoB;
  unsigned long i;

  /* compare element 0 */
  if (gth_sa_alphatype(saA) != gth_sa_alphatype(saB))
    return false;

  /* compare element 1 */
  if (gth_backtrace_path_length(saA->backtrace_path) !=
      gth_backtrace_path_length(saB->backtrace_path)) {
    return false;
  }
  for (i = 0; i < gth_backtrace_path_length(saA->backtrace_path); i++) {
    if (((Editoperation*) gth_backtrace_path_get(saA->backtrace_path))[i] !=
        ((Editoperation*) gth_backtrace_path_get(saB->backtrace_path))[i]) {
      return false;
    }
  }

  /* element 2 has been removed (indelcount) */

  /* compare element 3 */
  if (gth_sa_gen_dp_length(saA) != gth_sa_gen_dp_length(saB))
    return false;

  /* compare element 4 */
  if (saA->gen_total_length != saB->gen_total_length)
    return false;

  /* compare element 5 */
  if (saA->gen_offset != saB->gen_offset)
    return false;

  /* compare element 6 */
  if (gth_sa_ref_total_length(saA) != gth_sa_ref_total_length(saB))
    return false;

  /* compare element 7 */
  if (gth_sa_gen_dp_start(saA) != gth_sa_gen_dp_start(saB))
    return false;

  /* element 8 has been removed (gen_dp_end) */

  /* compare element 9 */
  if (saA->gen_file_num != saB->gen_file_num)
    return false;

  /* compare element 10 */
  if (saA->gen_seq_num != saB->gen_seq_num)
    return false;

  /* compare element 11 */
  if (saA->ref_file_num != saB->ref_file_num)
    return false;

  /* compare element 12 */
  if (saA->ref_seq_num != saB->ref_seq_num)
    return false;

  /* compare element 13 */
  if (gt_str_cmp(saA->gen_id, saB->gen_id))
    return false;

  /* compare element 14 */
  if (gt_str_cmp(saA->ref_id, saB->ref_id))
    return false;

  /* compare element 15 */
  if (saA->gen_strand_forward != saB->gen_strand_forward)
    return false;

  /* compare element 16 */
  if (saA->ref_strand_forward != saB->ref_strand_forward)
    return false;

  /* compare element 17 */
  if (gth_sa_genomiccutoff_start(saA) != gth_sa_genomiccutoff_start(saB))
    return false;
  if (gth_sa_referencecutoff_start(saA) != gth_sa_referencecutoff_start(saB))
    return false;
  if (gth_sa_eopcutoff_start(saA) != gth_sa_eopcutoff_start(saB))
    return false;
  if (gth_sa_genomiccutoff_end(saA) != gth_sa_genomiccutoff_end(saB))
    return false;
  if (gth_sa_referencecutoff_end(saA) != gth_sa_referencecutoff_end(saB))
    return false;
  if (gth_sa_eopcutoff_end(saA) != gth_sa_eopcutoff_end(saB))
    return false;

  /* compare element 18 */
  if (gt_array_size(saA->exons) != gt_array_size(saB->exons))
    return false;
  for (i = 0; i < gt_array_size(saA->exons); i++) {
    exoninfoA = (Exoninfo*) gt_array_get(saA->exons, i);
    exoninfoB = (Exoninfo*) gt_array_get(saB->exons, i);
    if (exoninfoA->leftgenomicexonborder != exoninfoB->leftgenomicexonborder)
      return false;
    if (exoninfoA->rightgenomicexonborder != exoninfoB->rightgenomicexonborder)
      return false;
    if (exoninfoA->leftreferenceexonborder !=
        exoninfoB->leftreferenceexonborder) {
      return false;
    }
    if (exoninfoA->rightreferenceexonborder !=
        exoninfoB->rightreferenceexonborder) {
      return false;
    }
    if (!gt_double_equals_double(exoninfoA->exonscore, exoninfoB->exonscore)) {
      return false;
    }
  }

  /* compare element 19 */
  if (gt_array_size(saA->introns) != gt_array_size(saB->introns))
    return false;
  for (i = 0; i < gt_array_size(saA->introns); i++) {
    introninfoA = (Introninfo*) gt_array_get(saA->introns, i);
    introninfoB = (Introninfo*) gt_array_get(saB->introns, i);
    if (!gt_double_equals_double(introninfoA->donorsiteprobability,
                                 introninfoB->donorsiteprobability)) {
      return false;
    }
    if (!gt_double_equals_double(introninfoA->acceptorsiteprobability,
                                 introninfoB->acceptorsiteprobability)) {
      return false;
    }
    if (!gt_double_equals_double(introninfoA->donorsitescore,
                                 introninfoB->donorsitescore)) {
      return false;
    }
    if (!gt_double_equals_double(introninfoA->acceptorsitescore,
                                 introninfoB->acceptorsitescore)) {
      return false;
    }
  }

  /* compare element 20 */
  if (saA->polyAtailpos.start != saB->polyAtailpos.start)
    return false;
  if (saA->polyAtailpos.end != saB->polyAtailpos.end)
    return false;

  /* compare element 21 */
  if (saA->alignmentscore != saB->alignmentscore)
    return false;

  /* compare element 22 */
  if (saA->coverage != saB->coverage)
    return false;

  /* compare element 23 */
  if (saA->genomic_cov_is_highest != saB->genomic_cov_is_highest)
    return false;

  /* compare element 24 */
  if (saA->cumlen_scored_exons != saB->cumlen_scored_exons)
    return false;

  return true;
}
