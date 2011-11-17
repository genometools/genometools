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

#include "core/ma_api.h"
#include "core/range.h"
#include "core/safearith.h"
#include "core/undef_api.h"
#include "core/unused_api.h"
#include "gth/gthoutput.h"
#include "gth/gthchain.h"

#define DPEXTENSION             300

/* An GthInvertedChain is the inversion of a GthChain.
   That is, its genomic ranges denote the potential introns.
   The variables <startpos> and <endpos> denote the leftmost genomic position of
   the orginal exons and the rightmost postions, respectively. They are needed
   to be able to convert a GthInvertedChain (back) to a GthChain. */
typedef struct {
  unsigned long gen_file_num,
                gen_seq_num,
                ref_file_num,
                ref_seq_num,
                startpos,
                endpos;
  GtArray *forwardranges;
} GthInvertedChain;

GthChain* gth_chain_new(void)
{
  GthChain *chain = gt_calloc(1, sizeof *chain);
  chain->gen_file_num  = GT_UNDEF_ULONG;
  chain->gen_seq_num   = GT_UNDEF_ULONG;
  chain->ref_file_num  = GT_UNDEF_ULONG;
  chain->ref_seq_num   = GT_UNDEF_ULONG;
  chain->forwardranges = gt_array_new(sizeof (GtRange));
  chain->reverseranges = gt_array_new(sizeof (GtRange));
  return chain;
}

void gth_chain_delete(GthChain *chain)
{
  if (!chain) return;
  if (chain->jump_table_delete) {
    chain->jump_table_delete(chain->reverse_jump_table);
    chain->jump_table_delete(chain->forward_jump_table);
  }
  gt_array_delete(chain->forwardranges);
  gt_array_delete(chain->reverseranges);
  gt_free(chain);
}

static void inverted_chain_init(GthInvertedChain *inverted_chain)
{
  inverted_chain->gen_file_num  = GT_UNDEF_ULONG;
  inverted_chain->gen_seq_num   = GT_UNDEF_ULONG;
  inverted_chain->ref_file_num  = GT_UNDEF_ULONG;
  inverted_chain->ref_seq_num   = GT_UNDEF_ULONG;
  inverted_chain->startpos      = GT_UNDEF_ULONG;
  inverted_chain->endpos        = GT_UNDEF_ULONG;
  inverted_chain->forwardranges = gt_array_new(sizeof (GtRange));
}

static void inverted_chain_free(GthInvertedChain *inverted_chain)
{
  if (!inverted_chain) return;
  gt_array_delete(inverted_chain->forwardranges);
}

static void chain_copy_core(GthChain *dest, const GthChain *src)
{
  /* destination chain is empty */
  gt_assert(!gt_array_size(dest->forwardranges) &&
            !gt_array_size(dest->reverseranges));
  /* source is not empty */
  gt_assert(gt_array_size(src->forwardranges) &&
            gt_array_size(src->reverseranges));
  /* number of forwardranges equals number of reverseranges */
  gt_assert(gt_array_size(src->forwardranges) ==
            gt_array_size(src->reverseranges));

  /* copy file and sequence numbers */
  dest->gen_file_num   = src->gen_file_num;
  dest->gen_seq_num    = src->gen_seq_num;
  dest->ref_file_num   = src->ref_file_num;
  dest->ref_seq_num    = src->ref_seq_num;
  dest->refseqcoverage = src->refseqcoverage;
}

void gth_chain_copy(GthChain *dest, const GthChain *src)
{
  gt_assert(dest&& src);
  chain_copy_core(dest, src);
  gt_array_add_array(dest->forwardranges, src->forwardranges);
  gt_array_add_array(dest->reverseranges, src->reverseranges);
}

#ifndef NDEBUG
static bool chain_is_filled_and_consistent(GthChain *chain,
                                           unsigned long gen_total_length,
                                           unsigned long gen_offset)
{
  GtArray *testranges;

  /* check of file sequence numbers are defined */
  if (chain->gen_file_num == GT_UNDEF_ULONG ||
      chain->gen_seq_num  == GT_UNDEF_ULONG ||
      chain->ref_file_num == GT_UNDEF_ULONG ||
      chain->ref_seq_num  == GT_UNDEF_ULONG) {
    return false;
  }

  if (!gt_ranges_are_consecutive(chain->forwardranges))
    return false;

  testranges = gt_array_new(sizeof (GtRange));
  gt_ranges_copy_to_opposite_strand(testranges, chain->reverseranges,
                                    gen_total_length, gen_offset);
  if (!gt_ranges_are_equal(testranges, chain->forwardranges)) {
    gt_array_delete(testranges);
    return false;
  }

  gt_array_delete(testranges);

  return true;
}
#endif

static void convert_chain_to_inverted_chain(GthInvertedChain *inverted_chain,
                                            GthChain *chain)
{
  unsigned long i, lastexonnum = gt_array_size(chain->forwardranges) - 1;
  GtRange range;

  /* inverted chain is empty */
  gt_assert(!gt_array_size(inverted_chain->forwardranges));
  /* chain is not empty */
  gt_assert(gt_array_size(chain->forwardranges));

  /* copy file and sequence numbers */
  inverted_chain->gen_file_num = chain->gen_file_num;
  inverted_chain->gen_seq_num  = chain->gen_seq_num;
  inverted_chain->ref_file_num = chain->ref_file_num;
  inverted_chain->ref_seq_num  = chain->ref_seq_num;

  /* save startpos */
  inverted_chain->startpos = ((GtRange*)
                              gt_array_get_first(chain->forwardranges))->start;

  /* save endpos */
  inverted_chain->endpos = ((GtRange*)
                             gt_array_get_last(chain->forwardranges))->end;

  /* convert (potential) exons to (potential) introns */
  for (i = 0; i < lastexonnum; i++) {
    range.start  = ((GtRange*) gt_array_get(chain->forwardranges, i))
                  ->end + 1;
    range.end = ((GtRange*) gt_array_get(chain->forwardranges, i+1))
                  ->start - 1;
    gt_array_add(inverted_chain->forwardranges, range);
  }
}

static void convert_inverted_chain_to_chain(GthChain *chain,
                                            GthInvertedChain *inverted_chain,
                                            unsigned long gen_total_length,
                                            unsigned long gen_offset)
{
  unsigned long i, numofintrons = gt_array_size(inverted_chain->forwardranges);
  GtRange range;

  /* chain is empty */
  gt_assert(!gt_array_size(chain->forwardranges) &&
            !gt_array_size(chain->reverseranges));

  /* copy file and sequence numbers */
  chain->gen_file_num = inverted_chain->gen_file_num;
  chain->gen_seq_num  = inverted_chain->gen_seq_num;
  chain->ref_file_num = inverted_chain->ref_file_num;
  chain->ref_seq_num  = inverted_chain->ref_seq_num;

  /* save forward ranges */
  range.start = inverted_chain->startpos;
  for (i = 0; i < numofintrons; i++) {
    range.end = ((GtRange*) gt_array_get(inverted_chain->forwardranges,
                                                i))->start - 1;
    gt_array_add(chain->forwardranges, range);
    range.start = ((GtRange*) gt_array_get(inverted_chain->forwardranges,
                                               i))->end + 1;
  }
  range.end = inverted_chain->endpos;
  gt_array_add(chain->forwardranges, range);

  /* save inverted forward ranges as reverse ranges */
  gt_ranges_copy_to_opposite_strand(chain->reverseranges, chain->forwardranges,
                                    gen_total_length, gen_offset);

  gt_assert(chain_is_filled_and_consistent(chain, gen_total_length,
                                           gen_offset));
}

#ifndef NDEBUG
static bool conversion_is_correct(GthChain *orig_chain,
                                  GthInvertedChain *inverted_chain,
                                  unsigned long gen_total_length,
                                  unsigned long gen_offset)
{
  GthChain *check_chain;
  unsigned long i;

  check_chain = gth_chain_new();
  convert_inverted_chain_to_chain(check_chain, inverted_chain, gen_total_length,
                                  gen_offset);

  /* compare number of (potential) exons */
  if ((gt_array_size(orig_chain->forwardranges) !=
       gt_array_size(check_chain->forwardranges)) ||
      (gt_array_size(orig_chain->reverseranges) !=
       gt_array_size(check_chain->reverseranges))) {
    gth_chain_delete(check_chain);
    return false;
  }

  /* compare positions of (potential) exon */
  for (i = 0; i < gt_array_size(orig_chain->forwardranges); i++) {
    if ((((GtRange*)gt_array_get(orig_chain->forwardranges, i))->start !=
         ((GtRange*)gt_array_get(check_chain->forwardranges, i))->start) ||
        (((GtRange*)gt_array_get(orig_chain->forwardranges, i))->end !=
         ((GtRange*)gt_array_get(check_chain->forwardranges, i))->end) ||
        (((GtRange*)gt_array_get(orig_chain->reverseranges, i))->start !=
         ((GtRange*)gt_array_get(check_chain->reverseranges, i))->start) ||
        (((GtRange*)gt_array_get(orig_chain->reverseranges, i))->end !=
         ((GtRange*)gt_array_get(check_chain->reverseranges, i))->end)) {
      gth_chain_delete(check_chain);
      return false;
    }
  }

  gth_chain_delete(check_chain);

  return true;
}
#endif

static void potentialintronspostpro(GtArray *intronstoprocess,
                                    unsigned long icdelta,
                                    unsigned long icminremintronlength)
{
  GtArray *originalintrons;
  GtRange potintron;
  unsigned long i, potintronlength,
       minintronlength = 2 * icdelta + icminremintronlength;

  originalintrons = gt_array_new(sizeof (GtRange));

  /* save all (potential) introns */
  gt_array_add_array(originalintrons, intronstoprocess);

  /* reset introns to process */
  gt_array_set_size(intronstoprocess, 0);

  /* store introns */
  for (i = 0; i < gt_array_size(originalintrons); i++) {
    potintron       = *(GtRange*) gt_array_get(originalintrons, i);
    potintronlength = potintron.end - potintron.start + 1;

    if (potintronlength >= minintronlength) {
      /* keep this intron (plus/minus intron deltas)
         that is, this intron is cut out later */
      potintron.start  += icdelta;
      potintron.end -= icdelta;
      gt_array_add(intronstoprocess, potintron);
    }
    /* else: skip this intron
       that is, this intron is not cut out later */
  }

  gt_array_delete(originalintrons);
}

/* XXX: change this function: add more sophisticated extension strategy */
void gth_chain_extend_borders(GthChain *chain, const GtRange *gen_seq_bounds,
                              const GtRange *gen_seq_bounds_rc,
                              GT_UNUSED unsigned long gen_total_length,
                              GT_UNUSED unsigned long gen_offset)
{
  long tmpborder;

  /* at least one range in chain */
  gt_assert(gt_array_size(chain->forwardranges));
  /* forward range borders are in considered genomic region */
  gt_assert(gt_ranges_borders_are_in_region(chain->forwardranges,
                                            gen_seq_bounds));
  /* reverse range borders are in considered genomic region */
  gt_assert(gt_ranges_borders_are_in_region(chain->reverseranges,
                                            gen_seq_bounds_rc));
  /* chain->forwardranges is forward and consecutive */
  gt_assert(gt_ranges_are_consecutive(chain->forwardranges));
  /* valid sequence bounds */
  gt_assert(gen_seq_bounds->start <= gen_seq_bounds->end);
  gt_assert(gen_seq_bounds_rc->start <= gen_seq_bounds_rc->end);

  /* set start border, forward strand */
  tmpborder = gt_safe_cast2long(((GtRange*)
                                 gt_array_get_first(chain->forwardranges))
                                 ->start);
  tmpborder -= DPEXTENSION;
  if (tmpborder < gt_safe_cast2long(gen_seq_bounds->start))
    tmpborder = gen_seq_bounds->start;
  ((GtRange*) gt_array_get_first(chain->forwardranges))->start =
    gt_safe_cast2ulong(tmpborder);

  /* set start border, reverse complement strand */
  tmpborder = gt_safe_cast2long(((GtRange*)
                                 gt_array_get_first(chain->reverseranges))
                                ->start);
  tmpborder -= DPEXTENSION;
  if (tmpborder < gt_safe_cast2long(gen_seq_bounds_rc->start))
    tmpborder = gen_seq_bounds_rc->start;
  ((GtRange*) gt_array_get_first(chain->reverseranges))->start =
    gt_safe_cast2ulong(tmpborder);

  /* set end border, forward strand */
  tmpborder = gt_safe_cast2long(((GtRange*)
                                gt_array_get_last(chain->forwardranges))
                                ->end);
  tmpborder += DPEXTENSION;
  if (tmpborder > gt_safe_cast2long(gen_seq_bounds->end))
    tmpborder = gen_seq_bounds->end;
  ((GtRange*) gt_array_get_last(chain->forwardranges))->end =
    gt_safe_cast2ulong(tmpborder);

  /* set end border, reverse complement strand */
  tmpborder = gt_safe_cast2long(((GtRange*)
                                gt_array_get_last(chain->reverseranges))
                                ->end);
  tmpborder += DPEXTENSION;
  if (tmpborder > gt_safe_cast2long(gen_seq_bounds_rc->end))
    tmpborder = gen_seq_bounds_rc->end;
  ((GtRange*) gt_array_get_last(chain->reverseranges))->end =
    gt_safe_cast2ulong(tmpborder);

  gt_assert(chain_is_filled_and_consistent(chain, gen_total_length,
                                           gen_offset));
}

void gth_chain_shorten_introns(GthChain *chain, unsigned long icdelta,
                               unsigned long icminremintronlength,
                               unsigned long gen_total_length,
                               unsigned long gen_offset, bool comments,
                               GtFile *outfp)
{
  GthInvertedChain inverted_chain;

  gt_assert(chain);

  /* init */
  inverted_chain_init(&inverted_chain);

  if (comments) {
    gt_file_xprintf(outfp, "%c forward DP ranges (before post processing of "
                       "potential introns):\n", COMMENTCHAR);
    gt_file_xprintf(outfp, "%c ", COMMENTCHAR);
    gt_ranges_show(chain->forwardranges, outfp);
  }

  /* chain -> inverted_chain */
  convert_chain_to_inverted_chain(&inverted_chain, chain);
  gt_assert(conversion_is_correct(chain, &inverted_chain, gen_total_length,
                                  gen_offset));

  /* post processing of potential introns */
  potentialintronspostpro(inverted_chain.forwardranges, icdelta,
                          icminremintronlength);

  /* reset chain */
  gt_array_set_size(chain->forwardranges, 0);
  gt_array_set_size(chain->reverseranges, 0);

  /* inverted_chain -> chain */
  convert_inverted_chain_to_chain(chain, &inverted_chain, gen_total_length,
                                  gen_offset);

  if (comments) {
    gt_file_xprintf(outfp,"%c forward DP ranges (after post processing of "
                       "potential introns):\n" , COMMENTCHAR);
    gt_file_xprintf(outfp, "%c ", COMMENTCHAR);
    gt_ranges_show(chain->forwardranges, outfp);
  }

  /* free space */
  inverted_chain_free(&inverted_chain);
}

static void showfragment(GtFragment *fragment, GtFile *outfp)
{
  gt_file_xprintf(outfp, "%c %lu %lu %lu %lu %ld\n", COMMENTCHAR,
                  fragment->startpos1, fragment->endpos1,
                  fragment->startpos2, fragment->endpos2, fragment->weight);
}

static unsigned long totallengthoffragments(GtChain *chain,
                                            GtFragment *fragments)
{
  GtRange currentrange, previousrange;
  unsigned long i, fragnum;
  long totallength = 0;

  previousrange.end = GT_UNDEF_ULONG;

  for (i = 0; i < gt_chain_size(chain); i++) {
    fragnum = gt_chain_get_fragnum(chain, i);
    currentrange.start  = fragments[fragnum].startpos1;
    currentrange.end = fragments[fragnum].endpos1;

    /* currentrange is forward */
    gt_assert(currentrange.start <= currentrange.end);

    totallength += currentrange.end - currentrange.start + 1;

    if (i > 0) {
      /* subtract overlaps if necessary */
      if (currentrange.start <= previousrange.end)
        totallength -= previousrange.end - currentrange.start + 1;
    }

    previousrange = currentrange;
  }

  gt_assert(totallength > 0);

  return totallength;
}

static bool globalchainislongenough(GtChain *chain, GtFragment *fragments,
                                    double *refseqcoverage,
                                    unsigned long gcfilterthreshold,
                                    unsigned long referencelength,
                                    GthStat *stat,
                                    bool comments,
                                    GtFile *outfp)
{
  unsigned long chain_length;

  chain_length = totallengthoffragments(chain, fragments);

  if (comments) {
    gt_file_xprintf(outfp, "%c chain_length=%lu\n", COMMENTCHAR,
                       chain_length);
    gt_file_xprintf(outfp, "%c referencelength=%lu\n", COMMENTCHAR,
                       referencelength);
  }

  *refseqcoverage = ((double) chain_length / (double) referencelength) * 100.0;

  gt_assert(*refseqcoverage >= 0.0 && *refseqcoverage <= 100.0);

  gth_stat_add_to_refseqcovdistri(stat, *refseqcoverage);
  if (comments) {
    gt_file_xprintf(outfp, "%c refseqcoverage=%f\n", COMMENTCHAR,
                    *refseqcoverage);
  }

  if (*refseqcoverage >= (double) gcfilterthreshold) {
    if (comments)
      gt_file_xprintf(outfp, "%c global chain long enough\n", COMMENTCHAR);
    return true;
  }

  return false;
}

#ifndef NDEBUG
static bool chain_is_colinear(GtChain *chain, GtFragment *fragments)
{
  GtFragment *firstfragment, *secondfragment;
  unsigned long i;

  if (gt_chain_size(chain) > 1) {
    for (i = 0; i < gt_chain_size(chain) - 1; i++) {
      firstfragment  = fragments + gt_chain_get_fragnum(chain, i);
      secondfragment = fragments + gt_chain_get_fragnum(chain, i+1);

      if (firstfragment->startpos1 >= secondfragment->startpos1 ||
          firstfragment->endpos1   >= secondfragment->endpos1   ||
          firstfragment->startpos2 >= secondfragment->startpos2 ||
          firstfragment->endpos2   >= secondfragment->endpos2) {
        return false;
      }
    }
  }
  return true;
}
#endif

static GtRange chain_get_genomicrange(GthChain *chain)
{
  GtRange range;
  gt_assert(chain);
  range.start = ((GtRange*) gt_array_get_first(chain->forwardranges))->start;
  range.end = ((GtRange*) gt_array_get_last(chain->forwardranges))->end;
  gt_assert(range.start <= range.end);
  return range;
}

static void enrich_chain(GthChain *chain, GtFragment *fragments,
                         unsigned long num_of_fragments, bool comments,
                         GtFile *outfp)
{
  GtRange genomicrange, fragmentrange;
  GtArray *enrichment;
  unsigned long i;
  gt_assert(chain && fragments && num_of_fragments);
  if (comments) {
    gt_file_xprintf(outfp, "%c enrich global chain with the following "
                              "forward ranges:\n",COMMENTCHAR);
    gt_file_xprintf(outfp, "%c ", COMMENTCHAR);
    gt_ranges_show(chain->forwardranges, outfp);
  }
  /* get genomic range of DP range */
  genomicrange = chain_get_genomicrange(chain);
  enrichment = gt_array_new(sizeof (GtRange));
  /* add each fragment which overlaps which DP range to the enrichment */
  for (i = 0; i < num_of_fragments; i++) {
    fragmentrange.start  = fragments[i].startpos2;
    fragmentrange.end = fragments[i].endpos2;
    if (gt_range_overlap(&genomicrange, &fragmentrange))
      gt_array_add(enrichment, fragmentrange);
  }
  gt_assert(gt_array_size(enrichment));
  /* sort the enrichment */
  qsort(gt_array_get_space(enrichment), gt_array_size(enrichment),
        sizeof (GtRange), (GtCompare) gt_range_compare);
  /* reset the current DP range array */
  gt_array_reset(chain->forwardranges);
  /* rebuild the DP range array which now includes the enrichment */
  genomicrange = *(GtRange*) gt_array_get_first(enrichment);
  gt_array_add(chain->forwardranges, genomicrange);
  for (i = 1; i < gt_array_size(enrichment); i++) {
    genomicrange = *(GtRange*) gt_array_get(enrichment, i);
    if (genomicrange.start <=
        ((GtRange*) gt_array_get_last(chain->forwardranges))->end) {
      /* overlap found -> modify last range, if necessary */
      if (((GtRange*) gt_array_get_last(chain->forwardranges))->end <
          genomicrange.end) {
        ((GtRange*) gt_array_get_last(chain->forwardranges))->end =
          genomicrange.end;
      }
    }
    else {
      /* save range */
      gt_array_add(chain->forwardranges, genomicrange);
    }
  }
  gt_array_delete(enrichment);
}

void gth_chain_contract(GthChain *dest, const GthChain *src)
{
  GtRange forwardrange, reverserange;

  gt_assert(gt_array_size(src->forwardranges) ==
            gt_array_size(src->reverseranges));

  /* copy core */
  chain_copy_core(dest, src);

  /* contract ranges */
  forwardrange.start  = ((GtRange*)
                        gt_array_get_first(src->forwardranges))->start;
  forwardrange.end = ((GtRange*)
                        gt_array_get_last(src->forwardranges))->end;
  reverserange.start  = ((GtRange*)
                        gt_array_get_first(src->reverseranges))->start;
  reverserange.end = ((GtRange*)
                        gt_array_get_last(src->reverseranges))->end;

  /* store contracted ranges */
  gt_array_add(dest->forwardranges, forwardrange);
  gt_array_add(dest->reverseranges, reverserange);
}

static GtArray* make_list_of_chain_fragments(GtChain *chain,
                                             GtFragment *fragments,
                                             unsigned long num_of_fragments,
                                             bool enrichchains,
                                             const GtRange *genomicrange)
{
  unsigned long i, fragnum;
  GtArray *chain_fragments;
  GthJTMatch match;
  gt_assert(chain && fragments && num_of_fragments);
  chain_fragments = gt_array_new(sizeof (GthJTMatch));
  if (!enrichchains) {
    /* no chain enrichment used -> store all fragments from chain */
    for (i = 0; i < gt_chain_size(chain); i++) {
      fragnum = gt_chain_get_fragnum(chain, i);
      match.gen_range.start = fragments[fragnum].startpos2;
      match.gen_range.end   = fragments[fragnum].endpos2;
      match.ref_range.start = fragments[fragnum].startpos1;
      match.ref_range.end   = fragments[fragnum].endpos1;
      gt_array_add(chain_fragments, match);
    }
  }
  else {
    GtRange fragmentrange;
    /* chain enrichment used -> store all fragments which overlap with genomic
       range of computed chain */
    for (i = 0; i < num_of_fragments; i++) {
      fragmentrange.start  = fragments[i].startpos2;
      fragmentrange.end = fragments[i].endpos2;
      if (gt_range_overlap(genomicrange, &fragmentrange)) {
        match.gen_range.start = fragments[i].startpos2;
        match.gen_range.end   = fragments[i].endpos2;
        match.ref_range.start = fragments[i].startpos1;
        match.ref_range.end   = fragments[i].endpos1;
        gt_array_add(chain_fragments, match);
      }
    }
  }
  return chain_fragments;
}

void gth_save_chain(GtChain *chain, GtFragment *fragments,
                    unsigned long num_of_fragments,
                    GT_UNUSED unsigned long max_gap_width,
                    void *data)
{
  GthSaveChainInfo *info = (GthSaveChainInfo*) data;
  GtRange range;
  GthChain *gth_chain;
  unsigned long i, fragnum;

  gt_assert(chain_is_colinear(chain, fragments));

  if (info->comments) {
    gt_file_xprintf(info->outfp, "%c process global chain with score %ld\n",
                       COMMENTCHAR, gt_chain_get_score(chain));
    gt_file_xprintf(info->outfp, "%c process global chain with the "
                       "following fragments\n", COMMENTCHAR);
    for (i = 0; i < gt_chain_size(chain); i++)
      showfragment(fragments + gt_chain_get_fragnum(chain, i), info->outfp);
  }

  /* init */
  gth_chain = gth_chain_new();
  gth_chain->gen_file_num = info->gen_file_num;
  gth_chain->gen_seq_num  = info->gen_seq_num;
  gth_chain->ref_file_num = info->ref_file_num;
  gth_chain->ref_seq_num  = info->ref_seq_num;

  /* chain has a minimum length of 1 */
  gt_assert(gt_chain_size(chain));

  /* global chain filter */
  if (globalchainislongenough(chain, fragments,
                              &gth_chain->refseqcoverage, info->gcmincoverage,
                              info->referencelength, info->stat, info->comments,
                              info->outfp)) {
    /* save all potential exons */
    for (i = 0; i < gt_chain_size(chain); i++) {
      fragnum = gt_chain_get_fragnum(chain, i);
      range.start = fragments[fragnum].startpos2;
      range.end = fragments[fragnum].endpos2;

      /* check for overlap */
      if (i > 0 &&
         range.start <=
         ((GtRange*) gt_array_get_last(gth_chain->forwardranges))->end) {
        /* overlap found -> modify last range */
        gt_assert(((GtRange*) gt_array_get_last(gth_chain->forwardranges))
                  ->end <= range.end);
        ((GtRange*) gt_array_get_last(gth_chain->forwardranges))->end =
          range.end;
      }
      else {
#ifndef NDEBUG
        if (i > 0) {
          /* gap width is smaller or equal than the maximum gap width */
          gt_assert((range.start - 1 -
                 ((GtRange*) gt_array_get_last(gth_chain->forwardranges))
                 ->end + 1 - 1) <= max_gap_width);
        }
#endif
        /* save range */
        gt_array_add(gth_chain->forwardranges, range);
      }
    }

    GtRange genomicrange = chain_get_genomicrange(gth_chain);

    if (info->enrichchains) {
      enrich_chain(gth_chain, fragments, num_of_fragments, info->comments,
                   info->outfp);
    }

    gt_assert(gt_ranges_are_consecutive(gth_chain->forwardranges));

    /* copy ranges to opposite strand */
    gt_ranges_copy_to_opposite_strand(gth_chain->reverseranges,
                                      gth_chain->forwardranges,
                                      info->gen_total_length,
                                      info->gen_offset);

    /* compute jump table if necessary */
    if (info->jump_table) {
      GthJumpTable *forward_jump_table, *reverse_jump_table;
      GtArray *chain_fragments;
      chain_fragments = make_list_of_chain_fragments(chain, fragments,
                                                     num_of_fragments,
                                                     info->enrichchains,
                                                     &genomicrange);
      forward_jump_table =
        info->jump_table_new(gt_array_get_space(chain_fragments),
                             gt_array_size(chain_fragments), info->jtdebug);
      reverse_jump_table =
        info->jump_table_new_reverse(forward_jump_table,
                                     info->gen_total_length, info->gen_offset,
                                     info->ref_total_length, info->ref_offset);
      gt_assert(!gth_chain->forward_jump_table);
      gth_chain->forward_jump_table = forward_jump_table;
      gt_assert(!gth_chain->reverse_jump_table);
      gth_chain->reverse_jump_table = reverse_jump_table;
      gt_array_delete(chain_fragments);
      gth_chain->jump_table_delete = info->jump_table_delete;
    }

    /* save array of potential exons */
    gth_chain_collection_add(info->chain_collection, gth_chain);
    if (info->comments) {
      gt_file_xprintf(info->outfp, "%c global chain with the following "
                                   "ranges has been saved\n",COMMENTCHAR);
      gt_file_xprintf(info->outfp, "%c forward ranges:\n", COMMENTCHAR);
      gt_file_xprintf(info->outfp, "%c ", COMMENTCHAR);
      gt_ranges_show(gth_chain->forwardranges, info->outfp);
      gt_file_xprintf(info->outfp, "%c reverse ranges:\n", COMMENTCHAR);
      gt_file_xprintf(info->outfp, "%c ", COMMENTCHAR);
      gt_ranges_show(gth_chain->reverseranges, info->outfp);
    }

    /* output stored chains here
       (Mohamed needed this to compare the chaining phase of gth with CHAINER)
     */
    if (info->stopafterchaining) {
      gt_file_xprintf(info->outfp,
                      "%c gl. chain with coverage=%.2f and score %ld "
                      "(genseq=%lu, str.=%c, refseq=%lu)\n", COMMENTCHAR,
                      gth_chain->refseqcoverage, gt_chain_get_score(chain),
                      gth_chain->gen_seq_num, SHOWSTRAND(info->directmatches),
                      gth_chain->ref_seq_num);

      for (i = 0; i < gt_chain_size(chain); i++)
        showfragment(fragments + gt_chain_get_fragnum(chain, i), info->outfp);
    }
  }
  else {
    /* for -paralogs this case is not supposed to occur */
    gt_assert(!info->paralogs);
    if (info->comments)
      gt_file_xprintf(info->outfp, "%c global chain discarded\n",
                         COMMENTCHAR);
    gth_chain_delete(gth_chain);
  }
}
