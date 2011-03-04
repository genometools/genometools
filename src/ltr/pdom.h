/*
  Copyright (c) 2008-2010 Sascha Steinbiss <steinbiss@zbh.uni-hamburg.de>
  Copyright (c) 2008-2010 Center for Bioinformatics, University of Hamburg

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

#ifndef PDOM_H
#define PDOM_H

#ifdef HAVE_HMMER

#include "core/array.h"
#include "core/error.h"
#include "core/hashmap.h"
#include "core/strand.h"
#include "core/str_array.h"
#include "ltr/ltrelement.h"

typedef enum {
  GT_PHMM_CUTOFF_TC,
  GT_PHMM_CUTOFF_GA,
  GT_PHMM_CUTOFF_NONE
} GtPdomCutoff;

typedef struct GtPdomOptions {
  double evalue_cutoff;
  GtStrArray *hmm_files;
  unsigned int chain_max_gap_length;
  bool write_alignments,
       write_aaseqs,
       output_all_chains;
  GtPdomCutoff cutoff;
} GtPdomOptions;

typedef struct GtPdomFinder GtPdomFinder;
typedef struct GtPdomModel GtPdomModel;
typedef struct GtPdomSingleHit GtPdomSingleHit;
typedef struct GtPdomModelHit GtPdomModelHit;
typedef struct GtPdomResults GtPdomResults;

typedef int (*GtPdomIteratorFunc)(GtPdomModel *model, GtPdomModelHit *hit,
                                  void *data, GtError*);

/* acts as a factory for GtPdomResults */
GtPdomFinder*    gt_pdom_finder_new(GtStrArray *hmmfiles,
                                    double eval_cutoff,
                                    unsigned int chain_max_gap_length,
                                    GtPdomCutoff cutoff,
                                    GtError*);
unsigned int     gt_pdom_finder_get_nof_threads(const GtPdomFinder*);
unsigned int     gt_pdom_finder_get_max_gap_length(const GtPdomFinder*);
double           gt_pdom_finder_get_eval_cutoff(const GtPdomFinder*);
GtStrArray*      gt_pdom_finder_get_model_filenames(const GtPdomFinder*);
unsigned long    gt_pdom_finder_get_number_of_models(const GtPdomFinder*);
GtPdomModel*     gt_pdom_finder_get_model(const GtPdomModel*, unsigned long);
GtPdomResults*   gt_pdom_finder_find(GtPdomFinder*, const char *seq,
                                     const char *rev_seq, GtLTRElement*,
                                     GtError *err);
void             gt_pdom_finder_delete(GtPdomFinder*);

/* holds results for all models */
int              gt_pdom_results_foreach_domain_hit(GtPdomResults*,
                                                    GtPdomIteratorFunc,
                                                    void*,
                                                    GtError*);
bool             gt_pdom_results_empty(GtPdomResults*);
double           gt_pdom_results_get_combined_evalue_fwd(GtPdomResults*);
double           gt_pdom_results_get_combined_evalue_rev(GtPdomResults*);
void             gt_pdom_results_delete(GtPdomResults*);

/* holds information about a single domain model */
const char*      gt_pdom_model_get_name(const GtPdomModel*);
const char*      gt_pdom_model_get_acc(const GtPdomModel*);

/* holds hits for a single domain model */
unsigned long    gt_pdom_model_hit_num_of_single_hits(const GtPdomModelHit*);
GtPdomSingleHit* gt_pdom_model_hit_single_hit(const GtPdomModelHit*,
                                                   unsigned long i);
GtStrand         gt_pdom_model_hit_get_best_strand(const GtPdomModelHit*);

/* holds information about an individual hit */
GtPhase          gt_pdom_single_hit_get_phase(const GtPdomSingleHit*);
GtRange          gt_pdom_single_hit_get_range(const GtPdomSingleHit*);
double           gt_pdom_single_hit_get_evalue(const GtPdomSingleHit*);
void             gt_pdom_single_hit_format_alignment(const GtPdomSingleHit*,
                                                     unsigned long width,
                                                     GtStr *dest);
void             gt_pdom_single_hit_get_aaseq(const GtPdomSingleHit*,
                                              GtStr *dest);
bool             gt_pdom_single_hit_is_chained(GtPdomSingleHit *singlehit);
GtArray*         gt_pdom_single_hit_get_chains(GtPdomSingleHit *singlehit);

#endif
#endif
