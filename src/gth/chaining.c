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

#include "core/unused_api.h"
#include "gth/chaining.h"

#define POLYATAILFILTERALPHASIZE                4
#define SHOW_CHAIN_CALCULATION_STATUS_BUF_SIZE  160

typedef struct {
  unsigned long seqnum1,  /* number of reference sequence */
                seqnum2,  /* number of genomic sequence */
                startpos,
                length;
} Bucket;

static int compare_matches(const void *dataA, const void *dataB)
{
  GthMatch *matchA = (GthMatch*) dataA;
  GthMatch *matchB = (GthMatch*) dataB;

  /* sort after genomic sequence number */
  if (matchA->Storeseqnumgenomic < matchB->Storeseqnumgenomic)
    return -1;
  if (matchA->Storeseqnumgenomic > matchB->Storeseqnumgenomic)
    return 1;

  /* matches have the same genomic sequence number,
     sort after reference sequence number */
  if (matchA->Storeseqnumreference < matchB->Storeseqnumreference)
    return -1;
  if (matchA->Storeseqnumreference > matchB->Storeseqnumreference)
    return 1;

  /* matches have the same reference sequence number
     sort after position reference */
  if (matchA->Storepositionreference < matchB->Storepositionreference)
    return -1;
  if (matchA->Storepositionreference > matchB->Storepositionreference)
    return 1;

  /* sort after length reference */
  if (matchA->Storelengthreference < matchB->Storelengthreference)
    return -1;
  if (matchA->Storelengthreference > matchB->Storelengthreference)
    return 1;

  /* sort after position genomic */
  if (matchA->Storepositiongenomic < matchB->Storepositiongenomic)
    return -1;
  if (matchA->Storepositiongenomic > matchB->Storepositiongenomic)
    return 1;

  /* sort after length genomic */
  if (matchA->Storelengthgenomic < matchB->Storelengthgenomic)
    return -1;
  if (matchA->Storelengthgenomic > matchB->Storelengthgenomic)
    return 1;

  /* matches are equal */
  return 0;
}

#ifndef NDEBUG
static bool sum_of_bucket_lengths_equals_num_of_matches(GtArray *buckets,
                                                        unsigned long
                                                        numofmatches)
{
  unsigned long i, sumofbucketlengths = 0;
  gt_assert(buckets);
  for (i = 0; i < gt_array_size(buckets); i++)
    sumofbucketlengths += ((Bucket*) gt_array_get(buckets, i))->length;
  if (sumofbucketlengths == numofmatches)
    return true;
  return false;
}
#endif

static void sort_matches_and_calc_buckets(GtArray *matches, GtArray *buckets,
                                          unsigned long *maxbucketlength)
{
  unsigned long i, currentstart = 0, currentend = 0;
  GthMatch *matchptr;
  Bucket bucket, *bucketptr;

  gt_assert(gt_array_size(matches));

  /* sort matches */
  qsort(gt_array_get_space(matches), gt_array_size(matches), sizeof (GthMatch),
        compare_matches);

  /* init first bucket */
  matchptr = gt_array_get_first(matches);
  bucket.seqnum1  = matchptr->Storeseqnumreference;
  bucket.seqnum2  = matchptr->Storeseqnumgenomic;
  bucket.startpos = 0;

  /* calc buckets */
  for (i = 1; i < gt_array_size(matches); i++) {
    matchptr = gt_array_get(matches, i);
    if (matchptr->Storeseqnumreference != bucket.seqnum1 ||
        matchptr->Storeseqnumgenomic != bucket.seqnum2) {
      /* save the current bucket */
      currentend    = i - 1;
      bucket.length = currentend - currentstart + 1;
      gt_array_add(buckets, bucket);

      /* create new bucket */
      currentstart    = i;
      bucket.seqnum1  = matchptr->Storeseqnumreference;
      bucket.seqnum2  = matchptr->Storeseqnumgenomic;
      bucket.startpos = i;
    }
  }

  /* save last bucket */
  currentend = i - 1;
  bucket.length = currentend - currentstart + 1;
  gt_array_add(buckets, bucket);

  /* compute maximum bucket length */
  *maxbucketlength = 0;
  for (i = 0; i < gt_array_size(buckets); i++) {
    bucketptr = gt_array_get(buckets, i);
    if (bucketptr->length > *maxbucketlength)
      *maxbucketlength = bucketptr->length;
  }

  gt_assert(sum_of_bucket_lengths_equals_num_of_matches(buckets,
                                                     gt_array_size(matches)));
}

static void outputbuckets(GtArray *buckets, GthMatch *sortedmatches,
                          GtFile *outfp)
{
  unsigned long i, startpos;
  Bucket *bucket;
  for (i = 0; i < gt_array_size(buckets); i++) {
    bucket = gt_array_get(buckets, i);
    gt_file_xprintf(outfp, "%c bucket %lu: seqnum1=%lu, startpos=%lu, "
                    "length=%lu\n", COMMENTCHAR, i, bucket->seqnum1,
                    bucket->startpos, bucket->length);
    gt_file_xprintf(outfp, "%c first match in bucket %lu:\n", COMMENTCHAR, i);
    startpos = bucket->startpos;
    gt_file_xprintf(outfp, "%c seqnum1=%lu, seqnum2=%lu, p1=%lu, l1=%lu, "
                    "p2=%lu, l2=%lu\n", COMMENTCHAR,
                    sortedmatches[startpos].Storeseqnumreference,
                    sortedmatches[startpos].Storeseqnumgenomic,
                    sortedmatches[startpos].Storepositionreference,
                    sortedmatches[startpos].Storelengthreference,
                    sortedmatches[startpos].Storepositiongenomic,
                    sortedmatches[startpos].Storelengthgenomic);
  }
}

/* The following function transforms the reference sequence positions to
   the opposite strand. This is necessary for proper chaining (only if
   palindromic matches have been calculated). */
static void transform_refseq_positions(GtArray *matches,
                                       GthSeqCon *ref_seq_con)
{
  unsigned long i, referencelength, referenceoffset;
  GtRange original, transformed;
  GthMatch *match;
  GtRange range;

  for (i = 0; i < gt_array_size(matches); i++) {
    match = gt_array_get(matches, i);
    /* get necessary data for transformation */
    range = gth_seq_con_get_range(ref_seq_con, match->Storeseqnumreference);
    referencelength = range.end - range.start + 1;
    referenceoffset = range.start;

    /* store original match range */
    original.start  = match->Storepositionreference;
    original.end = original.start + match->Storelengthreference - 1;
    gt_assert(original.end > original.start);

    /* transform match range */
    transformed.start = referencelength - 1
                        - (original.end - referenceoffset)
                        + referenceoffset;
    transformed.end = referencelength - 1
                      - (original.start - referenceoffset)
                      + referenceoffset;
    gt_assert(transformed.end > transformed.start);

    /* store transformed match range */
    match->Storepositionreference = transformed.start;
    match->Storelengthreference   = transformed.end - transformed.start + 1;
  }
}

/*
  The following function is a modified version of the function
  vmatchinitfragmentinfo from the file chainvm.c.
  A difference is that GthMatches are processed instead of StoreMatches
  and GTHDISTANCE2STORE is used instead of DISTANCE2SCORE.

  More importantly, consequtive matches which are the same are stored only once.
*/

static void gthinitfragments(GtFragment *fragments,
                             unsigned long *num_of_fragments,
                             GthMatch *storematchtab,
                             unsigned long numofmatches,
                             unsigned long rare,
                             double fragweightfactor)
{
  GthMatch *mptr;
  GtFragment *fragmentptr;
  long tmp, largestdim1 = 0, largestdim2 = 0;
  GtDiscDistri *startpointdistri = NULL;

  /* init number of fragments */
  *num_of_fragments = 0;
  if (rare)
    startpointdistri = gt_disc_distri_new();

  for (mptr = storematchtab; mptr < storematchtab + numofmatches; mptr++) {
    /* first dimension */
    tmp = mptr->Storepositionreference + mptr->Storelengthreference - 1;
    if (largestdim1 < tmp)
      largestdim1 = tmp;

    /* second dimension */
    tmp = mptr->Storepositiongenomic + mptr->Storelengthgenomic - 1;
    if (largestdim2 < tmp)
      largestdim2 = tmp;
  }

  for (mptr = storematchtab, fragmentptr = fragments;
       mptr < storematchtab + numofmatches;
       mptr++) {
    if (rare)
      gt_disc_distri_add(startpointdistri, mptr->Storepositionreference);
    if ((!rare ||
         gt_disc_distri_get(startpointdistri, mptr->Storepositionreference)
         <= rare) &&
        (mptr == storematchtab || /* is the first match */
         !gth_matches_are_equal(mptr, mptr-1))) { /* or is different from last
                                                     one */
      fragmentptr->weight     = (long) (fragweightfactor *
                                        (double) abs(mptr->Storescore));
      fragmentptr->startpos1  = mptr->Storepositionreference;
      fragmentptr->endpos1    = mptr->Storepositionreference
                                + mptr->Storelengthreference - 1;
      fragmentptr->startpos2  = mptr->Storepositiongenomic;
      fragmentptr->endpos2    = mptr->Storepositiongenomic
                                + mptr->Storelengthgenomic - 1;
      fragmentptr++;
      (*num_of_fragments)++;
    }
  }

  gt_disc_distri_delete(startpointdistri);

  gt_assert(*num_of_fragments <= numofmatches);
}

static void show_chain_calc_status(GthShowVerbose showverbose,
                                   unsigned long chainnum,
                                   unsigned long numofchains,
                                   unsigned long numofmatches,
                                   unsigned long currentgen_file_num,
                                   unsigned long numofgenomicfiles,
                                   unsigned long currentreffilenum,
                                   unsigned long numofreffiles,
                                   bool directmatches, bool verboseseqs,
                                   unsigned long genseqnum,
                                   unsigned long refseqnum)
{
  char buf[SHOW_CHAIN_CALCULATION_STATUS_BUF_SIZE];
  GT_UNUSED int rval;

  gt_assert(numofchains > 0);

  if (numofgenomicfiles == 1 && numofreffiles == 1) {
    rval = snprintf(buf, SHOW_CHAIN_CALCULATION_STATUS_BUF_SIZE,
                    "d=%c, compute chains for bucket %lu/%lu (matches in "
                    "bucket=%lu)", SHOWSTRAND(directmatches), chainnum,
                    numofchains, numofmatches);
  }
  else {
    rval = snprintf(buf, SHOW_CHAIN_CALCULATION_STATUS_BUF_SIZE,
                    "gf=%lu/%lu, d=%c, rf=%lu/%lu, compute chains for bucket "
                    "%lu/%lu (matches in bucket=%lu)",
                    currentgen_file_num + 1,
                    numofgenomicfiles, SHOWSTRAND(directmatches),
                    currentreffilenum + 1, numofreffiles, chainnum, numofchains,
                    numofmatches);
  }
  /* buf[SHOW_CHAIN_CALCULATION_STATUS_BUF_SIZE] is large enough */
  gt_assert(rval < SHOW_CHAIN_CALCULATION_STATUS_BUF_SIZE);
  showverbose(buf);

  if (verboseseqs) {
    rval = snprintf(buf, SHOW_CHAIN_CALCULATION_STATUS_BUF_SIZE,
                    "genseqnum=%lu, refseqnum=%lu", genseqnum, refseqnum);
    /* buf[SHOW_CHAIN_CALCULATION_STATUS_BUF_SIZE] is large enough */
    gt_assert(rval < SHOW_CHAIN_CALCULATION_STATUS_BUF_SIZE);
    showverbose(buf);
  }
}

static void chaining_info_init(GthChainingInfo *chaining_info,
                               bool directmatches,
                               bool refseqisdna,
                               GthCallInfo *call_info,
                               GthInput *input,
                               GthStat *stat,
                               unsigned long gen_file_num,
                               unsigned long ref_file_num)
{
  chaining_info->directmatches    = directmatches;
  chaining_info->refseqisindex    = call_info->simfilterparam.inverse ||
                                    !refseqisdna;
  chaining_info->jtdebug          = call_info->dp_options_core->jtdebug;
  chaining_info->call_info        = call_info;
  chaining_info->input            = input;
  chaining_info->stat             = stat;
  chaining_info->gen_file_num     = gen_file_num;
  chaining_info->ref_file_num     = ref_file_num;
  chaining_info->bucketnum        = 0;
  chaining_info->maxbucketnum     = ~0;
}

static void showmatches(GthMatch *matches, unsigned long numofmatches,
                        GtFile *outfp)
{
  unsigned long i;
  for (i = 0; i < numofmatches; i++) {
    gt_file_xprintf(outfp, "%c refseqnum=%lu, genseqnum=%lu, refpos=%lu, "
                           "reflen=%lu, genpos=%lu, genlen=%lu\n", COMMENTCHAR,
                    matches[i].Storeseqnumreference,
                    matches[i].Storeseqnumgenomic,
                    matches[i].Storepositionreference,
                    matches[i].Storelengthreference,
                    matches[i].Storepositiongenomic,
                    matches[i].Storelengthgenomic);
  }
}

static void calc_chains_from_matches(GthChainCollection *chain_collection,
                                     GtArray *matches,
                                     GthChainingInfo *chaining_info,
                                     GthSeqCon *gen_seq_con,
                                     GthSeqCon *ref_seq_con,
                                     unsigned long rare,
                                     double fragweightfactor,
                                     GthJumpTableNew jump_table_new,
                                     GthJumpTableNewReverse
                                     jump_table_new_reverse,
                                     GthJumpTableDelete jump_table_delete)
{
  unsigned long i, numofchains = 0, num_of_fragments, maxbucketlength = 0;
  GtRange range;
  GtFile *outfp = chaining_info->call_info->out->outfp;
  GtFragment *fragments;
  GthSaveChainInfo info;
  GtArray *buckets;
  Bucket *bucket;

  /* this is a random sample to check that no equal matches exist
     either one match to chain or if more than one the first two differ */
  gt_assert(gt_array_size(matches) == 1 ||
            (gt_array_size(matches) > 1 &&
             !gth_matches_are_equal(gt_array_get(matches, 0),
                                    gt_array_get(matches, 1))));

  /* init */
  buckets = gt_array_new(sizeof (Bucket));

  /* output unsorted matches */
  if (chaining_info->call_info->out->comments) {
    gt_file_xprintf(outfp, "%c output unsorted matches\n", COMMENTCHAR);
    showmatches(gt_array_get_space(matches), gt_array_size(matches), outfp);
  }

  /* transform reference sequence positions to opposite strand if necessary */
  if (!chaining_info->directmatches) {
    if (chaining_info->call_info->out->comments) {
      gt_file_xprintf(outfp, "%c\n", COMMENTCHAR);
      gt_file_xprintf(outfp, "%c transform reference sequence positions to "
                                "opposite strand\n", COMMENTCHAR);
      gt_file_xprintf(outfp, "%c\n", COMMENTCHAR);
    }

    transform_refseq_positions(matches, ref_seq_con);

    /* output transformed matches */
    if (chaining_info->call_info->out->comments) {
      gt_file_xprintf(outfp, "%c output transformed matches\n", COMMENTCHAR);
      showmatches(gt_array_get_space(matches), gt_array_size(matches), outfp);
    }
  }

  /* sort matches */
  sort_matches_and_calc_buckets(matches, buckets, &maxbucketlength);

  /* output sorted matches */
  if (chaining_info->call_info->out->comments) {
    gt_file_xprintf(outfp, "%c output sorted matches\n", COMMENTCHAR);
    showmatches(gt_array_get_space(matches), gt_array_size(matches), outfp);
  }

  /* output buckets */
  if (chaining_info->call_info->out->comments) {
    gt_file_xprintf(outfp, "%c output buckets\n", COMMENTCHAR);
    outputbuckets(buckets, gt_array_get_space(matches), outfp);
  }

  /* alloc space for fragments */
  fragments = gt_malloc(sizeof (GtFragment) * maxbucketlength);

  /* save data to process the chains with saveChainasDPrange; constant part */
  info.chain_collection  = chain_collection;
  info.gcmincoverage     = chaining_info->call_info->gcmincoverage;
  info.stat              = chaining_info->stat;
  info.comments          = chaining_info->call_info->out->comments;
  info.stopafterchaining = chaining_info->call_info->simfilterparam
                           .stopafterchaining;
  info.paralogs          = chaining_info->call_info->simfilterparam.paralogs;
  info.enrichchains      = chaining_info->call_info->simfilterparam
                           .enrichchains;
  info.jump_table        = chaining_info->call_info->simfilterparam.jump_table;
  info.jump_table_new    = jump_table_new;
  info.jump_table_new_reverse = jump_table_new_reverse;
  info.jump_table_delete = jump_table_delete;
  info.jtdebug           = chaining_info->jtdebug;
  info.directmatches     = chaining_info->directmatches;
  info.outfp             = outfp;
  info.gen_file_num      = chaining_info->gen_file_num;
  info.ref_file_num      = chaining_info->ref_file_num;

  /* for every bucket a chain and for every chain a DP call (later maybe more
     than one chain) */
  for (i = 0; i < gt_array_size(buckets); i++) {
    bucket = gt_array_get(buckets, i);
    if (chaining_info->call_info->out->showverbose) {
      if (chaining_info->refseqisindex &&
          !chaining_info->call_info->simfilterparam.online) {
        /* in this case the exact number of chains is known */
        numofchains = gt_array_size(buckets);
      }
      else {
        /* this expression gives an upper bound on the number of chains
           (because we do not know the exact number here) */
        numofchains = chaining_info->bucketnum +
                      gth_seq_con_num_of_seqs(gen_seq_con) *
                      (gth_seq_con_num_of_seqs(ref_seq_con) - bucket->seqnum1);

        if (numofchains > chaining_info->maxbucketnum)
          numofchains = chaining_info->maxbucketnum;
        else
          chaining_info->maxbucketnum = numofchains;
      }
    }

    /* compute a set of fragments from every bucket of matches */

    gthinitfragments(fragments, &num_of_fragments,
                     (GthMatch*) gt_array_get_space(matches) + bucket->startpos,
                     bucket->length, rare, fragweightfactor);

    if (chaining_info->call_info->out->showverbose) {
      show_chain_calc_status (chaining_info->call_info->out->showverbose,
                              ++chaining_info->bucketnum, numofchains,
                              num_of_fragments, chaining_info->gen_file_num,
                              gth_input_num_of_gen_files(chaining_info->input),
                              chaining_info->ref_file_num,
                              gth_input_num_of_ref_files(chaining_info->input),
                              chaining_info->directmatches,
                              chaining_info->call_info->out->verboseseqs,
                              bucket->seqnum2, bucket->seqnum1);
    }

    info.gen_seq_num = ((GthMatch*) gt_array_get(matches, bucket->startpos))
                       ->Storeseqnumgenomic;
    info.ref_seq_num = ((GthMatch*) gt_array_get(matches, bucket->startpos))
                       ->Storeseqnumreference;

    /* store genomic offset */
    range = gth_seq_con_get_range(gen_seq_con, info.gen_seq_num);
    info.gen_total_length = range.end - range.start + 1;
    info.gen_offset       = range.start;

    /* store length of reference sequence */
    range = gth_seq_con_get_range(ref_seq_con, info.ref_seq_num);
    info.ref_total_length = range.end - range.start + 1;
    info.ref_offset       = range.start;
    info.referencelength  = range.end - range.start + 1;

    /* set number of remaining buckets */
    info.numofremainingbuckets = gt_array_size(buckets) - i;

    if (chaining_info->call_info->simfilterparam.paralogs) {
      gt_globalchaining_coverage(fragments, num_of_fragments,
                                 chaining_info->call_info->gcmaxgapwidth,
                                 info.referencelength,
                                 ((double)
                                  chaining_info->call_info->gcmincoverage) /
                                  100.0, gth_save_chain, &info);
    }
    else {
      gt_globalchaining_max(fragments, num_of_fragments,
                            chaining_info->call_info->gcmaxgapwidth,
                            gth_save_chain, &info);
    }
  }

  /* free space */
  gt_array_delete(buckets);
  gt_free(fragments);
}

static void match_processor_info_init(GthMatchProcessorInfo
                                      *match_processor_info,
                                      GtArray *matches,
                                      GthChainCollection *chain_collection,
                                      bool directmatches, bool refseqisdna,
                                      bool online, bool inverse, GthStat *stat,
                                      GthChainingInfo *chaining_info,
                                      unsigned long maxnumofmatches,
                                      unsigned long rare,
                                      double fragweightfactor,
                                      GthJumpTableNew jump_table_new,
                                      GthJumpTableNewReverse
                                      jump_table_new_reverse,
                                      GthJumpTableDelete jump_table_delete)
{
  match_processor_info->chain_collection       = chain_collection;
  match_processor_info->matches                = matches;
  match_processor_info->directmatches          = directmatches;
  match_processor_info->refseqisdna            = refseqisdna;
  match_processor_info->online                 = online;
  match_processor_info->refseqisindex          = inverse || !refseqisdna;
  match_processor_info->stat                   = stat;
  match_processor_info->chaining_info          = chaining_info;
  match_processor_info->matchnumcounter        = NULL;
  match_processor_info->maxnumofmatches        = maxnumofmatches;
  match_processor_info->rare                   = rare;
  match_processor_info->lastrefseqnum          = 0;
  match_processor_info->fragweightfactor       = fragweightfactor;
  match_processor_info->gen_seq_con            = NULL;
  match_processor_info->ref_seq_con            = NULL;
  match_processor_info->jump_table_new         = jump_table_new;
  match_processor_info->jump_table_new_reverse = jump_table_new_reverse;
  match_processor_info->jump_table_delete      = jump_table_delete;
}

/*
  The following function expects <info> to point to an array of
  GthMatch-structures. It stores a match in the next free position of this
  array.
*/

int gth_match_processor(GthMatchProcessorInfo *info, GthSeqCon *gen_seq_con,
                        GthSeqCon *ref_seq_con, GthMatch *match)
{
  if (info->matchnumcounter) {
    info->matchnumcounter[match->Storeseqnumreference]++;

    if (info->maxnumofmatches > 0 &&
        info->matchnumcounter[match->Storeseqnumreference] >
        info->maxnumofmatches) {
      /* discard matchA */
      return 0;
    }
  }

  if (!(info->refseqisindex && !info->online) &&
      match->Storeseqnumreference != info->lastrefseqnum &&
      gt_array_size(info->matches)) {
    gt_assert(info->chain_collection && info->chaining_info);

    /* chain all current matches */
    calc_chains_from_matches(info->chain_collection, info->matches,
                             info->chaining_info, gen_seq_con, ref_seq_con,
                             info->rare, info->fragweightfactor,
                             info->jump_table_new, info->jump_table_new_reverse,
                             info->jump_table_delete);

    /* and remove them afterwards */
    gt_array_reset(info->matches);
  }

  /*...only if it does not equal the last one */
  if (gt_array_size(info->matches) &&
      gth_matches_are_equal(gt_array_get_last(info->matches), match)) {
    return 0;
  }
  gt_array_add_elem(info->matches, match, sizeof *match);

  /* update last reference sequence number */
  info->lastrefseqnum = match->Storeseqnumreference;

  return 0;
}

void gth_chaining(GthChainCollection *chain_collection,
                  unsigned long gen_file_num,
                  unsigned long ref_file_num,
                  GthCallInfo *call_info,
                  GthInput *input,
                  GthStat *stat,
                  bool directmatches,
                  const GthPlugins *plugins)
{
  unsigned long i, numofsequences = 0;
  GtArray *matches;
  GthChainingInfo chaining_info;
  void *matcher_arguments;
  GtFile *outfp = call_info->out->outfp;
  GthMatchProcessorInfo match_processor_info;
  bool refseqisdna = gth_input_ref_file_is_dna(input, ref_file_num);

  /* make sure matcher is defined */
  gt_assert(plugins);
  gt_assert(plugins->matcher_arguments_new);
  gt_assert(plugins->matcher_arguments_delete);
  gt_assert(plugins->matcher_runner);

  /* init */
  matches = gt_array_new(sizeof (GthMatch));

  chaining_info_init(&chaining_info, directmatches, refseqisdna, call_info,
                     input, stat, gen_file_num, ref_file_num);

  matcher_arguments =
    plugins->matcher_arguments_new(true,
                          input,
                          call_info->simfilterparam.inverse || !refseqisdna
                          ? gth_input_get_genomic_filename(input, gen_file_num)
                          : gth_input_get_reference_filename(input,
                                                             ref_file_num),
                          call_info->simfilterparam.inverse || !refseqisdna
                          ? gth_input_get_reference_filename(input,
                                                             ref_file_num)
                          : gth_input_get_genomic_filename(input, gen_file_num),
                          directmatches,
                          refseqisdna,
                          call_info->progname,
                          gt_str_get(gth_input_proteinsmap(input)),
                          call_info->simfilterparam.exact,
                          call_info->simfilterparam.edist,
                          false,
                          0,
                          call_info->simfilterparam.minmatchlength,
                          call_info->simfilterparam.seedlength,
                          call_info->simfilterparam.exdrop,
                          call_info->simfilterparam.prminmatchlen,
                          call_info->simfilterparam.prseedlength,
                          call_info->simfilterparam.prhdist,
                          call_info->translationtable,
                          call_info->simfilterparam.online,
                          call_info->simfilterparam.noautoindex,
                          call_info->simfilterparam.maskpolyAtails,
                          false);

  match_processor_info_init(&match_processor_info, matches, chain_collection,
                            directmatches, refseqisdna,
                            call_info->simfilterparam.online,
                            call_info->simfilterparam.inverse, stat,
                            &chaining_info,
                            call_info->simfilterparam.maxnumofmatches,
                            call_info->simfilterparam.rare,
                            call_info->fragweightfactor,
                            plugins->jump_table_new,
                            plugins->jump_table_new_reverse,
                            plugins->jump_table_delete);

  if (call_info->simfilterparam.maxnumofmatches > 0 ||
      gth_stat_get_matchnumdistri(stat)) {
    /* alloc space of match number counter */
    numofsequences = gth_input_num_of_ref_seqs(input, ref_file_num);
    match_processor_info.matchnumcounter = gt_malloc(sizeof (unsigned long) *
                                                     numofsequences);

    /* init match number counter to 0 */
    memset(match_processor_info.matchnumcounter, 0,
           (size_t) numofsequences * sizeof (unsigned long));
  }

  /* free input, which contains the virtual trees.
     because vmatch loads the virtual trees into memory, too.
     this prevents that the virtual trees are loaded twice. */
  gth_input_delete_current(input);

  /* call matcher */
  if (call_info->out->showverbose)
    call_info->out->showverbose("call vmatch to compute matches");

  plugins->matcher_runner(matcher_arguments, call_info->out->showverbose,
                          call_info->out->showverboseVM, &match_processor_info);

  /* free matcher stuff here, because otherwise the reference file is mapped
     twice below */
  plugins->matcher_arguments_delete(matcher_arguments);

  /* free sequence collections (if they have been filled by the matcher) */
  gth_seq_con_delete(match_processor_info.gen_seq_con);
  gth_seq_con_delete(match_processor_info.ref_seq_con);

  /* save match numbers of match number distribution, if necessary */
  if (gth_stat_get_matchnumdistri(stat)) {
    for (i = 0; i < numofsequences; i++) {
      if (match_processor_info.matchnumcounter[i] > 0) {
        gth_stat_add_to_matchnumdistri(stat,
                                      match_processor_info.matchnumcounter[i]);
      }
    }
  }

  /* free match number counter */
  gt_free(match_processor_info.matchnumcounter);

  /* return if no match has been found */
  if (!gt_array_size(matches)) {
    if (call_info->out->comments)
      gt_file_xprintf(outfp, "%c no match has been found\n", COMMENTCHAR);
    gt_array_delete(matches);
    return;
  }

  /* load genomic file back into memory */
  gth_input_load_genomic_file(input, gen_file_num, true);

  /* load reference file back into memory */
  gth_input_load_reference_file(input, ref_file_num, true);

  /* compute chains from matches */
  calc_chains_from_matches(chain_collection, matches, &chaining_info,
                           gth_input_current_gen_seq_con(input),
                           gth_input_current_ref_seq_con(input),
                           call_info->simfilterparam.rare,
                           call_info->fragweightfactor,
                           plugins->jump_table_new,
                           plugins->jump_table_new_reverse,
                           plugins->jump_table_delete);

  if (call_info->out->showverbose) {
    call_info->out->showverbose("sort global chains according to reference "
                                "sequence coverage");
  }

  /* sort chains */
  gth_chain_collection_sort(chain_collection);

  /* free */
  gt_array_delete(matches);
}
