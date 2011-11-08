/*
  Copyright (c) 2003-2011 Gordon Gremme <gremme@zbh.uni-hamburg.de>
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

#include "core/trans_table.h"
#include "core/undef_api.h"
#include "core/unused_api.h"
#include "core/warning_api.h"
#include "gth/chaining.h"
#include "gth/gtherror.h"
#include "gth/gthxml.h"
#include "gth/intermediate.h"
#include "gth/proc_sa_collection.h"
#include "gth/similarity_filter.h"

#define UNSUCCESSFULALIGNMENTSCORE      0.0

#define SHOW_MATRIX_CALCULATION_STATUS_BUF_SIZE 160
#define ICDELTASTRINGLENGTH                     40

#define SHOW_COMPUTE_MATCHES_STATUS_BUF_SIZE    160

typedef struct {
  unsigned long call_number;
  bool significant_match_found,
       max_call_number_reached,
       stop_amino_acid_warning;
} GthMatchInfo;

static void show_matrix_calculation_status(GthShowVerbose showverbose,
                                           bool gen_strand_forward,
                                           bool ref_strand_forward,
                                           bool introncutout,
                                           unsigned long chainctr,
                                           unsigned long num_of_chains,
                                           unsigned long icdelta,
                                           unsigned long gen_file_num,
                                           unsigned long num_of_gen_files,
                                           unsigned long ref_file_num,
                                           unsigned long num_of_ref_files,
                                           bool directmatches,
                                           bool verboseseqs,
                                           const char *gen_id,
                                           const char *ref_id)
{
  char buf[SHOW_MATRIX_CALCULATION_STATUS_BUF_SIZE],
       icdeltastring[ICDELTASTRINGLENGTH];
  GT_UNUSED int rval;

  if (introncutout) {
    rval =  snprintf(icdeltastring, ICDELTASTRINGLENGTH, ", icdelta=%lu",
                     icdelta);
    /* buffer icdeltastring[ICDELTASTRINGLENGTH] is large enough */
    gt_assert(rval <  ICDELTASTRINGLENGTH);
  }

  if (num_of_gen_files == 1 && num_of_ref_files == 1) {
    rval = snprintf(buf, SHOW_MATRIX_CALCULATION_STATUS_BUF_SIZE,
                    "d=%c, compute spliced alignment, genseq=%c, "
                    "chain=%lu/%lu, refseq=%c%s", SHOWSTRAND(directmatches),
                    SHOWSTRAND(gen_strand_forward),  chainctr + 1,
                    num_of_chains, SHOWSTRAND(ref_strand_forward),
                    introncutout ? icdeltastring : "");
  }
  else {
    rval = snprintf(buf, SHOW_MATRIX_CALCULATION_STATUS_BUF_SIZE,
                    "gf=%lu/%lu, d=%c, rf=%lu/%lu, compute spliced alignment, "
                    "gs=%c, chain=%lu/%lu, rs=%c%s", gen_file_num + 1,
                    num_of_gen_files, SHOWSTRAND(directmatches),
                    ref_file_num + 1, num_of_ref_files,
                    SHOWSTRAND(gen_strand_forward), chainctr + 1,
                    num_of_chains, SHOWSTRAND(ref_strand_forward),
                    introncutout ? icdeltastring : "");
  }
  /* buf[SHOW_MATRIX_CALCULATION_STATUS_BUF_SIZE] is large enough */
  gt_assert(rval <  SHOW_MATRIX_CALCULATION_STATUS_BUF_SIZE);
  showverbose(buf);

  if (verboseseqs) {
    rval = snprintf(buf, SHOW_MATRIX_CALCULATION_STATUS_BUF_SIZE,
                    "genomicid=%s, referenceid=%s", gen_id, ref_id);
    /* buf[SHOW_MATRIX_CALCULATION_STATUS_BUF_SIZE] is large enough */
    gt_assert(rval < SHOW_MATRIX_CALCULATION_STATUS_BUF_SIZE);
    showverbose(buf);
  }
}

static int callsahmt(bool call_dna_dp,
                     GthSA *sa,
                     bool forward,
                     unsigned long gen_file_num,
                     unsigned long ref_file_num,
                     GthChain *raw_chain,
                     unsigned long gen_total_length,
                     unsigned long gen_offset,
                     const GtRange *gen_seq_bounds,
                     const GtRange *gen_seq_bounds_rc,
                     const unsigned char *ref_seq_tran,
                     const unsigned char *ref_seq_orig,
                     unsigned long ref_total_length,
                     unsigned long ref_offset,
                     GthInput *input,
                     Introncutoutinfo *introncutoutinfo,
                     GthStat *stat,
                     unsigned long chainctr,
                     unsigned long num_of_chains,
                     unsigned long translationtable,
                     bool directmatches,
                     bool proteinexonpenal,
                     GthSpliceSiteModel *splice_site_model,
                     GthDPOptionsCore *dp_options_core,
                     GthDPOptionsEST *dp_options_est,
                     GthDPOptionsPostpro *dp_options_postpro,
                     GthDNACompletePathMatrixJT dna_complete_path_matrix_jt,
                     GthProteinCompletePathMatrixJT
                     protein_complete_path_matrix_jt,
                     GthOutput *out)
{
  int rval;
  GthChain *actual_chain, *contracted_chain, *used_chain;
  unsigned long icdelta = introncutoutinfo->icinitialdelta,
                iciterations = introncutoutinfo->iciterations;
  bool useintroncutout = introncutoutinfo->introncutout;
  /* initially useintron is set to the value of introncutoutinfo->introncutout,
     if the automatic intron cutotu technique is acitvated it can be set to
     true if an matrix allocation error (ERROR_MATRIX_ALLOCATION_FAILED) occurs
   */

  gt_assert(sa);

  actual_chain = gth_chain_new();
  contracted_chain = gth_chain_new();

  for (;;) {
    /* reset actualDPrange; */
    gt_array_set_size(actual_chain->forwardranges, 0);
    gt_array_set_size(actual_chain->reverseranges, 0);

    /* copy raw chain to actual chain */
    gth_chain_copy(actual_chain, raw_chain);

    /* shorten potential introns and compute spliced sequence, if the intron
       cutout technique is used */
    if (useintroncutout) {
      /* shorten potential introns */
      gth_chain_shorten_introns(actual_chain, icdelta,
                                introncutoutinfo->icminremintronlength,
                                gen_total_length, gen_offset, out->comments,
                                out->outfp);
    }
    else
      gth_chain_contract(contracted_chain, actual_chain);

    if (out->showverbose) {
      show_matrix_calculation_status(out->showverbose, forward,
                                     gth_sa_ref_strand_forward(sa),
                                     useintroncutout, chainctr, num_of_chains,
                                     icdelta, gen_file_num,
                                     gth_input_num_of_gen_files(input),
                                     ref_file_num,
                                     gth_input_num_of_ref_files(input),
                                     directmatches, out->verboseseqs,
                                     gth_sa_gen_id(sa), gth_sa_ref_id(sa));
    }

    /* allocate space for DP parameter */
    if (out->comments) {
      gt_file_xprintf(out->outfp, "%c alloc space for DP param "
                         "(genomicid=%s, referenceid=%s)\n", COMMENTCHAR,
                         gth_sa_gen_id(sa), gth_sa_ref_id(sa));
    }
    used_chain = useintroncutout ? actual_chain : contracted_chain;

    /* The variable 'forward' denotes the genomic strand on which the DP is
       applied. */
    if (forward) {
      if (call_dna_dp) {
        rval = gth_align_dna(sa, used_chain->forwardranges,
                             gth_input_current_gen_seq_tran(input),
                             gth_input_current_gen_seq_orig(input),
                             ref_seq_tran, ref_seq_orig, ref_total_length,
                             gth_input_current_gen_alphabet(input),
                             gth_input_current_ref_alphabet(input),
                             useintroncutout,
                             introncutoutinfo->autoicmaxmatrixsize,
                             out->showeops, out->comments, out->gs2out,
                             gen_seq_bounds, splice_site_model, dp_options_core,
                             dp_options_est, dp_options_postpro,
                             dna_complete_path_matrix_jt,
                             raw_chain->forward_jump_table, ref_offset, stat,
                             out->outfp);
      }
      else { /* call_protein_dp */
        rval = gth_align_protein(sa, used_chain->forwardranges,
                                 gth_input_current_gen_seq_tran(input),
                                 ref_seq_tran, ref_seq_orig, ref_total_length,
                                 gth_input_current_gen_alphabet(input),
                                 gth_input_current_ref_alphabet(input),
                                 input, useintroncutout,
                                 introncutoutinfo->autoicmaxmatrixsize,
                                 proteinexonpenal, out->showeops, out->comments,
                                 out->gs2out, translationtable, gen_seq_bounds,
                                 splice_site_model, dp_options_core,
                                 dp_options_postpro,
                                 protein_complete_path_matrix_jt,
                                 raw_chain->forward_jump_table, ref_offset,
                                 stat, out->outfp);
      }
    }
    else {
      /* the DP is called with the revers positions specifiers */
      if (call_dna_dp) {
        rval = gth_align_dna(sa, used_chain->reverseranges,
                             gth_input_current_gen_seq_tran_rc(input),
                             gth_input_current_gen_seq_orig_rc(input),
                             ref_seq_tran, ref_seq_orig, ref_total_length,
                             gth_input_current_gen_alphabet(input),
                             gth_input_current_ref_alphabet(input),
                             useintroncutout,
                             introncutoutinfo->autoicmaxmatrixsize,
                             out->showeops, out->comments, out->gs2out,
                             gen_seq_bounds_rc, splice_site_model,
                             dp_options_core, dp_options_est,
                             dp_options_postpro, dna_complete_path_matrix_jt,
                             raw_chain->reverse_jump_table, ref_offset, stat,
                             out->outfp);
      }
      else { /* call_protein_dp */
        rval = gth_align_protein(sa, used_chain->reverseranges,
                                 gth_input_current_gen_seq_tran_rc(input),
                                 ref_seq_tran, ref_seq_orig, ref_total_length,
                                 gth_input_current_gen_alphabet(input),
                                 gth_input_current_ref_alphabet(input),
                                 input, useintroncutout,
                                 introncutoutinfo->autoicmaxmatrixsize,
                                 proteinexonpenal, out->showeops, out->comments,
                                 out->gs2out, translationtable, gen_seq_bounds,
                                 splice_site_model, dp_options_core,
                                 dp_options_postpro,
                                 protein_complete_path_matrix_jt,
                                 raw_chain->reverse_jump_table, ref_offset,
                                 stat, out->outfp);
      }
    }

    if (rval == GTH_ERROR_DP_PARAMETER_ALLOCATION_FAILED)
      return GTH_ERROR_DP_PARAMETER_ALLOCATION_FAILED;

    /* handling of special error codes ERROR_CUTOUT_NOT_IN_INTRON and
       ERROR_MATRIX_ALLOCATION_FAILED from DP
       the only possible special error code given back by this function is
       ERROR_SA_COULD_NOT_BE_DETERMINED */
#ifndef NDEBUG
    if (!useintroncutout) gt_assert(rval != GTH_ERROR_CUTOUT_NOT_IN_INTRON);
#endif
    if (useintroncutout && rval == GTH_ERROR_CUTOUT_NOT_IN_INTRON) {
      /* the intron cutout technique failed -> increase counter */
      gth_stat_increment_numofunsuccessfulintroncutoutDPs(stat);
      if (--iciterations > 0) {
        /* if an iterations is left, increase icdelta, decrease the remaining
           iterations, and continue the while-loop */
        icdelta += introncutoutinfo->icdeltaincrease;
        continue;
      }
      else {
        /* no iteration left, discard SA */
        gth_stat_increment_numofundeterminedSAs(stat);
        gth_chain_delete(actual_chain);
        gth_chain_delete(contracted_chain);
        return GTH_ERROR_SA_COULD_NOT_BE_DETERMINED;
      }
    }
    else if (rval == GTH_ERROR_MATRIX_ALLOCATION_FAILED) {
      if (introncutoutinfo->autoicmaxmatrixsize > 0 && !useintroncutout) {
        /* if the automatic intron cutout technique is enabled and a ``normal''
           DP returned with the matrix allocation error, set useintroncutout,
           increase counter, and continue */
        if (out->showverbose) {
          out->showverbose("matrix allocation failed, use intron cutout "
                           "technique");
        }
        gth_stat_increment_numofautointroncutoutcalls(stat);
        useintroncutout = true;
        continue;
      }
      else {
        /* otherwise increase relevant statistics, free space and return with
           error */
        gth_stat_increment_numoffailedmatrixallocations(stat);
        gth_stat_increment_numofundeterminedSAs(stat);
        gth_chain_delete(actual_chain);
        gth_chain_delete(contracted_chain);
        return GTH_ERROR_SA_COULD_NOT_BE_DETERMINED;
      }
    }
    else if (rval) /* ``normal'' DP */
      return -1;
    break;
  }

#if 0
  if (out->comments) {
    gt_file_xprintf(out->outfp, "%c this SA has been computed:\n", COMMENTCHAR);
    gth_sa_show(sa, input, out->outfp);
  }
#endif

  /* free */
  gth_chain_delete(actual_chain);
  gth_chain_delete(contracted_chain);

  return 0;
}

/* the following function saves <sa> by inserting it into <sa_collection> and
   sets <significantmatchfound> to true, if the insertion was successful */
static void save_sa(GthSACollection *sa_collection, GthSA *sa,
                    GthSAFilter *sa_filter, GthMatchInfo *match_info,
                    GthStat *stat)
{
  if (!gth_sa_collection_insert_sa(sa_collection, sa, sa_filter, stat)) {
    /* unsuccessful insertion; discard sa */
    gth_sa_delete(sa);
    match_info->call_number--;
  }
  else {
    /* else successful insertion */
    match_info->significant_match_found = true;
  }
}

static bool isunsuccessfulalignment(GthSA *sa,
                                    bool comments,
                                    GtFile *outfp)
{
  if (gth_sa_score(sa) <= UNSUCCESSFULALIGNMENTSCORE) {
    if (comments)
      gt_file_xprintf(outfp, "%c discard alignment\n", COMMENTCHAR);
    return true;
  }
  return false;
}

static int call_dna_DP(bool directmatches, GthCallInfo *call_info,
                       GthInput *input, GthStat *stat,
                       GthSACollection *sa_collection, GthSA *saA,
                       unsigned long gen_file_num,
                       unsigned long ref_file_num,
                       unsigned long gen_total_length,
                       unsigned long gen_offset,
                       const GtRange *gen_seq_bounds,
                       const GtRange *gen_seq_bounds_rc,
                       unsigned long ref_total_length, unsigned long ref_offset,
                       unsigned long chainctr,
                       unsigned long num_of_chains, GthMatchInfo *match_info,
                       const unsigned char *ref_seq_tran,
                       const unsigned char *ref_seq_orig,
                       const unsigned char *ref_seq_tran_rc,
                       const unsigned char *ref_seq_orig_rc,
                       GthChain *chain,
                       GthDNACompletePathMatrixJT dna_complete_path_matrix_jt,
                       GthProteinCompletePathMatrixJT
                       protein_complete_path_matrix_jt)
{
  int rval;
  bool bothstrandsanalyzed, firstdp = true,
       GT_UNUSED gs2outdirectmatches = directmatches;
  GthSA *saB = NULL;
  GtFile *outfp = call_info->out->outfp;

  if (directmatches ? gth_input_forward(input)
                    : gth_input_reverse(input)) {
    /* calculate alignment */
    rval = callsahmt(true, saA, directmatches, gen_file_num, ref_file_num,
                     chain, gen_total_length, gen_offset, gen_seq_bounds,
                     gen_seq_bounds_rc,
                     ref_seq_tran, ref_seq_orig, ref_total_length, ref_offset,
                     input, &call_info->simfilterparam.introncutoutinfo, stat,
                     chainctr, num_of_chains, call_info->translationtable,
                     directmatches, call_info->proteinexonpenal,
                     call_info->splice_site_model, call_info->dp_options_core,
                     call_info->dp_options_est, call_info->dp_options_postpro,
                     dna_complete_path_matrix_jt,
                     protein_complete_path_matrix_jt, call_info->out);
    if (rval && rval != GTH_ERROR_SA_COULD_NOT_BE_DETERMINED) {
                     /* ^ this error is treated below */
      return rval;
    }

    firstdp = false;
    bothstrandsanalyzed = gth_input_both(input);

    if (rval == GTH_ERROR_SA_COULD_NOT_BE_DETERMINED ||
        isunsuccessfulalignment(saA, call_info->out->comments, outfp)) {
      match_info->call_number--;
      /* if the spliced alignment was unsuccessful, it is deleted and the
         next hit is considered. */
      gth_sa_delete(saA);
      return 0; /* continue */
    }

    /* if not both strands are analyzed, we can save this alignment now.
       Otherwise we have to calculate the alignment to the other strand
       first and then save the better one. */
    if (!bothstrandsanalyzed)
      save_sa(sa_collection, saA, call_info->sa_filter, match_info, stat);
  }

  if (directmatches ? gth_input_reverse(input)
                    : gth_input_forward(input)) {
    if ((firstdp || gth_sa_is_poor(saA, call_info->minaveragessp)) &&
        !call_info->cdnaforwardonly) {
      if (firstdp) {
        /* space for first alignment is already allocated, bu we have to
           change the direction of the genomic and the reference strand */
        gth_sa_set_gen_strand(saA, !directmatches);
        gth_sa_set_ref_strand(saA, false);
      }
      else {
        /* allocating space for second alignment */
        saB = gth_sa_new_and_set(!directmatches, false, input,
                                 chain->gen_file_num, chain->gen_seq_num,
                                 chain->ref_file_num, chain->ref_seq_num,
                                 match_info->call_number, gen_total_length,
                                 gen_offset, ref_total_length);
      }

      /* setting gs2outdirectmatches (for compatibility) */
      gs2outdirectmatches = (bool) !directmatches;

      /* calculate alignment */
      rval = callsahmt(true, firstdp ? saA : saB, !directmatches,
                       gen_file_num, ref_file_num, chain, gen_total_length,
                       gen_offset, gen_seq_bounds, gen_seq_bounds_rc,
                       ref_seq_tran_rc, ref_seq_orig_rc, ref_total_length,
                       ref_offset, input,
                       &call_info->simfilterparam.introncutoutinfo, stat,
                       chainctr, num_of_chains, call_info->translationtable,
                       directmatches, call_info->proteinexonpenal,
                       call_info->splice_site_model, call_info->dp_options_core,
                       call_info->dp_options_est, call_info->dp_options_postpro,
                       dna_complete_path_matrix_jt,
                       protein_complete_path_matrix_jt, call_info->out);
      if (rval && rval != GTH_ERROR_SA_COULD_NOT_BE_DETERMINED) {
                       /* ^ this error is treated below */
        return rval;
      }

      if (firstdp) {
        if (rval == GTH_ERROR_SA_COULD_NOT_BE_DETERMINED ||
            isunsuccessfulalignment(saA, call_info->out->comments, outfp)) {
          /* for compatibility with GS2 */
          /* XXX: makes no sense. Possibly only if -gs2out is used. */
          match_info->significant_match_found= true;

          /* if the spliced alignment was unsuccessful, it is deleted and
             the next hit is considered. */
          gth_sa_delete(saA);
          return 0; /* continue */
        }

        save_sa(sa_collection, saA, call_info->sa_filter, match_info, stat);
      }
      else /* !firstdp */
      {
        if (rval == GTH_ERROR_SA_COULD_NOT_BE_DETERMINED ||
            isunsuccessfulalignment(saB, call_info->out->comments, outfp) ||
            !gth_sa_B_is_better_than_A(saA, saB)) {
          /* insert first SA */
          save_sa(sa_collection, saA, call_info->sa_filter, match_info, stat);
          /* discard second SA */
          gth_sa_delete(saB);
        }
        else {
          /* insert second SA */
          save_sa(sa_collection, saB, call_info->sa_filter, match_info, stat);
          /* free first SA */
          gth_sa_delete(saA);
        }
      }
    }
    else
      save_sa(sa_collection, saA, call_info->sa_filter, match_info, stat);
  }

  return 0;
}

static int call_protein_DP(bool directmatches,
                           GthCallInfo *call_info,
                           GthInput *input,
                           GthStat *stat,
                           GthSACollection *sa_collection,
                           GthSA *saA,
                           unsigned long gen_file_num,
                           unsigned long ref_file_num,
                           unsigned long gen_total_length,
                           unsigned long gen_offset,
                           const GtRange *gen_seq_bounds,
                           const GtRange *gen_seq_bounds_rc,
                           unsigned long ref_total_length,
                           unsigned long ref_offset,
                           unsigned long chainctr,
                           unsigned long num_of_chains,
                           GthMatchInfo *match_info,
                           const unsigned char *ref_seq_tran,
                           const unsigned char *ref_seq_orig,
                           GthChain *chain,
                           GthDNACompletePathMatrixJT
                           dna_complete_path_matrix_jt,
                           GthProteinCompletePathMatrixJT
                           protein_complete_path_matrix_jt)
{
  GtFile *outfp = call_info->out->outfp;
  int rval;

#ifndef NDEBUG
  /* strand is in searchmode */
  if (directmatches)
    gt_assert(gth_input_forward(input));
  else
    gt_assert(gth_input_reverse(input));
#endif

  /* calculate alignment */
  rval = callsahmt(false, saA, directmatches, gen_file_num, ref_file_num,
                   chain, gen_total_length, gen_offset, gen_seq_bounds,
                   gen_seq_bounds_rc, ref_seq_tran, ref_seq_orig,
                   ref_total_length, ref_offset, input,
                   &call_info->simfilterparam.introncutoutinfo, stat, chainctr,
                   num_of_chains, call_info->translationtable, directmatches,
                   call_info->proteinexonpenal, call_info->splice_site_model,
                   call_info->dp_options_core, call_info->dp_options_est,
                   call_info->dp_options_postpro, dna_complete_path_matrix_jt,
                   protein_complete_path_matrix_jt, call_info->out);
  if (rval && rval != GTH_ERROR_SA_COULD_NOT_BE_DETERMINED) {
                   /* ^ this error is treated below */
    return rval;
  }

  if (rval == GTH_ERROR_SA_COULD_NOT_BE_DETERMINED ||
      isunsuccessfulalignment(saA, call_info->out->comments, outfp)) {
    match_info->call_number--;
    /* if the spliced alignment was unsuccessful, it is deleted and the
       next hit is considered. */
    gth_sa_delete(saA);
    /* continue */
    return 0;
  }

  /* we can save the alignment now */
  save_sa(sa_collection, saA, call_info->sa_filter, match_info, stat);

  return 0;
}

static void show_no_match_line(GthAlphatype overallalphatype, GtFile *outfp)
{
  gt_file_xprintf(outfp, "\nNo significant ");
  switch (overallalphatype)
  {
    case DNA_ALPHA:
      gt_file_xprintf(outfp, "EST");
      break;
    case PROTEIN_ALPHA:
      gt_file_xprintf(outfp, "protein");
      break;
    default: gt_assert(0);
  }
  gt_file_xprintf(outfp, " matches were found.\n");
}

static GthChainCollection* match_and_chain(GthCallInfo *call_info,
                                           GthInput *input,
                                           GthStat *stat,
                                           unsigned long gen_file_num,
                                           unsigned long ref_file_num,
                                           bool directmatches,
                                           GthMatchInfo *match_info,
                                           const GthPlugins *plugins)
{
  GtFile *outfp = call_info->out->outfp;
  GthChainCollection *chain_collection = gth_chain_collection_new();

  /* compute the chains */
  gth_chaining(chain_collection, gen_file_num, ref_file_num, call_info, input,
               stat, directmatches, plugins);

  /* update statistics */
  gth_stat_increase_numofchains(stat,
                                gth_chain_collection_size(chain_collection));

  /* stop after chaining phase */
  if (call_info->simfilterparam.stopafterchaining) {
    gth_chain_collection_delete(chain_collection);
    return NULL;
  }

  if (call_info->out->showverbose)
    call_info->out->showverbose("calculate spliced alignment for every chain");

  if (!gth_chain_collection_size(chain_collection)) {
    /* no matches found -> return */
    if (!call_info->out->xmlout && !call_info->out->gff3out && !directmatches &&
        !match_info->significant_match_found) {
      show_no_match_line(gth_input_get_alphatype(input, ref_file_num), outfp);
    }
    gth_chain_collection_delete(chain_collection);
    return NULL;
  }

  return chain_collection;
}

static int calc_spliced_alignments(GthSACollection *sa_collection,
                                   GthChainCollection *chain_collection,
                                   GthCallInfo *call_info,
                                   GthInput *input,
                                   GthStat *stat,
                                   unsigned long gen_file_num,
                                   unsigned long ref_file_num,
                                   bool directmatches,
                                   GthMatchInfo *match_info,
                                   GthDNACompletePathMatrixJT
                                   dna_complete_path_matrix_jt,
                                   GthProteinCompletePathMatrixJT
                                   protein_complete_path_matrix_jt)
{
  const unsigned char *ref_seq_tran, *ref_seq_orig, *ref_seq_tran_rc = NULL,
                      *ref_seq_orig_rc = NULL;
  unsigned long chainctr, gen_offset = GT_UNDEF_ULONG, gen_total_length,
                ref_total_length;
  GtFile *outfp = call_info->out->outfp;
  GtRange gen_seq_bounds, gen_seq_bounds_rc;
  bool refseqisdna;
  GthChain *chain;
  GtRange range;
  GthSA *saA;
  int rval;

  gt_assert(sa_collection && chain_collection);

  refseqisdna = gth_input_ref_file_is_dna(input, ref_file_num);

  for (chainctr = 0;
       chainctr < gth_chain_collection_size(chain_collection);
       chainctr++) {
       chain = gth_chain_collection_get(chain_collection, chainctr);
    if (++match_info->call_number > call_info->firstalshown &&
        call_info->firstalshown > 0) {
      if (!(call_info->out->xmlout || call_info->out->gff3out))
        gt_file_xfputc('\n', outfp);
      else if (call_info->out->xmlout)
        gt_file_xprintf(outfp, "<!--\n");

      if (!call_info->out->gff3out) {
        gt_file_xprintf(outfp, "Maximal matching %s count (%u) reached.\n",
                        refseqisdna ? "EST" : "protein",
                        call_info->firstalshown);
        gt_file_xprintf(outfp, "Only the first %u matches will be "
                           "displayed.\n", call_info->firstalshown);
      }

      if (!(call_info->out->xmlout || call_info->out->gff3out))
        gt_file_xfputc('\n', outfp);
      else if (call_info->out->xmlout)
        gt_file_xprintf(outfp, "-->\n");

      match_info->max_call_number_reached = true;
      break; /* break out of loop */
    }

    /* compute considered genomic regions if not set by -frompos */
    if (!gth_input_use_substring_spec(input)) {
      gen_seq_bounds = gth_input_get_genomic_range(input, chain->gen_file_num,
                                                   chain->gen_seq_num);
      gen_total_length      = gt_range_length(&gen_seq_bounds);
      gen_offset            = gen_seq_bounds.start;
      gen_seq_bounds_rc     = gen_seq_bounds;
    }
    else {
      /* genomic multiseq contains exactly one sequence */
      gt_assert(gth_input_num_of_gen_seqs(input, chain->gen_file_num) == 1);
      gen_total_length = gth_input_genomic_file_total_length(input,
                                                             chain
                                                             ->gen_file_num);
      gen_seq_bounds.start    = gth_input_genomic_substring_from(input);
      gen_seq_bounds.end      = gth_input_genomic_substring_to(input);
      gen_offset              = 0;
      gen_seq_bounds_rc.start = gen_total_length - 1 - gen_seq_bounds.end;
      gen_seq_bounds_rc.end   = gen_total_length - 1 - gen_seq_bounds.start;
    }

    /* "retrieving" the reference sequence */
    range = gth_input_get_reference_range(input, chain->ref_file_num,
                                          chain->ref_seq_num);
    ref_seq_tran = gth_input_current_ref_seq_tran(input) + range.start;
    ref_seq_orig = gth_input_current_ref_seq_orig(input) + range.start;
    if (refseqisdna) {
      ref_seq_tran_rc = gth_input_current_ref_seq_tran_rc(input) + range.start;
      ref_seq_orig_rc = gth_input_current_ref_seq_orig_rc(input) + range.start;
    }
    ref_total_length = range.end - range.start + 1;

    /* check if protein sequences have a stop amino acid */
    if (!refseqisdna && !match_info->stop_amino_acid_warning &&
       ref_seq_orig[ref_total_length - 1] != GT_STOP_AMINO) {
      GtStr *ref_id = gt_str_new();
      gth_input_save_ref_id(input, ref_id, chain->ref_file_num,
                            chain->ref_seq_num);
      gt_warning("protein sequence '%s' (#%lu in file %s) does not end with a "
                 "stop amino acid ('%c'). If it is not a protein fragment you "
                 "should add a stop amino acid to improve the prediction. "
                 "For example with `gt seqtransform -addstopaminos` (see "
                 "http://genometools.org for details).", gt_str_get(ref_id),
                 chain->ref_seq_num,
                 gth_input_get_reference_filename(input, chain->ref_file_num),
                 GT_STOP_AMINO);
      match_info->stop_amino_acid_warning = true;
      gt_str_delete(ref_id);
    }

    /* allocating space for alignment */
    saA = gth_sa_new_and_set(directmatches, true, input, chain->gen_file_num,
                             chain->gen_seq_num, chain->ref_file_num,
                             chain->ref_seq_num, match_info->call_number,
                             gen_total_length, gen_offset, ref_total_length);

    /* extend the DP borders to the left and to the right */
    gth_chain_extend_borders(chain, &gen_seq_bounds, &gen_seq_bounds_rc,
                             gen_total_length, gen_offset);

    /* From here on the dp positions always refer to the forward strand of the
       genomic DNA. */

    /* call the Dynamic Programming */
    if (refseqisdna) {
      rval = call_dna_DP(directmatches, call_info, input, stat,
                         sa_collection, saA, gen_file_num, ref_file_num,
                         gen_total_length, gen_offset, &gen_seq_bounds,
                         &gen_seq_bounds_rc, ref_total_length, range.start,
                         chainctr, gth_chain_collection_size(chain_collection),
                         match_info, ref_seq_tran, ref_seq_orig,
                         ref_seq_tran_rc, ref_seq_orig_rc, chain,
                         dna_complete_path_matrix_jt,
                         protein_complete_path_matrix_jt);
    }
    else {
      rval = call_protein_DP(directmatches, call_info, input,
                             stat, sa_collection, saA, gen_file_num,
                             ref_file_num, gen_total_length, gen_offset,
                             &gen_seq_bounds, &gen_seq_bounds_rc,
                             ref_total_length, range.start, chainctr,
                             gth_chain_collection_size(chain_collection),
                             match_info, ref_seq_tran, ref_seq_orig, chain,
                             dna_complete_path_matrix_jt,
                             protein_complete_path_matrix_jt);
    }
    /* check return value */
    if (rval == GTH_ERROR_DP_PARAMETER_ALLOCATION_FAILED) {
      /* statistics bookkeeping */
      gth_stat_increment_numoffailedDPparameterallocations(stat);
      gth_stat_increment_numofundeterminedSAs(stat);
      /* free space */
      gth_sa_delete(saA);
      match_info->call_number--;
      continue; /* continue with the next DP range */
    }
    else if (rval)
      return -1;
  }

  if (!call_info->out->xmlout && !call_info->out->gff3out && !directmatches &&
      !match_info->significant_match_found &&
      match_info->call_number <= call_info->firstalshown) {
    show_no_match_line(gth_input_get_alphatype(input, ref_file_num), outfp);
  }

  return 0;
}

static void show_compute_matches_status(bool direct, GthShowVerbose showverbose,
                                        unsigned long gen_file_num,
                                        unsigned long num_of_gen_files,
                                        unsigned long ref_file_num,
                                        unsigned long num_of_ref_files)
{
  char buf[SHOW_COMPUTE_MATCHES_STATUS_BUF_SIZE];
  GT_UNUSED int rval;
  gt_assert(num_of_gen_files && num_of_ref_files );
  if (num_of_gen_files == 1 && num_of_ref_files == 1) {
    if (direct)
      showverbose("compute direct matches");
    else
      showverbose("compute palindromic matches");
  }
  else {
    if (direct) {
      rval = snprintf(buf, SHOW_MATRIX_CALCULATION_STATUS_BUF_SIZE,
                      "compute direct matches for genomic file (gf) %lu/%lu "
                      "and reference file (rf) %lu/%lu",
                      gen_file_num + 1, num_of_gen_files,
                      ref_file_num + 1, num_of_ref_files);
    }
    else {
      rval = snprintf(buf, SHOW_MATRIX_CALCULATION_STATUS_BUF_SIZE,
                      "compute palindromic matches for genomic file (gf) "
                      "%lu/%lu and reference file (rf) %lu/%lu",
                      gen_file_num + 1,  num_of_gen_files,
                      ref_file_num + 1,  num_of_ref_files);
    }
    /* buf[SHOW_COMPUTE_MATCHES_STATUS_BUF_SIZE] is large enough */
    gt_assert(rval < SHOW_MATRIX_CALCULATION_STATUS_BUF_SIZE);
    showverbose(buf);
  }
}

static int compute_sa_collection(GthSACollection *sa_collection,
                                 GthCallInfo *call_info,
                                 GthInput *input,
                                 GthStat *stat,
                                 const GthPlugins *plugins)
{
  GthChainCollection *chain_collection;
  GthMatchInfo match_info;
  unsigned long g, r;
  int rval = 0;

  match_info.call_number = 0;
  match_info.significant_match_found = false;
  match_info.max_call_number_reached = false;
  match_info.stop_amino_acid_warning = false;

  for (g = 0; g < gth_input_num_of_gen_files(input); g++) {
    for (r = 0; r < gth_input_num_of_ref_files(input); r++) {
      if (gth_input_get_alphatype(input, r) == DNA_ALPHA ||
          gth_input_forward(input)) {
        if (call_info->out->showverbose) {
          show_compute_matches_status(true, call_info->out->showverbose, g,
                                      gth_input_num_of_gen_files(input), r,
                                      gth_input_num_of_ref_files(input));
        }
        /* compute direct matches */
        chain_collection = match_and_chain(call_info, input, stat, g, r, true,
                                           &match_info, plugins);
        if (chain_collection) {
          rval = calc_spliced_alignments(sa_collection, chain_collection,
                                         call_info, input, stat, g, r, true,
                                         &match_info,
                                         plugins->dna_complete_path_matrix_jt,
                                         plugins
                                         ->protein_complete_path_matrix_jt);
          gth_chain_collection_delete(chain_collection);
          if (rval)
            break;
        }
      }
      if (match_info.max_call_number_reached)
        break;

      if (gth_input_get_alphatype(input, r) == DNA_ALPHA ||
          gth_input_reverse(input)) {
        if (call_info->out->showverbose) {
          show_compute_matches_status(false, call_info->out->showverbose, g,
                                      gth_input_num_of_gen_files(input), r,
                                      gth_input_num_of_ref_files(input));
        }
        /* compute reverse complemented (palindromic) matches */
        chain_collection = match_and_chain(call_info, input, stat, g, r, false,
                                           &match_info, plugins);
        if (chain_collection) {
          rval = calc_spliced_alignments(sa_collection, chain_collection,
                                         call_info, input, stat, g, r, false,
                                         &match_info,
                                         plugins->dna_complete_path_matrix_jt,
                                         plugins
                                         ->protein_complete_path_matrix_jt);
          gth_chain_collection_delete(chain_collection);
          if (rval)
            break;
        }
        if (match_info.max_call_number_reached)
          break;
      }
    }
  }

  return rval;
}

int gth_similarity_filter(GthCallInfo *call_info, GthInput *input,
                          GthStat *stat, unsigned int indentlevel,
                          const GthPlugins *plugins, GT_UNUSED GtError *err)
{
  GthSACollection *sa_collection; /* stores the calculated spliced alignments */

  gt_error_check(err);

  /* initialization */
  sa_collection = gth_sa_collection_new(call_info->duplicate_check);

  /* compute the spliced alignments */
  if (compute_sa_collection(sa_collection, call_info, input, stat, plugins)) {
    gth_sa_collection_delete(sa_collection);
    return -1;
  }

  /* process the alignments */
  gth_proc_sa_collection(sa_collection, call_info, input, stat, indentlevel);

  /* show XML trailer */
  if (call_info->out->xmlout) {
    gth_xml_show_trailer(call_info->intermediate, call_info->out->outfp);
  }

  /* output statistics */
  if (!call_info->out->gff3out)
    gth_stat_show(stat, true, call_info->out->xmlout, call_info->out->outfp);

#ifndef NDEBUG
  if (call_info->intermediate && call_info->out->outputfile) {
    /* intermediate output equals tree of alignments */
    gt_assert(gth_intermediate_output_is_correct(call_info->out->outputfile,
                                                 sa_collection, input,
                                                 &call_info->out->outfp, err));
  }
#endif

  /* free spliced alignment collection */
  gth_sa_collection_delete(sa_collection);

  return 0;
}
