/*
  Copyright (c) 2015 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
  Copyright (c) 2015 Center for Bioinformatics, University of Hamburg

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
#ifndef SEED_EXTEND_H
#define SEED_EXTEND_H
#include "core/types_api.h"
#include "core/unused_api.h"
#include "querymatch.h"
#include "xdrop.h"
#include "ft-front-prune.h"

/* This header file describes the interface to two different
   methods for extending seeds, namely the xdrop-based method based on

   @ARTICLE{ZHA:SCHWA:WAG:MIL:2000,
   author = {Zhang, Z. and Schwartz, S. and Wagner, L. and Miller, W.},
   title = {{A Greedy Algorithm for Aligning DNA Sequences}},
   journal = JCB,
   year = {2000},
   volume = {{\PrintVol{7}}},
   pages = {{203-214}},
   number = {{1/2}}
   }

   and the greedy method of

   @inproceedings{MYE:2014,
   author    = {Gene Myers},
   title     = {Efficient Local Alignment Discovery amongst Noisy Long Reads},
   booktitle = {Algorithms in Bioinformatics - 14th International Workshop,
                {WABI} 2014, Wroclaw, Poland, September 8-10, 2014.
                Proceedings},
   year      = {2014},
   pages     = {52--67},
   url       = {http://dx.doi.org/10.1007/978-3-662-44753-6_5},
   doi       = {10.1007/978-3-662-44753-6_5},
   timestamp = {Tue, 30 Sep 2014 11:35:42 +0200},
   biburl    = {http://dblp.uni-trier.de/rec/bib/conf/wabi/Myers14},
   bibsource = {dblp computer science bibliography, http://dblp.org}
  }
*/

#define GT_DEFAULT_MATCHSCORE_BIAS 1.0  /* has no effect */

/* This is the minimum percentage value for extended seeds. */
#define GT_EXTEND_MIN_IDENTITY_PERCENTAGE 70

/* This is the type storing the relevant information for
   the xdrop-based seed extension method. */

typedef struct GtXdropmatchinfo GtXdropmatchinfo;

/* The constructor, which is called once before the first seed
   is to be extended. The parameter <userdefinedleastlength> is the minimum
   length of the the extension to both sides (including the seed itself).
   <errorpercentage> is the percentage of errors allowed in the
   extended seeds. <xdropbelowscore> is the parameter which influences the
   search space of the Xdrop-based extension. The larger this parameter,
   the larger the search space. If <xdropbelowscore> is 0,
   then a reasonable default value depending on the <errorpercentage>
   is chosen. */

typedef struct
{
  void *processinfo;
  GtQuerymatch *querymatchspaceptr;
} GtProcessinfo_and_querymatchspaceptr;;

GtXdropmatchinfo *gt_xdrop_matchinfo_new(GtUword userdefinedleastlength,
                                         GtUword errorpercentage,
                                         GtXdropscore xdropbelowscore,
                                         GtUword sensitivity);

/* reset seqabstract objects for a new run */

void gt_xdrop_matchinfo_reset_seqabstract(GtXdropmatchinfo *xdropmatchinfo);

/* The destructor-method. */

void gt_xdrop_matchinfo_delete(GtXdropmatchinfo *xdropmatchinfo);

/* The following function returns the optimal xdrop score depending
   on the error percentage. */

GtWord gt_optimalxdropbelowscore(GtUword errorpercentage,GtUword sensitivity);

/* Set the silent flag in the matchinfo object. */

void gt_xdrop_matchinfo_silent_set(GtXdropmatchinfo *xdropmatchinfo);

/* The following function is used for extending a seed obtained
   in a self comparison of the given <encseq>. The extension is performed
   using the xdrop strategy. The seed is specified
   by its length <len> and the two absolute positions <pos1> and <pos2> such
   that <pos1> is smaller than <pos2>.
   A <GtProcessinfo_and_querymatchspaceptr>-object is passed via the
   void pointer <info>.
   The extension is successful if it satisfies the considition specified
   in the <GtXdropmatchinfo>-object, passed as part of the
   <GtProcessinfo_and_querymatchspaceptr>-object. */

const GtQuerymatch *gt_xdrop_extend_selfmatch(void *info,
                                              const GtEncseq *encseq,
                                              GtUword len,
                                              GtUword pos1,
                                              GtUword pos2);

/*
   The following function performs an xdrop extension (using
   the previous function)
   and outputs the formatted match and possibly the alignment to stdout if
   the extension was successful.
   The function always returns 0, so the <GtError>-object <err> is not used.
*/

int gt_xdrop_extend_selfmatch_with_output(void *info,
                                          const GtEncseq *encseq,
                                          GtUword len,
                                          GtUword pos1,
                                          GtUword pos2,
                                          GT_UNUSED GtError *err);

/* The following function is identical to <gt_xdrop_extend_selfmatch>
   except that the positions of the seeds are defined by the number
   of the sequence they occur (dbseqnum for the first instance and querysenum
   for the second instance) and the relative position in that sequence
   (dbstart_relative for the first instance and querystart_relative
   for the second instance). */

const GtQuerymatch *gt_xdrop_extend_selfmatch_relative(void *info,
                                              const GtEncseq *encseq,
                                              GtUword dbseqnum,
                                              GtUword dbstart_relative,
                                              GtUword queryseqnum,
                                              GtUword querystart_relative,
                                              GtUword len,
                                              GtReadmode query_readmode);

/* The following function is used for extending a seed obtained
   in a comparison of the given sequence <query>
   against <encseq>. So here a byte sequence is compared against an
   encoded sequence and the seed is specified by <exactseed>.
   A <GtProcessinfo_and_querymatchspaceptr>-object is passed via the
   void pointer <info>.
   After the extension is performed and satisfies
   certain criteria, the resulting
   coordinates are returned as a reference to a <GtQuerymatch>-object. */

const GtQuerymatch *gt_xdrop_extend_querymatch(void *info,
                                               const GtEncseq *encseq,
                                               const GtQuerymatch *exactseed,
                                               const GtSeqorEncseq *query);

/*
   The following function performs an xdrop extension (as the previous function)
   and outputs the formatted match to stdout if
   the previous function returns a pointer different from NULL.
*/

void gt_xdrop_extend_querymatch_with_output(void *info,
                                            const GtEncseq *encseq,
                                            const GtQuerymatch *exactseed,
                                            const GtSeqorEncseq *query);

/* The following functions are used for the greedy extension. */

/* This is the type storing the relevant information for
   the greedy seed extension method. */

typedef struct GtGreedyextendmatchinfo GtGreedyextendmatchinfo;

/* The constructor, which is called once before the first seed
   is to be extended. <errorpercentage> is the percentage of errors
   allowed in the alignments reported.

   <maxalignedlendifference> is the maximum difference of the length
   of the aligned sequences for front-entries compared to the
   the arrow of the front. If <maxalignedlendifference> equals 0, then
   a reasonable default value depending on the the <errorpercentage>
   is automatically chosen.

   <history> is the size of the history. This is a value in the range
   from 1 to 64.

   <perc_mat_history> is the minimum percentage of the number of columns
   in the history are matches. This is a value in the range from
   1 to 100.

   <userdefinedleastlength> is the minimum
   length of the extension on both sides (including the seed itself).

   <extend_char_access> is the mode by which the characters are accessed
   in the encoded sequence.
   */

GtGreedyextendmatchinfo *gt_greedy_extend_matchinfo_new(
                                   GtUword errorpercentage,
                                   GtUword maxalignedlendifference,
                                   GtUword history,
                                   GtUword perc_mat_history,
                                   GtUword userdefinedleastlength,
                                   GtExtendCharAccess extend_char_access,
                                   GtUword sensitivity,
                                   const Polishing_info *pol_info);

/* the destructor-method for the gven object. */

void gt_greedy_extend_matchinfo_delete(GtGreedyextendmatchinfo *ggemi);

/* Set the check_extend_symmetry flag in the matchinfo object. */

void gt_greedy_extend_matchinfo_check_extend_symmetry_set(
                        GtGreedyextendmatchinfo *ggemi);

/* Set the silent flag in the matchinfo object. */

void gt_greedy_extend_matchinfo_silent_set(GtGreedyextendmatchinfo *ggemi);

/* Set the trimstat in the matchinfo object. */

void gt_greedy_extend_matchinfo_trimstat_set(GtGreedyextendmatchinfo *ggemi);

/* If <arg_maxalignedlendifference> and <arg_perc_mat_history> are 0, then
   an optimal value for the maximal alignment length difference and
   the percentage match history are determined, depending on the error
   percentage and sensitivity. The optimal values is determined by
   simulations with sequences of the given error percentage.
   The determined values are stored in the memory cells pointed
   to by <maxalignedlendifference> and <perc_mat_history>. If
   <arg_maxalignedlendifference> is not 0, then this value is
   used as the maximal alignment length difference and stored in the
   specified memory area. If <arg_perc_mat_history> is not 0, then this value is
   used as the percentage match history, and stored in the
   specified memory area.
*/

void gt_optimal_maxalilendiff_perc_mat_history(
                GtUword *maxalignedlendifference,
                GtUword *perc_mat_history,
                GtUword arg_maxalignedlendifference,
                GtUword arg_perc_mat_history,
                GtUword errorpercentage,
                GtUword sensitivity);

/* This function converts a string given as argument for option -cam
   and converts it to the given enum type <GtExtendCharAccess>. This
   option is used in the tool gt_repfind and gt_seedextend.
   In case of error, -1 is returned. */

GtExtendCharAccess gt_greedy_extend_char_access(const char *cam_string,
                                                GtError *err);

/* The following function returns a string specifying the possible arguments
   for the mentioned option -cam. */

const char *gt_cam_extendgreedy_comment(void);

/* The following function is used for extending a seed obtained
   in a self comparison of the given <encseq>. The extension is performed
   using the greedy strategy. The seed is specified
   by its length <len> and the two positions <pos1> and <pos2> such that
   <pos1> is smaller than <pos2>.
   A <GtProcessinfo_and_querymatchspaceptr>-object is passed via the
   void pointer <info>.
   After the extension is performed and the resulting alignment has an
   error rate below some threshold and a length above some threshold, then
   the cooordinates are delivered as a <GtQuerymatch>-object.
*/

const GtQuerymatch *gt_greedy_extend_selfmatch(void *info,
                                               const GtEncseq *encseq,
                                               GtUword len,
                                               GtUword pos1,
                                               GtUword pos2);

/*
   The following function performs a greedy extension (as the previous function)
   and outputs the formatted match and possibly the alignment to stdout.
   The function always returns 0, so the <GtError>-object <err> is not used.
*/

int gt_greedy_extend_selfmatch_with_output(void *info,
                                           const GtEncseq *encseq,
                                           GtUword len,
                                           GtUword pos1,
                                           GtUword pos2,
                                           GT_UNUSED GtError *err);

/* The following function is identical to <gt_greedy_extend_selfmatch>
   except that the positions of the seeds are defined by the number
   of the sequence they occur (dbseqnum for the first instance and querysenum
   for the second instance) and the relative position in that sequence
   (dbstart_relative for the first instance and querystart_relative
   for the second instance). */

const GtQuerymatch *gt_greedy_extend_selfmatch_relative(void *info,
                                              const GtEncseq *encseq,
                                              GtUword dbseqnum,
                                              GtUword dbstart_relative,
                                              GtUword queryseqnum,
                                              GtUword querystart_relative,
                                              GtUword len,
                                              GtReadmode query_readmode);

typedef const GtQuerymatch *(*GtExtendSelfmatchRelativeFunc)(void *,
                                                             const GtEncseq *,
                                                             GtUword,
                                                             GtUword,
                                                             GtUword,
                                                             GtUword,
                                                             GtUword,
                                                     GtReadmode query_readmode);

const GtQuerymatch* gt_xdrop_extend_querymatch_relative(
                                                  void *info,
                                                  const GtEncseq *dbencseq,
                                                  GtUword dbseqnum,
                                                  GtUword dbstart_relative,
                                                  const GtEncseq *queryencseq,
                                                  GtUword queryseqnum,
                                                  GtUword querystart_relative,
                                                  GtUword len,
                                                  GtReadmode query_readmode);

const GtQuerymatch* gt_greedy_extend_querymatch_relative(
                                                  void *info,
                                                  const GtEncseq *dbencseq,
                                                  GtUword dbseqnum,
                                                  GtUword dbstart_relative,
                                                  const GtEncseq *queryencseq,
                                                  GtUword queryseqnum,
                                                  GtUword querystart_relative,
                                                  GtUword len,
                                                  GtReadmode query_readmode);

typedef const GtQuerymatch *(*GtExtendQuerymatchRelativeFunc)(void *,
                                                             const GtEncseq *,
                                                             GtUword,
                                                             GtUword,
                                                             const GtEncseq *,
                                                             GtUword,
                                                             GtUword,
                                                             GtUword,
                                                             GtReadmode);

GtUword gt_align_front_prune_edist(bool rightextension,
                                   Polished_point *best_polished_point,
                                   GtFronttrace *front_trace,
                                   const GtEncseq *encseq,
                                   const GtSeqorEncseq *query,
                                   GtReadmode query_readmode,
                                   GtUword query_seqstartpos,
                                   GtUword query_totallength,
                                   GtGreedyextendmatchinfo *ggemi,
                                   bool greedyextension,
                                   GtUword seedlength,
                                   GtUword ustart,
                                   GtUword ulen,
                                   GtUword vstart,
                                   GtUword vlen);

GtUword gt_minidentity2errorpercentage(GtUword minidentity);

char *gt_seed_extend_params_keystring(bool use_greedy,
                                      bool forxdrop,
                                      unsigned int seedlength,
                                      unsigned int userdefinedleastlength,
                                      GtUword minidentity,
                                      GtUword maxalignedlendifference,
                                      GtUword perc_mat_history,
                                      GtUword extendgreedy,
                                      GtUword extendxdrop,
                                      GtUword xdropbelowscore);

void gt_greedy_extend_querymatch_with_output(void *info,
                                             const GtEncseq *dbencseq,
                                             const GtQuerymatch *exactseed,
                                             const GtSeqorEncseq *query);

double gt_greedy_dna_sequence_bias_get(const GtEncseq *encseq);

#endif
