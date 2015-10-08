/*
  Copyright (c) 2015 Annika <annika.seidel@studium.uni-hamburg.de>
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

#include <ctype.h>
#include <string.h>
#include "core/alphabet.h"
#include "core/chardef.h"
#include "core/cstr_api.h"
#include "core/fa.h"
#include "core/fasta_api.h"
#include "core/fasta_reader.h"
#include "core/fasta_reader_rec.h"
#include "core/ma.h"
#include "core/minmax.h"
#include "core/spacecalc.h"
#include "core/score_matrix.h"
#include "core/str.h"
#include "core/str_api.h"
#include "core/str_array.h"
#include "core/timer_api.h"
#include "core/types_api.h"
#include "core/unused_api.h"
#include "extended/diagonalbandalign.h"
#include "extended/diagonalbandalign_affinegapcost.h"
#include "extended/linearalign.h"
#include "extended/linearalign_affinegapcost.h"
#include "extended/linspaceManagement.h"
#include "extended/scorehandler.h"

#include "tools/gt_linspace_align.h"

#define LEFT_DIAGONAL_SHIFT(similarity, ulen, vlen) \
                                  -((1-similarity)*(MAX(ulen,vlen)) + 1 +\
                                   ((GtWord)ulen-(GtWord)vlen))

#define RIGHT_DIAGONAL_SHIFT(similarity, ulen, vlen) \
                                   ((1-similarity)*(MAX(ulen,vlen)) + 1 -\
                                   ((GtWord)ulen-(GtWord)vlen))

typedef struct{
  GtStr      *outputfile; /*default stdout*/
  GtStrArray *strings,
             *files,
             *linearcosts,
             *affinecosts,
             *diagonalbonds; /*left and right constant value shift of diagonal*/
  double     similarity; /*left and right shift of diagonal depends on
                                                                    similarity*/
  bool       global,
             local,
             diagonal, /* call diagonalband algorithm */
             dna,
             protein,
             costmatrix, /* special case of substituation matrix*/
             showscore,
             showsequences,
             scoreonly, /* dev option generate alignment, but do not show it*/
             wildcardshow, /* show symbol wildcards in output*/
             spacetime; /* write space peak and time overall on stdout*/
  GtUword timesquarefactor; /*factor to specified termination of recursion
                              and call 2dim algorithm */
} GtLinspaceArguments;

typedef struct {
  GtStr **seqarray;
  GtUword size, maxsize;
} GtSequences;

static void* gt_linspace_align_arguments_new(void)
{
  GtLinspaceArguments *arguments = gt_calloc((size_t) 1, sizeof *arguments);
  arguments->outputfile = gt_str_new();
  arguments->strings = gt_str_array_new();
  arguments->files = gt_str_array_new();
  arguments->linearcosts = gt_str_array_new();
  arguments->affinecosts = gt_str_array_new();
  arguments->diagonalbonds = gt_str_array_new();
  return arguments;
}

static void gt_linspace_align_arguments_delete(void *tool_arguments)
{
  GtLinspaceArguments *arguments = tool_arguments;
  if (arguments != NULL) {
    gt_str_delete(arguments->outputfile);
    gt_str_array_delete(arguments->strings);
    gt_str_array_delete(arguments->files);
    gt_str_array_delete(arguments->linearcosts);
    gt_str_array_delete(arguments->affinecosts);
    gt_str_array_delete(arguments->diagonalbonds);
    gt_free(arguments);
  }
}

static GtSequences* gt_sequences_new()
{
  GtSequences *sequences;

  sequences = gt_malloc(sizeof(*sequences));
  sequences->seqarray = gt_malloc(sizeof(*sequences->seqarray));
  sequences->size = 0;
  sequences->maxsize = 1;

  return sequences;
}

static void gt_sequences_delete(GtSequences *sequences)
{
  GtUword i;
  if (sequences != NULL) {
   for (i = 0; i < sequences->size ; i++)
    {
       gt_str_delete(sequences->seqarray[i]);
    }
    gt_free(sequences->seqarray);
    gt_free(sequences);
  }
}

static GtOptionParser* gt_linspace_align_option_parser_new(void *tool_arguments)
{
  GtLinspaceArguments *arguments = tool_arguments;
  GtOptionParser *op;
  GtOption *optionstrings, *optionfiles, *optionglobal, *optionlocal,
           *optiondna, *optionprotein, *optioncostmatrix, *optionlinearcosts,
           *optionaffinecosts, *optionoutputfile, *optionshowscore,
           *optionshowsequences, *optiondiagonal, *optiondiagonalbonds,
           *optionsimilarity, *optiontsfactor, *optionspacetime,
           *optionscoreonly, *optionwildcardsymbol;
           *optiontsfactor, *optionspacetime, *optionscoreonly,
           *optionwildcardsymbol;


  gt_assert(arguments);

  /* init */
  op = gt_option_parser_new("[ss|ff] sequence1 sequence2 [dna|protein] "
                            "[global|local] [a|l] costs/scores "
                            "[additional options]",
                            "Apply function to compute alignment.");
  gt_option_parser_set_mail_address(op,
                          "<annika.seidel@studium.uni-hamburg.de>");

  /* -bool */
  optionglobal = gt_option_new_bool("global", "global alignment",
                                    &arguments->global, false);
  gt_option_parser_add_option(op, optionglobal);

  optionlocal = gt_option_new_bool("local", "local alignment",
                                   &arguments->local, false);
  gt_option_parser_add_option(op, optionlocal);

  optiondiagonal = gt_option_new_bool("d", "diagonalband alignment",
                                      &arguments->diagonal, false);
  gt_option_parser_add_option(op, optiondiagonal);

  optiondna = gt_option_new_bool("dna", "type of sequences: dna",
                                 &arguments->dna, false);
  gt_option_parser_add_option(op, optiondna);

  optionprotein = gt_option_new_bool("protein", "type of sequences: protein",
                                      &arguments->protein, false);
  gt_option_parser_add_option(op, optionprotein);

  optionwildcardsymbol = gt_option_new_bool("wildcard", "show symbol used to "
                                            "represented wildcards in output",
                                            &arguments->wildcardshow, false);
  gt_option_parser_add_option(op, optionwildcardsymbol);

  /* special case, if given matrix includes cost values in place of scores */
  optioncostmatrix = gt_option_new_bool("costmatrix", "describes type of "
                                        "given substituation matrix",
                                        &arguments->costmatrix, false);
  gt_option_parser_add_option(op, optioncostmatrix);

  optionshowscore = gt_option_new_bool("showscore", "show score for alignment, "
                                       "please note it will calculate costs "
                                       "for global alignments and scores for "
                                       "local alignemnts always, independtly "
                                       "of input ",
                                       &arguments->showscore, false);
  gt_option_parser_add_option(op, optionshowscore);

  optionshowsequences = gt_option_new_bool("showsequences", "show sequences u "
                                           "and v in fornt of alignment",
                                       &arguments->showsequences, false);
  gt_option_parser_add_option(op, optionshowsequences);

  optionscoreonly = gt_option_new_bool("showonlyscore", "show only score for "
                                       "generated alignment to compare with "
                                       "other algorithmen",
                                       &arguments->scoreonly, false);
  gt_option_parser_add_option(op, optionscoreonly);

  optionspacetime = gt_option_new_bool("spacetime", "write space peak and time"
                                       " overall on stdout",
                                       &arguments->spacetime, false);
  gt_option_parser_add_option(op, optionspacetime);

  /* -str */
  optionstrings = gt_option_new_string_array("ss", "input, use two strings",
                                             arguments->strings);
  gt_option_parser_add_option(op, optionstrings);

  optionfiles = gt_option_new_filename_array("ff", "input, use two files",
                                             arguments->files);
  gt_option_parser_add_option(op, optionfiles);

  optionlinearcosts = gt_option_new_string_array("l", "lineargapcosts, "
                                                 "use match, mismatch and "
                                                 "gapcost, alternatively "
                                                 "substituationmatrix and "
                                                 "gapcost",
                                                 arguments->linearcosts);
  gt_option_parser_add_option(op, optionlinearcosts);

  optionaffinecosts = gt_option_new_string_array("a", "affinegapcosts, "
                                           "use match, mismatch, gap_extension "
                                           "and gap_opening, alternatively "
                                           "substituationmatrix, gap_extension "
                                           "and gap_opening",
                                           arguments->affinecosts);
  gt_option_parser_add_option(op, optionaffinecosts);

  optiondiagonalbonds = gt_option_new_string_array("lr", "specified left and "
                                                   "right shift of diagonal",
                                                   arguments->diagonalbonds);
  gt_option_parser_add_option(op, optiondiagonalbonds);

  optionoutputfile = gt_option_new_string("o", "print alignment, "
                                          "use outputfile",
                                          arguments->outputfile, "stdout");
  gt_option_parser_add_option(op, optionoutputfile);

  /* -ulong */
  optiontsfactor = gt_option_new_ulong("t", "timesquarefactor to organize "
                                       "time and space",
                                       &arguments->timesquarefactor,1);
  gt_option_parser_add_option(op, optiontsfactor);

  /* -double */
  optionsimilarity = gt_option_new_probability("similarity", "specified left "
                                               "and right shift of diagonal by "
                                               "similarity of sequences, "
                                               "0 <= similarty <= 1",
                                               &arguments->similarity, 0);
  gt_option_parser_add_option(op, optionsimilarity);

  /* dependencies */
  gt_option_is_mandatory_either(optionstrings, optionfiles);
  gt_option_is_mandatory_either(optiondna, optionprotein);
  gt_option_exclude(optionlocal, optionglobal);
  gt_option_exclude(optionlinearcosts, optionaffinecosts);
  gt_option_exclude(optiondna, optionprotein);
  gt_option_exclude(optionshowsequences, optionscoreonly);
  gt_option_exclude(optionsimilarity, optiondiagonalbonds);
  gt_option_imply_either_2(optionfiles, optionglobal, optionlocal);
  gt_option_imply_either_2(optiondna, optionstrings, optionfiles);
  gt_option_imply_either_2(optionstrings, optionglobal, optionlocal);
  gt_option_imply_either_2(optionprotein, optionstrings, optionfiles);
  gt_option_imply_either_2(optionlocal, optionlinearcosts, optionaffinecosts);
  gt_option_imply_either_2(optionglobal, optionlinearcosts, optionaffinecosts);
  gt_option_imply_either_2(optionshowscore,optionlinearcosts,optionaffinecosts);
  gt_option_imply_either_2(optionshowsequences, optionstrings, optionfiles);
  gt_option_imply_either_2(optionscoreonly,optionlinearcosts,optionaffinecosts);
  gt_option_imply(optiondiagonal, optionglobal);
  gt_option_imply(optiondiagonalbonds, optiondiagonal);
  gt_option_imply(optionsimilarity, optiondiagonal);
  gt_option_imply(optioncostmatrix, optionprotein);

  /* extended options */
  gt_option_is_extended_option(optiontsfactor);
  gt_option_is_extended_option(optionshowscore);
  gt_option_is_extended_option(optionwildcardsymbol);
  gt_option_is_extended_option(optioncostmatrix);

  /* development option(s) */
  gt_option_is_development_option(optionspacetime);
  gt_option_is_development_option(optionscoreonly);/*only useful to test*/

  return op;
}

static int gt_linspace_align_arguments_check(int rest_argc,
                                             void *tool_arguments,
                                             GtError *err)
{
  GtLinspaceArguments *arguments = tool_arguments;
  int had_err = 0;
  gt_error_check(err);
  gt_assert(arguments);

  if (rest_argc != 0)
  {
    gt_error_set(err,"superfluous arguments");
    return 1;
  }
  if ((gt_str_array_size(arguments->strings) > 0) &&
     (gt_str_array_size(arguments->strings) != 2UL))
  {
    gt_error_set(err, "option -ss requires two string arguments");
    had_err = 1;
  }
  if ((gt_str_array_size(arguments->files) > 0) &&
     (gt_str_array_size(arguments->files) != 2UL))
  {
    gt_error_set(err, "option -ff requires two file arguments");
    had_err = 1;
  }
  if (gt_str_array_size(arguments->linearcosts) > 0)
  {
    if (arguments->dna && gt_str_array_size(arguments->linearcosts) != 3UL)
    {
      gt_error_set(err, "option -l requires "
                        "match, mismatch, gap costs/scores when usign dna");
      had_err = 1;
    }
    if (arguments->protein && gt_str_array_size(arguments->linearcosts) != 2UL)
    {
      gt_error_set(err, "option -l requires  path of scorematrix "
                        "and gap costs/scores when usign protein");
      had_err = 1;
    }
  }
  if (gt_str_array_size(arguments->affinecosts) > 0)
  {
    if (arguments->dna && gt_str_array_size(arguments->affinecosts) != 4UL)
    {
      gt_error_set(err, "option -a requires match, mismatch, gap_opening, "
                        "gap_extending costs/scores when usign dna");
      had_err = 1;
    }
    if (arguments->protein && gt_str_array_size(arguments->affinecosts) != 3UL)
    {
      gt_error_set(err, "option -a requires path of scorematrix and "
                        "gap_opening, gap_extending costs/scores "
                        "when usign protein");
      had_err = 1;
    }
  }
  if ((gt_str_array_size(arguments->diagonalbonds) > 0) &&
     (gt_str_array_size(arguments->diagonalbonds) != 2UL))
  {
    gt_error_set(err, "option -lr requires left and right shift of diagonal");
    had_err = 1;
  }

  return had_err;
}

static void sequence_encode(GtUchar *seq, GtUword len, GtAlphabet *alphabet)
{
  GtUword idx;
  for (idx = 0; idx < len; idx++)
   seq[idx] = gt_alphabet_encode(alphabet, seq[idx]);
}

static void sequence_decode(GtUchar *seq, GtUword len, GtAlphabet *alphabet)
{
  GtUword idx;
  for (idx = 0; idx < len; idx++)
   seq[idx] = gt_alphabet_decode(alphabet, seq[idx]);
}

/*show sequence*/
static void print_sequence(const GtUchar *seq, GtUword len, FILE *fp)
{
  GtUword i = 0;
  fprintf(fp, "######\n");
  do{
    fprintf(fp, "%.80s\n",seq+i);
    i += 80;
  }while (i < len);

}

/*show sequences, alignment and score*/
static void alignment_show_with_sequences(GtUchar *useq, GtUword ulen,
                                          GtUchar *vseq, GtUword vlen,
                                          const GtAlignment *align,
                                          GtAlphabet *alphabet,
                                          GtUchar wildcardshow,
                                          bool showscore,
                                          bool showalign,
                                          bool showsequences,
                                          bool global,
                                          GtScoreHandler *scorehandler,
                                          FILE *fp)
{
  const GtUchar *characters;

  if (fp != NULL)
  {
    if (showsequences)
    {
      print_sequence(useq, ulen, fp);
      print_sequence(vseq, vlen, fp);
    }

    characters = gt_alphabet_characters(alphabet);
    sequence_encode(useq, ulen, alphabet);
    sequence_encode(vseq, vlen, alphabet);

    fprintf(fp, "######\n");
    if (showalign && gt_alignment_get_length(align) > 0)
    {
      gt_alignment_show_with_mapped_chars(align, characters,
                                          wildcardshow, fp, 80);
    }
    else if (showalign)
      fprintf(fp, "empty alignment\n");

    if (!showalign || showscore)
    {
      GtWord score = gt_scorehandler_eval_alignmentscore(scorehandler,
                                                         align, characters);

      fprintf(fp, "%s: "GT_WD"\n", global? "distance":"score", score);
    }

    sequence_decode(useq,ulen,alphabet);
    sequence_decode(vseq,vlen,alphabet);
  }
}

static GtWord select_digit_from_string(const char* str,
                                       bool global, GtError *err)
{
  bool haserr = false;
  GtWord evalue;

  if (sscanf(str ,GT_WD, &evalue) != 1)
  {
    gt_error_set(err, "found corrupt cost or score value");
    haserr = true;
  }

  if (global && evalue < 0)/*all costvalues have to be positive or 0*/
  {
    gt_error_set(err, "found invalid cost value "GT_WD, evalue);
    haserr = true;
  }

  if (haserr)
    return GT_WORD_MAX;

  return evalue;
}

static int save_fastasequence(const char *seqpart, GT_UNUSED GtUword length,
                              void *data, GT_UNUSED GtError* err)
{
  GtSequences *fasta_seqs = (GtSequences*) data;

  if (fasta_seqs->maxsize == fasta_seqs->size)
  {
    fasta_seqs->maxsize += 10;
    fasta_seqs->seqarray = gt_realloc(fasta_seqs->seqarray,
                                      fasta_seqs->maxsize*
                                      sizeof (*fasta_seqs->seqarray));
  }
  fasta_seqs->seqarray[fasta_seqs->size] = gt_str_new_cstr(seqpart);
  fasta_seqs->size++;

  return 0;
}

/*read fastasequences from file with GtFastaReader (-ff)*/
static int get_fastasequences(GtSequences *sequences, GtStr *filename,
                              GtError *err)
{
  int had_err = 0;
  GtFastaReader *reader;

  gt_assert(sequences != NULL);
  reader = gt_fasta_reader_rec_new (filename);
  had_err = gt_fasta_reader_run(reader, NULL, save_fastasequence,
                                NULL, sequences, err);
  gt_error_check(err);
  gt_fasta_reader_delete(reader);

  return had_err;
}

/* single sequences (-ss)*/
static int get_onesequence(GtSequences *sequences, const GtStrArray *strings,
                           GtUword idx, GtError *err)
{
  int had_err = 0;
  gt_assert(sequences != NULL && strings != NULL);

  if (gt_str_array_size(strings) <= idx)
  {
    gt_error_set(err, "out of range");
    return 1;
  }
  sequences->seqarray[0] = gt_str_new_cstr(gt_str_array_get(strings,idx));
  sequences->size++;

  return had_err;
}

#ifndef NDEBUG
/*checks if all character in <seq> are defined in <alphabet>. */
static inline int check_sequence(GtUchar *seq, GtUword len,
                                 GtAlphabet *alphabet, GtError *err)
{
  GtUword i;
  for (i = 0; i < len; i++)
  {
    if (!gt_alphabet_valid_input(alphabet, seq[i]))
    {
      gt_error_set(err, "found invalid character %c", seq[i]);
      return 1;
    }
    /*if (ISSPECIAL(gt_alphabet_encode(alphabet, seq[i])))
    {
      gt_error_set(err, "found wildcard character %c", seq[i]);
      return 1;
    }*/
  }
  return 0;
}
#endif

/*call function with linear gap costs for all given sequences */
static int alignment_with_linear_gap_costs(GtLinspaceArguments *arguments,
                                           GtError *err,
                                           GtScoreHandler *scorehandler,
                                           GtWord left_dist,
                                           GtWord right_dist,
                                           GtAlignment *align,
                                           LinspaceManagement *spacemanager,
                                           GtSequences *sequences1,
                                           GtSequences *sequences2,
                                           GtTimer *linspacetimer)
{
  GtUchar *useq, *vseq, wildcardshow;
  GtUword i, j, ulen, vlen;
  bool showsequences, showscore = false;
  int had_err = 0;

  gt_error_check(err);
  showsequences = arguments->showsequences;
  GtAlphabet *alphabet = gt_scorehandler_get_alphabet(scorehandler);
  wildcardshow = gt_alphabet_wildcard_show(alphabet);

  if (linspacetimer != NULL)
    gt_timer_start(linspacetimer);

  for (i = 0; i < sequences1->size; i++)
  {
      ulen = gt_str_length(sequences1->seqarray[i]);
      useq = (GtUchar*) gt_str_get(sequences1->seqarray[i]);
#ifndef NDEBUG
      had_err = check_sequence(useq, ulen, alphabet, err);
      if (had_err)
        return 1;
#endif

    for (j = 0; j< sequences2->size; j++)
    {
      vlen = gt_str_length(sequences2->seqarray[j]);
      vseq = (GtUchar*) gt_str_get(sequences2->seqarray[j]);
#ifndef NDEBUG
      had_err = check_sequence(vseq, vlen, alphabet, err);
      if (had_err)
        return 1;
#endif

      gt_alignment_reset(align);
      if (arguments->global)
      {
         if (arguments->diagonal)
         {
           if (gt_str_array_size(arguments->diagonalbonds) == 0)
           {
             left_dist =
                         LEFT_DIAGONAL_SHIFT(arguments->similarity, ulen, vlen);
             right_dist =
                        RIGHT_DIAGONAL_SHIFT(arguments->similarity, ulen, vlen);
           }

           if ((left_dist > MIN(0, (GtWord)vlen-(GtWord)ulen))||
               (right_dist < MAX(0, (GtWord)vlen-(GtWord)ulen)))
           {
              gt_error_set(err, "ERROR: invalid diagonalband for global "
                                "alignment (ulen: "GT_WU", vlen: "GT_WU")\n"
                                "left_dist <= MIN(0, vlen-ulen) and "
                                "right_dist >= MAX(0, vlen-ulen)", ulen, vlen);
              had_err = 1;
           }

           if (!had_err)
           {
             gt_computediagonalbandalign_generic(spacemanager, scorehandler,
                                                 align, useq, 0, ulen,
                                                 vseq, 0, vlen,
                                                 left_dist, right_dist);
           }
           else
             return 1;
         }
         else
         {
            gt_computelinearspace_generic(spacemanager, scorehandler, align,
                                          useq, 0, ulen, vseq, 0, vlen);
         }
      }
      else if (arguments->local)
      {
          gt_computelinearspace_local_generic(spacemanager, scorehandler, align,
                                              useq, 0, ulen, vseq, 0, vlen);
      }
      /* score */
      if (arguments->showscore)
      {
        showscore = true;
      }
      /* show alignment*/
      if (!had_err)
      {
        gt_assert(align != NULL);
        if (!strcmp(gt_str_get(arguments->outputfile),"stdout"))
        {
          alignment_show_with_sequences(useq, ulen, vseq, vlen, align, alphabet,
                                 wildcardshow, showscore, !arguments->scoreonly,
                                 showsequences, arguments->global,
                                 scorehandler, stdout);
        }
        else
        {
          FILE *fp = gt_fa_fopen_func(gt_str_get(arguments->outputfile),
                                                 "a", __FILE__,__LINE__,err);
          gt_error_check(err);
          alignment_show_with_sequences(useq, ulen, vseq, vlen, align, alphabet,
                                        wildcardshow, showscore,
                                        !arguments->scoreonly, showsequences,
                                        arguments->global, scorehandler,fp);
          gt_fa_fclose(fp);
        }
      }
    }
  }
  if (linspacetimer != NULL)
    gt_timer_stop(linspacetimer);

  if (!had_err && arguments->wildcardshow)
    printf("#wildcards are represented by %c\n", wildcardshow);

  return had_err;
}

/*call function with affine gap costs for all given sequences */
static int alignment_with_affine_gap_costs(GtLinspaceArguments *arguments,
                                           GtError *err,
                                           GtScoreHandler *scorehandler,
                                           GtWord left_dist,
                                           GtWord right_dist,
                                           GtAlignment *align,
                                           LinspaceManagement *spacemanager,
                                           GtSequences *sequences1,
                                           GtSequences *sequences2,
                                           GtTimer *linspacetimer)
{
  int had_err = 0;
  GtUchar *useq, *vseq, wildcardshow;
  GtUword i, j, ulen, vlen;
  bool showsequences, showscore = false;

  gt_error_check(err);
  showsequences = arguments->showsequences;
  GtAlphabet *alphabet = gt_scorehandler_get_alphabet(scorehandler);
  wildcardshow =  gt_alphabet_wildcard_show(alphabet);

  if (linspacetimer != NULL)
    gt_timer_start(linspacetimer);
  for (i = 0; i < sequences1->size; i++)
  {
    ulen = gt_str_length(sequences1->seqarray[i]);
    useq =  (GtUchar*) gt_str_get(sequences1->seqarray[i]);
#ifndef NDEBUG
    had_err = check_sequence(useq, ulen, alphabet, err);
    if (had_err)
      return 1;
#endif

    for (j = 0; j < sequences2->size; j++)
    {
      vlen = gt_str_length(sequences2->seqarray[j]);
      vseq = (GtUchar*) gt_str_get(sequences2->seqarray[j]);
#ifndef NDEBUG
      had_err = check_sequence(vseq, vlen, alphabet, err);
      if (had_err)
       return 1;
#endif

      gt_alignment_reset(align);
      if (arguments->global)
      {
        if (arguments->diagonal)
         {
           if (gt_str_array_size(arguments->diagonalbonds) == 0)
           {
             left_dist =
                         LEFT_DIAGONAL_SHIFT(arguments->similarity, ulen, vlen);
             right_dist =
                        RIGHT_DIAGONAL_SHIFT(arguments->similarity, ulen, vlen);
           }
           if ((left_dist > MIN(0, (GtWord)vlen-(GtWord)ulen))||
               (right_dist < MAX(0, (GtWord)vlen-(GtWord)ulen)))
           {
              gt_error_set(err, "ERROR: invalid diagonalband for global "
                                "alignment (ulen: "GT_WU", vlen: "GT_WU")\n"
                                "left_dist <= MIN(0, vlen-ulen) and "
                                "right_dist >= MAX(0, vlen-ulen)", ulen, vlen);
              had_err = 1;
           }

           if (!had_err)
           {
             gt_computediagonalbandaffinealign_generic(spacemanager,
                                                       scorehandler, align,
                                                       useq, 0, ulen,
                                                       vseq, 0, vlen,
                                                       left_dist, right_dist);
           }
           else
             return 1;
        }
        else
        {
          gt_computeaffinelinearspace_generic(spacemanager, scorehandler, align,
                                              useq, 0, ulen, vseq, 0, vlen);
        }
      }
      else if (arguments->local)
      {
        gt_computeaffinelinearspace_local_generic(spacemanager, scorehandler,
                                          align, useq, 0, ulen, vseq, 0, vlen);
      }

      /* score */
      if (arguments->showscore)
      {
        showscore = true;
      }
      /* show */
      if (!had_err)
      {
        gt_assert(align != NULL);
        if (!strcmp(gt_str_get(arguments->outputfile),"stdout"))
        {
          alignment_show_with_sequences(useq, ulen, vseq, vlen, align, alphabet,
                                 wildcardshow, showscore,
                                 !arguments->scoreonly, showsequences,
                                 arguments->global, scorehandler, stdout);
        }
        else
        {
          FILE *fp = gt_fa_fopen_func(gt_str_get(arguments->outputfile),
                                                 "a", __FILE__,__LINE__,err);
          gt_error_check(err);
          alignment_show_with_sequences(useq, ulen, vseq, vlen, align, alphabet,
                                        wildcardshow, showscore,
                                        !arguments->scoreonly, showsequences,
                                        arguments->global, scorehandler,fp);
          gt_fa_fclose(fp);
        }
      }
    }
  }
  if (linspacetimer != NULL)
    gt_timer_stop(linspacetimer);

  if (!had_err && arguments->wildcardshow)
    printf("wildcards are represented by %c\n", wildcardshow);

  return had_err;
}

/* handle score and cost values */
static int fill_scorehandler(GtScoreHandler **scorehandler,
                             GtLinspaceArguments *arguments,
                             GtError *err)
{
  int had_err = 0;
  GtWord matchscore, mismatchscore, gap_open, gap_extension;
  GtScoreMatrix *sm = NULL;

  gt_error_check(err);
  if (arguments->dna)
  {
    if (gt_str_array_size(arguments->linearcosts) > 0)
    {
      matchscore = select_digit_from_string(
                                     gt_str_array_get(arguments->linearcosts,0),
                                     arguments->global,err);
      mismatchscore = select_digit_from_string(
                                     gt_str_array_get(arguments->linearcosts,1),
                                     arguments->global,err);
      gap_open = 0;
      gap_extension = select_digit_from_string(
                                     gt_str_array_get(arguments->linearcosts,2),
                                     arguments->global,err);
    }
    else /*if (gt_str_array_size(arguments->affinecosts) > 0)*/
    {
      matchscore = select_digit_from_string(
                                     gt_str_array_get(arguments->affinecosts,0),
                                     arguments->global,err);
      mismatchscore = select_digit_from_string(
                                     gt_str_array_get(arguments->affinecosts,1),
                                     arguments->global,err);
      gap_open = select_digit_from_string(
                                     gt_str_array_get(arguments->affinecosts,2),
                                     arguments->global,err);
      gap_extension = select_digit_from_string(
                                     gt_str_array_get(arguments->affinecosts,3),
                                     arguments->global,err);
    }
    if (gt_error_is_set(err))
      return 1;
    *scorehandler = gt_scorehandler_new_DNA(matchscore, mismatchscore,
                                            gap_open, gap_extension);
  }
  else if (arguments->protein)
  {
    if (gt_str_array_size(arguments->linearcosts) > 0)
    {
      sm = gt_score_matrix_new_read_protein(
                               gt_str_array_get(arguments->linearcosts,0), err);
      if (gt_error_is_set(err))
        return 1;
      gap_open = 0;
      gap_extension = select_digit_from_string(
                          gt_str_array_get(arguments->linearcosts,1),false,err);

    }
    else /*if (gt_str_array_size(arguments->affinecosts) > 0)*/
    {
      sm = gt_score_matrix_new_read_protein(
                               gt_str_array_get(arguments->affinecosts,0), err);
      if (gt_error_is_set(err))
        return 1;
      gap_open = select_digit_from_string(
                          gt_str_array_get(arguments->affinecosts,1),false,err);
      gap_extension = select_digit_from_string(
                          gt_str_array_get(arguments->affinecosts,2),false,err);
    }

    if (gt_error_is_set(err))
      return 1;
    *scorehandler = gt_scorehandler_new_Protein(sm, gap_open, gap_extension);

    /* change score to costs if necessary */
    if (arguments->global && !arguments->costmatrix)
      gt_scorehandler_change_score_to_cost_without_costhandler(*scorehandler);
  }
  return had_err;
}

static int gt_linspace_align_runner(GT_UNUSED int argc,
                                    GT_UNUSED const char **argv,
                                    GT_UNUSED int parsed_args,
                                    void *tool_arguments,
                                    GtError *err)
{
  GtLinspaceArguments *arguments = tool_arguments;
  int had_err = 0;
  GtAlignment *align;
  GtWord left_dist = 0, right_dist = 0;
  GtSequences *sequences1, *sequences2;
  LinspaceManagement *spacemanager;
  GtScoreHandler *scorehandler = NULL;
  GtTimer *linspacetimer = NULL;

  gt_error_check(err);
  gt_assert(arguments);

  sequences1 = gt_sequences_new();
  sequences2 = gt_sequences_new();
  align = gt_alignment_new();
  spacemanager = gt_linspaceManagement_new();
  gt_linspaceManagement_set_TSfactor(spacemanager,
                                             arguments->timesquarefactor);

  /* get sequences */
  if (gt_str_array_size(arguments->strings) > 0)
  {
    get_onesequence(sequences1, arguments->strings, 0, err);
    gt_error_check(err);
    get_onesequence(sequences2, arguments->strings, 1, err);
    gt_error_check(err);
  }
  else if (gt_str_array_size(arguments->files) > 0)
  {
    get_fastasequences(sequences1,gt_str_array_get_str(arguments->files,0),err);
    gt_error_check(err);
    get_fastasequences(sequences2,gt_str_array_get_str(arguments->files,1),err);
    gt_error_check(err);
  }

  /* get costs/scores */
  had_err = fill_scorehandler(&scorehandler, arguments, err);

  /* get diagonal band */
  if (!had_err && arguments->diagonal)
  {
    if (gt_str_array_size(arguments->diagonalbonds) > 0)
    {
       left_dist = select_digit_from_string(
                      gt_str_array_get(arguments->diagonalbonds,0), false, err);
       right_dist = select_digit_from_string(
                      gt_str_array_get(arguments->diagonalbonds,1), false, err);
       gt_error_check(err);
    }
  }

  if (arguments->spacetime)
  {
    linspacetimer = gt_timer_new();
  }

  /* alignment functions with linear gap costs */
  if (!had_err && gt_str_array_size(arguments->linearcosts) > 0)
  {
    had_err  = alignment_with_linear_gap_costs(arguments, err, scorehandler,
                                               left_dist, right_dist,
                                               align, spacemanager,
                                               sequences1, sequences2,
                                               linspacetimer);
  }/* alignment functions with affine gap costs */
  else if (!had_err && gt_str_array_size(arguments->affinecosts) > 0)
  {
    had_err = alignment_with_affine_gap_costs(arguments, err, scorehandler,
                                              left_dist, right_dist,
                                              align, spacemanager,
                                              sequences1, sequences2,
                                              linspacetimer);
  }

  gt_sequences_delete(sequences1);
  gt_sequences_delete(sequences2);
  gt_alignment_delete(align);
  gt_scorehandler_delete(scorehandler);

  /*spacetime option*/
  if (!had_err && arguments->spacetime)
  {
    printf("# combined space peak in kilobytes: %f\n",
           GT_KILOBYTES(gt_linspaceManagement_get_spacepeak(spacemanager)));
    printf("# TIME");

    gt_timer_show_formatted(linspacetimer," overall " GT_WD ".%02ld\n",stdout);
    gt_timer_delete(linspacetimer);
  }
  gt_linspaceManagement_delete(spacemanager);
  return had_err;
}

GtTool* gt_linspace_align(void)
{
  return gt_tool_new(gt_linspace_align_arguments_new,
                     gt_linspace_align_arguments_delete,
                     gt_linspace_align_option_parser_new,
                     gt_linspace_align_arguments_check,
                     gt_linspace_align_runner);
}
