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
#include "core/cstr_api.h"
#include "core/fa.h"
#include "core/fasta_api.h"
#include "core/fasta_reader.h"
#include "core/fasta_reader_rec.h"
#include "core/ma.h"
#include "core/str.h"
#include "core/str_api.h"
#include "core/str_array.h"
#include "core/types_api.h"
#include "core/unused_api.h"
#include "extended/diagonalbandalign.h"
#include "extended/diagonalbandalign_affinegapcost.h"
#include "extended/linearalign.h"
#include "extended/linearalign_affinegapcost.h"
#include "extended/linspaceManagement.h"

#include "tools/gt_linspace_align.h"
#define CHECK_COST(cost) { if (cost < 0) {\
                        fprintf(stderr, "invalid cost value " GT_WD"\n", cost);\
                        exit(GT_EXIT_PROGRAMMING_ERROR);}}

typedef struct{
  GtStr      *outputfile;
  GtStrArray *strings,
             *files,
             *linearcosts,
             *affinecosts,
             *diagonalbonds;
  bool       global,
             local,
             diagonal,
             showscore,
             spacetime;
  GtUword timesquarefactor;
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
           *optionlinearcosts, *optionaffinecosts, *optionoutputfile,
           *optionshowscore, *optiondiagonal, *optiondiagonalbonds,
           *optiontsfactor, *optionspacetime;
  gt_assert(arguments);

  /* init */
  op = gt_option_parser_new("[ss|ff] sequence1 sequence2 [global|local] "
                            "[a|l] [additional options]",
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

  optionshowscore = gt_option_new_bool("showscore", "show score for alignment",
                                      &arguments->showscore, false);
  gt_option_parser_add_option(op, optionshowscore);

  optionspacetime = gt_option_new_bool("spacetime", "write space peak and time"
                                       " overall on stdout",
                                       &arguments->spacetime, false);
  gt_option_parser_add_option(op, optionspacetime);

  /* -str */
  optionstrings = gt_option_new_string_array("ss", "use two strings",
                                             arguments->strings);
  gt_option_parser_add_option(op, optionstrings);

  optionfiles = gt_option_new_filename_array("ff", "use two files",
                                             arguments->files);
  gt_option_parser_add_option(op, optionfiles);

  optionlinearcosts = gt_option_new_string_array("l", "lineargapcosts, "
                                                 "use three values",
                                                 arguments->linearcosts);
  gt_option_parser_add_option(op, optionlinearcosts);

  optionaffinecosts = gt_option_new_string_array("a", "affinegapcosts, "
                                                 "use four values",
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

  /* dependencies */
  gt_option_is_mandatory_either(optionstrings, optionfiles);
  gt_option_exclude(optionlocal, optionglobal);
  gt_option_exclude(optionlinearcosts, optionaffinecosts);
  gt_option_imply_either_2(optionstrings, optionglobal, optionlocal);
  gt_option_imply_either_2(optionfiles, optionglobal, optionlocal);
  gt_option_imply_either_2(optionlocal, optionlinearcosts, optionaffinecosts);
  gt_option_imply_either_2(optionglobal, optionlinearcosts, optionaffinecosts);
  gt_option_imply_either_2(optionshowscore,optionlinearcosts,optionaffinecosts);
  gt_option_imply(optiondiagonal, optionglobal);
  gt_option_imply(optiondiagonalbonds, optiondiagonal);
  /* development option(s) */
  gt_option_is_development_option(optionspacetime);

  return op;
}

static int gt_linspace_align_arguments_check(GT_UNUSED int rest_argc,
                                             void *tool_arguments,
                                             GtError *err)
{
  GtLinspaceArguments *arguments = tool_arguments;
  int had_err = 0;
  gt_error_check(err);
  gt_assert(arguments);

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
  if ((gt_str_array_size(arguments->linearcosts) > 0) &&
     (gt_str_array_size(arguments->linearcosts) != 3UL))
  {
    gt_error_set(err, "option -l requires "
                      "match, mismatch, gap costs/scores");
    had_err = 1;
  }
  if ((gt_str_array_size(arguments->affinecosts) > 0) &&
     (gt_str_array_size(arguments->affinecosts) != 4UL))
  {
    gt_error_set(err, "option -a requires match, mismatch, "
                      "gap_opening, gap_extending costs/scores");
    had_err = 1;
  }
  if ((gt_str_array_size(arguments->diagonalbonds) > 0) &&
     (gt_str_array_size(arguments->diagonalbonds) != 2UL))
  {
    gt_error_set(err, "option -lr requires left and right shift of diagonal");
    had_err = 1;
  }

  return had_err;
}

static void print_sequence(const GtUchar *seq, GtUword len, FILE *fp)
{
  GtUword i = 0;
  fprintf(fp, "######\n");
  do{
    fprintf(fp, "%.80s\n",seq+i);
    i += 80;
  }while (i < len);
}

static void alignment_with_seqs_show(const GtUchar *useq, GtUword ulen,
                                     const GtUchar *vseq, GtUword vlen,
                                     const GtAlignment *align,
                                     GtWord score, FILE *fp)
{
  if (fp != NULL)
  {
    print_sequence(useq, ulen, fp);
    print_sequence(vseq, vlen, fp);
    fprintf(fp, "######\n");

    if (gt_alignment_get_length(align) > 0)
      gt_alignment_show(align, fp, 80);
    else
      fprintf(fp, "empty alignment\n");
    if (score != GT_WORD_MAX)
      fprintf(fp, "score: "GT_WD"\n", score);
  }
}

static GtWord* select_digit_from_string(const GtStrArray *arr,
                                        bool global, GtError *err)
{
  bool haserr = false;
  GtWord *evalues;
  GtUword size, i;

  gt_assert(arr != NULL);
  size = gt_str_array_size(arr);
  evalues = gt_malloc(sizeof(*evalues)*size);
  for (i = 0; !haserr && i < size; i++)
  {
    if (sscanf(gt_str_array_get(arr,i),GT_WD, &evalues[i]) != 1)
    {
      gt_error_set(err, "found invalid cost or score");
      haserr = true;
    }
    if (global)
      CHECK_COST(evalues[i]);
  }
  if (haserr)
  {
    gt_free(evalues);
    return NULL;
  }
  return evalues;
}

static int save_fastasequence(const char *seqpart, GT_UNUSED GtUword length,
                              void *data, GT_UNUSED GtError* err)
{
  GtSequences *fasta_seqs = (GtSequences*) data;

  if (fasta_seqs->maxsize == fasta_seqs->size)
  {
    fasta_seqs->maxsize += 5;
    fasta_seqs->seqarray = gt_realloc(fasta_seqs->seqarray,
                                      fasta_seqs->maxsize*
                                      sizeof (*fasta_seqs->seqarray));
  }
  fasta_seqs->seqarray[fasta_seqs->size] = gt_str_new_cstr(seqpart);
  fasta_seqs->size++;

  return 0;
}

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

static inline void sequence_to_lower(GtUchar *seq, GtUword len)
{
  GtUword i;

  for (i = 0; i < len; i++)
    seq[i] = tolower((int)seq[i]);
}
/*call function with linear gap costs for all given sequences */
static void alignment_with_linear_gap_costs(GtLinspaceArguments *arguments,
                                            GtError *err, GtWord *linearcosts,
                                            GtWord *diagonalbonds,
                                            GtAlignment *align,
                                            LinspaceManagement *spacemanager,
                                            GtSequences *sequences1,
                                            GtSequences *sequences2)
{
  GtUchar *useq, *vseq;
  GtUword i, j, ulen, vlen;
  GtWord score = GT_WORD_MAX;

  for (i = 0; i < sequences1->size; i++)
  {
    for (j = 0; j< sequences2->size; j++)
    {
      ulen = gt_str_length(sequences1->seqarray[i]);
      useq = (GtUchar*) gt_str_get(sequences1->seqarray[i]);
      sequence_to_lower(useq, ulen);
      vlen = gt_str_length(sequences2->seqarray[j]);
      vseq = (GtUchar*) gt_str_get(sequences2->seqarray[j]);
      sequence_to_lower(vseq, vlen);

      gt_alignment_reset(align);
      if (arguments->global)
      {
         if (arguments->diagonal)
         {
           if (!gt_str_array_size(arguments->diagonalbonds) > 0)
           {
             gt_assert (diagonalbonds != NULL);
             diagonalbonds[0] = -ulen;
             diagonalbonds[1] = vlen;
           }
           gt_computediagonalbandalign(spacemanager,align,
                                       useq, 0, ulen, vseq, 0, vlen,
                                       diagonalbonds[0], diagonalbonds[1],
                                       linearcosts[0], linearcosts[1],
                                       linearcosts[2]);
         }
         else
         {
           gt_computelinearspace(spacemanager, align,
                                 useq, 0, ulen, vseq, 0, vlen,
                                 linearcosts[0], linearcosts[1],
                                 linearcosts[2]);
         }
      }
      else if (arguments->local)
      {
          gt_computelinearspace_local(spacemanager, align,
                                      useq, 0, ulen, vseq, 0, vlen,
                                      linearcosts[0], linearcosts[1],
                                      linearcosts[2]);
      }
      /* score */
      if (arguments->showscore)
      {
        score = gt_alignment_eval_generic_with_score(false, align,
                                                     linearcosts[0],
                                                     linearcosts[1],
                                                     linearcosts[2]);
      }
      /* show alignment*/
      gt_assert(align != NULL);
      if (!strcmp(gt_str_get(arguments->outputfile),"stdout"))
        alignment_with_seqs_show(useq, ulen, vseq, vlen, align, score, stdout);
      else
      {
        FILE *fp = gt_fa_fopen_func(gt_str_get(arguments->outputfile),
                                               "a", __FILE__,__LINE__,err);
        gt_error_check(err);
        alignment_with_seqs_show(useq, ulen, vseq, vlen, align, score, fp);
        gt_fa_fclose(fp);
      }
    }
  }
  if (diagonalbonds != NULL)
    gt_free(diagonalbonds);
}

/*call function with affine gap costs for all given sequences */
static void alignment_with_affine_gap_costs(GtLinspaceArguments *arguments,
                                            GtError *err,
                                            GtWord *affinecosts,
                                            GtWord *diagonalbonds,
                                            GtAlignment *align,
                                            LinspaceManagement *spacemanager,
                                            GtSequences *sequences1,
                                            GtSequences *sequences2)
{
  GtUchar *useq, *vseq;
  GtUword i, j, ulen, vlen;
  GtWord score = GT_WORD_MAX;

  for (i = 0; i < sequences1->size; i++)
  {
    for (j = 0; j< sequences2->size; j++)
    {
      ulen = gt_str_length(sequences1->seqarray[i]);
      useq =  (GtUchar*) gt_str_get(sequences1->seqarray[i]);
      sequence_to_lower(useq, ulen);
      vlen = gt_str_length(sequences2->seqarray[j]);
      vseq = (GtUchar*) gt_str_get(sequences2->seqarray[j]);
      sequence_to_lower(vseq, vlen);

      gt_alignment_reset(align);
      if (arguments->global)
      {
        if (arguments->diagonal)
         {
           if (!gt_str_array_size(arguments->diagonalbonds) > 0)
           {
             gt_assert (diagonalbonds != NULL);
             diagonalbonds[0] = -ulen;
             diagonalbonds[1] = vlen;
           }
           gt_computediagonalbandaffinealign(spacemanager, align,
                                             useq, 0, ulen, vseq, 0, vlen,
                                             diagonalbonds[0], diagonalbonds[1],
                                             affinecosts[0], affinecosts[1],
                                             affinecosts[2], affinecosts[3]);
        }
        else
        {
          gt_computeaffinelinearspace(spacemanager, align,
                                      useq, 0, ulen, vseq, 0, vlen,
                                      affinecosts[0], affinecosts[1],
                                      affinecosts[2], affinecosts[3]);
        }
      }
      else if (arguments->local)
      {
        gt_computeaffinelinearspace_local(spacemanager, align,
                                          useq, 0, ulen, vseq, 0,
                                          vlen, affinecosts[0],
                                          affinecosts[1],
                                          affinecosts[2],
                                          affinecosts[3]);
      }

      /* score */
      if (arguments->showscore)
      {
        score = gt_alignment_eval_generic_with_affine_score(false, align,
                                                            affinecosts[0],
                                                            affinecosts[1],
                                                            affinecosts[2],
                                                            affinecosts[3]);
      }
      /* show */
      gt_assert(align != NULL);
      if (!strcmp(gt_str_get(arguments->outputfile),"stdout"))
        alignment_with_seqs_show(useq, ulen, vseq, vlen, align, score, stdout);
      else
      {
        FILE *fp = gt_fa_fopen_func(gt_str_get(arguments->outputfile),
                                               "a", __FILE__,__LINE__,err);
        gt_error_check(err);
        alignment_with_seqs_show(useq, ulen, vseq, vlen, align, score, fp);
        gt_fa_fclose(fp);
      }
    }
  }
  if (diagonalbonds != NULL)
    gt_free(diagonalbonds);
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
  GtWord *diagonalbonds = NULL;
  GtSequences *sequences1, *sequences2;
  LinspaceManagement *spacemanager;

  gt_error_check(err);
  gt_assert(arguments);

  sequences1 = gt_sequences_new();
  sequences2 = gt_sequences_new();
  align = gt_alignment_new();
  spacemanager = gt_linspaceManagement_new();
  gt_linspaceManagement_set_TSfactor(spacemanager,
                                             arguments->timesquarefactor);

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

  /* call alignment functions */

  if (arguments->diagonal)
  {
    if (gt_str_array_size(arguments->diagonalbonds) > 0)
    {
       diagonalbonds = select_digit_from_string(arguments->diagonalbonds,
                                                false, err);
       gt_error_check(err);
    }
    else
      diagonalbonds = gt_malloc(sizeof (*diagonalbonds)*2);
  }
  /* alignment functions with linear gap costs */
  if (gt_str_array_size(arguments->linearcosts) > 0)
  {
    GtWord *linearcosts;
    gt_assert(gt_str_array_size(arguments->linearcosts) == 3UL);
    linearcosts = select_digit_from_string(arguments->linearcosts,
                                           arguments->global, err);
    gt_error_check(err);
    if (linearcosts == NULL)
      return 1;

    alignment_with_linear_gap_costs(arguments, err,
                                    linearcosts, diagonalbonds,
                                    align, spacemanager,
                                    sequences1, sequences2);
    gt_free(linearcosts);
  }/* alignment functions with affine gap costs */
  else if (gt_str_array_size(arguments->affinecosts) > 0)
  {
    GtWord *affinecosts;

    gt_assert(gt_str_array_size(arguments->affinecosts) == 4UL);
    affinecosts = select_digit_from_string(arguments->affinecosts,
                                           arguments->global, err);
    gt_error_check(err);
    if (affinecosts == NULL)
      return 1;

    alignment_with_affine_gap_costs(arguments, err,
                                    affinecosts, diagonalbonds,
                                    align, spacemanager,
                                    sequences1, sequences2);
     gt_free(affinecosts);
  }

  gt_sequences_delete(sequences1);
  gt_sequences_delete(sequences2);
  gt_alignment_delete(align);

  if (arguments->spacetime)
  {
    printf("# space peak: "GT_ZU"Bytes\n",
                             gt_linspaceManagement_get_spacepeak(spacemanager));
    /*TODO: add time*/
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
