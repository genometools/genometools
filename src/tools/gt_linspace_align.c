/*
  Copyright (c) 2015 Annika Seidel <annika.seidel@studium.uni-hamburg.de>
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
#include "extended/linspace_management.h"
#include "extended/scorehandler.h"
#include "tools/gt_linspace_align.h"

#define LEFT_DIAGONAL_SHIFT(similarity, ulen, vlen) \
                            -((1-similarity)*(MAX(ulen,vlen)) +\
                            MIN(((GtWord)ulen-(GtWord)vlen),0))

#define RIGHT_DIAGONAL_SHIFT(similarity, ulen, vlen) \
                             ((1-similarity)*(MAX(ulen,vlen)) -\
                             MAX(((GtWord)ulen-(GtWord)vlen),0))

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
             has_costmatrix, /* special case of substituation matrix*/
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
} GtSequenceTable;

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

static GtSequenceTable* gt_sequence_table_new()
{
  GtSequenceTable *sequence_table;

  sequence_table = gt_malloc(sizeof(*sequence_table));
  sequence_table->seqarray = gt_malloc(sizeof(*sequence_table->seqarray));
  sequence_table->size = 0;
  sequence_table->maxsize = 1;
  return sequence_table;
}

static void gt_sequence_table_delete(GtSequenceTable *sequence_table)
{
  GtUword i;
  if (sequence_table != NULL)
  {
    for (i = 0; i < sequence_table->size ; i++)
    {
      gt_str_delete(sequence_table->seqarray[i]);
    }
    gt_free(sequence_table->seqarray);
    gt_free(sequence_table);
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

  optiondna = gt_option_new_bool("dna", "type of sequences: DNA",
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
                                        &arguments->has_costmatrix, false);
  gt_option_parser_add_option(op, optioncostmatrix);

  optionshowscore = gt_option_new_bool("showscore", "show score for alignment, "
                                       "please note it will calculate costs "
                                       "for global alignments and scores for "
                                       "local alignments always, independently "
                                       "of input ",
                                       &arguments->showscore, false);
  gt_option_parser_add_option(op, optionshowscore);

  optionshowsequences = gt_option_new_bool("showsequences", "show sequences u "
                                           "and v in front of alignment",
                                       &arguments->showsequences, false);
  gt_option_parser_add_option(op, optionshowsequences);

  optionscoreonly = gt_option_new_bool("showonlyscore", "show only score for "
                                       "generated alignment to compare with "
                                       "other algorithms",
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
  gt_assert(arguments != NULL);
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

/*show sequence*/
static void gt_linspace_print_sequence(const GtUchar *characters,
                                       GtUchar wildcardshow,
                                       const GtUchar *seq, GtUword len,
                                       FILE *fp)
{
  const GtUword linewidth = 80;
  GtUword idx;

  fprintf(fp, "######\n");
  for (idx = 0; idx < len; idx++)
  {
    GtUchar cc;
    if (characters == NULL)
    {
      cc = seq[idx];
    } else
    {
      if (ISSPECIAL(seq[idx]))
      {
        cc = wildcardshow;
      } else
      {
        cc = characters[seq[idx]];
      }
    }
    fputc(cc,fp);
    if ((idx + 1) % linewidth == 0)
    {
      fputc('\n',fp);
    }
  }
  if (idx % linewidth != 0)
  {
    fputc('\n',fp);
  }
}

/*show sequences, alignment and score*/
static void alignment_show_with_sequences(const GtUchar *useq, GtUword ulen,
                                          const GtUchar *vseq, GtUword vlen,
                                          const GtAlignment *align,
                                          const GtUchar *characters,
                                          GtUchar wildcardshow,
                                          bool showscore,
                                          bool showalign,
                                          bool showsequences,
                                          bool global,
                                          const GtScoreHandler *scorehandler,
                                          FILE *fp)
{
  if (fp != NULL)
  {
    if (showsequences)
    {
      gt_linspace_print_sequence(characters,wildcardshow,useq, ulen, fp);
      gt_linspace_print_sequence(characters,wildcardshow,vseq, vlen, fp);
    }
    fprintf(fp, "######\n");
    if (showalign && gt_alignment_get_length(align) > 0)
    {
      gt_alignment_show_with_mapped_chars(align, characters,
                                          wildcardshow, fp, 80);
    } else
    {
      if (showalign)
      {
        fprintf(fp, "empty alignment\n");
      }
    }

    if (!showalign || showscore)
    {
      GtWord score = gt_scorehandler_eval_alignmentscore(scorehandler,
                                                         align, characters);

      fprintf(fp, "%s: "GT_WD"\n", global? "distance" : "score", score);
    }
  }
}

static int gt_parse_score_value(int line,GtWord *value,const char* str,
                                GT_UNUSED bool non_negative, GtError *err)
{
  if (sscanf(str, GT_WD, value) != 1 || (non_negative && *value < 0))
  {
    gt_error_set(err,"line %d: invalid %s value \"%s\"",line,
                 non_negative ? "cost" : "score",str);
    return -1;
  }
  return 0;
}

static int save_fastasequence(const char *seqpart, GT_UNUSED GtUword length,
                              void *data, GT_UNUSED GtError* err)
{
  GtSequenceTable *fasta_seqs = (GtSequenceTable*) data;

  if (fasta_seqs->maxsize == fasta_seqs->size)
  {
    fasta_seqs->maxsize += 10;
    fasta_seqs->seqarray = gt_realloc(fasta_seqs->seqarray,
                                      fasta_seqs->maxsize*
                                      sizeof (*fasta_seqs->seqarray));
  }
  fasta_seqs->seqarray[fasta_seqs->size++] = gt_str_new_cstr(seqpart);
  return 0;
}

/*read fastasequences from file with GtFastaReader (-ff)*/
static int get_fastasequences(GtSequenceTable *sequence_table,
                              const GtStr *filename,
                              GtError *err)
{
  int had_err = 0;
  GtFastaReader *reader;

  gt_assert(sequence_table != NULL);
  reader = gt_fasta_reader_rec_new (filename);
  had_err = gt_fasta_reader_run(reader, NULL, save_fastasequence,
                                NULL, sequence_table, err);
  gt_fasta_reader_delete(reader);
  return had_err;
}

/* single sequences (-ss)*/
static void get_onesequence(const GtSequenceTable *sequence_table,
                            const GtStrArray *strings,
                            GtUword idx)
{
  gt_assert(sequence_table != NULL && strings != NULL &&
            idx < gt_str_array_size(strings));

  sequence_table->seqarray[0] = gt_str_new_cstr(gt_str_array_get(strings,idx));
}

/*call function with linear gap costs for all given sequences */
static int gt_all_against_all_alignment_check(bool affine,
                                        GtAlignment *align,
                                        const GtLinspaceArguments *arguments,
                                        GtLinspaceManagement *spacemanager,
                                        const GtScoreHandler *scorehandler,
                                        const GtUchar *characters,
                                        GtUchar wildcardshow,
                                        const GtSequenceTable *sequence_table1,
                                        const GtSequenceTable *sequence_table2,
                                        GtWord left_dist,
                                        GtWord right_dist,
                                        GtTimer *linspacetimer,
                                        GtError *err)
{
  int had_err = 0;
  const GtUchar *useq, *vseq;
  GtUword i, j, ulen, vlen;

  gt_error_check(err);
  if (linspacetimer != NULL)
  {
    gt_timer_start(linspacetimer);
  }
  for (i = 0; !had_err && i < sequence_table1->size; i++)
  {
    ulen = gt_str_length(sequence_table1->seqarray[i]);
    useq = (const GtUchar*) gt_str_get(sequence_table1->seqarray[i]);
    for (j = 0; j< sequence_table2->size; j++)
    {
      vlen = gt_str_length(sequence_table2->seqarray[j]);
      vseq = (const GtUchar*) gt_str_get(sequence_table2->seqarray[j]);
      gt_alignment_reset(align);
      if (arguments->global)
      {
        if (arguments->diagonal)
        {
          if (gt_str_array_size(arguments->diagonalbonds) == 0)
          {
            left_dist = LEFT_DIAGONAL_SHIFT(arguments->similarity, ulen, vlen);
            right_dist = RIGHT_DIAGONAL_SHIFT(arguments->similarity, ulen,
                                              vlen);
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
            (affine ? gt_diagonalbandalign_affinegapcost_compute_generic
                    : gt_diagonalbandalign_compute_generic)
                       (spacemanager, scorehandler, align,
                        useq, 0, ulen, vseq, 0, vlen,
                        left_dist, right_dist);
          }
        } else
        {
          (affine ? gt_linearalign_affinegapcost_compute_generic
                  : gt_linearalign_compute_generic)
                             (spacemanager, scorehandler, align,
                              useq, 0, ulen, vseq, 0, vlen);
        }
      }
      else if (arguments->local)
      {
        (affine ? gt_linearalign_affinegapcost_compute_local_generic
                : gt_linearalign_compute_local_generic)
                    (spacemanager, scorehandler, align,
                     useq, 0, ulen, vseq, 0, vlen);
      }
      /* show alignment*/
      if (!had_err)
      {
        gt_assert(align != NULL);
        if (!strcmp(gt_str_get(arguments->outputfile),"stdout"))
        {
          alignment_show_with_sequences(useq, ulen, vseq, vlen, align,
                                        characters,
                                        wildcardshow, arguments->showscore,
                                        !arguments->scoreonly,
                                        arguments->showsequences,
                                        arguments->global,
                                        scorehandler, stdout);
        } else
        {
          FILE *fp = gt_fa_fopen_func(gt_str_get(arguments->outputfile),
                                                 "a", __FILE__,__LINE__,err);
          if (fp == NULL)
          {
            had_err = -1;
          } else
          {
            alignment_show_with_sequences(useq, ulen, vseq, vlen, align,
                                          characters, wildcardshow,
                                          arguments->showscore,
                                          !arguments->scoreonly,
                                          arguments->showsequences,
                                          arguments->global, scorehandler,fp);
            gt_fa_fclose(fp);
          }
        }
      }
    }
  }
  if (linspacetimer != NULL)
  {
    gt_timer_stop(linspacetimer);
  }
  if (!had_err && arguments->wildcardshow)
  {
    printf("# wildcards are represented by %c\n", wildcardshow);
  }
  return had_err;
}

/* handle score and cost values */
static GtScoreHandler *gt_arguments2scorehandler(
                             const GtLinspaceArguments *arguments,
                             GtError *err)
{
  GtWord matchscore, mismatchscore, gap_open, gap_extension;
  GtScoreHandler *scorehandler = NULL;
  GtScoreMatrix *scorematrix = NULL;
  int had_err = 0;

  gt_error_check(err);
  if (gt_str_array_size(arguments->linearcosts) > 0)
  {
    GtUword wordindex = 0;

    if (arguments->protein)
    {
      scorematrix
        = gt_score_matrix_new_read_protein(
                          gt_str_array_get(arguments->linearcosts,wordindex++),
                          err);
      if (scorematrix == NULL)
      {
        had_err = -1;
      }
      matchscore = 0;
      mismatchscore = 0;
    } else
    {
      had_err = gt_parse_score_value(__LINE__,&matchscore,
                                     gt_str_array_get(arguments->linearcosts,
                                                      wordindex++),
                                     arguments->global,err);
      if (!had_err)
      {
        had_err = gt_parse_score_value(__LINE__,&mismatchscore,
                                       gt_str_array_get(arguments->linearcosts,
                                                        wordindex++),
                                       arguments->global,err);
      }
    }
    if (!had_err)
    {
      gap_open = 0;
      had_err = gt_parse_score_value(__LINE__,&gap_extension,
                                     gt_str_array_get(arguments->linearcosts,
                                                      wordindex++),
                                     false,err);
    }
  } else /*if (gt_str_array_size(arguments->affinecosts) > 0)*/
  {
    GtUword wordindex = 0;

    if (arguments->protein)
    {
      scorematrix = gt_score_matrix_new_read_protein(
                               gt_str_array_get(arguments->affinecosts,
                                                wordindex++), err);
      if (scorematrix == NULL)
      {
        had_err = -1;
      }
      matchscore = mismatchscore = 0;
    } else
    {
      had_err = gt_parse_score_value(__LINE__,&matchscore,
                                     gt_str_array_get(arguments->affinecosts,
                                                      wordindex++),
                                     arguments->global,err);
      if (!had_err)
      {
        had_err = gt_parse_score_value(__LINE__,&mismatchscore,
                                       gt_str_array_get(arguments->affinecosts,
                                                        wordindex++),
                                       arguments->global,err);
      }
    }
    if (!had_err)
    {
      had_err = gt_parse_score_value(__LINE__,&gap_open,
                                     gt_str_array_get(arguments->affinecosts,
                                                      wordindex++),
                                     false,err);
    }
    if (!had_err)
    {
      had_err = gt_parse_score_value(__LINE__,&gap_extension,
                                     gt_str_array_get(arguments->affinecosts,
                                                      wordindex),
                                     false,err);
    }
  }
  if (!had_err)
  {
    scorehandler = gt_scorehandler_new(matchscore, mismatchscore,
                                       gap_open, gap_extension);
    if (scorematrix != NULL)
    {
      gt_scorehandler_add_scorematrix(scorehandler,scorematrix);
    }
  }
  return scorehandler;
}

static void gt_encode_sequence_table(const GtAlphabet *alphabet,
                                     GtSequenceTable *sequence_table)
{
  GtUword idx;

  for (idx = 0; idx < sequence_table->size; idx++)
  {
    GtStr *gt_str = sequence_table->seqarray[idx];
    GtUchar *seq = (GtUchar*) gt_str_get(gt_str);

    gt_alphabet_encode_seq(alphabet, seq,(const char*) seq,
                           gt_str_length(gt_str));
  }
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
  GtSequenceTable *sequence_table1, *sequence_table2;
  GtLinspaceManagement *spacemanager;
  GtScoreHandler *scorehandler = NULL;
  GtTimer *linspacetimer = NULL;
  GtAlphabet *alphabet = NULL;

  gt_error_check(err);
  gt_assert(arguments);
  sequence_table1 = gt_sequence_table_new();
  sequence_table2 = gt_sequence_table_new();
  align = gt_alignment_new();
  spacemanager = gt_linspace_management_new();
  gt_linspace_management_set_TSfactor(spacemanager,arguments->timesquarefactor);

  /* get sequences */
  if (gt_str_array_size(arguments->strings) > 0)
  {
    get_onesequence(sequence_table1, arguments->strings, 0);
    sequence_table1->size++;
    get_onesequence(sequence_table2, arguments->strings, 1);
    sequence_table2->size++;
  }
  else if (gt_str_array_size(arguments->files) > 0)
  {
    had_err = get_fastasequences(sequence_table1,
                                 gt_str_array_get_str(arguments->files,0),err);
    if (!had_err)
    {
      had_err = get_fastasequences(sequence_table2,
                                  gt_str_array_get_str(arguments->files,1),err);
    }
  }
  if (arguments->dna)
  {
    alphabet = gt_alphabet_new_dna();
  } else
  {
    gt_assert(arguments->protein);
    alphabet = gt_alphabet_new_protein();
  }
  gt_encode_sequence_table(alphabet,sequence_table1);
  gt_encode_sequence_table(alphabet,sequence_table2);
  if (!had_err)
  {
    scorehandler = gt_arguments2scorehandler(arguments,err);
    if (scorehandler == NULL)
    {
      had_err = -1;
    } else
    {
      if (arguments->global && arguments->protein && !arguments->has_costmatrix)
      {
        GtScoreHandler *costhandler = gt_scorehandler2costhandler(scorehandler);
        gt_scorehandler_delete(scorehandler);
        scorehandler = costhandler;
      }
    }
  }
  /* get diagonal band */
  if (!had_err && arguments->diagonal)
  {
    if (gt_str_array_size(arguments->diagonalbonds) > 0)
    {
      had_err = gt_parse_score_value(__LINE__,&left_dist,
                                  gt_str_array_get(arguments->diagonalbonds,0),
                                  false, err);
      if (!had_err)
      {
        had_err = gt_parse_score_value(__LINE__,&right_dist,
                                  gt_str_array_get(arguments->diagonalbonds,1),
                                  false, err);
      }
    }
  }
  if (!had_err && arguments->spacetime)
  {
    linspacetimer = gt_timer_new();
  }

  /* alignment functions with linear gap costs */
  if (!had_err)
  {
    bool affine;

    if (gt_str_array_size(arguments->linearcosts) > 0)
    {
      affine = false;
    } else
    {
      gt_assert(gt_str_array_size(arguments->affinecosts) > 0);
      affine = true;
    }
    had_err = gt_all_against_all_alignment_check (
                            affine, align, arguments,
                            spacemanager,
                            scorehandler,
                            gt_alphabet_characters(alphabet),
                            gt_alphabet_wildcard_show(alphabet),
                            sequence_table1,
                            sequence_table2,
                            left_dist,
                            right_dist,
                            linspacetimer,err);
  }
  /*spacetime option*/
  if (!had_err && arguments->spacetime)
  {
    printf("# combined space peak in kilobytes: %f\n",
           GT_KILOBYTES(gt_linspace_management_get_spacepeak(spacemanager)));
    gt_timer_show_formatted(linspacetimer,"# TIME overall " GT_WD ".%02ld\n",
                            stdout);
  }
  gt_timer_delete(linspacetimer);
  gt_linspace_management_delete(spacemanager);
  gt_sequence_table_delete(sequence_table1);
  gt_sequence_table_delete(sequence_table2);
  gt_alignment_delete(align);
  gt_alphabet_delete(alphabet);
  gt_scorehandler_delete(scorehandler);
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
