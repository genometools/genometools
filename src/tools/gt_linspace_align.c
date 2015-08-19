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
#include "extended/linearalign_affinegapcost.h"
#include "extended/linearalign.h"
#include "tools/gt_linspace_align.h"

typedef struct {
  GtStr *outputfile;
  GtStrArray *strings,
             *files,
             *linearcosts,
             *affinecosts;
  bool global,
       local;
} GtLinspaceArguments;

typedef struct {
  GtUchar* seq;
  GtUword len;
} Fastaentry;

typedef struct {
  Fastaentry *seqarray;
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
      gt_free(sequences->seqarray[i].seq);
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
  *optionlinearcosts, *optionaffinecosts, *optionoutputfile;
  gt_assert(arguments);

  /* init */
  op = gt_option_parser_new("options",
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

  optionoutputfile = gt_option_new_string("o", "use outputfile",
                                          arguments->outputfile, "stdout");
  gt_option_parser_add_option(op, optionoutputfile);

  /* dependencies*/
  gt_option_is_mandatory_either(optionstrings, optionfiles);
  gt_option_exclude(optionlocal, optionglobal);
  gt_option_exclude(optionlinearcosts, optionaffinecosts);
  gt_option_imply_either_2(optionstrings, optionglobal, optionlocal);
  gt_option_imply_either_2(optionfiles, optionglobal, optionlocal);
  gt_option_imply_either_2(optionlocal, optionlinearcosts, optionaffinecosts);
  gt_option_imply_either_2(optionglobal, optionlinearcosts, optionaffinecosts);

  return op;
}

static int gt_linspace_align_arguments_check(GT_UNUSED int rest_argc,
                                       void *tool_arguments,
                                       GT_UNUSED GtError *err)
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

  return had_err;
}

static void print_sequence(const GtUchar *seq, const GtUword len, FILE *fp)
{
  GtUword i = 0;
  fprintf(fp, "######\n");
  do{
    fprintf(fp, "%.80s\n",seq+i);
    i += 80;
  }while (i < len);
}

static void show(const GtUchar *useq, const GtUword ulen,
                 const GtUchar *vseq, const GtUword vlen,
                 const GtAlignment *align, FILE *fp)
{
  if (fp != NULL)
  {
    /*fprintf(fp,"# two sequences \"%s\" \"%s\"\n", useq, vseq);*/
    print_sequence(useq, ulen, fp);
    print_sequence(vseq, vlen, fp);
    fprintf(fp, "######\n");

    gt_alignment_show(align, fp, true, 80);
    gt_alignment_eval_with_score(align, 0,1,1);
  }

}

static GtWord* select_costs(const GtStrArray *arr,
                            GT_UNUSED GtError *err)
{
  GtWord *evalues, check;
  GtUword size, i;

  size = gt_str_array_size(arr);
  evalues = gt_malloc(sizeof(*evalues)*size);
  for (i = 0; i < size; i++)
  {
    check = sscanf(gt_str_array_get(arr,i),GT_WD, &evalues[i]);
    if (check != 1)
    {
      gt_error_set(err, "find invalid cost or score");
      gt_free(evalues);
      return NULL;
    }
  }

  return evalues;
}

static int save_fastasequence(const char *seqpart, GtUword length,
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
  fasta_seqs->seqarray[fasta_seqs->size].seq
    = gt_malloc(sizeof(char)*(length+1));
  memcpy(fasta_seqs->seqarray[fasta_seqs->size].seq, seqpart, length+1);
  fasta_seqs->seqarray[fasta_seqs->size].len = length;
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

static int gt_linspace_align_runner(GT_UNUSED int argc,
                                 GT_UNUSED const char **argv,
                                 GT_UNUSED int parsed_args,
                                 void *tool_arguments,
                                 GtError *err)
{
  GtLinspaceArguments *arguments = tool_arguments;
  int had_err = 0;
  const GtUchar *useq, *vseq;
  GtUword i, j, ulen, vlen;
  GtWord *linearcosts, *affinecosts;
  GtAlignment *align = NULL;
  FILE *fp;

  GtSequences *sequences1, *sequences2;
  sequences1 = gt_sequences_new();
  sequences2 = gt_sequences_new();

  gt_error_check(err);
  gt_assert(arguments);

  if (gt_str_array_size(arguments->strings) > 0)
  {
    sequences1->seqarray[0].seq =
                    (GtUchar *) gt_str_array_get(arguments->strings,0);
    sequences1->seqarray[0].len =
               (GtUword) strlen(gt_str_array_get(arguments->strings,0));
    sequences1->size++;

    sequences2->seqarray[0].seq =
                    (GtUchar *) gt_str_array_get(arguments->strings,1UL);
    sequences2->seqarray[0].len =
               (GtUword) strlen(gt_str_array_get(arguments->strings,1UL));
    sequences2->size++;
  }
  else if (gt_str_array_size(arguments->files) > 0)
  {
    had_err = get_fastasequences(sequences1,
                                 gt_str_array_get_str(arguments->files,0),
                                 err);
    if (had_err)
    {
      gt_sequences_delete(sequences1);
      gt_sequences_delete(sequences2);
      return had_err;
    }
    had_err = get_fastasequences(sequences2,
                                 gt_str_array_get_str(arguments->files,1),
                                 err);
    if (had_err)
    {
      gt_sequences_delete(sequences1);
      gt_sequences_delete(sequences2);
      return had_err;
    }
  }

  /* call functions */
  for (i = 0; i < sequences1->size; i++) {
    for (j = 0; j< sequences2->size; j++) {
      useq = sequences1->seqarray[i].seq;
      ulen = sequences1->seqarray[i].len;
      vseq = sequences2->seqarray[j].seq;
      vlen = sequences2->seqarray[j].len;

      /* linear gap costs */
      if (gt_str_array_size(arguments->linearcosts) > 0)
      {
        gt_assert(gt_str_array_size(arguments->linearcosts) == 3UL);
        linearcosts = select_costs(arguments->linearcosts, err);
        if (linearcosts == NULL)
          return 1;

        if (arguments->global)
        {
          align = gt_computelinearspace(useq, 0, ulen, vseq, 0, vlen,
                  linearcosts[0],linearcosts[1],linearcosts[2]);
        }
        else if (arguments->local)
        {
          align = gt_computelinearspace_local(useq, 0, ulen, vseq, 0, vlen,
                  linearcosts[0],linearcosts[1],linearcosts[2]);
        }

        gt_free(linearcosts);
      }/* affine gap costs */
      else if (gt_str_array_size(arguments->affinecosts) > 0)
      {
        gt_assert(gt_str_array_size(arguments->affinecosts) == 4UL);
        affinecosts = select_costs(arguments->affinecosts, err);
        if (affinecosts == NULL)
          return 1;

        if (arguments->global)
        {
          align = gt_computeaffinelinearspace(useq, 0, ulen, vseq, 0, vlen,
                  affinecosts[0],affinecosts[1],affinecosts[2], affinecosts[3]);
        }
        else if (arguments->local)
        {
          align = gt_computeaffinelinearspace_local(useq, 0, ulen, vseq, 0,
                                                    vlen, affinecosts[0],
                                                    affinecosts[1],
                                                    affinecosts[2],
                                                    affinecosts[3]);
        }
         gt_free(affinecosts);
      }

      /* show */
      gt_assert(align != NULL);
      if (!strcmp(gt_str_get(arguments->outputfile),"stdout"))
        show(useq, ulen, vseq, vlen, align, stdout);
      else
      {
        fp = gt_fa_fopen_func(gt_str_get(arguments->outputfile),
                              "a", __FILE__,__LINE__,err);
        gt_error_check(err);
        show(useq, ulen, vseq, vlen, align, fp);
        gt_fa_fclose(fp);
      }
      gt_alignment_delete(align);
    }
  }
  gt_sequences_delete(sequences1);
  gt_sequences_delete(sequences2);

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
