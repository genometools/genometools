/*
  Copyright (c) 2008 Sascha Steinbiss <ssteinbiss@zbh.uni-hamburg.de>
  Copyright (c) 2008 Center for Bioinformatics, University of Hamburg

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
#include <ctype.h>
#include <math.h>
#include "libgtcore/array.h"
#include "libgtcore/error.h"
#include "libgtcore/option.h"
#include "libgtcore/versionfunc.h"
#include "libgtcore/fileutils.h"
#include "libgtcore/warning.h"
#include "libgtcore/ma.h"
#include "libgtcore/bioseq.h"
#include "libgtext/reverse.h"
#include "libgtext/swalign.h"
#include "libgtltr/repeattypes.h"
#include "libgtltr/ltrharvest-opt.h"

#define ALILEN 30
#define ALIMINLEN 12

typedef struct {
  unsigned int alilen;
} PBSOptions;

static OPrval parse_options(int *parsed_args, PBSOptions *opts, int argc,
                            const char **argv, Error *err)
{
  OptionParser *op;
/*  Option *option; */
  OPrval oprval;
  error_check(err);
  op = option_parser_new("[option ...] trna_file fasta_file tabular_file",
                         "Find and score primer tRNA binding sites.");
  option_parser_set_mailaddress(op, "<ssteinbiss@stud.zbh.uni-hamburg.de>");
  oprval = option_parser_parse_min_max_args(op, parsed_args, argc, argv,
                                            versionfunc, 3, 3, err);

  if (oprval == OPTIONPARSER_OK && !file_exists(argv[*parsed_args]))
  {
    error_set(err, "could not open file %s", argv[*parsed_args]);
    oprval = OPTIONPARSER_ERROR;
  }
  if (oprval == OPTIONPARSER_OK && !file_exists(argv[*parsed_args+1]))
  {
    error_set(err, "could not open file %s", argv[*parsed_args+1]);
    oprval = OPTIONPARSER_ERROR;
  }
  if (oprval == OPTIONPARSER_OK && !file_exists(argv[*parsed_args+2]))
  {
    error_set(err, "could not open file %s", argv[*parsed_args+2]);
    oprval = OPTIONPARSER_ERROR;
  }
  option_parser_delete(op);
  return oprval;
}

static void read_ltrboundaries_from_file(const char *filename, Array *ltrboundaries)
{
  FILE *fp = fopen(filename,"r");
  char sline[BUFSIZ];
  while (fgets(sline, BUFSIZ, fp) != NULL)
  {
    LTRboundaries *ltr = ma_malloc(sizeof (LTRboundaries));
    unsigned long leftLTR_5=0, leftLTR_3, rightLTR_3, rightLTR_5;
    float similarity;
    sscanf(sline,"%*d %*d %*d %lu %lu %*d %lu %lu %*d %f %*d",
                                           &leftLTR_5,
                                           &leftLTR_3,
                                           &rightLTR_5,
                                           &rightLTR_3,
                                           &similarity);
   ltr->leftLTR_3 = (Seqpos) leftLTR_3;
   ltr->leftLTR_5 = (Seqpos) leftLTR_5;
   ltr->rightLTR_3 = (Seqpos) rightLTR_3;
   ltr->rightLTR_5 = (Seqpos) rightLTR_5;
   ltr->similarity = (double) similarity;
   if (leftLTR_5>0)
     array_add(ltrboundaries, ltr);
   else
     ma_free(ltr);
  }
  fclose(fp);
}

void reverse_seq(char *seq, unsigned long seqlen)
{
  char *front_char, *back_char, tmp_char;
  for (front_char = seq, back_char = seq + seqlen - 1;
       front_char <= back_char;
       front_char++, back_char--) {
    tmp_char = *back_char;
    *back_char = *front_char;
    *front_char = tmp_char;
  }
}

static ScoreFunction* dna_scorefunc_new(Alpha *a, int match, int mismatch,
                                        int insertion, int deletion)
{
  ScoreMatrix *sm = score_matrix_new(a);
  ScoreFunction *sf = scorefunction_new(sm, insertion, deletion);
  unsigned int m,n;

  for(m=0;m<5;m++)
  {
    for(n=0;n<5;n++)
    {
      score_matrix_set_score(sm, m,n,(n==m ? match : mismatch));
    }
  }
  return sf;
}

static LTRboundaries* find_element(Seq *seq, Array *annos)
{
  unsigned long ltrstart, i;
  LTRboundaries *line;
  /* identify element, quick and dirty */
  char *desc = (char*)seq_get_description(seq);
  /* search for beginning of interval */
  while (*desc != '[')
    desc++;
  sscanf(desc,"[%lu,%*d]", &ltrstart);

  /* find annotation line for element */
  for (i=0;i<array_size(annos);i++)
  {
    line  = *(LTRboundaries**) array_get(annos, i);
    if (line->leftLTR_5 == ltrstart)
      break;
  }
  return line;
}

int gt_findpbs(int argc, const char **argv, Error *err)
{
  unsigned long i,j;
  int parsed_args, had_err = 0;
  Bioseq *trnas, *fastas;
  PBSOptions opts;
  Array *annos;
  Alignment *ali;
  long dist;
  Alpha *a = (Alpha*) alpha_new_dna();
  ScoreFunction *sf = dna_scorefunc_new(a, 5, -4, -10, -10);

  /* option parsing */
  switch (parse_options(&parsed_args, &opts, argc, argv, err)) {
    case OPTIONPARSER_OK: break;
    case OPTIONPARSER_ERROR: return -1;
    case OPTIONPARSER_REQUESTS_EXIT: return 0;
  }

  /* get sequences from FASTA file */
  trnas = bioseq_new(argv[parsed_args], err);
  fastas = bioseq_new(argv[parsed_args+1], err);
  annos = array_new(sizeof(LTRboundaries*));
  read_ltrboundaries_from_file(argv[parsed_args+2], annos);

  for(i=0;i<bioseq_number_of_sequences(fastas);i++)
  {
    Seq *seq, *seq_forward, *seq_rev;
    unsigned long seqlen;
    char *seq_rev_full;
    LTRboundaries *line;
    seq = bioseq_get_seq(fastas, i);
    seqlen = seq_length(seq);

    /* get reverse complement */
    seq_rev_full = ma_malloc(sizeof(char)*seqlen);
    memcpy(seq_rev_full, seq_get_orig(seq), sizeof(char)*seqlen);
    reverse_complement(seq_rev_full, seqlen, err);

    /* get element boundaries */
    line = find_element(seq, annos);
    printf("%s\n", seq_get_description(seq));
    seq_forward = seq_new(seq_get_orig(seq)+
                            line->leftLTR_3-line->leftLTR_5-ALILEN+1,
                          2*ALILEN,
                          a);

    seq_rev = seq_new_own(seq_rev_full+
                            line->rightLTR_3-line->rightLTR_5-ALILEN+1,
                          2*ALILEN,
                          a);


    for(j=0;j<bioseq_number_of_sequences(trnas);j++)
    {
      Seq *trna_seq, *trna_rev;
      char *trna_rev_full;
      unsigned long trna_seqlen;
      Range urange, vrange;

      trna_seq = bioseq_get_seq(trnas, j);
      trna_seqlen = seq_length(trna_seq);
      trna_rev_full = ma_malloc(sizeof (char)*trna_seqlen);
      memcpy(trna_rev_full, seq_get_orig(trna_seq), sizeof(char)*trna_seqlen);
      reverse_complement(trna_rev_full, trna_seqlen, err);
      trna_rev = seq_new_own(trna_rev_full, trna_seqlen, a);

      ali = swalign(seq_forward, trna_seq, sf);
      dist = alignment_eval(ali);
      urange = alignment_get_urange(ali);
      vrange = alignment_get_vrange(ali);
      if (dist < 2 && urange.end-urange.start+1 >= ALIMINLEN)
      {
        printf("%s\n", seq_get_description(trna_seq));
        printf("%lu-%lu (%d)\n", urange.start, urange.end, abs(ALILEN-urange.start));
        alignment_show(ali, stdout);
        printf("%lu-%lu\n", vrange.start, vrange.end);
        printf("on fwd: %lu\n\n", dist);
      }
      ali = swalign(seq_rev, trna_seq, sf);
      dist = alignment_eval(ali);
      urange = alignment_get_urange(ali);
      vrange = alignment_get_vrange(ali);
      if (dist < 2 && urange.end-urange.start+1 >= ALIMINLEN)
      {
        printf("%s\n", seq_get_description(trna_seq));
        printf("%lu-%lu (%d)\n", urange.start, urange.end, abs(ALILEN-urange.start));
        alignment_show(ali, stdout);
        printf("%lu-%lu\n", vrange.start, vrange.end);
        printf("on rev: %lu\n\n", dist);
      }
    }
  }


  bioseq_delete(trnas);
  bioseq_delete(fastas);
  for (i=0;i<array_size(annos);i++)
  {
    LTRboundaries *line = *(LTRboundaries**) array_get(annos, i);
    ma_free(line);
  }
  array_delete(annos);
  return had_err;
}
