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
#include "libgtltr/ppt.h"

/* number of nucleotides to print around PPT hits in CSV output */
#define CSV_SEQ_SPACING 10

typedef struct {
  unsigned int ppt_minlen,
               ubox_minlen,
               radius;
  bool csv;
} HMMTestOptions;

static OPrval parse_options(int *parsed_args, HMMTestOptions *opts, int argc,
                            const char **argv, Error *err)
{
  OptionParser *op;
  Option *option;
  OPrval oprval;
  error_check(err);
  op = option_parser_new("[option ...] fasta_file tabular_file",
                         "Find and score polypurine tracts in LTRharvest "
                         "predictions.");
  option = option_new_uint("pptminlen", "minimum length for PPTs",
                          &opts->ppt_minlen, 6);
  option_parser_add_option(op, option);
  option = option_new_uint("uboxminlen", "minimum length for U-boxes",
                          &opts->ubox_minlen, 3);
  option_parser_add_option(op, option);
  option = option_new_uint("radius", "radius around beginning of 3' LTR "
                                    "to search in",
                          &opts->radius, 40);
  option_parser_add_option(op, option);
  option = option_new_bool("csv", "print output in CSV format",
                          &opts->csv, 0);
  option_parser_add_option(op, option);
  option_parser_set_mailaddress(op, "<ssteinbiss@stud.zbh.uni-hamburg.de>");
  oprval = option_parser_parse_min_max_args(op, parsed_args, argc, argv,
                                            versionfunc, 2, 2, err);

  if (oprval == OPTIONPARSER_OK && !file_exists(argv[*parsed_args]))
  {
    error_set(err, "could not open FASTA file %s", argv[*parsed_args]);
    oprval = OPTIONPARSER_ERROR;
  }
  if (oprval == OPTIONPARSER_OK && !file_exists(argv[*parsed_args+1]))
  {
    error_set(err, "could not open annotation file %s", argv[*parsed_args+1]);
    oprval = OPTIONPARSER_ERROR;
  }
  option_parser_delete(op);
  return oprval;
}

static void print_str_upper(char* str)
{
  register int t;

  for (t=0; str[t]; ++t)  {
    putchar(toupper(str[t]));
  }
}

static void output_tabular(Array* results, LTRboundaries *b,
                       unsigned long startpos,
                       unsigned int index, const char* seqoffset,
                       LTRharvestoptions options)
{
  unsigned long i;
  /* print reference sequence */
  char *prseq = ma_malloc(sizeof (char) * 2*options.ppt_radius+1);
  strncpy(prseq,
          seqoffset-1,
          2*options.ppt_radius);
  prseq[2*options.ppt_radius] = '\0';
  printf("%s\n", prseq);

  for (i=0;i<options.ppt_radius;i++)
  {
    printf(" ");
  }
  printf("|\n");

  for (i=0;i<array_size(results);i++)
  {
    unsigned long j;
    PPT_Hit *hit = *(PPT_Hit**) array_get(results, i);
    for (j=0;j<hit->end-hit->start+1;j++)
    {
      switch (hit->state)
      {
        case PPT_OUT:
          printf("-");
          break;
        case PPT_IN:
          printf("P");
          break;
        case PPT_UBOX:
          printf("U");
          break;
        case PPT_NOF_STATES:
        default:
          break;
      }
    }
  }
  ma_free(prseq);
  printf("\n");
}

static void output_csv(PPT_Hit *showhit, LTRboundaries *b,
                       unsigned long startpos,
                       const char* seqoffset)
{
  if (showhit)
  {
    unsigned long s, e;
    char *seq, tmp[CSV_SEQ_SPACING+1];

    tmp[CSV_SEQ_SPACING]='\0';
    printf("%lu;%lu;", (unsigned long) b->leftLTR_5,
                       (unsigned long) b->rightLTR_3);
    switch (showhit->strand)
    {
      case STRAND_FORWARD:
        printf("+;"); break;
      case STRAND_REVERSE:
        printf("-;"); break;
      case STRAND_UNKNOWN:
      case STRAND_BOTH:
      case NUM_OF_STRAND_TYPES:
        return;
    }
    printf("%lu;", startpos+showhit->start);
    printf("%lu;", showhit->end-showhit->start+1);
    s = showhit->start;
    e = showhit->end;
    if (showhit->ubox)
    {
      printf("%lu;", startpos+showhit->ubox->start);
      printf("%lu;",
            showhit->ubox->end-showhit->ubox->start+1);
      s = showhit->ubox->start;
    }
    else
      printf(";0;");

    strncpy(tmp, seqoffset+s-1-CSV_SEQ_SPACING, CSV_SEQ_SPACING);
    printf("%s", tmp);

    seq = ma_malloc(sizeof (char) * e-s+2);
    strncpy(seq, seqoffset+s-1, e-s+1);
    seq[e-s+1]='\0';
    print_str_upper(seq);

    strncpy(tmp, seqoffset+e, CSV_SEQ_SPACING);
    printf("%s;%f;\n", tmp,showhit->score);

    ma_free(seq);
  };
}

int gt_findppt(int argc, const char **argv, Error *err)
{
  unsigned int i;
  unsigned long seqindex;
  int parsed_args, had_err = 0;
  Bioseq *bs;
  HMMTestOptions opts;
  long nofseq;
  Array *ltrboundaries;
  LTRharvestoptions options;
  FILE *fp=NULL;
  char sline[BUFSIZ];

  /* option parsing */
  switch (parse_options(&parsed_args, &opts, argc, argv, err)) {
    case OPTIONPARSER_OK: break;
    case OPTIONPARSER_ERROR: return -1;
    case OPTIONPARSER_REQUESTS_EXIT: return 0;
  }

  options.ppt_minlen = opts.ppt_minlen;
  options.ubox_minlen = opts.ubox_minlen;
  options.ppt_radius = opts.radius;

  /* get sequences from FASTA file */
  bs = bioseq_new(argv[parsed_args], err);

  if (!had_err)
  {
    nofseq = bioseq_number_of_sequences(bs);
    fprintf(stderr,"%li sequences loaded.\n", nofseq);

    /* read annotations from file */
    ltrboundaries = array_new(sizeof (LTRboundaries*));
    fp = fopen(argv[parsed_args+1],"r");

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

    fprintf(stderr, "%lu annotations loaded.\n", array_size(ltrboundaries));
    for (seqindex=0UL;seqindex<nofseq;seqindex++)
    {
      Seq *seq;
      unsigned long seqlen;
      const char *seqc, *seqoffset;
      char *seqc_m, *seqoffset_m;
      unsigned long ltrstart=0UL, idx_p=0UL, idx_m=0UL;
      LTRboundaries *line = NULL;
      const Alpha *alpha = alpha_new_dna();
      char *desc = NULL;
      Array *results_p, *results_m;

      /* get element */
      seq = bioseq_get_seq(bs, seqindex);
      seqlen = seq_length(seq);

      /* identify element, quick and dirty */
      desc = (char*)seq_get_description(seq);
      /* search for beginning of interval */
      while (*desc != '[')
        desc++;
      sscanf(desc,"[%lu,%*d]", &ltrstart);

      /* find annotation line for element */
      for (i=0;i<array_size(ltrboundaries);i++)
      {
        line  = *(LTRboundaries**) array_get(ltrboundaries, i);
        if (line->leftLTR_5 == ltrstart)
          break;
      }
      /* if no annotation could be found, skip this prediction */
      if (line)
      {
        if (!opts.csv)
          printf("\nsequence #%lu [%lu,%lu] (%lu/%lubp):\n", seqindex,
                              (unsigned long) line->leftLTR_5,
                              (unsigned long) line->rightLTR_3,
                              (unsigned long) line->rightLTR_3-
                                 (unsigned long) line->leftLTR_5+1,
                              seqlen);

        /* prepare sequences */
        seqc = seq_get_orig(seq);
        seqc_m = ma_malloc(sizeof (char) * seqlen);
        memcpy(seqc_m, seqc, (sizeof (char) * seqlen));
        reverse_complement(seqc_m, seqlen, err);

        /* precalculate search offsets in sequences */
        seqoffset = seqc+(seqlen-(line->rightLTR_3-line->rightLTR_5)
                     -options.ppt_radius);
        seqoffset_m = seqc_m+(seqlen-(line->leftLTR_3-line->leftLTR_5)
                     -options.ppt_radius);

        /* run PPT identification on both strands */
        results_p = array_new(sizeof (PPT_Hit*));
        results_m = array_new(sizeof (PPT_Hit*));
        idx_p = ppt_find(seqc,
                         seqlen,
                         line->rightLTR_3-line->rightLTR_5+1,
                         results_p,
                         &options,
                         ppt_score,
                         STRAND_FORWARD);
        idx_m = ppt_find(seqc_m,
                         seqlen,
                         line->leftLTR_3-line->leftLTR_5+1,
                         results_m,
                         &options,
                         ppt_score,
                         STRAND_REVERSE);
        double pscore = 0.0, mscore = 0.0;
        PPT_Hit *maxhit_p = NULL, *maxhit_m = NULL;
        if (idx_p != UNDEF_ULONG)
        {
          maxhit_p = *(PPT_Hit**) array_get(results_p, idx_p);
          pscore = maxhit_p->score;
        }
        if (idx_m != UNDEF_ULONG)
        {
          maxhit_m = *(PPT_Hit**) array_get(results_m, idx_m);
          mscore = maxhit_m->score;
        }
        if (pscore > mscore)
        {
          if (opts.csv)
          output_csv(maxhit_p,
                  line,
                  line->leftLTR_5+seqlen-(line->rightLTR_3-line->rightLTR_5)
                    -options.ppt_radius,
                  seqoffset);
          else
          output_tabular(results_p,
                  line,
                  line->leftLTR_5+seqlen-(line->rightLTR_3-line->rightLTR_5)
                    -options.ppt_radius,
                  idx_p,
                  seqoffset,
                  options);
        }
        else
        {
          if (opts.csv)
          output_csv(maxhit_m,
                  line,
                  line->leftLTR_5+(line->leftLTR_3-line->leftLTR_5),
                  seqoffset_m);
          else
          output_tabular(results_m,
                  line,
                  line->leftLTR_5+seqlen-(line->rightLTR_3-line->rightLTR_5)
                    -options.ppt_radius,
                  idx_m,
                  seqoffset_m,
                  options);
        }

        /* free stuff */
        for (i=0;i<array_size(results_p);i++)
        {
          ma_free(*(PPT_Hit**) array_get(results_p,i));
        }
        for (i=0;i<array_size(results_m);i++)
        {
          ma_free(*(PPT_Hit**) array_get(results_m,i));
        }
        alpha_delete((Alpha*)alpha);
        array_delete(results_p);
        array_delete(results_m);
        ma_free(seqc_m);
      }
      else
        warning("LTR transposon found with no matching annotation!");
    }
    /* free annotations */
    for (i=0;i<array_size(ltrboundaries);i++)
    {
      LTRboundaries *line = *(LTRboundaries**) array_get(ltrboundaries, i);
      ma_free(line);
    }

  }

  array_delete(ltrboundaries);
  bioseq_delete(bs);
  return had_err;
}
