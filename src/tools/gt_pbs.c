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
#include "libgtcore/strand.h"
#include "libgtcore/warning.h"
#include "libgtcore/ma.h"
#include "libgtltr/pbs.h"
#include "libgtltr/repeattypes.h"
#include "libgtltr/ltrharvest-opt.h"

typedef struct {
  unsigned int radius,
               ali_min_len,
               max_offset,
               max_offset_trna,
               max_edist;
} PBSOptions;

static OPrval parse_options(int *parsed_args, PBSOptions *opts, int argc,
                            const char **argv, Error *err)
{
  OptionParser *op;
  Option *option;
  OPrval oprval;
  error_check(err);
  op = option_parser_new("[option ...] trna_file fasta_file tabular_file",
                         "Find and score primer tRNA binding sites.");

  option = option_new_uint("radius", "radius around 5' LTR to search in",
                          &opts->radius, 30);
  option_parser_add_option(op, option);
  option = option_new_uint("aliminlen", "minimum PBS local alignment length",
                          &opts->ali_min_len, 11);
  option_parser_add_option(op, option);
  option = option_new_uint("maxoffsetltr", "maximal allowed offset from LTR boundary",
                          &opts->max_offset, 4);
  option_parser_add_option(op, option);
  option = option_new_uint("maxoffsettrna", "maximal allowed offset from tRNA 3' end",
                          &opts->max_offset_trna, 10);
  option_parser_add_option(op, option);
  option = option_new_uint("maxedist", "maximal allowed alignment edit distance",
                          &opts->max_edist, 1);
  option_parser_add_option(op, option);

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

static LTRboundaries* find_element(Seq *seq, Array *annos)
{
  unsigned long ltrstart, i;
  LTRboundaries *line = NULL;
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
  unsigned long i;
  int parsed_args, had_err = 0;
  Bioseq *trnas, *fastas;
  PBSOptions opts;
  LTRharvestoptions lo;
  Array *annos;

  /* option parsing */
  switch (parse_options(&parsed_args, &opts, argc, argv, err)) {
    case OPTIONPARSER_OK: break;
    case OPTIONPARSER_ERROR: return -1;
    case OPTIONPARSER_REQUESTS_EXIT: return 0;
  }

  lo.pbs_aliminlen = opts.ali_min_len;
  lo.pbs_maxedist = opts.max_edist;
  lo.pbs_radius = opts.radius;
  lo.pbs_maxoffset_5_ltr = opts.max_offset;
  lo.pbs_maxoffset_trna = opts.max_offset_trna;
  lo.pbs_ali_score_match = 5;
  lo.pbs_ali_score_mismatch = -10;
  lo.pbs_ali_score_insertion = -20;
  lo.pbs_ali_score_deletion = -20;

  /* get sequences from FASTA file */
  trnas = bioseq_new(argv[parsed_args], err);
  fastas = bioseq_new(argv[parsed_args+1], err);
  annos = array_new(sizeof(LTRboundaries*));
  read_ltrboundaries_from_file(argv[parsed_args+2], annos);

  for(i=0;i<bioseq_number_of_sequences(fastas);i++)
  {
    PBS_Hit *hit;
    Seq *seq = bioseq_get_seq(fastas, i);
    LTRboundaries *line = find_element(seq, annos);

    if (!line)
      continue;

    hit = pbs_find(seq_get_orig(seq),
                   line,
                   seq_length(seq),
                   trnas,
                   &lo,
                   err);

    if(hit)
    {
      printf("%lu %lu %lu %lu %lu %lu %s %lu %c\n",
                (unsigned long) line->leftLTR_5,
                (unsigned long) line->rightLTR_3,
                hit->offset,
                hit->edist,
                hit->tstart,
                hit->alilen,
                hit->trna,
                hit->start,
                STRANDCHARS[hit->strand]);
      ma_free(hit);
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
