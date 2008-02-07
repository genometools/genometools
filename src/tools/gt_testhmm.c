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
#include <math.h>
#include "libgtcore/array.h"
#include "libgtcore/error.h"
#include "libgtcore/option.h"
#include "libgtcore/ma.h"
#include "libgtcore/bioseq.h"
#include "libgtltr/ppt.h"

int gt_testhmm(int argc, const char **argv, Error *err)
{
  unsigned int i;
  unsigned long seqindex;
  Bioseq *bs;
  long nofseq;
  Array *ltrboundaries;
  LTRharvestoptions options;

  /* set options */
  options.ppt_minlen = 6;
  options.ubox_minlen = 3;
  options.ppt_radius = 50;

  /* get sequences from FASTA file */
  bs = bioseq_new(argv[1], err);
  nofseq = bioseq_number_of_sequences(bs);
  printf("%li sequences loaded.\n", nofseq);

  /* read annotations from file */
  ltrboundaries = array_new(sizeof(LTRboundaries*));
  FILE *fp;
  char sline[BUFSIZ];
  fp = fopen(argv[2],"r");
  while(fgets(sline, BUFSIZ, fp) != NULL)
  {
    LTRboundaries *ltr = ma_malloc(sizeof(LTRboundaries));
    unsigned long leftLTR_5, leftLTR_3, rightLTR_3, rightLTR_5;
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
   array_add(ltrboundaries, ltr);
  }
  fclose(fp);
  printf("%lu annotations loaded.\n", array_size(ltrboundaries));

  for(seqindex=0UL;seqindex<nofseq;seqindex++)
  {
    Seq *seq;
    unsigned long seqlen;
    char *seqc, *encoded_seq_c;
    unsigned long ltrstart=0UL;
    LTRboundaries *line = NULL;
    const Alpha *alpha = alpha_new_dna();
    char out[BUFSIZ] = "\0" ;  /* TODO: make this overflow-safe, IF NEEDED*/
    char *desc = NULL;

    /* get element */
    seq = bioseq_get_seq(bs, seqindex);
    encoded_seq_c = (char*) seq_get_encoded(seq);
    seqlen = seq_length(seq);

    /* identify element, quick and dirty */
    desc = (char*)seq_get_description(seq);
    /* search for beginning of interval */
    while(*desc != '[')
      desc++;
    sscanf(desc,"[%lu,%*d]", &ltrstart);

    /* find annotation line for element */
    for(i=0;i<array_size(ltrboundaries);i++)
    {
      line  = *(LTRboundaries**) array_get(ltrboundaries, i);
      if(line->leftLTR_5 == ltrstart)
        break;
    }
    printf("\nsequence %lu [%lu,%lu] (%lu/%lubp):\n", seqindex,
                        (unsigned long) line->leftLTR_5,
                        (unsigned long) line->rightLTR_3,
                        (unsigned long) line->rightLTR_3-
                           (unsigned long) line->leftLTR_5+1,
                        seqlen);

    seqc = ma_malloc(sizeof(char) * seqlen+1);
    alpha_decode_seq(alpha, seqc, encoded_seq_c, seqlen);
    seqc[seqlen]='\0';
    /* print reference sequence */
    char *prseq = ma_malloc(sizeof(char) * 2*options.ppt_radius+1);
    strncpy(prseq,
            seqc+(seqlen-(line->rightLTR_3-line->rightLTR_5)+1-options.ppt_radius),
            2*options.ppt_radius);
    prseq[2*options.ppt_radius] = '\0';
    printf("%s\n", prseq);

    for(i=0;i<options.ppt_radius;i++)
    {
      printf(" ");
    }
    printf("|\n");

    Array *results = array_new(sizeof(PPT_Hit*));

    ppt_find(line, results, options.ppt_radius, seqc, &options, ppt_score);

    for(i=0;i<array_size(results);i++)
    {
      unsigned long j;
      PPT_Hit *hit = *(PPT_Hit**) array_get(results, i);
      for(j=0;j<hit->end-hit->start+1;j++)
      {
        switch(hit->state)
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
      if (hit->state == PPT_IN)
      {
        if (i>0)
        {
          PPT_Hit *uhit = *(PPT_Hit**) array_get(results, i-1);
          if (uhit->state == PPT_UBOX)
          {
            snprintf(out, BUFSIZ, "hit: PPT length %lu"
                                  ", U-box length %lu score %f",
                                hit->end-hit->start+1,
                                uhit->end-uhit->start+1,
                                hit->score);
          } else {
            snprintf(out, BUFSIZ, "hit: PPT length %lu, score %f",
                                  hit->end-hit->start+1,
                                  hit->score);
          }
        } else {
          snprintf(out, BUFSIZ, "hit: PPT length %lu, score %f",
                                hit->end-hit->start+1,
                                hit->score);
        }
      }
    }
    printf("\n%s\n\n", out);

    /* free stuff */

    for(i=0;i<array_size(results);i++)
    {
      ma_free(*(PPT_Hit**) array_get(results,i));
    }
    ma_free(seqc);
    ma_free(prseq);
    alpha_delete((Alpha*)alpha);
    array_delete(results);
  }

  /* free annotations */
  for(i=0;i<array_size(ltrboundaries);i++)
  {
    LTRboundaries *line = *(LTRboundaries**) array_get(ltrboundaries, i);
    ma_free(line);
  }

  array_delete(ltrboundaries);
  bioseq_delete(bs);

  return 0;
}
