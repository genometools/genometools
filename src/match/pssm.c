/*
  Copyright (c) 2008 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
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
#include <stdarg.h>
#include "core/assert_api.h"
#include "core/symboldef.h"
#include "core/chardef.h"
#include "core/ma_api.h"
#include "core/arraydef.h"
#include "seqpos-def.h"
#include "spacedef.h"
#include "encodedsequence.h"
#include "absdfstrans-imp.h"

typedef int ProfScore;
#define SHOWProfScore(FP,VAL) fprintf(FP,"%d",VAL)

#define GETPROFSCORE(S,A,D)   (S)->scoretab[A][D]

typedef struct
{
  unsigned long dimension,
                numofcharacters;
  ProfScore *itmthreshold,  /* intermediate thresholds */
            maxtotalscore,  /* the maximum score */
            mintotalscore,  /* the minimum score */
            **scoretab;     /* the scores */
} Profilematrix;

void showProfilematrix(const Profilematrix *prof,
                              const GtUchar *characters)
{
  unsigned long d, a;

  printf("# %lu x %lu matrix\n",prof->numofcharacters,
                                prof->dimension);
  printf("# mintotalscore=");
  SHOWProfScore(stdout,prof->mintotalscore);
  printf("\n");
  printf("# maxtotalscore=");
  SHOWProfScore(stdout,prof->maxtotalscore);
  printf("\n");
  printf("   ");
  for (a = 0; a < prof->numofcharacters; a++)
  {
    printf("%c",(int) characters[a]);
    printf("%s",(a < prof->numofcharacters - 1) ? "   " : "\n");
  }
  for (d = 0; d < prof->dimension; d++)
  {
    for (a = 0; a < prof->numofcharacters; a++)
    {
      SHOWProfScore(stdout,GETPROFSCORE(prof,a,d));
      printf("%s",(a < prof->numofcharacters - 1) ? " " : " \n");
    }
  }
}

void makeitmthresholds(Profilematrix *prof,
                              ProfScore minscore)
{
  unsigned long d, a;
  long ddown;
  ProfScore partsum,
            score,
            *maxscore;

  ALLOCASSIGNSPACE(maxscore,NULL,ProfScore,prof->dimension);
  for (d=0; d<prof->dimension; d++)
  {
    for (a=0; a<prof->numofcharacters; a++)
    {
      score = GETPROFSCORE(prof,a,d);
      if (a == 0 || maxscore[d] < score)
      {
        maxscore[d] = score;
      }
    }
  }
  partsum = (ProfScore) 0;
  gt_assert(prof->itmthreshold != NULL);
  for (ddown = (long) (prof->dimension-1); ddown>=0; ddown--)
  {
    prof->itmthreshold[ddown] = minscore - partsum;
    partsum += maxscore[ddown];
  }
  FREESPACE(maxscore);
}

void lookaheadsearchPSSM(const GtEncodedsequence *encseq,
                                const Profilematrix *prof)
{
  unsigned long firstpos, bufsize;
  GtUchar currentchar;
  Seqpos pos;
  GtEncodedsequenceScanstate *esr;
  Seqpos totallength = gt_encodedsequence_total_length(encseq);
  GtUchar *buffer;

  esr = gt_encodedsequence_scanstate_new();
  gt_encodedsequence_scanstate_init(esr,encseq,Forwardmode,0);
  ALLOCASSIGNSPACE(buffer,NULL,GtUchar,prof->dimension);
  firstpos = bufsize = 0;
  for (pos=0; pos < totallength; pos++)
  {
    currentchar = gt_encodedsequence_sequentialgetencodedchar(encseq,esr,pos,
                                                              Forwardmode);
    if (ISSPECIAL(currentchar))
    {
      bufsize = firstpos = 0;
    } else
    {
      if (bufsize < prof->dimension)
      {
        buffer[bufsize++] = currentchar;
      } else
      {
        buffer[firstpos++] = currentchar;
        if (firstpos == prof->dimension)
        {
          firstpos = 0;
        }
      }
    }
  }
  gt_encodedsequence_scanstate_delete(esr);
  FREESPACE(buffer);
}
