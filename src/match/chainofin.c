/*
  Copyright (c) 2009 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
  Copyright (c) 2009 Center for Bioinformatics, University of Hamburg

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

#ifndef S_SPLINT_S
#include <stdbool.h>
#endif
#include "core/fa.h"
#include "core/error_api.h"
#include "chain2dim.h"

#define READNUMS 5

#define CANNOTPARSELINE(S)\
        gt_error_set(err,"matchfile \"%s\", line %lu, column %lu: %s",\
                     matchfile,linenum+1,countcolumns+1,S)

static int numberoflinesinfile(unsigned long *linenum,
                               const char *filename,
                               GtError *err)
{
  FILE *fp;
  int cc;

  fp = gt_fa_fopen(filename,"r",err);
  *linenum = 0;
  if (fp == NULL)
  {
    return -1;
  }
  while ((cc = getc(fp)) != EOF)
  {
    if (cc == '\n')
    {
      (*linenum)++;
    }
  }
  gt_fa_fclose(fp);
  return 0;
}

GtChain2Dimmatchtable *gt_chain_analyzeopenformatfile(double weightfactor,
                                                  const char *matchfile,
                                                  GtError *err)
{
  GtChain2Dimmatchtable *matchtable;
  unsigned long linenum;
  long storeinteger[READNUMS];
  FILE *matchfp;
  bool haserr = false;
  GtChain2Dimmatchvalues fragment;

  if (numberoflinesinfile(&linenum,matchfile,err) != 0)
  {
    return NULL;
  }
  matchfp = gt_fa_fopen(matchfile,"r",err);
  if (matchfp == NULL)
  {
    return NULL;
  }
  matchtable = gt_chain_matchtable_new(linenum);
  for (linenum = 0; fscanf(matchfp,"%ld %ld %ld %ld %ld\n",
                           &storeinteger[0],
                           &storeinteger[1],
                           &storeinteger[2],
                           &storeinteger[3],
                           &storeinteger[4]) == READNUMS; linenum++)
  {
    unsigned long countcolumns;

    for (countcolumns = 0; countcolumns < (unsigned long) (READNUMS-1);
         countcolumns++)
    {
      if (storeinteger[countcolumns] < 0)
      {
        CANNOTPARSELINE("non-negative integer expected");
        haserr = true;
      }
    }
    if (storeinteger[0] > storeinteger[1])
    {
      CANNOTPARSELINE("startpos1 <= endpos1 expected");
      haserr = true;
      break;
    }
    if (storeinteger[2] > storeinteger[3])
    {
      CANNOTPARSELINE("startpos2 <= endpos2 expected");
      haserr = true;
      break;
    }
    fragment.startpos[0] = (GtChain2Dimpostype) storeinteger[0];
    fragment.endpos[0] = (GtChain2Dimpostype) storeinteger[1];
    fragment.startpos[1] = (GtChain2Dimpostype) storeinteger[2];
    fragment.endpos[1] = (GtChain2Dimpostype) storeinteger[3];
    fragment.weight
      = (GtChain2Dimscoretype) (weightfactor * (double) storeinteger[4]);
    gt_chain_matchtable_add(matchtable,&fragment);
    /*gt_chain_printchainelem(stdout,&fragment); */
  }
  gt_fa_fclose(matchfp);
  if (haserr)
  {
    gt_chain_matchtable_delete(matchtable);
    return NULL;
  }
  gt_chain_fillthegapvalues(matchtable);
  return matchtable;
}
