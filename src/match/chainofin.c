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

#include <stdbool.h>
#include "chaindef.h"
#include "core/fa.h"
#include "core/error_api.h"
#include "core/str_api.h"
#include "safecast-gen.h"

#define READNUMS 5

typedef struct
{
  bool silent;
  GtFragmentinfotable *fragementinfotable;
  unsigned long chaincounter;
} Ofchainoutinfo;

#define CANNOTPARSELINE(S)\
        gt_error_set("matchfile \"%s\", line %lu, column %lu: %s",\
                     matchfile,linenum,countcolumns,S)

DECLARESAFECASTFUNCTION(size_t,size_t,unsigned long,unsigned_long)

static size_t numberoflinesinfile(const char *filename,GtError *err)
{
  size_t numberoflinesinfile = 0;
  FILE *fp;
  int cc;

  fp = gt_fa_fopen(filename,"r",err);
  if (fp == NULL)
  {
    return -1;
  }
  while((cc = getc(fp)) != EOF)
  {
    if (cc == '\n')
    {
      numberoflinesinfile++;
    }
  }
  return CALLCASTFUNC(size_t,unsigned_long,numberoflinesinfile)\
}

static GtFragmentinfotable *analyzeopenformatfile(double weightfactor,
                                                  const char *matchfile,
                                                  GtError *err)
{
  GtFragmentinfotable *finfotab;
  GtStr *currentline;
  char *matchline;
  unsigned long linenum = 0, countcolumns;
  GtChainpostype storeinteger[READNUMS];
  FILE *matchfp;
  long readint;
  bool haserr = false;
  Chainscoretype weight;

  matchfp = gt_fa_fopen(matchfile,"r",err);
  if (matchfp == NULL)
  {
    return NULL;
  }
  finfotab = fragmentinfotable_new(numberoflinesinfile(matchfile));
  currentline = gt_str_new();
  for (linenum = 0; gt_str_read_next_line(currentline, matchfp) != EOF;
       linenum++)
  {
    unsigned long idx = 0;
    matchline = gt_str_get(currenline);
    matchlinelength = gt_str_length(currentline);
    gt_assert(matchline != NULL);
    idx = 0;
    if(matchline[idx] != '#')
    {
      countcolumns = 0;
      while(!haserr && countcolumns < (unsigned long) READNUMS &&
            idx < matchlinelength)
      {
        while(isspace((int) matchline[idx]))
        {
          idx++;
        }
        if (sscanf(matchline,"%ld",&readint) == 1 || readint < 0)
        {
          storeinteger[countcolumns] = (GtChainpostype) readint;
        } else
        {
          CANNOTPARSE("cannot read positive integer");
          haserr = true;
        }
        if (!haserr)
        {
          while(isdigit((Ctypeargumenttype) *matchline))
          {
            matchline++;
          }
          countcolumns++;
        }
      }
      if (haserr)
      {
        break;
      }
      if(countcolumns != READNUMS)
      {
        CANNOTPARSE("not enough integers: there must be exactly five integers");
        haserr = true;
        break;
      }
      if(storeinteger[0] > storeinteger[1])
      {
        CANNOTPARSELINE("startpos1 > endpos1");
        haserr = true;
        break;
      }
      if(storeinteger[2] > storeinteger[3])
      {
        CANNOTPARSELINE("startpos1 > endpos1");
        haserr = true;
        break;
      }
      weight = (Chainscoretype) (weightfactor *
                                 (double) storeinteger[READNUMS-1]);
      fragmentinfotable_add(fragmentinfotable,
                            storeinteger[0],
                            storeinteger[1],
                            storeinteger[2],
                            storeinteger[3],
                            weight);
    }
    gt_str_reset(currentline);
  }
  gt_str_delete(currentline);
  gt_fa_fclose(matchfp);
  return haserr ? -1 : 0;
}
