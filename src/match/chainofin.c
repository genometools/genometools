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
#include <ctype.h>
#include "core/fa.h"
#include "core/error_api.h"
#include "core/str_api.h"
#include "chain2dim.h"

#define READNUMS 5

#define CANNOTPARSELINE(S)\
        gt_error_set(err,"matchfile \"%s\", line %lu, column %lu: %s",\
                     matchfile,linenum+1,countcolumns+1,S)

static int numberoflinesinfile(unsigned long *linenum,
                               const char *filename,GtError *err)
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

GtFragmentinfotable *gt_chain_analyzeopenformatfile(double weightfactor,
                                                    const char *matchfile,
                                                    GtError *err)
{
  GtFragmentinfotable *fragmentinfotable;
  GtStr *currentline;
  unsigned long linenum;
  GtChainpostype storeinteger[READNUMS];
  FILE *matchfp;
  long readint;
  bool haserr = false;
  GtChainscoretype weight;

  if (numberoflinesinfile(&linenum,matchfile,err) != 0)
  {
    return NULL;
  }
  matchfp = gt_fa_fopen(matchfile,"r",err);
  if (matchfp == NULL)
  {
    return NULL;
  }
  fragmentinfotable = gt_chain_fragmentinfotable_new(linenum);
  currentline = gt_str_new();
  for (linenum = 0; gt_str_read_next_line(currentline, matchfp) != EOF;
       linenum++)
  {
    const char *matchline = (const char *) gt_str_get(currentline);
    gt_assert(matchline != NULL);
    if (matchline[0] != '#')
    {
      unsigned long idx = 0,
                    countcolumns = 0,
                    matchlinelength = gt_str_length(currentline);
      while (!haserr && countcolumns < (unsigned long) READNUMS &&
             idx < matchlinelength)
      {
        while (isspace((int) matchline[idx]))
        {
          idx++;
        }
        if (sscanf(matchline+idx,"%ld",&readint) == 1 && readint >= 0)
        {
          storeinteger[countcolumns] = (GtChainpostype) readint;
        } else
        {
          CANNOTPARSELINE("non-negative integer expected");
          haserr = true;
        }
        if (!haserr)
        {
          while (isdigit((int) matchline[idx]))
          {
            idx++;
          }
          countcolumns++;
        }
      }
      if (haserr)
      {
        break;
      }
      if (countcolumns != (unsigned long) READNUMS)
      {
        CANNOTPARSELINE("five integers per line expected");
        haserr = true;
        break;
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
      weight = (GtChainscoretype) (weightfactor *
                                  (double) storeinteger[READNUMS-1]);
      gt_chain_fragmentinfotable_add(fragmentinfotable,
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
  if (haserr)
  {
    gt_chain_fragmentinfotable_delete(fragmentinfotable);
    return NULL;
  }
  gt_chain_fillthegapvalues(fragmentinfotable);
  return fragmentinfotable;
}
