/*
  Copyright (c) 2016 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
  Copyright (c) 2016 Center for Bioinformatics, University of Hamburg

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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <limits.h>
#include "core/arraydef_api.h"
#include "core/csvline_reader.h"

struct GtCsvlineReader
{
  GtArraychar line;
  GtArrayGtUword columnoffset;
  struct
  {
    GtArraychar alphabet;
    GtUword *charcount;
  } dist;
  bool empty;
  char separator;
};

GtCsvlineReader *gt_csvline_reader_new(void)
{
  GtCsvlineReader *csvline_reader = gt_malloc(sizeof *csvline_reader);
  gt_assert(csvline_reader != NULL);
  GT_INITARRAY(&csvline_reader->line,char);
  GT_INITARRAY(&csvline_reader->dist.alphabet,char);
  GT_INITARRAY(&csvline_reader->columnoffset,GtUword);
  csvline_reader->empty = true;
  csvline_reader->separator = 0;
  csvline_reader->dist.charcount
    = gt_calloc(UCHAR_MAX+1,sizeof *csvline_reader->dist.charcount);
  return csvline_reader;
}

static void csvline_reader_append_char(GtCsvlineReader *csvline_reader,char cc)
{
  GT_STOREINARRAY(&csvline_reader->line,char,1024UL,cc);
  if (csvline_reader->dist.charcount[(int) cc] == 0)
  {
    GT_STOREINARRAY(&csvline_reader->dist.alphabet,char,256UL,cc);
  }
  csvline_reader->dist.charcount[(int) cc]++;
}

void gt_csvline_reader_delete(GtCsvlineReader *csvline_reader)
{
  if (csvline_reader != NULL)
  {
    GT_FREEARRAY(&csvline_reader->columnoffset,GtUword);
    GT_FREEARRAY(&csvline_reader->line,char);
    GT_FREEARRAY(&csvline_reader->dist.alphabet,char);
    gt_free(csvline_reader->dist.charcount);
    gt_free(csvline_reader);
  }
}

void gt_csvline_reader_clear(GtCsvlineReader *csvline_reader)
{
  GtUword idx;

  csvline_reader->line.nextfreechar = 0;
  for (idx = 0; idx < csvline_reader->dist.alphabet.nextfreechar; idx++)
  {
    char cc = csvline_reader->dist.alphabet.spacechar[idx];
    csvline_reader->dist.charcount[(int) cc] = 0;
  }
  csvline_reader->dist.alphabet.nextfreechar = 0;
  for (idx = 0; idx <= UCHAR_MAX; idx++)
  {
    gt_assert(csvline_reader->dist.charcount[idx] == 0);
  }
  csvline_reader->columnoffset.nextfreeGtUword = 0;
  csvline_reader->empty = true;
}

bool gt_csvline_reader_next(GtCsvlineReader *csvline_reader,
                            FILE *inputfp,char separator)
{
  gt_assert(csvline_reader != NULL);
  gt_csvline_reader_clear(csvline_reader);
  csvline_reader->separator = separator;
  while (true)
  {
    int cc = fgetc(inputfp);

    if (cc == EOF)
    {
      break;
    }
    if (cc == '\n')
    {
      csvline_reader_append_char(csvline_reader,'\0');
      return true;
    }
    if (cc == csvline_reader->separator)
    {
      csvline_reader_append_char(csvline_reader,cc);
      GT_CHECKARRAYSPACE(&csvline_reader->columnoffset,GtUword,1024);
      csvline_reader->columnoffset.spaceGtUword[csvline_reader->
                                                columnoffset.nextfreeGtUword++]
        = csvline_reader->line.nextfreechar;
    } else
    {
      if (csvline_reader->empty && !isspace(cc))
      {
        csvline_reader->empty = false;
      }
      csvline_reader_append_char(csvline_reader,cc);
    }
  }
  return false;
}

bool gt_csvline_reader_white_space_line(const GtCsvlineReader *csvline_reader)
{
  return csvline_reader->empty;
}

GtCsvcolumn gt_csvline_reader_column(const GtCsvlineReader *csvline_reader,
                                     GtUword colnum)
{
  GtCsvcolumn col;

  gt_assert(csvline_reader != NULL &&
            !gt_csvline_reader_white_space_line(csvline_reader));
  if (colnum == 0)
  {
    gt_assert(csvline_reader->line.nextfreechar >= 2UL);
    col.content = csvline_reader->line.spacechar;
  } else
  {
    gt_assert(colnum - 1 < csvline_reader->columnoffset.nextfreeGtUword);
    col.content = csvline_reader->line.spacechar +
                  csvline_reader->columnoffset.spaceGtUword[colnum-1];
  }
  if (colnum == 0)
  {
    col.width = csvline_reader->columnoffset.spaceGtUword[colnum] - 1;
  } else
  {
    if (colnum == csvline_reader->columnoffset.nextfreeGtUword)
    {
      col.width = csvline_reader->line.nextfreechar -
                  csvline_reader->columnoffset.spaceGtUword[colnum-1] - 1;
    } else
    {
      col.width = csvline_reader->columnoffset.spaceGtUword[colnum] -
                  csvline_reader->columnoffset.spaceGtUword[colnum-1] - 1;
    }
  }
  return col;
}

void gt_csvline_reader_check(const GtCsvlineReader *csvline_reader)
{
  GtUword idx, bfdist[UCHAR_MAX+1] = {0};

  for (idx = 0; idx < csvline_reader->line.nextfreechar; idx++)
  {
    bfdist[(int) csvline_reader->line.spacechar[idx]]++;
  }
  for (idx = 0; idx <= UCHAR_MAX; idx++)
  {
    if (bfdist[idx] != csvline_reader->dist.charcount[idx])
    {
      fprintf(stderr,"%s\nidx=" GT_WU ",bfdist=" GT_WU " != "
                                                 GT_WU " = chardist\n",
                csvline_reader->line.spacechar,idx,bfdist[idx],
                csvline_reader->dist.charcount[idx]);
      exit(GT_EXIT_PROGRAMMING_ERROR);
    }
  }
}

GtUword gt_csvline_reader_column_number(const GtCsvlineReader *csvline_reader)
{
  if (gt_csvline_reader_white_space_line(csvline_reader))
  {
    return 0;
  }
  return csvline_reader->columnoffset.nextfreeGtUword + 1;
}

void gt_csvline_reader_dist_only_for_column(GtCsvlineReader *csvline_reader,
                                            GtUword colnum)
{
  GtUword idx, write_idx,
          numofcols = gt_csvline_reader_column_number(csvline_reader);

  gt_assert(numofcols > 0 && colnum < numofcols &&
            csvline_reader->dist.charcount[(int) csvline_reader->separator]
            == numofcols - 1);
  csvline_reader->dist.charcount[(int) csvline_reader->separator] = 0;
  gt_assert(csvline_reader->dist.charcount[(int) '\0'] == 1);
  csvline_reader->dist.charcount[(int) '\0'] = 0;
  for (idx = 0; idx < numofcols; idx++)
  {
    if (idx != colnum)
    {
      GtCsvcolumn col = gt_csvline_reader_column(csvline_reader,idx);
      GtUword j;

      for (j = 0; j < col.width; j++)
      {
        char currentcc = col.content[j];
        gt_assert(csvline_reader->dist.charcount[(int) currentcc] > 0);
        csvline_reader->dist.charcount[(int) currentcc]--;
      }
    }
  }
  for (idx = 0, write_idx = 0;
       idx < csvline_reader->dist.alphabet.nextfreechar; idx++)
  {
    char currentcc = csvline_reader->dist.alphabet.spacechar[idx];
    if (csvline_reader->dist.charcount[(int) currentcc] > 0)
    {
      if (write_idx < idx)
      {
        csvline_reader->dist.alphabet.spacechar[write_idx] = currentcc;
      }
      write_idx++;
    }
  }
  csvline_reader->dist.alphabet.nextfreechar = write_idx;
}

#ifndef NDEBUG
void gt_csvline_reader_dist_check(const GtCsvlineReader *csvline_reader,
                                  const char *string,GtUword len)
{
  GtUword idx, numchars = 0, dist_tab[UCHAR_MAX+1] = {0};

  for (idx = 0; idx < len; idx++)
  {
    dist_tab[(int) string[idx]]++;
  }
  for (idx = 0; idx <= UCHAR_MAX; idx++)
  {
    gt_assert(dist_tab[idx] == csvline_reader->dist.charcount[idx]);
    if (csvline_reader->dist.charcount[idx] > 0)
    {
      numchars++;
    }
  }
  gt_assert(numchars == csvline_reader->dist.alphabet.nextfreechar);
  for (idx = 0; idx < numchars; idx++)
  {
    char cc = csvline_reader->dist.alphabet.spacechar[idx];
    gt_assert(csvline_reader->dist.charcount[(int) cc] > 0);
  }
}
#endif

void gt_csvline_reader_dist_show(const GtCsvlineReader *csvline_reader)
{
  GtUword idx;

  for (idx = 0; idx < csvline_reader->dist.alphabet.nextfreechar; idx++)
  {
    char cc = csvline_reader->dist.alphabet.spacechar[idx];
    printf("%c/" GT_WU,cc,csvline_reader->dist.charcount[(int) cc]);
    if (idx + 1 < csvline_reader->dist.alphabet.nextfreechar)
    {
      printf("%c",csvline_reader->separator);
    } else
    {
      printf("\n");
    }
  }
}

void gt_csvline_reader_dist_accumulate(GtUword *charcount,
                                       GT_UNUSED GtUword max_character,
                                       const GtCsvlineReader *csvline_reader)
{
  GtUword idx;

  for (idx = 0; idx < csvline_reader->dist.alphabet.nextfreechar; idx++)
  {
    char cc = csvline_reader->dist.alphabet.spacechar[idx];
    gt_assert(cc <= max_character);
    charcount[(int) cc] += csvline_reader->dist.charcount[(int) cc];
  }
}
