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

#ifndef CSVLINE_READER_H
#define CSVLINE_READER_H
#include <stdbool.h>
#include "core/types_api.h"
#include "core/unused_api.h"

/* Objects of the <GtCsvlineReader> represents the machinery to
   read files consisting of columns separated by a single separator
   character, aka comma separated values.
   The object does not save the lines, but uses the same
   internal structures for each line. So the user is responsible to
   save the contents appropriately. A <GtCsvlineReader>-reader reads each
   character of a line exactly once and maintains the character
   distribution for the entire line. */
typedef struct GtCsvlineReader GtCsvlineReader;
/* returns an empty reader. We do not pass a file pointer here
   since we want to be able to read different files one after the
   other, so that only the next-method gets the file pointer. */
GtCsvlineReader *gt_csvline_reader_new(void);
/* delete the reader */
void gt_csvline_reader_delete(GtCsvlineReader *csvline_reader);
/* clear the internal state, i.e. everything is set as if nothing was
   read yet. */
void gt_csvline_reader_clear(GtCsvlineReader *csvline_reader);
/* try to read the next line via the <FILE>-pointer <inputfp> and use the
   <separator> to separate the columns. Return <true> on success
   and store the line in <csvline_reader>. If there is no more
   line left, return <false>. */
bool gt_csvline_reader_next(GtCsvlineReader *csvline_reader,
                            FILE *inputfp,char separator);
/* Each column of a line is represented by the following struct.
   The content of the return column is NOT a \0-terminated string. */
typedef struct
{
  const char *content; /* content of the column */
  GtUword width;       /* number of character column consists of */
} GtCsvcolumn;

/* Return the number of columns of the current line read by
   <csvline_reader>. */
GtUword gt_csvline_reader_column_number(const GtCsvlineReader *csvline_reader);
/* Return the column with the given number <colnum>. Columns are counted
   from 0. A <colnum> not in the proper range leads to a failed assertion. */
GtCsvcolumn gt_csvline_reader_column(const GtCsvlineReader *csvline_reader,
                                     GtUword colnum);
/* Return <true> iff the current line only contains white spaces. */
bool gt_csvline_reader_white_space_line(const GtCsvlineReader *csvline_reader);
/* While the <GtCsvlineReader>-reader maintains the distribution of the
   characters for the entire line, it is often required to obtain the
   character distribution for a single column, whose number is not known
   when reading the line. In such a situation the following method can
   be used. It subtracts the contribution of the characters in all columns
   except for the column number <colnum>, from the character distribution.
   In this way, the remaining distribution is the one for the
   specified column. The original distribution is not saved.
   This method is efficient for cases where the
   chosen column is much longer than all the other columns of the line. */
void gt_csvline_reader_dist_only_for_column(GtCsvlineReader *csvline_reader,
                                            GtUword colnum);
/* We currently only provide a single method to access the character
   distribution. The following method adds the character distribution
   of the current line the to array pointed to by <charcount>. This must
   be an array one longer than the <max_character> parameter.
   So <UCHAR_MAX+1> is an appropriate size of the array. */
void gt_csvline_reader_dist_accumulate(GtUword *charcount,
                                       GT_UNUSED GtUword max_character,
                                       const GtCsvlineReader *csvline_reader);
#endif
