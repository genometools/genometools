/*
  Copyright (c) 2007 David Ellinghaus <d.ellinghaus@ikmb.uni-kiel.de>
  Copyright (c) 2007 Center for Bioinformatics, University of Hamburg

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

#include "repeattypes.h"

/*
 The following function removes exact duplicates from the array of
 predicted LTR elements. Exact duplicates occur, when
 different seeds are extended to same boundary coordinates.
 */
void removeduplicates(ArrayLTRboundaries *arrayLTRboundaries)
{
  unsigned long i, j;
  Seqpos startpos_i, endpos_i, startpos_j, endpos_j;
  LTRboundaries *boundaries_i,
                *boundaries_j;

  for (i = 0; i < arrayLTRboundaries->nextfreeLTRboundaries; i++)
  {
    boundaries_i = &(arrayLTRboundaries->spaceLTRboundaries[i]);
    if ( boundaries_i->skipped )
    {
      continue;
    }
    startpos_i = boundaries_i->leftLTR_5;
    endpos_i   = boundaries_i->rightLTR_3;

    for (j = i + 1; j < arrayLTRboundaries->nextfreeLTRboundaries; j++)
    {
      boundaries_j = &(arrayLTRboundaries->spaceLTRboundaries[j]);
      if ( boundaries_j->skipped )
      {
        continue;
      }
      startpos_j = boundaries_j->leftLTR_5;
      endpos_j   = boundaries_j->rightLTR_3;

      if (startpos_i == startpos_j && endpos_i == endpos_j)
      {
        boundaries_j->skipped = true;
      }
    }
  }
}

/* The following function removes overlaps and deletes the prediction with
   a lower similarity value. If "nooverlapallowed" is set, then all
   overlapping predicitions are deleted completely.
 */
void removeoverlapswithlowersimilarity(
  ArrayLTRboundaries *arrayLTRboundaries,
  bool nooverlapallowed
  )
{
  unsigned long i, j;
  Seqpos startpos_i, endpos_i, startpos_j, endpos_j;
  LTRboundaries *boundaries_i,
                *boundaries_j;

  for (i = 0; i < arrayLTRboundaries->nextfreeLTRboundaries; i++)
  {
    boundaries_i = &(arrayLTRboundaries->spaceLTRboundaries[i]);
    if ( boundaries_i->skipped )
    {
      continue;
    }
    startpos_i = boundaries_i->leftLTR_5;
    endpos_i   = boundaries_i->rightLTR_3;

    for (j = i + 1; j < arrayLTRboundaries->nextfreeLTRboundaries; j++)
    {
      boundaries_j = &(arrayLTRboundaries->spaceLTRboundaries[j]);
      if ( boundaries_j->skipped )
      {
        continue;
      }
      startpos_j = boundaries_j->leftLTR_5;
      endpos_j   = boundaries_j->rightLTR_3;

      /* if overlap */
      if ( !((endpos_i < startpos_j) || (endpos_j < startpos_i)) )
      {
        if (nooverlapallowed)
        {
          /* All predictions in a cluster will be deleted. */

          /* take min(startpos_i, startpos_j) */
          if (startpos_j < startpos_i)
          {
            startpos_i = startpos_j;
          }
          /* take max(endpos_i, endpos_j) */
          if (endpos_i < endpos_j)
          {
            endpos_i = endpos_j;
          }
          /* delete both predictions */
          boundaries_i->skipped = true;
          boundaries_j->skipped = true;
        }
        else
        {
          /* take prediction with higher similarity */
          if ( boundaries_i->similarity >= boundaries_j->similarity )
          {
            boundaries_j->skipped = true;
          }
          else
          {
            boundaries_i->skipped = true;
            break;
          }
        }
      }
    }
  }
}
