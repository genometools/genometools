/*
  Copyright (c) 2007 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
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

#include <ctype.h>
#include <stdbool.h>
#include <inttypes.h>
#include "core/error.h"
#include "core/sequence_buffer_fasta.h"
#include "core/str.h"
#include "stamp.h"

int guessifproteinsequencestream(const GtStrArray *filenametab,GtError *err)
{
  unsigned int countnonbases = 0,
               currentposition;
  GtUchar currentchar;
  GtSequenceBuffer *fb;
  int retval;

  gt_error_check(err);
  fb = gt_sequence_buffer_new_guess_type((GtStrArray*) filenametab, err);
  if (!fb) return -1;

  for (currentposition = 0; currentposition < 1000U;
       currentposition++)
  {
    retval = gt_sequence_buffer_next(fb,&currentchar,err);
    if (retval < 0)
    {
      gt_sequence_buffer_delete(fb);
      return -1;
    }
    if (retval == 0)
    {
      break;
    }
    switch (currentchar)
    {
      case 'L':
      case 'I':
      case 'F':
      case 'E':
      case 'Q':
      case 'P':
      case 'X':
      case 'Z': countnonbases++;
                break;
      default:  break;
    }
    if (countnonbases > 0)
    {
      break;
    }
  }
  gt_sequence_buffer_delete(fb);
  if (countnonbases > 0)
  {
    return 1; /* guess it is a protein sequence */
  }
  return 0; /* guess it is a dna sequence */
}
