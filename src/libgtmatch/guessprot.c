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
#include "libgtcore/error.h"
#include "libgtcore/fastabuffer.h"
#include "libgtcore/str.h"
#include "stamp.h"

int guessifproteinsequencestream(const StrArray *filenametab,Error *err)
{
  unsigned int countnonbases = 0,
               currentposition;
  Uchar currentchar;
  FastaBuffer *fb;
  int retval;

  error_check(err);
  fb = fastabuffer_new(filenametab,
                       NULL,
                       false,
                       NULL,
                       NULL,
                       NULL);
  for (currentposition = 0; currentposition < 1000U;
       currentposition++)
  {
    retval = fastabuffer_next(fb,&currentchar,err);
    if (retval < 0)
    {
      fastabuffer_delete(fb);
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
  }
  fastabuffer_delete(fb);
  if (countnonbases > 0 && countnonbases >= currentposition/10)
  {
    return 1;
  }
  return 0;
}
