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

#include <string.h>
#include "libgtcore/error.h"
#include "readmode-def.h"

static char *readmodes[] = {"fwd",
                            "rev",
                            "cpl",
                            "rcl"};

const char *showreadmode(Readmode readmode)
{
  return readmodes[(int) readmode];
}

int parsereadmode(const char *dirargstring,Error *err)
{
  size_t i;

  error_check(err);
  for (i=0; i<sizeof (readmodes)/sizeof (readmodes[0]); i++)
  {
    if (strcmp(dirargstring,readmodes[i]) == 0)
    {
      return (int) i;
    }
  }
  error_set(err,"argument to option -dir must be fwd or rev or cpl or rcl");
  return -1;
}
