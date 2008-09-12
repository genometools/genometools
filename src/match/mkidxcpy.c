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
#include <inttypes.h>
#include "core/error.h"
#include "core/fa.h"
#include "core/str.h"

#include "opensfxfile.pr"

int makeindexfilecopy(const GtStr *destindex,
                      const GtStr *sourceindex,
                      const char *suffix,
                      uint64_t maxlength,
                      GtError *err)
{
  FILE *fpdest = NULL, *fpsource = NULL;
  int cc;
  bool haserr = false;

  gt_error_check(err);
  fpdest = opensfxfile(destindex,suffix,"wb",err);
  if (fpdest == NULL)
  {
    haserr = true;
  }
  if (!haserr)
  {
    fpsource = opensfxfile(sourceindex,suffix,"rb",err);
    if (fpsource == NULL)
    {
      haserr = true;
    }
  }
  printf("# cp %s%s %s%s\n",
           gt_str_get(sourceindex),suffix,gt_str_get(destindex),suffix);
  if (!haserr)
  {
    if (maxlength == 0)
    {
      while ((cc = fgetc(fpsource)) != EOF)
      {
        (void) putc(cc,fpdest);
      }
    } else
    {
      uint64_t pos;

      for (pos = 0; pos < maxlength; pos++)
      {
        if ((cc = fgetc(fpsource)) == EOF)
        {
          break;
        }
        (void) putc(cc,fpdest);
      }
    }
  }
  gt_xfclose(fpdest);
  gt_xfclose(fpsource);
  return haserr ? -1 : 0;
}
