/*
  Copyright (c) 2007      Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
  Copyright (c) 2007-2008 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2007-2008 Center for Bioinformatics, University of Hamburg

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

#include <sys/types.h>
#include <sys/stat.h>
#include "libgtcore/fileutils.h"
#include "libgtcore/ma.h"
#include "libgtcore/option.h"
#include "libgtcore/seqiterator.h"
#include "libgtcore/versionfunc.h"
#include "libgtcore/xposix.h"
#include "libgtcore/progressbar.h"
#include "format64.h"

static const char *desc2ginumber(unsigned long *ginumlen,const char *desc,
                                 Error *err)
{
  unsigned long i, firstpipe = 0, secondpipe = 0;

  for (i=0; desc[i] != '\0'; i++)
  {
    if (desc[i] == '|')
    {
      if (firstpipe > 0)
      {
        assert(i>0);
        secondpipe = i;
        break;
      }
      assert(i>0);
      firstpipe = i;
    }
  }
  if (firstpipe == 0 || secondpipe == 0)
  {
    error_set(err,"Cannot find gi-number in description \"%s\"\n",desc);
    return NULL;
  }
  assert(firstpipe < secondpipe);
  *ginumlen = firstpipe - secondpipe - 1;
  return desc + firstpipe + 1;
}

int extractginumbers(bool verbose,
                     /*@unused@*/ const Str *ginumberfile,
                     int argc,const char **argv,Error *err)
{
  StrArray *files;
  SeqIterator *seqit;
  const Uchar *sequence;
  char *desc;
  unsigned long len, ginumlen;
  int i, had_err;
  off_t totalsize;

  error_check(err);

  files = strarray_new();
  for (i = 0; i < argc; i++)
  {
    printf("add file %s\n",argv[i]);
    strarray_add_cstr(files, argv[i]);
  }
  totalsize = files_estimate_total_size(files);
  printf("# estimated total size is " Formatuint64_t "\n",
            PRINTuint64_tcast(totalsize));
  seqit = seqiterator_new(files, NULL, true);
  if (verbose)
  {
    progressbar_start(seqiterator_getcurrentcounter(seqit, (unsigned long long)
                                                           totalsize),
                                                           (unsigned long long)
                                                           totalsize);
  }
  while (true)
  {
    had_err = seqiterator_next(seqit, &sequence, &len, &desc, err);
    if (had_err != 1)
    {
      break;
    }
    if (desc2ginumber(&ginumlen,desc,err) == NULL)
    {
      had_err = -1;
    }
    ma_free(desc);
  }
  if (verbose)
  {
    progressbar_stop();
  }
  seqiterator_delete(seqit);
  strarray_delete(files);
  return had_err;
}
