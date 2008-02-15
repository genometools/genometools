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

#include <errno.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include "libgtcore/error.h"
#include "libgtcore/fa.h"
#include "libgtcore/str.h"

/*@null@*/ FILE *opensfxfile(const Str *indexname,
                             const char *suffix,
                             const char *mode,
                             Error *err)
{
  Str *tmpfilename;
  FILE *fp;

  error_check(err);
  tmpfilename = str_clone(indexname);
  str_append_cstr(tmpfilename,suffix);
  fp = fa_fopen(str_get(tmpfilename),mode);
  if (fp == NULL)
  {
    error_set(err,"fa_fopen: cannot open file \"%s\": %s",
                  str_get(tmpfilename),
                  strerror(errno));
  }
  str_delete(tmpfilename);
  return fp;
}

bool indexfilealreadyexists(const Str *indexname,const char *suffix)
{
  struct stat statbuf;
  Str *tmpfilename;

  tmpfilename = str_clone(indexname);
  str_append_cstr(tmpfilename,suffix);

  if (stat(str_get(tmpfilename),&statbuf) == 0)
  {
    str_delete(tmpfilename);
    return true;
  }
  str_delete(tmpfilename);
  return false;
}
