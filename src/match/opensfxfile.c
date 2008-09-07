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
#include "core/error.h"
#include "core/fa.h"
#include "core/str.h"

/*@null@*/ FILE *opensfxfile(const GT_Str *indexname,
                             const char *suffix,
                             const char *mode,
                             GT_Error *err)
{
  GT_Str *tmpfilename;
  FILE *fp;

  gt_error_check(err);
  tmpfilename = gt_str_clone(indexname);
  gt_str_append_cstr(tmpfilename,suffix);
  fp = fa_fopen(gt_str_get(tmpfilename),mode,err);
  gt_str_delete(tmpfilename);
  return fp;
}

bool indexfilealreadyexists(const GT_Str *indexname,const char *suffix)
{
  struct stat statbuf;
  GT_Str *tmpfilename;

  tmpfilename = gt_str_clone(indexname);
  gt_str_append_cstr(tmpfilename,suffix);

  if (stat(gt_str_get(tmpfilename),&statbuf) == 0)
  {
    gt_str_delete(tmpfilename);
    return true;
  }
  gt_str_delete(tmpfilename);
  return false;
}
