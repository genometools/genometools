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
#include "opensfxfile.h"

/*@null@*/ FILE *opensfxfile(const GtStr *indexname,
                             const char *suffix,
                             const char *mode,
                             GtError *err)
{
  GtStr *tmpfilename;
  FILE *fp;

  gt_error_check(err);
  tmpfilename = gt_str_clone(indexname);
  gt_str_append_cstr(tmpfilename,suffix);
  fp = gt_fa_fopen(gt_str_get(tmpfilename),mode,err);
  gt_str_delete(tmpfilename);
  return fp;
}

bool indexfilealreadyexists(const GtStr *indexname,const char *suffix)
{
  struct stat statbuf;
  GtStr *tmpfilename;

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

void *genericmaponlytable(const GtStr *indexname,const char *suffix,
                          size_t *numofbytes,GtError *err)
{
  GtStr *tmpfilename;
  void *ptr;
  bool haserr = false;

  gt_error_check(err);
  tmpfilename = gt_str_clone(indexname);
  gt_str_append_cstr(tmpfilename,suffix);
  ptr = gt_fa_mmap_read(gt_str_get(tmpfilename),numofbytes);
  if (ptr == NULL)
  {
    gt_error_set(err,"cannot map file \"%s\": %s",gt_str_get(tmpfilename),
                  strerror(errno));
    haserr = true;
  }
  gt_str_delete(tmpfilename);
  return haserr ? NULL : ptr;
}

static int checkmappedfilesize(const GtStr *indexname,
                               const char *suffix,
                               size_t numofbytes,
                               unsigned long expectedunits,
                               size_t sizeofunit,
                               GtError *err)
{
  gt_error_check(err);
  if (expectedunits != (unsigned long) (numofbytes/sizeofunit))
  {
    gt_error_set(err,"mapping file %s%s: number of mapped units (of size %u) "
                     " = %lu != %lu = expected number of mapped units",
                      gt_str_get(indexname),
                      suffix,
                      (unsigned int) sizeofunit,
                      (unsigned long) (numofbytes/sizeofunit),
                      expectedunits);
    return -1;
  }
  return 0;
}

void *genericmaptable(const GtStr *indexname,
                      const char *suffix,
                      unsigned long expectedunits,
                      size_t sizeofunit,
                      GtError *err)
{
  size_t numofbytes;

  void *ptr = genericmaponlytable(indexname,suffix,&numofbytes,err);
  if (ptr == NULL)
  {
    return NULL;
  }
  if (checkmappedfilesize(indexname,suffix,
                          numofbytes,expectedunits,sizeofunit,err) != 0)
  {
    gt_fa_xmunmap(ptr);
    return NULL;
  }
  return ptr;
}
