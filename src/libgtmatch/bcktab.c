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
#include <errno.h>
#include "libgtcore/error.h"
#include "libgtcore/str.h"
#include "libgtcore/fa.h"
#include "esafileend.h"
#include "spacedef.h"
#include "bckbound.h"

#include "initbasepower.pr"

/* DELETE THE FOLLOWING THREE FUNCTIONS */

static void *genericmaponlytable(const Str *indexname,const char *suffix,
                                 size_t *numofbytes,Error *err)
{
  Str *tmpfilename;
  void *ptr;
  bool haserr = false;

  error_check(err);
  tmpfilename = str_clone(indexname);
  str_append_cstr(tmpfilename,suffix);
  ptr = fa_mmap_read(str_get(tmpfilename),numofbytes);
  if (ptr == NULL)
  {
    error_set(err,"cannot map file \"%s\": %s",str_get(tmpfilename),
                  strerror(errno));
    haserr = true;
  }
  str_delete(tmpfilename);
  return haserr ? NULL : ptr;
}

static int checkmappedfilesize(size_t numofbytes,Seqpos expectedunits,
                               size_t sizeofunit,Error *err)
{
  error_check(err);
  if (expectedunits != (Seqpos) (numofbytes/sizeofunit))
  {
    error_set(err,"number of mapped units = %lu != " FormatSeqpos
                      " = expected number of integers",
                      (unsigned long) (numofbytes/sizeofunit),
                      PRINTSeqposcast(expectedunits));
    return -1;
  }
  return 0;
}

static void *genericmaptable(const Str *indexname,
                             const char *suffix,
                             Seqpos expectedunits,size_t sizeofunit,
                             Error *err)
{
  size_t numofbytes;

  void *ptr = genericmaponlytable(indexname,suffix,&numofbytes,err);
  if (ptr == NULL)
  {
    return NULL;
  }
  if (checkmappedfilesize(numofbytes,expectedunits,sizeofunit,err) != 0)
  {
    fa_xmunmap(ptr);
    return NULL;
  }
  return ptr;
}

int mapbcktab(Bcktab *bcktab,
              const Str *indexname,
              unsigned int numofchars,
              unsigned int prefixlength,
              Error *err)
{
  Codetype numofspecialcodes;

  bcktab->basepower = initbasepower(numofchars,prefixlength);
  bcktab->numofallcodes = bcktab->basepower[prefixlength];
  numofspecialcodes = bcktab->basepower[prefixlength-1];
  bcktab->filltable = initfilltable(bcktab->basepower,prefixlength);
  bcktab->multimappower = initmultimappower(numofchars,prefixlength);
  bcktab->leftborder
    = genericmaptable(indexname,
                      BCKTABSUFFIX,
                      (Seqpos) (bcktab->numofallcodes + 1 + numofspecialcodes),
                      sizeof (Seqpos),
                      err);
  if (bcktab->leftborder == NULL)
  {
    bcktab->countspecialcodes = NULL;
    return -1;
  }
  bcktab->countspecialcodes = bcktab->leftborder + bcktab->numofallcodes + 1;
  return 0;
}

void freebcktab(Bcktab *bcktab)
{
  fa_xmunmap((void *) bcktab->leftborder);
  bcktab->leftborder = NULL;
  bcktab->countspecialcodes = NULL;
  multimappowerfree(&bcktab->multimappower);
  FREESPACE(bcktab->filltable);
  FREESPACE(bcktab->basepower);
}

void initbcktabwithNULL(Bcktab *bcktab)
{
  bcktab->leftborder = NULL;
  bcktab->countspecialcodes = NULL;
  bcktab->multimappower = NULL;
  bcktab->filltable = NULL;
  bcktab->basepower = NULL;
}
