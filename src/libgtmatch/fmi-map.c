/*
  Copyright (c) 2007 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
  Copyright (c) 2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#include <stdio.h>
#include <errno.h>
#include <stdbool.h>
#include "libgtcore/env.h"
#include "libgtcore/array.h"
#include "libgtcore/str.h"
#include "types.h"
#include "arraydef.h"
#include "sfx-ri-def.h"
#include "fmindex.h"

#include "readnextline.pr"
#include "opensfxfile.pr"
#include "sfx-readint.pr"

int scanfmprjfileviafileptr(Fmindex *fmindex,
                                   Seqpos *totallength,
                                   const Str *indexname,FILE *fpin,Env *env)
{
  ArrayUchar linebuffer;
  bool haserr = false;
  Array *riktab;
  Seqpos bwtlength;
  uint32_t linenum,
           storeindexpos,
           log2blocksize,
           log2markdist,
           suffixlength;

  env_error_check(env);
  riktab = array_new(sizeofReadintkeys(),env);
  SETREADINTKEYS("bwtlength",&bwtlength,NULL);
  SETREADINTKEYS("longest",&fmindex->longestsuffixpos,NULL);
  SETREADINTKEYS("storeindexpos",&storeindexpos,NULL);
  SETREADINTKEYS("log2blocksize",&log2blocksize,NULL);
  SETREADINTKEYS("log2markdist",&log2markdist,NULL);
  SETREADINTKEYS("specialcharacters",
                 &fmindex->specialcharinfo.specialcharacters,NULL);
  SETREADINTKEYS("specialranges",
                 &fmindex->specialcharinfo.specialranges,NULL);
  SETREADINTKEYS("lengthofspecialprefix",
                 &fmindex->specialcharinfo.lengthofspecialprefix,NULL);
  SETREADINTKEYS("lengthofspecialsuffix",
                 &fmindex->specialcharinfo.lengthofspecialsuffix,NULL);
  SETREADINTKEYS("suffixlength",&suffixlength,NULL);

  INITARRAY(&linebuffer,Uchar);
  for (linenum = 0; /* Nothing */; linenum++)
  {
    linebuffer.nextfreeUchar = 0;
    if (readnextline(fpin,&linebuffer,env) == EOF)
    {
      break;
    }
    if(analyzeuintline(indexname,
                       linenum, 
                       linebuffer.spaceUchar,
                       linebuffer.nextfreeUchar,
                       riktab,
                       env) != 0)
    {
      haserr = true;
      break;
    }
  }
  if (!haserr && allkeysdefined(indexname,riktab,env) != 0)
  {
    haserr = true;
  }
  FREEARRAY(&linebuffer,Uchar);
  array_delete(riktab,env);
  return haserr ? -1 : 0;
}
