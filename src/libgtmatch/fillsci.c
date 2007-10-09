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

#include <stdbool.h>
#include <ctype.h>
#include "libgtcore/chardef.h"
#include "libgtcore/discdistri.h"
#include "libgtcore/env.h"
#include "libgtcore/fastabuffer.h"
#include "libgtcore/strarray.h"
#include "seqpos-def.h"
#include "safecast-gen.h"

static unsigned long currentrangevalue(unsigned long i,unsigned long distvalue)
{
  if (i <= UCHAR_MAX)
  {
    return distvalue;
  }
  if (i % UCHAR_MAX == 0)
  {
    return i/UCHAR_MAX * distvalue;
  }
  return (((unsigned long) 1) + i/UCHAR_MAX) * distvalue;
}

static void updatesumranges(unsigned long key, unsigned long long value,
                            void *data)
{
  unsigned long distvalue,
                *specialrangesptr = (unsigned long *) data;

  distvalue = (unsigned long) value; /* XXX: is this cast always OK? */
  (*specialrangesptr) += currentrangevalue(key,distvalue);
  /* printf("# specialranges of length %lu: %lu\n",keyvalue,distvalue);
     XXX integrate later */
}

 DECLARESAFECASTFUNCTION(Seqpos,Seqpos,unsigned long,unsigned_long)

int fasta2sequencekeyvalues(
        unsigned long *numofsequences,
        Seqpos *totallength,
        Specialcharinfo *specialcharinfo,
        const StrArray *filenametab,
        Filelengthvalues **filelengthtab,
        const Uchar *symbolmap,
        bool plainformat,
        unsigned long *characterdistribution,
        Env *env)
{
  FastaBuffer *fb;
  Uchar charcode;
  Seqpos pos;
  int retval;
  bool specialprefix = true;
  Seqpos lastspeciallength = 0;
  DiscDistri *specialrangelengths;
  unsigned long idx;
  bool haserr = false;

  env_error_check(env);
  *numofsequences = 0;
  specialcharinfo->specialcharacters = 0;
  specialcharinfo->lengthofspecialprefix = 0;
  specialcharinfo->lengthofspecialsuffix = 0;

  fb = fastabuffer_new(filenametab,
                       symbolmap,
                       plainformat,
                       filelengthtab,
                       NULL,
                       characterdistribution,
                       env);
  specialrangelengths = discdistri_new(env);
  for (pos = 0; /* Nothing */; pos++)
  {
    retval = fastabuffer_next(fb,&charcode,env);
    if (retval < 0)
    {
      haserr = true;
      break;
    }
    if (retval == 0)
    {
      if (lastspeciallength > 0)
      {
        idx = CALLCASTFUNC(Seqpos,unsigned_long,lastspeciallength);
        discdistri_add(specialrangelengths,idx,env);
      }
      break;
    }
    if (ISSPECIAL(charcode))
    {
      if (specialprefix)
      {
        specialcharinfo->lengthofspecialprefix++;
      }
      specialcharinfo->specialcharacters++;
      if (lastspeciallength == 0)
      {
        lastspeciallength = (Seqpos) 1;
      } else
      {
        lastspeciallength++;
      }
      if (charcode == (Uchar) SEPARATOR)
      {
        (*numofsequences)++;
      }
    } else
    {
      if (specialprefix)
      {
        specialprefix = false;
      }
      if (lastspeciallength > 0)
      {
        idx = CALLCASTFUNC(Seqpos,unsigned_long,lastspeciallength);
        discdistri_add(specialrangelengths,idx,env);
        lastspeciallength = 0;
      }
    }
  }
  if (!haserr)
  {
    specialcharinfo->specialranges = 0;
    discdistri_foreach(specialrangelengths,updatesumranges,
                       &specialcharinfo->specialranges,env);
    specialcharinfo->lengthofspecialsuffix = lastspeciallength;
    (*numofsequences)++;
    *totallength = pos;
  }
  discdistri_delete(specialrangelengths,env);
  fastabuffer_delete(fb, env);
  return haserr ? -1 : 0;
}

void sequence2specialcharinfo(Specialcharinfo *specialcharinfo,
                              const Uchar *seq,
                              const Seqpos len,
                              Env *env)
{
  Uchar charcode;
  Seqpos pos;
  bool specialprefix = true;
  Seqpos lastspeciallength = 0;
  DiscDistri *specialrangelengths;
  unsigned long idx;

  env_error_check(env);
  specialcharinfo->specialcharacters = 0;
  specialcharinfo->lengthofspecialprefix = 0;
  specialcharinfo->lengthofspecialsuffix = 0;
  specialrangelengths = discdistri_new(env);
  for (pos = 0; pos < len; pos++)
  {
    charcode = seq[pos];
    if (ISSPECIAL(charcode))
    {
      if (specialprefix)
      {
        specialcharinfo->lengthofspecialprefix++;
      }
      specialcharinfo->specialcharacters++;
      if (lastspeciallength == 0)
      {
        lastspeciallength = (Seqpos) 1;
      } else
      {
        lastspeciallength++;
      }
    } else
    {
      if (specialprefix)
      {
        specialprefix = false;
      }
      if (lastspeciallength > 0)
      {
        idx = CALLCASTFUNC(Seqpos,unsigned_long,lastspeciallength);
        discdistri_add(specialrangelengths,idx,env);
        lastspeciallength = 0;
      }
    }
  }
  if (lastspeciallength > 0)
  {
    idx = CALLCASTFUNC(Seqpos,unsigned_long,lastspeciallength);
    discdistri_add(specialrangelengths,idx,env);
  }
  specialcharinfo->specialranges = 0;
  discdistri_foreach(specialrangelengths,updatesumranges,
                     &specialcharinfo->specialranges,env);
  specialcharinfo->lengthofspecialsuffix = lastspeciallength;
  discdistri_delete(specialrangelengths,env);
}
