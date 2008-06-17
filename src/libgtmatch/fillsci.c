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
#include "libgtcore/disc_distri.h"
#include "libgtcore/error.h"
#include "libgtcore/fa.h"
#include "libgtcore/fastabuffer.h"
#include "libgtcore/strarray.h"
#include "esafileend.h"
#include "seqpos-def.h"
#include "verbose-def.h"
#include "spacedef.h"
#include "safecast-gen.h"

#include "opensfxfile.pr"

static unsigned long currentrangevalue(unsigned long len,
                                       unsigned long occcount)
{
  if (len <= (unsigned long) (UCHAR_MAX+1))
  {
    return occcount;
  }
  if (len % (UCHAR_MAX+1) == 0)
  {
    return len/(UCHAR_MAX+1) * occcount;
  }
  return (1UL + len/(UCHAR_MAX+1)) * occcount;
}

typedef struct
{
  Verboseinfo *verboseinfo;
  Seqpos *specialrangesptr,
         *realspecialrangesptr;
} Updatesumrangeinfo;

static void updatesumranges(unsigned long key, unsigned long long value,
                            void *data)
{
  unsigned long distvalue;
  Updatesumrangeinfo *updatesumrangeinfo = (Updatesumrangeinfo *) data;

  distvalue = (unsigned long) value; /* XXX: is this cast always OK? */
  (*updatesumrangeinfo->specialrangesptr) += currentrangevalue(key,distvalue);
  (*updatesumrangeinfo->realspecialrangesptr) += (unsigned long) distvalue;
  showverbose(updatesumrangeinfo->verboseinfo,
              "specialranges of length %lu: %lu",key,distvalue);
}

 DECLARESAFECASTFUNCTION(Seqpos,Seqpos,unsigned long,unsigned_long)

int fasta2sequencekeyvalues(
        const Str *indexname,
        unsigned long *numofsequences,
        Seqpos *totallength,
        Specialcharinfo *specialcharinfo,
        const StrArray *filenametab,
        Filelengthvalues **filelengthtab,
        const Uchar *symbolmap,
        bool plainformat,
        bool withdestab,
        unsigned long *characterdistribution,
        Verboseinfo *verboseinfo,
        Error *err)
{
  FastaBuffer *fb = NULL;
  Uchar charcode;
  Seqpos pos = 0;
  int retval;
  bool specialprefix = true;
  Seqpos lastspeciallength = 0;
  DiscDistri *distspralen = NULL;
  unsigned long idx;
  bool haserr = false;
  Queue *descqueue = NULL;
  char *desc;
  FILE *desfp = NULL;

  error_check(err);
  *numofsequences = 0;
  specialcharinfo->specialcharacters = 0;
  specialcharinfo->lengthofspecialprefix = 0;
  specialcharinfo->lengthofspecialsuffix = 0;

  if (withdestab)
  {
    descqueue = queue_new();
    desfp = opensfxfile(indexname,DESTABSUFFIX,"wb",err);
    if (desfp == NULL)
    {
      haserr = true;
    }
  }
  if (!haserr)
  {
    fb = fastabuffer_new(filenametab,
                         symbolmap,
                         plainformat,
                         filelengthtab,
                         descqueue,
                         characterdistribution);
    distspralen = disc_distri_new();
    for (pos = 0; /* Nothing */; pos++)
    {
      retval = fastabuffer_next(fb,&charcode,err);
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
          disc_distri_add(distspralen,idx);
        }
        break;
      }
      if (ISSPECIAL(charcode))
      {
        if (desfp != NULL && charcode == (Uchar) SEPARATOR)
        {
          desc = queue_get(descqueue);
          if (fputs(desc,desfp) == EOF)
          {
            error_set(err,"cannot write description to file %s.%s",
                              str_get(indexname),DESTABSUFFIX);
            haserr = true;
            break;
          }
          (void) putc((int) '\n',desfp);
          FREESPACE(desc);
        }
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
          disc_distri_add(distspralen,idx);
          lastspeciallength = 0;
        }
      }
    }
  }
  if (!haserr)
  {
    Updatesumrangeinfo updatesumrangeinfo;

    if (desfp != NULL)
    {
      desc = queue_get(descqueue);
      if (fputs(desc,desfp) == EOF)
      {
        error_set(err,"cannot write description to file %s.%s",
                          str_get(indexname),DESTABSUFFIX);
        haserr = true;
      }
      (void) putc((int) '\n',desfp);
      FREESPACE(desc);
    }
    specialcharinfo->specialranges = 0;
    specialcharinfo->realspecialranges = 0;
    updatesumrangeinfo.specialrangesptr = &specialcharinfo->specialranges;
    updatesumrangeinfo.realspecialrangesptr
      = &specialcharinfo->realspecialranges;
    updatesumrangeinfo.verboseinfo = verboseinfo;
    disc_distri_foreach(distspralen,updatesumranges,
                       &updatesumrangeinfo);
    specialcharinfo->lengthofspecialsuffix = lastspeciallength;
    (*numofsequences)++;
    *totallength = pos;
  }
  fa_xfclose(desfp);
  disc_distri_delete(distspralen);
  fastabuffer_delete(fb);
  queue_delete_with_contents(descqueue);
  return haserr ? -1 : 0;
}

void sequence2specialcharinfo(Specialcharinfo *specialcharinfo,
                              const Uchar *seq,
                              const Seqpos len,
                              Verboseinfo *verboseinfo)
{
  Uchar charcode;
  Seqpos pos;
  bool specialprefix = true;
  Seqpos lastspeciallength = 0;
  DiscDistri *distspralen;
  unsigned long idx;
  Updatesumrangeinfo updatesumrangeinfo;

  specialcharinfo->specialcharacters = 0;
  specialcharinfo->lengthofspecialprefix = 0;
  specialcharinfo->lengthofspecialsuffix = 0;
  distspralen = disc_distri_new();
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
        disc_distri_add(distspralen,idx);
        lastspeciallength = 0;
      }
    }
  }
  if (lastspeciallength > 0)
  {
    idx = CALLCASTFUNC(Seqpos,unsigned_long,lastspeciallength);
    disc_distri_add(distspralen,idx);
  }
  specialcharinfo->specialranges = 0;
  specialcharinfo->realspecialranges = 0;
  updatesumrangeinfo.specialrangesptr = &specialcharinfo->specialranges;
  updatesumrangeinfo.realspecialrangesptr = &specialcharinfo->realspecialranges;
  updatesumrangeinfo.verboseinfo = verboseinfo;
  disc_distri_foreach(distspralen,updatesumranges,&updatesumrangeinfo);
  specialcharinfo->lengthofspecialsuffix = lastspeciallength;
  disc_distri_delete(distspralen);
}
