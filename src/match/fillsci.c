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
#include "core/chardef.h"
#include "core/disc_distri.h"
#include "core/error.h"
#include "core/fa.h"
#include "core/fastabuffer.h"
#include "core/strarray.h"
#include "esa-fileend.h"
#include "seqpos-def.h"
#include "verbose-def.h"
#include "spacedef.h"
#include "safecast-gen.h"
#include "encseq-def.h"
#include "stamp.h"

#include "opensfxfile.pr"

static unsigned long currentspecialrangevalue(
                             unsigned long len,
                             unsigned long occcount,
                             unsigned long maxspecialtype)
{
/*
  printf("len=%lu,occcount=%lu,maxspecialtype=%lu\n",
           len,occcount,maxspecialtype);
*/
  if (maxspecialtype == UINT32_MAX)
  {
    assert(len - 1 <= UINT32_MAX);
    return occcount;
  }
  if (len <= maxspecialtype+1)
  {
    return occcount;
  }
  if (len % (maxspecialtype+1) == 0)
  {
    return len/(maxspecialtype+1) * occcount;
  }
  return (1UL + len/(maxspecialtype+1)) * occcount;
}

typedef struct
{
  Verboseinfo *verboseinfo;
  Seqpos specialrangesUchar,
         specialrangesUshort,
         specialrangesUint32,
         realspecialranges;
} Updatesumrangeinfo;

static void updatesumranges(unsigned long key, unsigned long long value,
                            void *data)
{
  unsigned long distvalue;
  Updatesumrangeinfo *updatesumrangeinfo = (Updatesumrangeinfo *) data;

  distvalue = (unsigned long) value; /* XXX: is this cast always OK? */
  updatesumrangeinfo->specialrangesUchar
     += currentspecialrangevalue(key,distvalue,(unsigned long) UCHAR_MAX);
  updatesumrangeinfo->specialrangesUshort
     += currentspecialrangevalue(key,distvalue,(unsigned long) USHRT_MAX);
  updatesumrangeinfo->specialrangesUint32
     += currentspecialrangevalue(key,distvalue,(unsigned long) UINT32_MAX);
  updatesumrangeinfo->realspecialranges += distvalue;
  showverbose(updatesumrangeinfo->verboseinfo,
              "specialranges of length %lu: %lu",key,distvalue);
}

 DECLARESAFECASTFUNCTION(Seqpos,Seqpos,unsigned long,unsigned_long)

static void doupdatesumranges(Specialcharinfo *specialcharinfo,
                              Seqpos totallength,
                              unsigned int mapsize,
                              DiscDistri *distspralen,
                              Verboseinfo *verboseinfo)
{
  Updatesumrangeinfo updatesumrangeinfo;
  uint64_t smallestsize, tmp;

  updatesumrangeinfo.specialrangesUchar = 0;
  updatesumrangeinfo.specialrangesUshort = 0;
  updatesumrangeinfo.specialrangesUint32 = 0;
  updatesumrangeinfo.realspecialranges = 0;
  updatesumrangeinfo.verboseinfo = verboseinfo;
  disc_distri_foreach(distspralen,updatesumranges,&updatesumrangeinfo);
  smallestsize = detsizeencseq(0,totallength,
                               updatesumrangeinfo.specialrangesUchar,mapsize);
  specialcharinfo->specialranges = updatesumrangeinfo.specialrangesUchar;
  tmp = detsizeencseq(1,totallength,updatesumrangeinfo.specialrangesUshort,
                      mapsize);
  if (tmp < smallestsize)
  {
    smallestsize = tmp;
    specialcharinfo->specialranges = updatesumrangeinfo.specialrangesUshort;
  }
  tmp = detsizeencseq(2,totallength,
                      updatesumrangeinfo.specialrangesUint32,mapsize);
  if (tmp < smallestsize)
  {
    specialcharinfo->specialranges = updatesumrangeinfo.specialrangesUint32;
  }
  /*
  printf("specialrangesUchar=%lu\n",
         (unsigned long) updatesumrangeinfo.specialrangesUchar);
  printf("specialrangesUshort=%lu\n",
         (unsigned long) updatesumrangeinfo.specialrangesUshort);
  printf("specialrangesUint32=%lu\n",
         (unsigned long) updatesumrangeinfo.specialrangesUint32);
  printf("specialranges%lu\n",
         (unsigned long) specialcharinfo->specialranges);
  */
  specialcharinfo->realspecialranges = updatesumrangeinfo.realspecialranges;
}

int fasta2sequencekeyvalues(
        const GtStr *indexname,
        unsigned long *numofsequences,
        Seqpos *totallength,
        Specialcharinfo *specialcharinfo,
        const GtStrArray *filenametab,
        Filelengthvalues **filelengthtab,
        const Alphabet *alpha,
        bool plainformat,
        bool withdestab,
        unsigned long *characterdistribution,
        Verboseinfo *verboseinfo,
        GtError *err)
{
  GT_FastaBuffer *fb = NULL;
  Uchar charcode;
  Seqpos pos = 0;
  int retval;
  bool specialprefix = true;
  Seqpos lastspeciallength = 0;
  DiscDistri *distspralen = NULL;
  unsigned long idx;
  bool haserr = false;
  GtQueue *descqueue = NULL;
  char *desc;
  FILE *desfp = NULL;

  gt_error_check(err);
  *numofsequences = 0;
  specialcharinfo->specialcharacters = 0;
  specialcharinfo->lengthofspecialprefix = 0;
  specialcharinfo->lengthofspecialsuffix = 0;

  if (withdestab)
  {
    descqueue = gt_queue_new();
    desfp = opensfxfile(indexname,DESTABSUFFIX,"wb",err);
    if (desfp == NULL)
    {
      haserr = true;
    }
  }
  if (!haserr)
  {
    fb = gt_fastabuffer_new(filenametab,
                         getsymbolmapAlphabet(alpha),
                         plainformat,
                         filelengthtab,
                         descqueue,
                         characterdistribution);
    distspralen = disc_distri_new();
    for (pos = 0; /* Nothing */; pos++)
    {
      retval = gt_fastabuffer_next(fb,&charcode,err);
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
          desc = gt_queue_get(descqueue);
          if (fputs(desc,desfp) == EOF)
          {
            gt_error_set(err,"cannot write description to file %s.%s",
                              gt_str_get(indexname),DESTABSUFFIX);
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
    if (desfp != NULL)
    {
      desc = gt_queue_get(descqueue);
      if (fputs(desc,desfp) == EOF)
      {
        gt_error_set(err,"cannot write description to file %s.%s",
                          gt_str_get(indexname),DESTABSUFFIX);
        haserr = true;
      }
      (void) putc((int) '\n',desfp);
      FREESPACE(desc);
    }
    *totallength = pos;
    specialcharinfo->lengthofspecialsuffix = lastspeciallength;
    doupdatesumranges(specialcharinfo,pos,
                      getmapsizeAlphabet(alpha),distspralen,verboseinfo);
    (*numofsequences)++;
  }
  gt_xfclose(desfp);
  disc_distri_delete(distspralen);
  gt_fastabuffer_delete(fb);
  gt_queue_delete_with_contents(descqueue);
  return haserr ? -1 : 0;
}

void sequence2specialcharinfo(Specialcharinfo *specialcharinfo,
                              const Uchar *seq,
                              const Seqpos len,
                              unsigned int mapsize,
                              Verboseinfo *verboseinfo)
{
  Uchar charcode;
  Seqpos pos;
  bool specialprefix = true;
  Seqpos lastspeciallength = 0;
  DiscDistri *distspralen;
  unsigned long idx;

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
  specialcharinfo->lengthofspecialsuffix = lastspeciallength;
  doupdatesumranges(specialcharinfo,len,mapsize,distspralen,verboseinfo);
  disc_distri_delete(distspralen);
}
