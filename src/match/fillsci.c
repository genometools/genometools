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
#include "core/sequence_buffer_fasta.h"
#include "core/sequence_buffer_plain.h"
#include "core/str_array.h"
#include "esa-fileend.h"
#include "seqpos-def.h"
#include "verbose-def.h"
#include "spacedef.h"
#include "safecast-gen.h"
#include "alphadef.h"
#include "encseq-def.h"
#include "stamp.h"
#include "opensfxfile.h"
#include "fillsci.h"

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
    gt_assert(len - 1 <= UINT32_MAX);
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

  gt_assert(value <= (unsigned long long) ULONG_MAX);
  distvalue = (unsigned long) value;
  updatesumrangeinfo->specialrangesUchar
     += currentspecialrangevalue(key,distvalue,(unsigned long) UCHAR_MAX);
  updatesumrangeinfo->specialrangesUshort
     += currentspecialrangevalue(key,distvalue,(unsigned long) USHRT_MAX);
  updatesumrangeinfo->specialrangesUint32
     += currentspecialrangevalue(key,distvalue,(unsigned long) UINT32_MAX);
  updatesumrangeinfo->realspecialranges += distvalue;
  showverbose(updatesumrangeinfo->verboseinfo,
              "specialranges of length %lu=%lu",key,distvalue);
}

 DECLARESAFECASTFUNCTION(Seqpos,Seqpos,unsigned long,unsigned_long)

static Seqpos calcspecialranges(Seqpos *specialrangestab,
                                GtDiscDistri *distspralen,
                                Verboseinfo *verboseinfo)
{
  Updatesumrangeinfo updatesumrangeinfo;

  updatesumrangeinfo.specialrangesUchar = 0;
  updatesumrangeinfo.specialrangesUshort = 0;
  updatesumrangeinfo.specialrangesUint32 = 0;
  updatesumrangeinfo.realspecialranges = 0;
  updatesumrangeinfo.verboseinfo = verboseinfo;
  gt_disc_distri_foreach(distspralen,updatesumranges,&updatesumrangeinfo);
  if (specialrangestab != NULL)
  {
    specialrangestab[0] = updatesumrangeinfo.specialrangesUchar;
    specialrangestab[1] = updatesumrangeinfo.specialrangesUshort;
    specialrangestab[2] = updatesumrangeinfo.specialrangesUint32;
  }
  return updatesumrangeinfo.realspecialranges;
}

static void doupdatesumranges(Specialcharinfo *specialcharinfo,
                              unsigned int forcetable,
                              Seqpos *specialrangestab,
                              Seqpos totallength,
                              unsigned int numofchars,
                              GtDiscDistri *distspralen,
                              Verboseinfo *verboseinfo)
{
  uint64_t smallestsize = 0, tmp;
  bool smallestdefined = false;
  int c;

  specialcharinfo->realspecialranges
    = calcspecialranges(specialrangestab,distspralen,verboseinfo);
  assert (forcetable <= 3U);
  for (c = 0; c<3; c++)
  {
    if (forcetable == 3U || c == (int) forcetable)
    {
      tmp = detencseqofsatviatables(c,totallength,specialrangestab[c],
                                    numofchars);
      if (!smallestdefined || tmp < smallestsize)
      {
        smallestdefined = true;
        smallestsize = tmp;
        specialcharinfo->specialranges = specialrangestab[c];
      }
    }
  }
}

int fasta2sequencekeyvalues(
        const GtStr *indexname,
        Seqpos *totallength,
        Specialcharinfo *specialcharinfo,
        unsigned int forcetable,
        Seqpos *specialrangestab,
        const GtStrArray *filenametab,
        Filelengthvalues **filelengthtab,
        const SfxAlphabet *alpha,
        bool plainformat,
        bool withdestab,
        unsigned long *characterdistribution,
        bool withssptab,
        ArraySeqpos *sequenceseppos,
        Verboseinfo *verboseinfo,
        GtError *err)
{
  GtSequenceBuffer *fb = NULL;
  Uchar charcode;
  Seqpos currentpos = 0;
  int retval;
  bool specialprefix = true;
  Seqpos lastspeciallength = 0;
  GtDiscDistri *distspralen = NULL;
  unsigned long idx;
  bool haserr = false;
  GtQueue *descqueue = NULL;
  char *desc;
  FILE *desfp = NULL;

  gt_error_check(err);
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
    if (plainformat) {
      fb = gt_sequence_buffer_plain_new(filenametab);
    } else {
      fb = gt_sequence_buffer_new_guess_type((GtStrArray*) filenametab, err);
    }
    if (!fb)
      haserr = true;
    if (!haserr) {
      gt_sequence_buffer_set_symbolmap(fb, getsymbolmapAlphabet(alpha));
      *filelengthtab = gt_calloc((size_t) gt_str_array_size(filenametab),
                                 sizeof (Filelengthvalues));
      gt_sequence_buffer_set_filelengthtab(fb, *filelengthtab);
      if (descqueue != NULL)
        gt_sequence_buffer_set_desc_queue(fb, descqueue);
      gt_sequence_buffer_set_chardisttab(fb, characterdistribution);

      distspralen = gt_disc_distri_new();
      for (currentpos = 0; /* Nothing */; currentpos++)
      {
        retval = gt_sequence_buffer_next(fb,&charcode,err);
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
            gt_disc_distri_add(distspralen,idx);
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
            if (withssptab)
            {
              STOREINARRAY(sequenceseppos,Seqpos,128,currentpos);
            } else
            {
              sequenceseppos->nextfreeSeqpos++;
            }
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
            gt_disc_distri_add(distspralen,idx);
            lastspeciallength = 0;
          }
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
    *totallength = currentpos;
    specialcharinfo->lengthofspecialsuffix = lastspeciallength;
    doupdatesumranges(specialcharinfo,forcetable,specialrangestab,currentpos,
                      getnumofcharsAlphabet(alpha),distspralen,verboseinfo);
  }
  gt_fa_xfclose(desfp);
  gt_disc_distri_delete(distspralen);
  gt_sequence_buffer_delete(fb);
  gt_queue_delete_with_contents(descqueue);
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
  GtDiscDistri *distspralen;
  unsigned long idx;

  specialcharinfo->specialcharacters = 0;
  specialcharinfo->lengthofspecialprefix = 0;
  specialcharinfo->lengthofspecialsuffix = 0;
  distspralen = gt_disc_distri_new();
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
        gt_disc_distri_add(distspralen,idx);
        lastspeciallength = 0;
      }
    }
  }
  if (lastspeciallength > 0)
  {
    idx = CALLCASTFUNC(Seqpos,unsigned_long,lastspeciallength);
    gt_disc_distri_add(distspralen,idx);
  }
  specialcharinfo->lengthofspecialsuffix = lastspeciallength;
  specialcharinfo->realspecialranges
    = calcspecialranges(NULL,distspralen,verboseinfo);
  specialcharinfo->specialranges = specialcharinfo->realspecialranges;
  gt_disc_distri_delete(distspralen);
}
