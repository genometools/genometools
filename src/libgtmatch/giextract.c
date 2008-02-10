/*
  Copyright (c) 2007      Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
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

#include <assert.h>
#include <errno.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "libgtcore/fileutils.h"
#include "libgtcore/error.h"
#include "libgtcore/ma.h"
#include "libgtcore/fa.h"
#include "libgtcore/seqiterator.h"
#include "libgtcore/progressbar.h"
#include "libgtcore/fasta.h"
#include "spacedef.h"
#include "divmodmul.h"
#include "format64.h"

#define EXTRABUF 128

#define CHECKPOSITIVE(VAL,FORMAT,WHICH)\
        if ((VAL) <= 0)\
        {\
          error_set(err,"file \"%s\", line %lu: illegal format: %s element "\
                        " = " FORMAT " is not a positive integer",\
                        str_get(ginumberfile),\
                        linenum+1,\
                        WHICH,\
                        VAL);\
          haserr = true;\
          break;\
        }

typedef struct
{
  uint64_t ginumber;
  unsigned long frompos, topos;
  bool markhit;
} Ginumberwithrange;

static int compareginumbers(const void *a,const void *b)
{
  if (((Ginumberwithrange *) a)->ginumber <
      ((Ginumberwithrange *) b)->ginumber)
  {
    return -1;
  }
  if (((Ginumberwithrange *) a)->ginumber >
      ((Ginumberwithrange *) b)->ginumber)
  {
    return 1;
  }
  if (((Ginumberwithrange *) a)->frompos <
      ((Ginumberwithrange *) b)->frompos)
  {
    return -1;
  }
  if (((Ginumberwithrange *) a)->frompos >
      ((Ginumberwithrange *) b)->frompos)
  {
    return 1;
  }
  if (((Ginumberwithrange *) a)->topos <
      ((Ginumberwithrange *) b)->topos)
  {
    return -1;
  }
  if (((Ginumberwithrange *) a)->topos >
      ((Ginumberwithrange *) b)->topos)
  {
    return 1;
  }
  return 0;
}

static Ginumberwithrange *readginumberfile(bool verbose,
                                           unsigned long *numofentries,
                                           const Str *ginumberfile,
                                           Error *err)
{
  FILE *fp;
  bool haserr = false;
  unsigned long linenum;
  int64_t readint64;
  long readlongfrom, readlongto;
  Ginumberwithrange *ginumbertable;
#ifdef DEBUG
  unsigned long i;
#endif

  error_check(err);
  *numofentries = file_number_of_lines(str_get(ginumberfile));
  if (*numofentries == 0)
  {
    error_set(err,"empty file \"%s\" not allowed",str_get(ginumberfile));
    return NULL;
  }
  fp = fa_fopen(str_get(ginumberfile),"r");
  if (fp == NULL)
  {
    error_set(err,"fa_fopen: cannot open file \"%s\": %s",
                  str_get(ginumberfile),
                  strerror(errno));
    return NULL;
  }
  if (verbose)
  {
    printf("# opened gi-queryfile \"%s\"\n",str_get(ginumberfile));
  }
  ALLOCASSIGNSPACE(ginumbertable,NULL,Ginumberwithrange,*numofentries);
  for (linenum = 0; !feof(fp); linenum++)
  {
    if (fscanf(fp,FormatScanint64_t " %ld %ld\n",
               SCANint64_tcast(&readint64),&readlongfrom,&readlongto) != 3)
    {
      error_set(err,"file \"%s\", line %lu: illegal format",
                  str_get(ginumberfile),
                  linenum+1);
      haserr = true;
      break;
    }
    CHECKPOSITIVE(readint64,FormatScanint64_t,"first");
    ginumbertable[linenum].ginumber = (uint64_t) readint64;
    CHECKPOSITIVE(readlongfrom,"%ld","second");
    ginumbertable[linenum].frompos = (unsigned long) readlongfrom;
    CHECKPOSITIVE(readlongfrom,"%ld","third");
    ginumbertable[linenum].topos = (unsigned long) readlongto;
    ginumbertable[linenum].markhit = false;
    if (ginumbertable[linenum].frompos > ginumbertable[linenum].topos)
    {
      error_set(err,"file \"%s\", line %lu: illegal format: second value %lu "
                    "is larger than third value %lu",
                  str_get(ginumberfile),
                  linenum+1,
                  ginumbertable[linenum].frompos,
                  ginumbertable[linenum].topos);
      haserr = true;
      break;
    }
  }
  fa_fclose(fp);
  if (haserr)
  {
    FREESPACE(ginumbertable);
    return NULL;
  }
  qsort(ginumbertable,(size_t) *numofentries,sizeof(Ginumberwithrange),
        compareginumbers);
  if (verbose)
  {
    printf("# %lu gi-queries successfully parsed and sorted\n",*numofentries);
  }
#ifdef DEBUG
  for (i=0; i<*numofentries; i++)
  {
    printf("%lu %lu\n",i,ginumbertable[i].ginumber);
  }
#endif
  return ginumbertable;
}

static unsigned long findginumber(uint64_t ginumber,
                                  const Ginumberwithrange *ginumbertable,
                                  unsigned long numofentries)
{
  const Ginumberwithrange *leftptr, *rightptr, *midptr;

  leftptr = ginumbertable;
  rightptr = ginumbertable + numofentries - 1;
  while (leftptr <= rightptr)
  {
    midptr = leftptr + DIV2((unsigned long) (rightptr-leftptr));
    if (midptr->ginumber == ginumber)
    {
      if (midptr > ginumbertable && (midptr-1)->ginumber == ginumber)
      {
        rightptr = midptr - 1;
      } else
      {
        return (unsigned long) (midptr - ginumbertable);
      }
    } else
    {
      if (ginumber < midptr->ginumber)
      {
        rightptr = midptr-1;
      } else
      {
        leftptr = midptr + 1;
      }
    }
  }
  return numofentries;
}

static void outputnonmarked(const Ginumberwithrange *ginumbertable,
                            unsigned long numofentries)
{
  unsigned long i, countmissing = 0;

  for (i=0; i<numofentries; i++)
  {
    if (!ginumbertable[i].markhit)
    {
      printf("unsatisfied " Formatuint64_t " %lu %lu\n",
              PRINTuint64_tcast(ginumbertable[i].ginumber),
              ginumbertable[i].frompos,
              ginumbertable[i].topos);
      countmissing++;
    }
  }
  printf("# number of unsatified gi queries: %lu\n",countmissing);
}

static const char *desc2ginumber(unsigned long *ginumlen,const char *desc,
                                 Error *err)
{
  unsigned long i, firstpipe = 0, secondpipe = 0;

  error_check(err);
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
                     GenFile *outfp,
                     unsigned long width,
                     const Str *ginumberfile,
                     StrArray *referencefiletab,
                     Error *err)
{
  SeqIterator *seqit;
  const Uchar *sequence;
  char *desc, *headerbufferspace = NULL;
  const char *ginumberasstring;
  uint64_t referenceginumber;
  unsigned long len, ginumlen, numofentries, ginumberhit, countmarkhit = 0;
  int had_err = 0;
  int64_t readint64;
  off_t totalsize;
  Ginumberwithrange *ginumbertable;
  size_t headerbuffersize = 0, headerlength;

  error_check(err);
  ginumbertable = readginumberfile(verbose,&numofentries,ginumberfile,err);
  if (ginumbertable == NULL)
  {
    return -1;
  }
  totalsize = files_estimate_total_size(referencefiletab);
  printf("# estimated total size is " Formatuint64_t "\n",
            PRINTuint64_tcast(totalsize));
  seqit = seqiterator_new(referencefiletab, NULL, true);
  if (verbose)
  {
    progressbar_start(seqiterator_getcurrentcounter(seqit, (unsigned long long)
                                                           totalsize),
                                                           (unsigned long long)
                                                           totalsize);
  }
  while (had_err != -1 && countmarkhit < numofentries)
  {
    had_err = seqiterator_next(seqit, &sequence, &len, &desc, err);
    if (had_err != 1)
    {
      break;
    }
    ginumberasstring = desc2ginumber(&ginumlen,desc,err);
    if (ginumberasstring == NULL)
    {
      had_err = -1;
    } else
    {
      if (sscanf(ginumberasstring,FormatScanint64_t "|",
                 SCANint64_tcast(&readint64)) != 1)
      {
        error_set(err,"cannot parse ginumber(integer) in \"%s\"",
                    ginumberasstring);
        had_err = -1;
      }
      if (had_err != -1 && readint64 <= 0)
      {
        error_set(err,"gi number " Formatuint64_t "must be positive integer",
                      readint64);
        had_err = -1;
      }
      referenceginumber = (uint64_t) readint64;
      ginumberhit = findginumber(referenceginumber,ginumbertable,numofentries);
      if (ginumberhit < numofentries)
      {
        while (ginumberhit < numofentries &&
               ginumbertable[ginumberhit].ginumber == referenceginumber)
        {
          if (ginumbertable[ginumberhit].markhit)
          {
            fprintf(stderr,"ginumber " Formatuint64_t
                           " was already found before\n",
                     PRINTuint64_tcast(ginumbertable[ginumberhit].ginumber));
            exit(EXIT_FAILURE); /* programming error */
          }
          headerlength = strlen(desc);
          if (headerbuffersize < headerlength + EXTRABUF + 1)
          {
            headerbuffersize = headerlength + EXTRABUF + 1;
            headerbufferspace = ma_realloc(headerbufferspace,
                                           sizeof (char) * headerbuffersize);
          }
          (void) snprintf(headerbufferspace,headerbuffersize,Formatuint64_t
                          " %lu %lu %s",
                          PRINTuint64_tcast(referenceginumber),
                          ginumbertable[ginumberhit].frompos,
                          ginumbertable[ginumberhit].topos,
                          desc);
          fasta_show_entry_generic(headerbufferspace,
                                   (const char *) (sequence +
                                                   ginumbertable[ginumberhit].
                                                   frompos - 1),
                                   ginumbertable[ginumberhit].topos -
                                   ginumbertable[ginumberhit].frompos+1,
                                   width, outfp);
          ginumbertable[ginumberhit].markhit = true;
          countmarkhit++;
          ginumberhit++;
        }
      }
//#ifdef DEBUG
      printf(Formatuint64_t " 1 %lu\n",PRINTuint64_tcast(referenceginumber),
             len);
//#endif
    }
    ma_free(desc);
  }
  ma_free(headerbufferspace);
  if (verbose)
  {
    progressbar_stop();
  }
  outputnonmarked(ginumbertable,numofentries);
  FREESPACE(ginumbertable);
  seqiterator_delete(seqit);
  return had_err;
}
