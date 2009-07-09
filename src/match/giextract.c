/*
  Copyright (c) 2007-2008 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
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

#include <errno.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <ctype.h>
#include "core/assert_api.h"
#include "core/fileutils.h"
#include "core/error.h"
#include "core/ma.h"
#include "core/fa.h"
#include "core/seqiterator.h"
#include "core/progressbar.h"
#include "core/fasta.h"
#include "divmodmul.h"
#include "format64.h"

#define COMPLETE(VALUE)\
        ((VALUE).frompos == 1UL && (VALUE).topos == 0)

#define EXTRABUF 128

#define CHECKPOSITIVE(VAL,FORMAT,WHICH)\
        if ((VAL) <= 0)\
        {\
          gt_error_set(err,"file \"%s\", line %lu: illegal format: %s element "\
                        " = " FORMAT " is not a positive integer",\
                        gt_str_get(ginumberfile),\
                        linenum+1,\
                        WHICH,\
                        VAL);\
          haserr = true;\
          break;\
        }

typedef struct
{
  char *ginumber;
  unsigned long frompos, topos;
  bool markhit;
} Giquery;

static int compareginumbers(const void *a,const void *b)
{
  int cmp = strcmp(((Giquery *) a)->ginumber,((Giquery *) b)->ginumber);

  if (cmp < 0)
  {
    return -1;
  }
  if (cmp > 0)
  {
    return 1;
  }
  if (((Giquery *) a)->frompos < ((Giquery *) b)->frompos)
  {
    return -1;
  }
  if (((Giquery *) a)->frompos > ((Giquery *) b)->frompos)
  {
    return 1;
  }
  if (((Giquery *) a)->topos < ((Giquery *) b)->topos)
  {
    return -1;
  }
  if (((Giquery *) a)->topos > ((Giquery *) b)->topos)
  {
    return 1;
  }
  return 0;
}

static unsigned long remdupsgiqueries(Giquery *giqueries,
                                      unsigned long numofqueries)
{
  if (numofqueries == 0)
  {
    return 0;
  } else
  {
    Giquery *storeptr, *readptr;
    unsigned long newnumofqueries;

    for (storeptr = giqueries, readptr = giqueries+1;
         readptr < giqueries + numofqueries;
         readptr++)
    {
      if (storeptr->ginumber != readptr->ginumber ||
          storeptr->frompos != readptr->frompos ||
          storeptr->topos != readptr->topos)
      {
        storeptr++;
        if (storeptr != readptr)
        {
          *storeptr = *readptr;
        }
      }
    }
    newnumofqueries = (unsigned long) (storeptr - giqueries + 1);
    if (newnumofqueries < numofqueries)
    {
      printf("# removed %lu duplicate gi-queries\n",
              numofqueries - newnumofqueries);
    }
    return newnumofqueries;
  }
}

static Giquery *readginumberfile(bool verbose,
                                 unsigned long *numofqueries,
                                 const GtStr *ginumberfile,
                                 GtError *err)
{
  FILE *fp;
  GtStr *currentline;
  bool haserr = false;
  unsigned long linenum;
  long readlongfrompos, readlongtopos;
  Giquery *giqueries;
#undef SKDEBUG
#ifdef SKDEBUG
  unsigned long i;
#endif

  gt_error_check(err);
  *numofqueries = gt_file_number_of_lines(gt_str_get(ginumberfile));
  if (*numofqueries == 0)
  {
    gt_error_set(err,"empty file \"%s\" not allowed",gt_str_get(ginumberfile));
    return NULL;
  }
  fp = gt_fa_fopen(gt_str_get(ginumberfile),"r",err);
  if (fp == NULL)
  {
    return NULL;
  }
  if (verbose)
  {
    printf("# opened gi-queryfile \"%s\"\n",gt_str_get(ginumberfile));
  }
  giqueries = gt_malloc(sizeof (*giqueries) * (*numofqueries));
  currentline = gt_str_new();
  for (linenum = 0; gt_str_read_next_line(currentline, fp) != EOF; linenum++)
  {
    char *lineptr = gt_str_get(currentline);
    size_t idx;

    for (idx = 0; !isspace(lineptr[idx]); idx++)
      /* nothing */ ;
    if (sscanf(lineptr+idx,"%ld %ld\n",&readlongfrompos,&readlongtopos) != 2)
    {
      gt_error_set(err,"file \"%s\", line %lu: illegal format",
                  gt_str_get(ginumberfile),
                  linenum+1);
      haserr = true;
      break;
    }
    giqueries[linenum].ginumber = gt_malloc(sizeof(char) * (idx+1));
    strncpy(giqueries[linenum].ginumber,lineptr,idx);
    giqueries[linenum].ginumber[idx] = '\0';
    CHECKPOSITIVE(readlongfrompos,"%ld","second");
    giqueries[linenum].frompos = (unsigned long) readlongfrompos;
    if (readlongfrompos != 1L || readlongtopos != 0)
    {
      CHECKPOSITIVE(readlongtopos,"%ld","third");
    }
    giqueries[linenum].topos = (unsigned long) readlongtopos;
    giqueries[linenum].markhit = false;
    if (!COMPLETE(giqueries[linenum]) &&
        giqueries[linenum].frompos > giqueries[linenum].topos)
    {
      gt_error_set(err, "file \"%s\", line %lu: illegal format: second value "
                   "%lu is larger than third value %lu",
                   gt_str_get(ginumberfile),
                   linenum+1,
                   giqueries[linenum].frompos,
                   giqueries[linenum].topos);
      haserr = true;
      break;
    }
    gt_str_reset(currentline);
  }
  gt_str_delete(currentline);
  gt_fa_fclose(fp);
  if (haserr)
  {
    gt_free(giqueries);
    return NULL;
  }
  qsort(giqueries,(size_t) *numofqueries,sizeof (*giqueries),
        compareginumbers);
  if (verbose)
  {
    printf("# %lu gi-queries successfully parsed and sorted\n",*numofqueries);
  }
  *numofqueries = remdupsgiqueries(giqueries,*numofqueries);
#ifdef SKDEBUG
  for (i=0; i<*numofqueries; i++)
  {
    printf("%lu %s\n",i,giqueries[i].ginumber);
  }
#endif
  return giqueries;
}

static unsigned long findginumber(const char *ginumber,
                                  unsigned long ginumlen,
                                  const Giquery *giqueries,
                                  unsigned long numofqueries)
{
  const Giquery *leftptr, *rightptr, *midptr;
  int cmp;

  leftptr = giqueries;
  rightptr = giqueries + numofqueries - 1;
  while (leftptr <= rightptr)
  {
    midptr = leftptr + DIV2((unsigned long) (rightptr-leftptr));
    cmp = strncmp(ginumber,midptr->ginumber,(size_t) ginumlen);
    if (cmp == 0)
    {
      if (midptr > giqueries &&
          strncmp(ginumber,(midptr-1)->ginumber,(size_t) ginumlen) == 0)
      {
        rightptr = midptr - 1;
      } else
      {
        return (unsigned long) (midptr - giqueries);
      }
    } else
    {
      if (cmp < 0)
      {
        rightptr = midptr-1;
      } else
      {
        leftptr = midptr + 1;
      }
    }
  }
  return numofqueries;
}

static void outputnonmarked(const Giquery *giqueries,
                            unsigned long numofqueries)
{
  unsigned long idx, countmissing = 0;

  for (idx=0; idx<numofqueries; idx++)
  {
    if (!giqueries[idx].markhit)
    {
      printf("unsatisfied %s",giqueries[idx].ginumber);
      if (COMPLETE(giqueries[idx]))
      {
        printf(" complete\n");
      } else
      {
        printf(" %lu %lu\n",giqueries[idx].frompos,giqueries[idx].topos);
      }
      countmissing++;
    }
  }
  printf("# number of unsatified gi queries: %lu\n",countmissing);
}

static void giqueries_delete(Giquery *giqueries,unsigned long numofqueries)
{
  unsigned long idx;

  for (idx=0; idx<numofqueries; idx++)
  {
    gt_free(giqueries[idx].ginumber);
  }
  gt_free(giqueries);
}

static const char *desc2ginumber(unsigned long *ginumlen,const char *desc,
                                 GtError *err)
{
  unsigned long i, firstpipe = 0, secondpipe = 0;

  gt_error_check(err);
  for (i=0; desc[i] != '\0'; i++)
  {
    if (desc[i] == '|')
    {
      if (firstpipe > 0)
      {
        gt_assert(i>0);
        secondpipe = i;
        break;
      }
      gt_assert(i>0);
      firstpipe = i;
    }
  }
  if (firstpipe == 0 || secondpipe == 0)
  {
    gt_error_set(err,"Cannot find gi-number in description \"%s\"\n",desc);
    return NULL;
  }
  gt_assert(firstpipe < secondpipe);
  *ginumlen = secondpipe - firstpipe - 1;
  return desc + firstpipe + 1;
}

int extractginumbers(bool verbose,
                     GtGenFile *outfp,
                     unsigned long width,
                     const GtStr *ginumberfile,
                     GtStrArray *referencefiletab,
                     GtError *err)
{
  GtSeqIterator *seqit;
  const GtUchar *sequence;
  char *desc, *headerbufferspace = NULL;
  const char *ginumberasstring;
  unsigned long len, ginumlen, numofqueries, ginumberhit, countmarkhit = 0;
  int had_err = 0;
  off_t totalsize;
  Giquery *giqueries;
  size_t headerbuffersize = 0, headerlength;

  gt_error_check(err);
  giqueries = readginumberfile(verbose,&numofqueries,ginumberfile,err);
  if (giqueries == NULL)
  {
    return -1;
  }
  totalsize = gt_files_estimate_total_size(referencefiletab);
  printf("# estimated total size is " Formatuint64_t "\n",
            PRINTuint64_tcast(totalsize));
  seqit = gt_seqiterator_new(referencefiletab, err);
  if (!seqit)
    had_err = -1;
  if (!had_err && verbose)
  {
    gt_progressbar_start(gt_seqiterator_getcurrentcounter(seqit,
                                                          (unsigned long long)
                                                          totalsize),
                                                          (unsigned long long)
                                                          totalsize);
  }
  while (had_err != -1 && countmarkhit < numofqueries)
  {
    had_err = gt_seqiterator_next(seqit, &sequence, &len, &desc, err);
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
      ginumberhit = findginumber(ginumberasstring,ginumlen,giqueries,
                                 numofqueries);
      if (ginumberhit < numofqueries)
      {
        while (ginumberhit < numofqueries &&
               strncmp(giqueries[ginumberhit].ginumber,ginumberasstring,
                       (size_t) ginumlen) == 0)
        {
#ifndef NDEBUG
          if (giqueries[ginumberhit].markhit)
          {
            fprintf(stderr,"ginumber %s was already found before\n",
                     giqueries[ginumberhit].ginumber);
            exit(GT_EXIT_PROGRAMMING_ERROR);
          }
#endif
          headerlength = strlen(desc);
          if (headerbuffersize < headerlength + EXTRABUF + 1)
          {
            headerbuffersize = headerlength + EXTRABUF + 1;
            headerbufferspace = gt_realloc(headerbufferspace,
                                           sizeof (*headerbufferspace)
                                           * headerbuffersize);
          }
          if (COMPLETE(giqueries[ginumberhit]))
          {
            /*
            (void) snprintf(headerbufferspace,headerbuffersize,
                            "%*.*s complete %s",
                            (int) ginumlen,(int) ginumlen,ginumberasstring,
                            desc);
            */
            gt_fasta_show_entry_generic(desc,
                                        (const char *) sequence,
                                        len, width, outfp);
          } else
          {
            (void) snprintf(headerbufferspace,headerbuffersize,
                            "%*.*s %lu %lu %s",
                            (int) ginumlen,(int) ginumlen,ginumberasstring,
                            giqueries[ginumberhit].frompos,
                            giqueries[ginumberhit].topos,
                            desc);
            gt_fasta_show_entry_generic(headerbufferspace,
                                        (const char *) (sequence +
                                                        giqueries[ginumberhit].
                                                        frompos - 1),
                                        giqueries[ginumberhit].topos -
                                        giqueries[ginumberhit].frompos+1,
                                        width, outfp);
          }
          giqueries[ginumberhit].markhit = true;
          countmarkhit++;
          ginumberhit++;
        }
      }
#ifdef SKDEBUG
      printf("%*.*s 1 %lu\n",ginumlen,ginumlen,ginumberasstring, len);
#endif
    }
    gt_free(desc);
  }
  gt_free(headerbufferspace);
  if (verbose)
  {
    gt_progressbar_stop();
  }
  outputnonmarked(giqueries,numofqueries);
  giqueries_delete(giqueries,numofqueries);
  gt_seqiterator_delete(seqit);
  return had_err;
}
