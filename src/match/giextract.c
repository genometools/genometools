/*
  Copyright (c) 2007-2009 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
  Copyright (c) 2007-2009 Center for Bioinformatics, University of Hamburg

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
#include "core/divmodmul.h"
#include "core/fileutils_api.h"
#include "core/error.h"
#include "core/ma.h"
#include "core/fa.h"
#include "core/seqiterator.h"
#include "core/progressbar.h"
#include "core/fasta.h"
#include "giextract.h"
#include "format64.h"
#include "opensfxfile.h"
#include "esa-fileend.h"
#include "encseq-def.h"
#include "echoseq.h"

#define COMPLETE(VALUE)\
        ((VALUE).frompos == 1UL && (VALUE).topos == 0)

#define EXTRABUF 128

#define CHECKPOSITIVE(VAL,FORMAT,WHICH)\
        if ((VAL) <= 0)\
        {\
          gt_error_set(err,"file \"%s\", line %lu: illegal format: %s element "\
                        " = " FORMAT " is not a positive integer",\
                        gt_str_get(fileofkeystoextract),\
                        linenum+1,\
                        WHICH,\
                        VAL);\
          haserr = true;\
          break;\
        }

typedef struct
{
  char *fastakey;
  unsigned long frompos, topos;
  bool markhit;
} Fastakeyquery;

static int comparefastakeys(const void *a,const void *b)
{
  int cmp = strcmp(((Fastakeyquery *) a)->fastakey,
                   ((Fastakeyquery *) b)->fastakey);

  if (cmp < 0)
  {
    return -1;
  }
  if (cmp > 0)
  {
    return 1;
  }
  if (((Fastakeyquery *) a)->frompos < ((Fastakeyquery *) b)->frompos)
  {
    return -1;
  }
  if (((Fastakeyquery *) a)->frompos > ((Fastakeyquery *) b)->frompos)
  {
    return 1;
  }
  if (((Fastakeyquery *) a)->topos < ((Fastakeyquery *) b)->topos)
  {
    return -1;
  }
  if (((Fastakeyquery *) a)->topos > ((Fastakeyquery *) b)->topos)
  {
    return 1;
  }
  return 0;
}

static unsigned long remdupsfastakeyqueries(Fastakeyquery *fastakeyqueries,
                                            unsigned long numofqueries,
                                            bool verbose)
{
  if (numofqueries == 0)
  {
    return 0;
  } else
  {
    Fastakeyquery *storeptr, *readptr;
    unsigned long newnumofqueries;

    for (storeptr = fastakeyqueries, readptr = fastakeyqueries+1;
         readptr < fastakeyqueries + numofqueries;
         readptr++)
    {
      if (strcmp(storeptr->fastakey,readptr->fastakey) != 0 ||
          storeptr->frompos != readptr->frompos ||
          storeptr->topos != readptr->topos)
      {
        storeptr++;
        if (storeptr != readptr)
        {
          size_t len;

          storeptr->frompos = readptr->frompos;
          storeptr->topos = readptr->topos;
          len = strlen(readptr->fastakey);
          storeptr->fastakey = gt_realloc(storeptr->fastakey,
                                          sizeof (char) * (len+1));
          strcpy(storeptr->fastakey,readptr->fastakey);
        }
      }
    }
    newnumofqueries = (unsigned long) (storeptr - fastakeyqueries + 1);
    if (newnumofqueries < numofqueries)
    {
      if (verbose)
      {
        printf("# removed %lu duplicate queries\n",
                numofqueries - newnumofqueries);
      }
      for (storeptr = fastakeyqueries + newnumofqueries;
           storeptr < fastakeyqueries + numofqueries;
           storeptr++)
      {
        gt_free(storeptr->fastakey);
        storeptr->fastakey = NULL;
      }
    }
    return newnumofqueries;
  }
}

static void fastakeyqueries_delete(Fastakeyquery *fastakeyqueries,
                                   unsigned long numofqueries)
{
  if (fastakeyqueries != NULL)
  {
    unsigned long idx;

    for (idx=0; idx<numofqueries; idx++)
    {
      gt_free(fastakeyqueries[idx].fastakey);
    }
    gt_free(fastakeyqueries);
  }
}

static Fastakeyquery *readfileofkeystoextract(bool verbose,
                                              unsigned long *numofqueries,
                                              const GtStr *fileofkeystoextract,
                                              GtError *err)
{
  FILE *fp;
  GtStr *currentline;
  bool haserr = false;
  unsigned long linenum;
  long readlongfrompos, readlongtopos;
  Fastakeyquery *fastakeyqueries;
#undef SKDEBUG
#ifdef SKDEBUG
  unsigned long i;
#endif

  gt_error_check(err);
  *numofqueries = gt_file_number_of_lines(gt_str_get(fileofkeystoextract));
  if (*numofqueries == 0)
  {
    gt_error_set(err,"empty file \"%s\" not allowed",
                 gt_str_get(fileofkeystoextract));
    return NULL;
  }
  fp = gt_fa_fopen(gt_str_get(fileofkeystoextract),"r",err);
  if (fp == NULL)
  {
    return NULL;
  }
  if (verbose)
  {
    printf("# opened keyfile \"%s\"\n",gt_str_get(fileofkeystoextract));
  }
  fastakeyqueries = gt_malloc(sizeof (*fastakeyqueries) * (*numofqueries));
  currentline = gt_str_new();
  for (linenum = 0; gt_str_read_next_line(currentline, fp) != EOF; linenum++)
  {
    char *lineptr = gt_str_get(currentline);
    size_t idx;

    for (idx = 0; lineptr[idx] != '\0' && !isspace(lineptr[idx]); idx++)
      /* nothing */ ;
    fastakeyqueries[linenum].fastakey = gt_malloc(sizeof(char) * (idx+1));
    strncpy(fastakeyqueries[linenum].fastakey,lineptr,idx);
    fastakeyqueries[linenum].fastakey[idx] = '\0';
    if (sscanf(lineptr+idx,"%ld %ld\n",&readlongfrompos,&readlongtopos) == 2)
    {
      CHECKPOSITIVE(readlongfrompos,"%ld","second");
      fastakeyqueries[linenum].frompos = (unsigned long) readlongfrompos;
      CHECKPOSITIVE(readlongtopos,"%ld","third");
      fastakeyqueries[linenum].topos = (unsigned long) readlongtopos;
    } else
    {
      fastakeyqueries[linenum].frompos = 1UL;
      fastakeyqueries[linenum].topos = 0;
    }
    fastakeyqueries[linenum].markhit = false;
    if (!COMPLETE(fastakeyqueries[linenum]) &&
        fastakeyqueries[linenum].frompos > fastakeyqueries[linenum].topos)
    {
      gt_error_set(err, "file \"%s\", line %lu: illegal format: second value "
                   "%lu is larger than third value %lu",
                   gt_str_get(fileofkeystoextract),
                   linenum+1,
                   fastakeyqueries[linenum].frompos,
                   fastakeyqueries[linenum].topos);
      haserr = true;
      break;
    }
    gt_str_reset(currentline);
  }
  gt_str_delete(currentline);
  gt_fa_fclose(fp);
  if (haserr)
  {
    fastakeyqueries_delete(fastakeyqueries,*numofqueries);
    return NULL;
  }
  qsort(fastakeyqueries,(size_t) *numofqueries,sizeof (*fastakeyqueries),
        comparefastakeys);
  if (verbose)
  {
    printf("# %lu fastakey-queries successfully parsed and sorted\n",
            *numofqueries);
  }
  *numofqueries = remdupsfastakeyqueries(fastakeyqueries,*numofqueries,verbose);
#ifdef SKDEBUG
  for (i=0; i<*numofqueries; i++)
  {
    printf("%lu %s\n",i,fastakeyqueries[i].fastakey);
  }
#endif
  return fastakeyqueries;
}

static unsigned long searchdesinfastakeyqueries(const char *extractkey,
                                                const Fastakeyquery
                                                  *fastakeyqueries,
                                                unsigned long numofqueries)
{
  const Fastakeyquery *leftptr, *rightptr, *midptr;
  int cmp;

  leftptr = fastakeyqueries;
  rightptr = fastakeyqueries + numofqueries - 1;
  while (leftptr <= rightptr)
  {
    midptr = leftptr + GT_DIV2((unsigned long) (rightptr-leftptr));
    cmp = strcmp(extractkey,midptr->fastakey);
    if (cmp == 0)
    {
      if (midptr > fastakeyqueries &&
          strcmp(extractkey,(midptr-1)->fastakey) == 0)
      {
        rightptr = midptr - 1;
      } else
      {
        return (unsigned long) (midptr - fastakeyqueries);
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

static void outputnonmarked(const Fastakeyquery *fastakeyqueries,
                            unsigned long numofqueries)
{
  unsigned long idx, countmissing = 0;

  for (idx=0; idx<numofqueries; idx++)
  {
    if (!fastakeyqueries[idx].markhit)
    {
      printf("unsatisfied %s",fastakeyqueries[idx].fastakey);
      if (COMPLETE(fastakeyqueries[idx]))
      {
        printf(" complete\n");
      } else
      {
        printf(" %lu %lu\n",fastakeyqueries[idx].frompos,
                            fastakeyqueries[idx].topos);
      }
      countmissing++;
    }
  }
  printf("# number of unsatified fastakey-queries: %lu\n",countmissing);
}

static const char *desc2key(unsigned long *keylen,const char *desc,
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
    gt_error_set(err,"Cannot find key in description \"%s\"",desc);
    return NULL;
  }
  gt_assert(firstpipe < secondpipe);
  *keylen = secondpipe - firstpipe - 1;
  return desc + firstpipe + 1;
}

#define KEYSTABSUFFIX ".kys"

int gt_extractkeysfromdesfile(const GtStr *indexname, GtError *err)
{
  FILE *fpin, *fpout;
  GtStr *line = NULL;
  uint64_t linenum;
  const char *keyptr;
  unsigned long keylen, constantkeylen = 0;
  bool haserr = false, firstdesc = true;

  fpin = opensfxfile(indexname,DESTABSUFFIX,"rb",err);
  if (fpin == NULL)
  {
    return -1;
  }
  fpout = opensfxfile(indexname,KEYSTABSUFFIX,"wb",err);
  if (fpout == NULL)
  {
    haserr = true;
  }
  if (!haserr)
  {
    line = gt_str_new();
  }
  for (linenum = 0; !haserr && gt_str_read_next_line(line, fpin) != EOF;
       linenum++)
  {
    keyptr = desc2key(&keylen,gt_str_get(line),err);
    if (keyptr == NULL)
    {
      haserr = true;
      break;
    }
    if (keylen == 0)
    {
      gt_error_set(err,"key of length 0 in \"%s\" not expected",
                   gt_str_get(line));
      haserr = true;
      break;
    }
    if (firstdesc)
    {
      constantkeylen = keylen;
      if (constantkeylen > (unsigned long) CHAR_MAX)
      {
        gt_error_set(err,"key of length %lu not allowed; must not be larger "
                         "than %d",constantkeylen,CHAR_MAX);
        haserr = true;
        break;
      }
      (void) fputc((char) constantkeylen,fpout);
      firstdesc = false;
    } else
    {
      if (constantkeylen != keylen)
      {
        gt_error_set(err,"key of length %lu: keys must of length %lu",
                         keylen,constantkeylen);
        haserr = true;
        break;
      }
    }
    if (fwrite(keyptr,sizeof *keyptr,(size_t) keylen,fpout) != (size_t) keylen)
    {
      gt_error_set(err,"cannot write key of length %lu",keylen);
      haserr = true;
      break;
    }
    (void) fputc('\0',fpout);
    gt_str_reset(line);
  }
  if (!haserr)
  {
    printf("number of keys of length %lu = " Formatuint64_t ", ",
            constantkeylen,PRINTuint64_tcast(linenum));
  }
  gt_str_delete(line);
  gt_fa_fclose(fpin);
  gt_fa_fclose(fpout);
  return haserr ? -1 : 0;
}

bool gt_deskeysfileexists(const GtStr *indexname)
{
  return indexfilealreadyexists(indexname,KEYSTABSUFFIX);
}

static unsigned long searchfastaqueryindes(const char *extractkey,
                                           const char *keytab,
                                           unsigned long numofkeys,
                                           unsigned long keysize)
{
  unsigned long left = 0, right = numofkeys - 1, mid;
  int cmp;

  while (left <= right)
  {
    mid = left + GT_DIV2((unsigned long) (right-left));
    cmp = strcmp(extractkey,keytab + 1UL + mid * (keysize+1));
    if (cmp < 0)
    {
      gt_assert(mid > 0);
      right = mid-1;
    } else
    {
      if (cmp > 0)
      {
        left = mid+1;
      } else
      {
        gt_assert(mid < numofkeys);
        return mid;
      }
    }
  }
  return numofkeys;
}

static int itersearchoverallkeys(const Encodedsequence *encseq,
                                 const char *keytab,
                                 unsigned long numofkeys,
                                 unsigned long keysize,
                                 const GtStr *fileofkeystoextract,
                                 unsigned long linewidth,
                                 GtError *err)
{
  FILE *fp;
  GtStr *currentline;
  uint64_t linenum;
  char *extractkey;
  unsigned long seqnum, desclen, countmissing = 0;
  bool haserr = false;
  const char *desc;
  Seqinfo seqinfo;

  if (linewidth == 0)
  {
    gt_error_set(err,"use option width to specify line width for formatting");
    return -1;
  }
  fp = gt_fa_fopen(gt_str_get(fileofkeystoextract),"r",err);
  if (fp == NULL)
  {
    return -1;
  }
  currentline = gt_str_new();
  extractkey = gt_malloc(sizeof(char) * (keysize+1));
  for (linenum = 0; gt_str_read_next_line(currentline, fp) != EOF; linenum++)
  {
    char *lineptr = gt_str_get(currentline);
    size_t idx;

    for (idx = 0; lineptr[idx] != '\0' && !isspace(lineptr[idx]); idx++)
      /* nothing */ ;
    if (idx != (size_t) keysize)
    {
      gt_error_set(err,"key \"%*.*s\" is not of size %lu",(int) idx,
                   (int) idx,lineptr,keysize);
      haserr = true;
      break;
    }
    strncpy(extractkey,lineptr,idx);
    extractkey[idx] = '\0';
    seqnum = searchfastaqueryindes(extractkey,keytab,numofkeys,keysize);
    if (seqnum < numofkeys)
    {
      desc = retrievesequencedescription(&desclen, encseq, seqnum);
      (void) fputc('>',stdout);
      if (fwrite(desc,sizeof *desc,(size_t) desclen,stdout) != (size_t) desclen)
      {
        gt_error_set(err,"cannot write header of length %lu",desclen);
        haserr = true;
        break;
      }
      (void) fputc('\n',stdout);
      getencseqSeqinfo(&seqinfo,encseq,seqnum);
      encseq2symbolstring(stdout,
                          encseq,
                          Forwardmode,
                          seqinfo.seqstartpos,
                          seqinfo.seqlength,
                          linewidth);
    } else
    {
      countmissing++;
    }
    gt_str_reset(currentline);
  }
  if (countmissing > 0)
  {
    printf("# number of unsatified fastakey-queries: %lu\n",countmissing);
  }
  gt_str_delete(currentline);
  gt_fa_fclose(fp);
  gt_free(extractkey);
  return haserr ? - 1 : 0;
}

static int readkeysize(const GtStr *indexname,GtError *err)
{
  FILE *fp;
  bool haserr = false;
  char cc;

  gt_error_check(err);
  fp = opensfxfile(indexname,KEYSTABSUFFIX,"rb",err);
  if (fp == NULL)
  {
    haserr = true;
  }
  if (!haserr)
  {
    size_t ret;

    ret = fread(&cc,sizeof cc, (size_t) 1, fp);
    if (ferror(fp))
    {
      gt_error_set(err,"error when trying to read first byte of file %s%s: %s",
                   gt_str_get(indexname),KEYSTABSUFFIX,strerror(errno));
      haserr = true;
    }
  }
  gt_assert(cc >= 0);
  gt_fa_xfclose(fp);
  return haserr ? -1 : (int) cc;
}

int gt_extractkeysfromfastaindex(const GtStr *indexname,
                                 const GtStr *fileofkeystoextract,
                                 unsigned long linewidth,GtError *err)
{
  Encodedsequence *encseq = NULL;
  bool haserr = false;
  unsigned long numofdbsequences = 0, keysize = 0;

  encseq = mapencodedsequence(false,
                              indexname,
                              true,
                              true,
                              true,
                              true,
                              NULL,
                              err);
  if (encseq == NULL)
  {
    haserr = true;
  }
  if (!haserr)
  {
    int retval;

    numofdbsequences = getencseqnumofdbsequences(encseq);
    retval = readkeysize(indexname,err);
    if (retval < 0)
    {
      haserr = true;
    }
    keysize = (unsigned long) retval;
  }
  if (!haserr)
  {
    char *keytab;
    unsigned long keytablength;

    keytablength = 1UL + numofdbsequences * (keysize+1);
    keytab = genericmaptable(indexname,
                             KEYSTABSUFFIX,
                             keytablength,
                             sizeof (GtUchar),
                             err);
    if (keytab == NULL)
    {
      haserr = true;
    } else
    {
      if (itersearchoverallkeys(encseq,keytab,numofdbsequences,
                                keysize,fileofkeystoextract,
                                linewidth,err) != 0)
      {
        haserr = true;
      }
    }
    gt_fa_xmunmap(keytab);
  }
  if (encseq != NULL)
  {
    encodedsequence_free(&encseq);
  }
  return haserr ? -1 : 0;
}

int gt_extractkeysfromfastafile(bool verbose,
                                GtFile *outfp,
                                unsigned long width,
                                const GtStr *fileofkeystoextract,
                                GtStrArray *referencefiletab,
                                GtError *err)
{
  GtSeqIterator *seqit;
  const GtUchar *sequence;
  char *desc, *headerbufferspace = NULL;
  const char *keyptr;
  char *keyspace = NULL;
  unsigned long allockeyspace = 0, len, keylen, numofqueries, keyposition,
                countmarkhit = 0;
  int had_err = 0;
  off_t totalsize;
  Fastakeyquery *fastakeyqueries;
  size_t headerbuffersize = 0, headerlength;

  gt_error_check(err);
  fastakeyqueries = readfileofkeystoextract(verbose,&numofqueries,
                                            fileofkeystoextract,err);
  if (fastakeyqueries == NULL)
  {
    return -1;
  }
  totalsize = gt_files_estimate_total_size(referencefiletab);
  if (verbose)
  {
    printf("# estimated total size is " Formatuint64_t "\n",
            PRINTuint64_tcast(totalsize));
  }
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
    keyptr = desc2key(&keylen,desc,err);
    if (keyptr == NULL)
    {
      had_err = -1;
    } else
    {
      if (allockeyspace < keylen)
      {
        keyspace = gt_realloc(keyspace,sizeof(*keyspace) * (keylen+1));
        allockeyspace = keylen;
      }
      gt_assert(keyspace != NULL);
      strncpy(keyspace,keyptr,(size_t) keylen);
      keyspace[keylen] = '\0';
      keyposition = searchdesinfastakeyqueries(keyspace,fastakeyqueries,
                                               numofqueries);
      if (keyposition < numofqueries)
      {
        while (keyposition < numofqueries &&
               strcmp(fastakeyqueries[keyposition].fastakey,keyspace) == 0)
        {
#ifndef NDEBUG
          if (fastakeyqueries[keyposition].markhit)
          {
            fprintf(stderr,"key %s was already found before\n",
                     fastakeyqueries[keyposition].fastakey);
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
          if (COMPLETE(fastakeyqueries[keyposition]))
          {
            /*
            (void) snprintf(headerbufferspace,headerbuffersize,
                            "%*.*s complete %s",
                            (int) keylen,(int) keylen,keyspace,
                            desc);
            */
            gt_fasta_show_entry_generic(desc,
                                        (const char *) sequence,
                                        len, width, outfp);
          } else
          {
            (void) snprintf(headerbufferspace,headerbuffersize,
                            "%*.*s %lu %lu %s",
                            (int) keylen,(int) keylen,keyspace,
                            fastakeyqueries[keyposition].frompos,
                            fastakeyqueries[keyposition].topos,
                            desc);
            gt_fasta_show_entry_generic(headerbufferspace,
                                        (const char *)
                                        (sequence+fastakeyqueries[keyposition].
                                                                  frompos - 1),
                                        fastakeyqueries[keyposition].topos -
                                        fastakeyqueries[keyposition].frompos+1,
                                        width, outfp);
          }
          fastakeyqueries[keyposition].markhit = true;
          countmarkhit++;
          keyposition++;
        }
      }
#ifdef SKDEBUG
      printf("%s 1 %lu\n",keyspace, len);
#endif
    }
    gt_free(desc);
  }
  gt_free(headerbufferspace);
  gt_free(keyspace);
  if (verbose)
  {
    gt_progressbar_stop();
  }
  if (verbose)
  {
    outputnonmarked(fastakeyqueries,numofqueries);
  }
  fastakeyqueries_delete(fastakeyqueries,numofqueries);
  gt_seqiterator_delete(seqit);
  return had_err;
}
