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
#ifndef S_SPLINT_S
#include <ctype.h>
#endif
#include "core/fileutils.h"
#include "core/fa.h"
#include "core/assert_api.h"
#include "core/divmodmul.h"
#include "core/error.h"
#include "core/fasta.h"
#include "core/fileutils.h"
#include "core/format64.h"
#include "core/ma_api.h"
#include "core/progressbar.h"
#include "core/seq_iterator_sequence_buffer_api.h"
#include "core/xansi_api.h"
#include "giextract.h"
#include "echoseq.h"

#define COMPLETE(VALUE)\
        ((VALUE)->frompos == 1UL && (VALUE)->topos == 0)

#define EXTRABUF 128

#define CHECKPOSITIVE(VAL,FORMAT,WHICH)\
        if ((VAL) <= 0)\
        {\
          gt_error_set(err,"file \"%s\", line " Formatuint64_t \
                           ": illegal format: %s element = " FORMAT \
                           " is not a positive integer",\
                           gt_str_get(fileofkeystoextract),\
                           PRINTuint64_tcast(linenum+1),\
                           WHICH,\
                           VAL);\
          haserr = true;\
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

static int extractkeyfromcurrentline(Fastakeyquery *fastakeyptr,
                                     unsigned long keysize,
                                     const GtStr *currentline,
                                     uint64_t linenum,
                                     const GtStr *fileofkeystoextract,
                                     GtError *err)
{
  char *lineptr = gt_str_get(currentline);
  long readlongfrompos, readlongtopos;
  size_t idx;
  bool haserr = false;

  for (idx = 0; lineptr[idx] != '\0' && !isspace(lineptr[idx]); idx++)
    /* nothing */ ;
  if (keysize == 0)
  {
    fastakeyptr->fastakey = gt_malloc(sizeof (char) * (idx+1));
  } else
  {
    if (idx != (size_t) keysize)
    {
      gt_error_set(err,"key \"%*.*s\" is not of size %lu",(int) idx,
                   (int) idx,lineptr,keysize);
      haserr = true;
    }
  }
  if (!haserr)
  {
    strncpy(fastakeyptr->fastakey,lineptr,idx);
    fastakeyptr->fastakey[idx] = '\0';
    fastakeyptr->frompos = 1UL;
    fastakeyptr->topos = 0;
    if (sscanf(lineptr+idx,"%ld %ld\n",&readlongfrompos,&readlongtopos) == 2)
    {
      CHECKPOSITIVE(readlongfrompos,"%ld","second");
      if (!haserr)
      {
        fastakeyptr->frompos = (unsigned long) readlongfrompos;
        CHECKPOSITIVE(readlongtopos,"%ld","third");
      }
      if (!haserr)
      {
        fastakeyptr->topos = (unsigned long) readlongtopos;
      }
    }
  }
  if (!haserr)
  {
    fastakeyptr->markhit = false;
    if (!COMPLETE(fastakeyptr) && fastakeyptr->frompos > fastakeyptr->topos)
    {
      gt_error_set(err, "file \"%s\", line " Formatuint64_t
                        "illegal format: second value "
                        "%lu is larger than third value %lu",
                        gt_str_get(fileofkeystoextract),
                        PRINTuint64_tcast(linenum+1),
                        fastakeyptr->frompos,
                        fastakeyptr->topos);
      haserr = true;
    }
  }
  return haserr ? -1 : 0;
}

static Fastakeyquery *readfileofkeystoextract(bool verbose,
                                              unsigned long *numofqueries,
                                              const GtStr *fileofkeystoextract,
                                              GtError *err)
{
  FILE *fp;
  GtStr *currentline;
  bool haserr = false;
  uint64_t linenum;
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
    if (extractkeyfromcurrentline(fastakeyqueries + linenum,
                                  0,
                                  currentline,
                                  linenum,
                                  fileofkeystoextract,
                                  err) != 0)
    {
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
      if (COMPLETE(fastakeyqueries + idx))
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

static int giextract_encodedseq2fasta(FILE *fpout,
                                      const GtEncseq *encseq,
                                      unsigned long seqnum,
                                      const Fastakeyquery *fastakeyquery,
                                      unsigned long linewidth,
                                      GT_UNUSED GtError *err)
{
  const char *desc;
  unsigned long desclen;
  bool haserr = false;

  desc = gt_encseq_description(encseq, &desclen, seqnum);
  gt_xfputc('>',fpout);
  if (fastakeyquery != NULL && !COMPLETE(fastakeyquery))
  {
    printf("%s %lu %lu ",fastakeyquery->fastakey,
                         fastakeyquery->frompos,
                         fastakeyquery->topos);
  }
  gt_xfwrite(desc,sizeof *desc,(size_t) desclen,fpout);
  if (!haserr)
  {
    unsigned long frompos, topos, seqstartpos, seqlength ;

    gt_xfputc('\n',fpout);
    seqstartpos = gt_encseq_seqstartpos(encseq, seqnum);
    seqlength = gt_encseq_seqlength(encseq, seqnum);
    if (fastakeyquery != NULL && !COMPLETE(fastakeyquery))
    {
      frompos = fastakeyquery->frompos-1;
      topos = fastakeyquery->topos - fastakeyquery->frompos + 1;
    } else
    {
      frompos = 0;
      topos = seqlength;
    }
    gt_encseq2symbolstring(fpout,
                           encseq,
                           GT_READMODE_FORWARD,
                           seqstartpos + frompos,
                           topos,
                           linewidth);
  }
  return haserr ? -1 : 0;
}

#define MAXFIXEDKEYSIZE 11

typedef struct
{
  char key[MAXFIXEDKEYSIZE+1];
  unsigned long seqnum;
} Fixedsizekey;

static int compareFixedkeys(const void *a,const void *b)
{
  const Fixedsizekey *fa = a;
  const Fixedsizekey *fb = b;
  return strcmp(((const Fixedsizekey *) fa)->key,
                ((const Fixedsizekey *) fb)->key);
}

#define GT_KEYSTABFILESUFFIX ".kys"

int gt_extractkeysfromdesfile(const char *indexname,
                              bool sortkeys,
                              GtLogger *logger,
                              GtError *err)
{
  FILE *fpin, *fpout = NULL;
  GtStr *line = NULL;
  const char *keyptr;
  unsigned long keylen, constantkeylen = 0, linenum;/* incorrectorder = 0;*/
  bool haserr = false, firstdesc = true;
  char *previouskey = NULL;
  Fixedsizekey *keytab = NULL, *keytabptr = NULL;
  GtEncseq *encseq = NULL;
  unsigned long numofentries = 0;
  const unsigned long linewidth = 60UL;

  fpin = gt_fa_fopen_with_suffix(indexname,GT_DESTABFILESUFFIX,"rb",err);
  if (fpin == NULL)
  {
    return -1;
  }
  if (!sortkeys)
  {
    fpout = gt_fa_fopen_with_suffix(indexname,GT_KEYSTABFILESUFFIX,"wb",err);
    if (fpout == NULL)
    {
      haserr = true;
    }
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
      if (keylen > (unsigned long) CHAR_MAX)
      {
        gt_error_set(err,"key \"%*.*s\" of length %lu not allowed; "
                         "no key must be larger than %d",
                          (int) keylen,(int) keylen,keyptr,keylen,CHAR_MAX);
        haserr = true;
        break;
      }
      constantkeylen = keylen;
      previouskey = gt_malloc(sizeof (char) * (constantkeylen+1));
      firstdesc = false;
      if (!sortkeys)
      {
        gt_xfputc((char) constantkeylen,fpout);
      } else
      {
        GtEncseqLoader *el;
        if (constantkeylen > (unsigned long) MAXFIXEDKEYSIZE)
        {
          gt_error_set(err,"key \"%*.*s\" of length %lu not allowed; "
                           "no key must be larger than %d",
                            (int) keylen,(int) keylen,keyptr,keylen,
                            MAXFIXEDKEYSIZE);
          haserr = true;
          break;
        }
        el = gt_encseq_loader_new();
        gt_encseq_loader_set_logger(el, logger);
        encseq = gt_encseq_loader_load(el, indexname, err);
        gt_encseq_loader_delete(el);
        if (encseq == NULL)
        {
          haserr = true;
          break;
        }
        numofentries = gt_encseq_num_of_sequences(encseq);
        gt_assert(numofentries > 0);
        keytab = gt_malloc(sizeof (*keytab) * numofentries);
        keytabptr = keytab;
      }
    } else
    {
      if (constantkeylen != keylen)
      {
        gt_error_set(err,"key \"%*.*s\" of length %lu: all keys must be of "
                         "the same length which for all previously seen "
                         "headers is %lu",
                         (int) keylen,(int) keylen,keyptr,keylen,
                         constantkeylen);
        haserr = true;
        break;
      }
      gt_assert(previouskey != NULL);
      if (!sortkeys && strncmp(previouskey,keyptr,(size_t) constantkeylen) >= 0)
      {
        gt_error_set(err,"previous key \"%s\" is not lexicographically smaller "
                         "than current key \"%*.*s\"",
                         previouskey,(int) keylen,(int) keylen,keyptr);
        haserr = true;
        break;
        /*
        printf("previous key \"%s\" (no %lu) is lexicographically larger "
               "than current key \"%*.*s\"\n",
               previouskey,linenum,(int) keylen,(int) keylen,keyptr);
        incorrectorder++;
        */
      }
    }
    if (!sortkeys)
    {
      gt_xfwrite(keyptr,sizeof *keyptr,(size_t) keylen,fpout);
      gt_xfputc('\0',fpout);
    } else
    {
      gt_assert(keytabptr != NULL);
      strncpy(keytabptr->key,keyptr,(size_t) constantkeylen);
      keytabptr->key[constantkeylen] = '\0';
      keytabptr->seqnum = linenum;
      keytabptr++;
    }
    strncpy(previouskey,keyptr,(size_t) constantkeylen);
    previouskey[constantkeylen] = '\0';
    gt_str_reset(line);
  }
  if (!haserr)
  {
    gt_logger_log(logger,"number of keys of length %lu = %lu",
                constantkeylen,linenum);
    /*
    gt_logger_log(logger,"number of incorrectly ordered keys = %lu",
                incorrectorder);
    */
  }
  gt_str_delete(line);
  gt_fa_fclose(fpin);
  gt_fa_fclose(fpout);
  gt_free(previouskey);
  if (!haserr && sortkeys)
  {
    gt_assert(keytabptr != NULL);
    gt_assert(numofentries > 0);
    gt_assert(keytabptr == keytab + numofentries);
    qsort(keytab,(size_t) numofentries,sizeof (*keytab),compareFixedkeys);
    gt_assert(keytabptr != NULL);
    for (keytabptr = keytab; !haserr && keytabptr < keytab + numofentries;
         keytabptr++)
    {
      if (giextract_encodedseq2fasta(stdout,
                                     encseq,
                                     keytabptr->seqnum,
                                     NULL,
                                     linewidth,
                                     err) != 0)
      {
        haserr = true;
        break;
      }
    }
  }
  if (encseq != NULL)
  {
    gt_encseq_delete(encseq);
    encseq = NULL;
  }
  gt_free(keytab);
  return haserr ? -1 : 0;
}

bool gt_deskeysfileexists(const char *indexname)
{
  return gt_file_exists_with_suffix(indexname,GT_KEYSTABFILESUFFIX);
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

static int itersearchoverallkeys(const GtEncseq *encseq,
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
  unsigned long seqnum, countmissing = 0;
  bool haserr = false;
  Fastakeyquery fastakeyquery;

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
  fastakeyquery.fastakey = gt_malloc(sizeof (char) * (keysize+1));
  for (linenum = 0; gt_str_read_next_line(currentline, fp) != EOF; linenum++)
  {
    if (extractkeyfromcurrentline(&fastakeyquery,
                                  keysize,
                                  currentline,
                                  linenum,
                                  fileofkeystoextract,
                                  err) != 0)
    {
      haserr = true;
      break;
    }
    seqnum = searchfastaqueryindes(fastakeyquery.fastakey,keytab,numofkeys,
                                   keysize);
    if (seqnum < numofkeys)
    {
      if (giextract_encodedseq2fasta(stdout,
                                     encseq,
                                     seqnum,
                                     &fastakeyquery,
                                     linewidth,
                                     err) != 0)
      {
        haserr = true;
        break;
      }
    } else
    {
      countmissing++;
    }
    gt_str_reset(currentline);
  }
  if (!haserr && countmissing > 0)
  {
    printf("# number of unsatified fastakey-queries: %lu\n",countmissing);
  }
  gt_str_delete(currentline);
  gt_fa_fclose(fp);
  gt_free(fastakeyquery.fastakey);
  return haserr ? - 1 : 0;
}

static int readkeysize(const char *indexname,GtError *err)
{
  FILE *fp;
  bool haserr = false;
  char cc;

  gt_error_check(err);
  fp = gt_fa_fopen_with_suffix(indexname,GT_KEYSTABFILESUFFIX,"rb",err);
  if (fp == NULL)
  {
    haserr = true;
  }
  if (!haserr)
  {
    GT_UNUSED size_t ret;

    ret = fread(&cc,sizeof cc, (size_t) 1, fp);
    if (ferror(fp))
    {
      gt_error_set(err,"error when trying to read first byte of file %s%s: %s",
                   indexname,GT_KEYSTABFILESUFFIX,strerror(errno));
      haserr = true;
    }
  }
  gt_assert(cc >= 0);
  gt_fa_xfclose(fp);
  return haserr ? -1 : (int) cc;
}

int gt_extractkeysfromfastaindex(const char *indexname,
                                 const GtStr *fileofkeystoextract,
                                 unsigned long linewidth,GtError *err)
{
  GtEncseq *encseq = NULL;
  GtEncseqLoader *el = NULL;
  bool haserr = false;
  unsigned long numofdbsequences = 0, keysize = 0;

  el = gt_encseq_loader_new();
  encseq = gt_encseq_loader_load(el, indexname, err);
  gt_encseq_loader_delete(el);
  if (encseq == NULL)
  {
    haserr = true;
  }
  if (!haserr)
  {
    int retval;

    numofdbsequences = gt_encseq_num_of_sequences(encseq);
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
    keytab = gt_fa_mmap_check_size_with_suffix(indexname,
                                               GT_KEYSTABFILESUFFIX,
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
    gt_encseq_delete(encseq);
    encseq = NULL;
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
  char *desc, *headerbufferspace = NULL, *keyspace = NULL;
  const char *keyptr;
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
  seqit = gt_seq_iterator_sequence_buffer_new(referencefiletab, err);
  if (!seqit)
  {
    had_err = -1;
  }
  if (!had_err && verbose)
  {
    gt_progressbar_start(gt_seq_iterator_getcurrentcounter(seqit,
                                                          (unsigned long long)
                                                          totalsize),
                                                          (unsigned long long)
                                                          totalsize);
  }
  while (had_err != -1 && countmarkhit < numofqueries)
  {
    had_err = gt_seq_iterator_next(seqit, &sequence, &len, &desc, err);
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
        keyspace = gt_realloc(keyspace,sizeof (*keyspace) * (keylen+1));
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
          if (COMPLETE(fastakeyqueries + keyposition))
          {
            /*
            (void) snprintf(headerbufferspace,headerbuffersize,
                            "%*.*s complete %s",
                            (int) keylen,(int) keylen,keyspace,
                            desc);
            */
            gt_fasta_show_entry(desc, (const char *) sequence, len, width,
                                outfp);
          } else
          {
            (void) snprintf(headerbufferspace,headerbuffersize,
                            "%*.*s %lu %lu %s",
                            (int) keylen,(int) keylen,keyspace,
                            fastakeyqueries[keyposition].frompos,
                            fastakeyqueries[keyposition].topos,
                            desc);
            gt_fasta_show_entry(headerbufferspace,
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
  gt_seq_iterator_delete(seqit);
  return had_err;
}
