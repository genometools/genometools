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

#include <stdio.h>
#include <errno.h>
#include <stdbool.h>
#include <limits.h>
#include <string.h>
#include "libgtcore/array.h"
#include "libgtcore/endianess.h"
#include "libgtcore/env.h"
#include "libgtcore/fa.h"
#include "format64.h"
#include "sarr-def.h"
#include "encseq-def.h"
#include "esafileend.h"
#include "sfx-ri-def.h"
#include "spacedef.h"
#include "stamp.h"

#include "opensfxfile.pr"

#define DBFILEKEY "dbfile="

#define INITBufferedfile(INDEXNAME,STREAM,SUFFIX)\
        (STREAM)->fp = opensfxfile(INDEXNAME,SUFFIX,"rb",env);\
        if ((STREAM)->fp == NULL)\
        {\
          haserr = true;\
        } else\
        {\
          (STREAM)->nextread = 0;\
          (STREAM)->nextfree = 0;\
        }

static int scanprjfileviafileptr(Suffixarray *suffixarray,
                                 Seqpos *totallength,
                                 const Str *indexname,
                                 Verboseinfo *verboseinfo,
                                 FILE *fpin,
                                 Env *env)
{
  uint32_t integersize, littleendian, readmodeint;
  unsigned int linenum;
  unsigned long numofsequences, numofquerysequences,
                numoffiles = 0, numofallocatedfiles = 0, currentlinelength;
  DefinedSeqpos maxbranchdepth;
  size_t dbfilelen = strlen(DBFILEKEY);
  bool haserr = false;
  Array *riktab;
  Str *currentline;

  env_error_check(env);
  riktab = array_new(sizeofReadintkeys(),env);
  SETREADINTKEYS("totallength",totallength,NULL);
  SETREADINTKEYS("specialcharacters",
                 &suffixarray->specialcharinfo.specialcharacters,NULL);
  SETREADINTKEYS("specialranges",
                 &suffixarray->specialcharinfo.specialranges,NULL);
  SETREADINTKEYS("lengthofspecialprefix",
                 &suffixarray->specialcharinfo.lengthofspecialprefix,NULL);
  SETREADINTKEYS("lengthofspecialsuffix",
                 &suffixarray->specialcharinfo.lengthofspecialsuffix,NULL);
  SETREADINTKEYS("numofsequences",&numofsequences,NULL);
  SETREADINTKEYS("numofdbsequences",&suffixarray->numofdbsequences,NULL);
  setreadintkeys(riktab,"numofquerysequences",&numofquerysequences,
                 0,NULL,env);
  SETREADINTKEYS("longest",&suffixarray->longest.valueseqpos,
                           &suffixarray->longest.defined);
  SETREADINTKEYS("prefixlength",&suffixarray->prefixlength,NULL);
  SETREADINTKEYS("largelcpvalues",&suffixarray->numoflargelcpvalues.valueseqpos,
                 &suffixarray->numoflargelcpvalues.defined);
  SETREADINTKEYS("maxbranchdepth",&maxbranchdepth.valueseqpos,
                 &maxbranchdepth.defined);
  SETREADINTKEYS("integersize",&integersize,NULL);
  SETREADINTKEYS("littleendian",&littleendian,NULL);
  SETREADINTKEYS("readmode",&readmodeint,NULL);
  suffixarray->filenametab = strarray_new(env);
  suffixarray->filelengthtab = NULL;
  currentline = str_new(env);
  for (linenum = 0; str_read_next_line(currentline, fpin, env) != EOF;
       linenum++)
  {
    currentlinelength = str_length(currentline);
    if (dbfilelen <= (size_t) currentlinelength &&
       memcmp(DBFILEKEY,str_get(currentline),dbfilelen) == 0)
    {
      char *tmpfilename;
      int64_t readint1, readint2;

      if (numoffiles >= numofallocatedfiles)
      {
        numofallocatedfiles += 2;
        ALLOCASSIGNSPACE(suffixarray->filelengthtab,suffixarray->filelengthtab,
                         Filelengthvalues,numofallocatedfiles);
      }
      assert(suffixarray->filelengthtab != NULL);
      ALLOCASSIGNSPACE(tmpfilename,NULL,char,currentlinelength);
      if (sscanf((const char *) str_get(currentline),
                  "dbfile=%s " FormatScanint64_t " " FormatScanint64_t "\n",
                   tmpfilename,
                   Scanuint64_tcast(&readint1),
                   Scanuint64_tcast(&readint2)) != 3)
      {
        env_error_set(env,"cannot parse line %*.*s",
                          (int) currentlinelength,
                          (int) currentlinelength,
                          (const char *) str_get(currentline));
        FREESPACE(tmpfilename);
        FREESPACE(suffixarray->filelengthtab);
        haserr = true;
        break;
      }
      if (readint1 < (int64_t) 1 || readint2 < (int64_t) 1)
      {
        env_error_set(env,"need positive integers in line %*.*s",
                          (int) currentlinelength,
                          (int) currentlinelength,
                          (const char *) str_get(currentline));
        FREESPACE(tmpfilename);
        FREESPACE(suffixarray->filelengthtab);
        haserr = true;
        break;
      }
      if (!haserr)
      {
        strarray_add_cstr(suffixarray->filenametab,tmpfilename,env);
        FREESPACE(tmpfilename);
        assert(suffixarray->filelengthtab != NULL);
        suffixarray->filelengthtab[numoffiles].length = (Seqpos) readint1;
        suffixarray->filelengthtab[numoffiles].effectivelength
                                               = (Seqpos) readint2;
        showverbose(verboseinfo,
                    "%s%s " Formatuint64_t " " Formatuint64_t,
                    DBFILEKEY,
                    strarray_get(suffixarray->filenametab,numoffiles),
                    PRINTuint64_tcast(suffixarray->filelengthtab[numoffiles].
                                      length),
                    PRINTuint64_tcast(suffixarray->filelengthtab[numoffiles].
                                      effectivelength));
        numoffiles++;
      }
    } else
    {
      if (analyzeuintline(indexname,
                         PROJECTFILESUFFIX,
                         linenum,
                         str_get(currentline),
                         currentlinelength,
                         riktab,
                         env) != 0)
      {
        haserr = true;
        break;
      }
    }
    str_reset(currentline);
  }
  str_delete(currentline,env);
  if (!haserr && allkeysdefined(indexname,PROJECTFILESUFFIX,riktab,
                                verboseinfo,env) != 0)
  {
    haserr = true;
  }
  if (!haserr &&
      integersize != (uint32_t) 32 &&
      integersize != (uint32_t) 64)
  {
    env_error_set(env,"%s%s contains illegal line defining the integer size",
                  str_get(indexname),PROJECTFILESUFFIX);
    haserr = true;
  }
  if (!haserr && integersize != (uint32_t) (sizeof (Seqpos) * CHAR_BIT))
  {
    env_error_set(env,"index was generated for %u-bit integers while "
                      "this program uses %u-bit integers",
                      (unsigned int) integersize,
                      (unsigned int) (sizeof (Seqpos) * CHAR_BIT));
    haserr = true;
  }
  if (!haserr)
  {
    if (is_little_endian())
    {
      if (littleendian != (uint32_t) 1)
      {
        env_error_set(env,"computer has little endian byte order, while index "
                          "was build on computer with big endian byte order");
        haserr = true;
      }
    } else
    {
      if (littleendian == (uint32_t) 1)
      {
        env_error_set(env,"computer has big endian byte order, while index "
                          "was build on computer with little endian byte "
                          "order");
        haserr = true;
      }
    }
  }
  if (!haserr)
  {
    if (readmodeint > (uint32_t) 3)
    {
      env_error_set(env,"illegal readmode %u",(unsigned int) readmodeint);
      haserr = true;
    }
    suffixarray->readmode = (Readmode) readmodeint;
  }
  array_delete(riktab,env);
  return haserr ? -1 : 0;
}

static void *genericmaponlytable(const Str *indexname,const char *suffix,
                                 size_t *numofbytes,Env *env)
{
  Str *tmpfilename;
  void *ptr;
  bool haserr = false;

  env_error_check(env);
  tmpfilename = str_clone(indexname,env);
  str_append_cstr(tmpfilename,suffix,env);
  ptr = fa_mmap_read(str_get(tmpfilename),numofbytes);
  if (ptr == NULL)
  {
    env_error_set(env,"cannot map file \"%s\": %s",str_get(tmpfilename),
                  strerror(errno));
    haserr = true;
  }
  str_delete(tmpfilename,env);
  return haserr ? NULL : ptr;
}

static int checkmappedfilesize(size_t numofbytes,Seqpos expectedunits,
                               size_t sizeofunit,Env *env)
{
  env_error_check(env);
  if (expectedunits != (Seqpos) (numofbytes/sizeofunit))
  {
    env_error_set(env,"number of mapped units = %lu != " FormatSeqpos
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
                             Env *env)
{
  size_t numofbytes;

  void *ptr = genericmaponlytable(indexname,suffix,&numofbytes,env);
  if (ptr == NULL)
  {
    return NULL;
  }
  if (checkmappedfilesize(numofbytes,expectedunits,sizeofunit,env) != 0)
  {
    fa_xmunmap(ptr);
    return NULL;
  }
  return ptr;
}

static void initsuffixarray(Suffixarray *suffixarray)
{
  suffixarray->filelengthtab = NULL;
  suffixarray->filenametab = NULL;
  suffixarray->encseq = NULL;
  suffixarray->suftab = NULL;
  suffixarray->lcptab = NULL;
  suffixarray->llvtab = NULL;
  suffixarray->bwttab = NULL;
  suffixarray->destab = NULL;
  suffixarray->alpha = NULL;
  suffixarray->encseq = NULL;
  suffixarray->bwttabstream.fp = NULL;
  suffixarray->suftabstream.fp = NULL;
  suffixarray->llvtabstream.fp = NULL;
  suffixarray->lcptabstream.fp = NULL;
  suffixarray->destablength = 0;
}

static bool scanprjfile(Suffixarray *suffixarray,Seqpos *totallength,
                        const Str *indexname,Verboseinfo *verboseinfo,
                        Env *env)
{
  bool haserr = false;
  FILE *fp;

  env_error_check(env);
  fp = opensfxfile(indexname,PROJECTFILESUFFIX,"rb",env);
  if (fp == NULL)
  {
    haserr = true;
  }
  if (!haserr && scanprjfileviafileptr(suffixarray,totallength,
                                       indexname,verboseinfo,
                                       fp,env) != 0)
  {
    haserr = true;
  }
  fa_xfclose(fp);
  return haserr;
}

static bool scanal1file(Suffixarray *suffixarray,const Str *indexname,Env *env)
{
  Str *tmpfilename;
  bool haserr = false;

  env_error_check(env);
  tmpfilename = str_clone(indexname,env);
  str_append_cstr(tmpfilename,ALPHABETFILESUFFIX,env);
  suffixarray->alpha = assigninputalphabet(false,
                                           false,
                                           tmpfilename,
                                           NULL,
                                           env);
  if (suffixarray->alpha == NULL)
  {
    haserr = true;
  }
  str_delete(tmpfilename,env);
  return haserr;
}

void freesuffixarray(Suffixarray *suffixarray,Env *env)
{
  fa_xmunmap((void *) suffixarray->suftab);
  suffixarray->suftab = NULL;
  fa_xmunmap((void *) suffixarray->lcptab);
  suffixarray->lcptab = NULL;
  fa_xmunmap((void *) suffixarray->llvtab);
  suffixarray->llvtab = NULL;
  fa_xmunmap((void *) suffixarray->bwttab);
  suffixarray->bwttab = NULL;
  fa_xmunmap((void *) suffixarray->destab);
  suffixarray->destab = NULL;
  fa_xfclose(suffixarray->suftabstream.fp);
  suffixarray->suftabstream.fp = NULL;
  fa_xfclose(suffixarray->lcptabstream.fp);
  suffixarray->lcptabstream.fp = NULL;
  fa_xfclose(suffixarray->llvtabstream.fp);
  suffixarray->llvtabstream.fp = NULL;
  fa_xfclose(suffixarray->bwttabstream.fp);
  suffixarray->bwttabstream.fp = NULL;
  if (suffixarray->alpha != NULL)
  {
    freeAlphabet(&suffixarray->alpha,env);
  }
  freeEncodedsequence(&suffixarray->encseq,env);
  strarray_delete(suffixarray->filenametab,env);
  suffixarray->filenametab = NULL;
  FREESPACE(suffixarray->filelengthtab);
}

static int inputsuffixarray(bool map,
                            Suffixarray *suffixarray,
                            Seqpos *totallength,
                            unsigned int demand,
                            const Str *indexname,
                            Verboseinfo *verboseinfo,
                            Env *env)
{
  bool haserr = false;

  env_error_check(env);
  initsuffixarray(suffixarray);
  haserr = scanprjfile(suffixarray,totallength,indexname,verboseinfo,env);
  if (!haserr)
  {
    haserr = scanal1file(suffixarray,indexname,env);
  }
  if (!haserr && (demand & SARR_ESQTAB))
  {
    suffixarray->encseq = mapencodedsequence(true,
                                             indexname,
                                             *totallength,
                                             suffixarray->specialcharinfo.
                                                          specialranges,
                                             getmapsizeAlphabet(suffixarray->
                                                                alpha),
                                             verboseinfo,
                                             env);
    if (suffixarray->encseq == NULL)
    {
      haserr = true;
    }
  }
  if (!haserr && (demand & SARR_DESTAB))
  {
    size_t numofbytes;

    suffixarray->destab = genericmaponlytable(indexname,
                                              DESTABSUFFIX,
                                              &numofbytes,
                                              env);
    suffixarray->destablength = (unsigned long) numofbytes;
    if (suffixarray->destab == NULL)
    {
      haserr = true;
    }
  }
  if (!haserr && (demand & SARR_SUFTAB))
  {
    if (map)
    {
      suffixarray->suftab = genericmaptable(indexname,
                                            SUFTABSUFFIX,
                                            (*totallength)+1,
                                            sizeof (Seqpos),
                                            env);
      if (suffixarray->suftab == NULL)
      {
        haserr = true;
      }
    } else
    {
      INITBufferedfile(indexname,&suffixarray->suftabstream,SUFTABSUFFIX);
    }
    if (!haserr && !suffixarray->longest.defined)
    {
      env_error_set(env,"longest not defined");
      haserr = true;
    }
  }
  if (!haserr && (demand & SARR_LCPTAB))
  {
    if (map)
    {
      suffixarray->lcptab = genericmaptable(indexname,
                                            LCPTABSUFFIX,
                                            (*totallength)+1,
                                            sizeof (Uchar),
                                            env);
      if (suffixarray->lcptab == NULL)
      {
        haserr = true;
      }
    } else
    {
      INITBufferedfile(indexname,&suffixarray->lcptabstream,LCPTABSUFFIX);
      if (!haserr &&
          fseek(suffixarray->lcptabstream.fp,(long) sizeof (Uchar),SEEK_SET))
      {
        env_error_set(env,"fseek(esastream) failed: %s",strerror(errno));
        haserr = true;
      }
    }
    if (!haserr && !suffixarray->numoflargelcpvalues.defined)
    {
      env_error_set(env,"numoflargelcpvalues not defined");
      haserr = true;
    }
    if (!haserr && suffixarray->numoflargelcpvalues.valueseqpos > 0)
    {
      if (map)
      {
        suffixarray->llvtab = genericmaptable(indexname,
                                              LARGELCPTABSUFFIX,
                                              suffixarray->numoflargelcpvalues.
                                              valueseqpos,
                                              sizeof (Largelcpvalue),
                                              env);
        if (suffixarray->llvtab == NULL)
        {
          haserr = true;
        }
      } else
      {
        INITBufferedfile(indexname,&suffixarray->llvtabstream,
                         LARGELCPTABSUFFIX);
      }
    }
  }
  if (!haserr && (demand & SARR_BWTTAB))
  {
    if (map)
    {
      suffixarray->bwttab = genericmaptable(indexname,
                                            BWTTABSUFFIX,
                                            (*totallength)+1,
                                            sizeof (Uchar),
                                            env);
      if (suffixarray->bwttab == NULL)
      {
        haserr = true;
      }
    } else
    {
      INITBufferedfile(indexname,&suffixarray->bwttabstream,BWTTABSUFFIX);
    }
  }
  if (haserr)
  {
    freesuffixarray(suffixarray,env);
  }
  return haserr ? -1 : 0;
}

int streamsuffixarray(Suffixarray *suffixarray,
                      Seqpos *totallength,
                      unsigned int demand,
                      const Str *indexname,
                      Verboseinfo *verboseinfo,
                      Env *env)
{
  env_error_check(env);
  return inputsuffixarray(false,
                          suffixarray,
                          totallength,
                          demand,
                          indexname,
                          verboseinfo,
                          env);
}

int mapsuffixarray(Suffixarray *suffixarray,
                   Seqpos *totallength,
                   unsigned int demand,
                   const Str *indexname,
                   Verboseinfo *verboseinfo,
                   Env *env)
{
  env_error_check(env);
  return inputsuffixarray(true,
                          suffixarray,
                          totallength,
                          demand,
                          indexname,
                          verboseinfo,
                          env);
}
