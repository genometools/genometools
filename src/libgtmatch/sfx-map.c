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
#include "types.h"
#include "arraydef.h"
#include "sarr-def.h"
#include "encseq-def.h"
#include "esafileend.h"
#include "sfx-ri-def.h"
#include "stamp.h"

#include "readnextline.pr"
#include "endianess.pr"
#include "alphabet.pr"
#include "opensfxfile.pr"
#include "sfx-readint.pr"

#define DBFILEKEY "dbfile="

#define INITBufferedfile(INDEXNAME,STREAM,SUFFIX)\
        (STREAM)->fp = opensfxfile(INDEXNAME,SUFFIX,"rb",env);\
        if((STREAM)->fp == NULL)\
        {\
          return -1;\
        }\
        (STREAM)->nextread = 0;\
        (STREAM)->nextfree = 0

static int scanprjfileviafileptr(Suffixarray *suffixarray,
                                 Seqpos *totallength,
                                 const Str *indexname,FILE *fpin,Env *env)
{
  ArrayUchar linebuffer;
  uint32_t integersize, littleendian, linenum;
  unsigned long numofsequences, numoffiles = 0, numofallocatedfiles = 0;
  DefinedSeqpos maxbranchdepth;
  size_t dbfilelen = strlen(DBFILEKEY);
  bool haserr = false;
  Array *riktab;

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
  setreadintkeys(riktab,"numofquerysequences",&suffixarray->numofdbsequences,
                 0,NULL,env);
  SETREADINTKEYS("longest",&suffixarray->longest.value,
                           &suffixarray->longest.defined);
  SETREADINTKEYS("prefixlength",&suffixarray->prefixlength,NULL);
  SETREADINTKEYS("largelcpvalues",&suffixarray->numoflargelcpvalues.value,
                 &suffixarray->numoflargelcpvalues.defined);
  SETREADINTKEYS("maxbranchdepth",&maxbranchdepth.value,
                 &maxbranchdepth.defined);
  SETREADINTKEYS("integersize",&integersize,NULL);
  SETREADINTKEYS("littleendian",&littleendian,NULL);
  INITARRAY(&linebuffer,Uchar);
  suffixarray->filenametab = strarray_new(env);
  suffixarray->filelengthtab = NULL;
  for (linenum = 0; /* Nothing */; linenum++)
  {
    linebuffer.nextfreeUchar = 0;
    if (readnextline(fpin,&linebuffer,env) == EOF)
    {
      break;
    }
    if (dbfilelen <= (size_t) linebuffer.nextfreeUchar &&
       memcmp(DBFILEKEY,linebuffer.spaceUchar,dbfilelen) == 0)
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
      ALLOCASSIGNSPACE(tmpfilename,NULL,char,linebuffer.nextfreeUchar);
      if (sscanf((const char *) linebuffer.spaceUchar,
                  "dbfile=%s " FormatScanint64_t " " FormatScanint64_t "\n",
                   tmpfilename,
                   Scanuint64_tcast(&readint1),
                   Scanuint64_tcast(&readint2)) != 3)
      {
        env_error_set(env,"cannot parse line %*.*s",
                          (int) linebuffer.nextfreeUchar,
                          (int) linebuffer.nextfreeUchar,
                          (const char *) linebuffer.spaceUchar);
        FREESPACE(tmpfilename);
        FREESPACE(suffixarray->filelengthtab);
        haserr = true;
        break;
      }
      if (readint1 < (int64_t) 1 || readint2 < (int64_t) 1)
      {
        env_error_set(env,"need positive integers in line %*.*s",
                          (int) linebuffer.nextfreeUchar,
                          (int) linebuffer.nextfreeUchar,
                          (const char *) linebuffer.spaceUchar);
        FREESPACE(tmpfilename);
        FREESPACE(suffixarray->filelengthtab);
        haserr = true;
      }
      if (!haserr)
      {
        strarray_add_cstr(suffixarray->filenametab,tmpfilename,env);
        FREESPACE(tmpfilename);
        assert(suffixarray->filelengthtab != NULL);
        suffixarray->filelengthtab[numoffiles].length = (Seqpos) readint1;
        suffixarray->filelengthtab[numoffiles].effectivelength 
                                               = (Seqpos) readint2;
        printf("%s%s " FormatSeqpos " " FormatSeqpos "\n",
                DBFILEKEY,
                strarray_get(suffixarray->filenametab,numoffiles),
                PRINTSeqposcast(suffixarray->filelengthtab[numoffiles].length),
                PRINTSeqposcast(suffixarray->filelengthtab[numoffiles].
                                effectivelength));
        numoffiles++;
      }
    } else
    {
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
  }
  if (!haserr && allkeysdefined(indexname,riktab,env) != 0)
  {
    haserr = true;
  }
  if (!haserr && 
      integersize != (uint32_t) 32 && 
      integersize != (uint32_t) 64)
  {
    env_error_set(env,"%s%s contains illegal line defining the integer size",
                  str_get(indexname),".prj");
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
    if (islittleendian())
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
  FREEARRAY(&linebuffer,Uchar);
  array_delete(riktab,env);
  return haserr ? -1 : 0;
}

static void *genericmaptable(const Str *indexname,const char *suffix,
                             Seqpos expectedunits,size_t sizeofunit,
                             Env *env)
{
  size_t numofbytes;
  Str *tmpfilename;
  void *ptr;
  bool haserr = false;

  tmpfilename = str_clone(indexname,env);
  str_append_cstr(tmpfilename,suffix,env);
  ptr = env_fa_mmap_read(env,str_get(tmpfilename),&numofbytes);
  if (ptr == NULL)
  {
    env_error_set(env,"cannot map file \"%s\": %s",str_get(tmpfilename),
                  strerror(errno));
    haserr = true;
  }
  if (!haserr && expectedunits != (Seqpos) (numofbytes/sizeofunit))
  {
    env_error_set(env,"number of mapped integers = %lu != " FormatSeqpos
                      " = expected number of integers",
                         (unsigned long) (numofbytes/sizeofunit),
                         PRINTSeqposcast(expectedunits));
    haserr = true;
  }
  str_delete(tmpfilename,env);
  return haserr ? NULL : ptr;
}

int initUcharBufferedfile(UcharBufferedfile *stream,
                          const Str *indexname,
                          const char *suffix,Env *env)
{
  INITBufferedfile(indexname,stream,suffix);
  return 0;
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
  suffixarray->alpha = NULL;
  suffixarray->bwttabstream.fp = NULL;
  suffixarray->suftabstream.fp = NULL;
  suffixarray->llvtabstream.fp = NULL;
  suffixarray->lcptabstream.fp = NULL;
}

static bool scanprjfile(Suffixarray *suffixarray,Seqpos *totallength,
                        const Str *indexname,Env *env)
{
  bool haserr = false;
  FILE *fp;

  fp = opensfxfile(indexname,".prj","rb",env);
  if (fp == NULL)
  {
    haserr = true;
  }
  if (!haserr && scanprjfileviafileptr(suffixarray,totallength,
                                       indexname,fp,env) != 0)
  {
    haserr = true;
  }
  env_fa_xfclose(fp,env);
  return haserr;
}

static bool scanal1file(Suffixarray *suffixarray,const Str *indexname,Env *env)
{
  Str *tmpfilename;
  bool haserr = false;

  tmpfilename = str_clone(indexname,env);
  str_append_cstr(tmpfilename,ALPHATABSUFFIX,env);
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

int mapsuffixarray(Suffixarray *suffixarray,
                   unsigned int demand,
                   const Str *indexname,
                   Env *env)
{
  bool haserr = false;
  Seqpos totallength;

  env_error_check(env);
  initsuffixarray(suffixarray);
  haserr = scanprjfile(suffixarray,&totallength,indexname,env);
  if (!haserr)
  {
    haserr = scanal1file(suffixarray,indexname,env);
  }
  if (!haserr && (demand & SARR_ESQTAB))
  {
    suffixarray->encseq = initencodedseq(true,
					 NULL,
					 indexname,
					 totallength,
					 &suffixarray->specialcharinfo,
					 suffixarray->alpha,
                                         false,  /* not relevant, since
                                                 indexname != NULL*/
                                         NULL,
					 env);
    if (suffixarray->encseq == NULL)
    {
      haserr = true;
    }
  }
  if (!haserr && (demand & SARR_SUFTAB))
  {
    suffixarray->suftab = genericmaptable(indexname,SUFTABSUFFIX,
                                          totallength+1,
                                          sizeof(Seqpos),
                                          env);
    if (suffixarray->suftab == NULL)
    {
      haserr = true;
    }
  }
  if(!haserr && (demand & SARR_LCPTAB))
  {
    suffixarray->lcptab = genericmaptable(indexname,LCPTABSUFFIX,
                                          totallength+1,
                                          sizeof(Uchar),
                                          env);
    if (suffixarray->lcptab == NULL)
    {
      haserr = true;
    }
    if(!haserr && !suffixarray->numoflargelcpvalues.defined)
    {
      env_error_set(env,"numoflargelcpvalues not defined");
      haserr = true;
    }
    if(!haserr && suffixarray->numoflargelcpvalues.value > 0)
    {
      suffixarray->llvtab = genericmaptable(indexname,LARGELCPTABSUFFIX,
                                            suffixarray->numoflargelcpvalues.
                                            value,
                                            sizeof(Largelcpvalue),
                                            env);
      if (suffixarray->llvtab == NULL)
      {
        haserr = true;
      }
    }
  }
  if(!haserr && (demand & SARR_BWTTAB))
  {
    suffixarray->bwttab = genericmaptable(indexname,BWTTABSUFFIX,
                                          totallength+1,
                                          sizeof(Uchar),
                                          env);
    if (suffixarray->bwttab == NULL)
    {
      haserr = true;
    }
  }
  return haserr ? -1 : 0;
}

int streamsuffixarray(Suffixarray *suffixarray,
                      unsigned int demand,
                      const Str *indexname,
                      Env *env)
{
  bool haserr = false;
  Seqpos totallength;

  env_error_check(env);
  initsuffixarray(suffixarray);
  haserr = scanprjfile(suffixarray,&totallength,indexname,env);
  if (!haserr)
  {
    haserr = scanal1file(suffixarray,indexname,env);
  }
  if (!haserr && (demand & SARR_ESQTAB))
  {
    suffixarray->encseq = initencodedseq(true,
					 NULL,
					 indexname,
					 totallength,
					 &suffixarray->specialcharinfo,
					 suffixarray->alpha,
                                         false, /* not relevant since 
                                                   indexname != NULL */
                                         NULL,
					 env);
    if (suffixarray->encseq == NULL)
    {
      haserr = true;
    }
  }
  if (!haserr && (demand & SARR_SUFTAB))
  {
    INITBufferedfile(indexname,&suffixarray->suftabstream,SUFTABSUFFIX);
    if(!suffixarray->longest.defined)
    {
      env_error_set(env,"longest not defined");
      haserr = true;
    }
  }
  if(!haserr && (demand & SARR_BWTTAB))
  {
    INITBufferedfile(indexname,&suffixarray->bwttabstream,BWTTABSUFFIX);
  }
  if(!haserr && (demand & SARR_LCPTAB))
  {
    INITBufferedfile(indexname,&suffixarray->lcptabstream,LCPTABSUFFIX);
    if(fseek(suffixarray->lcptabstream.fp,(long) sizeof(Uchar),SEEK_SET) != 0)
    {
      env_error_set(env,"fseek(esastream) failed: %s",strerror(errno));
      haserr = true;
    }
    if(!haserr && !suffixarray->numoflargelcpvalues.defined)
    {
      env_error_set(env,"numoflargelcpvalues not defined");
      haserr = true;
    }
    if(!haserr && suffixarray->numoflargelcpvalues.value > 0)
    {
      INITBufferedfile(indexname,&suffixarray->llvtabstream,LARGELCPTABSUFFIX);
    } 
  }
  return haserr ? -1 : 0;
}

void freesuffixarray(Suffixarray *suffixarray,Env *env)
{
  env_fa_xmunmap((void *) suffixarray->suftab,env);
  env_fa_xmunmap((void *) suffixarray->lcptab,env);
  env_fa_xmunmap((void *) suffixarray->llvtab,env);
  env_fa_xmunmap((void *) suffixarray->bwttab,env);
  env_fa_xfclose(suffixarray->suftabstream.fp,env);
  env_fa_xfclose(suffixarray->lcptabstream.fp,env);
  env_fa_xfclose(suffixarray->llvtabstream.fp,env);
  env_fa_xfclose(suffixarray->bwttabstream.fp,env);
  freeAlphabet(&suffixarray->alpha,env);
  freeEncodedsequence(&suffixarray->encseq,env);
  strarray_delete(suffixarray->filenametab,env);
  FREESPACE(suffixarray->filelengthtab);
}
