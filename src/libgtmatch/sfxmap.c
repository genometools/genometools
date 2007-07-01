/*
  Copyright (c) 2007 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
  Copyright (c) 2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#include <stdio.h>
#include <errno.h>
#include <stdbool.h>
#include "libgtcore/env.h"
#include "types.h"
#include "arraydef.h"
#include "sarr-def.h"
#include "encseq-def.h"

#include "readnextline.pr"
#include "endianess.pr"
#include "alphabet.pr"

#define DBFILEKEY "dbfile="

#define SETREADINTKEYS(VALNAME,VAL)\
        setreadintkeys(&rik,VALNAME,VAL,sizeof(*(VAL)),env)

#ifdef S_SPLINT_S
#define FormatScanint64_t "%lu"
#else
#define FormatScanint64_t "%020" SCNd64
#endif

typedef union
{
  uint32_t smallvalue;
  uint64_t bigvalue;
} Smallorbigint;

typedef struct
{
  const char *keystring;
  uint32_t *smallvalueptr;
  uint64_t *bigvalueptr;
  bool ptrdefined, found;
} Readintkeys;

DECLAREARRAYSTRUCT(Readintkeys);

static void setreadintkeys(ArrayReadintkeys *riktab,
                           const char *keystring,
                           void *valueptr,
                           size_t sizeval,
                           Env *env)
{
  Readintkeys *riktabptr;

  GETNEXTFREEINARRAY(riktabptr,riktab,Readintkeys,1);
  riktabptr->keystring = keystring;
  switch(sizeval)
  {
    case 0: riktabptr->smallvalueptr = NULL;
            riktabptr->bigvalueptr = NULL;
            riktabptr->ptrdefined = false;
            break;
    case 4: assert(sizeof(uint32_t) == (size_t) 4);
            riktabptr->smallvalueptr = valueptr;
            riktabptr->bigvalueptr = NULL;
            riktabptr->ptrdefined = true;
            break;
    case 8: assert(sizeof(uint64_t) == (size_t) 8);
            riktabptr->bigvalueptr = valueptr;
            riktabptr->smallvalueptr = NULL;
            riktabptr->ptrdefined = true;
            break;
    default: fprintf(stderr,"sizeval must be 0 or 4 or 8\n");
             exit(EXIT_FAILURE);
  }
  riktabptr->found = false;
}

static int scanuintintline(uint32_t *lengthofkey,
                           Smallorbigint *smallorbigint,
                           const ArrayUchar *linebuffer,
                           Env *env)
{
  int64_t readint;
  unsigned long i;
  bool found = false;
  int retval = 0;

  env_error_check(env);
  for (i=0; i<linebuffer->nextfreeUchar; i++)
  {
    if (linebuffer->spaceUchar[i] == '=')
    {
      *lengthofkey = (uint32_t) i;
      found = true;

      if (sscanf((const char *) (linebuffer->spaceUchar + i + 1),
                 FormatScanint64_t,
                 Scanuint64_tcast(&readint)) != 1 ||
         readint < (int64_t) 0)
      {
        env_error_set(env,"cannot find non-negative integer in \"%*.*s\"",
                           (Fieldwidthtype) (linebuffer->nextfreeUchar - (i+1)),
                           (Fieldwidthtype) (linebuffer->nextfreeUchar - (i+1)),
                           (const char *) (linebuffer->spaceUchar + i + 1));
        return -1;
      }
      if(readint <= (int64_t) UINT32_MAX)
      {
        smallorbigint->smallvalue = (uint32_t) readint;
        retval = 0;
      } else
      {
        smallorbigint->bigvalue = (uint64_t) readint;
        retval = 1;
      }
      break;
    }
  }
  if (!found)
  {
    env_error_set(env,"missing equality symbol in \"%*.*s\"",
                       (Fieldwidthtype) linebuffer->nextfreeUchar,
                       (Fieldwidthtype) linebuffer->nextfreeUchar,
                       (const char *) linebuffer->spaceUchar);
    return -2;
  }
  return retval;
}

static int allkeysdefined(const Str *prjfile,const ArrayReadintkeys *rik,
                          Env *env)
{
  unsigned long i;

  env_error_check(env);
  for (i=0; i<rik->nextfreeReadintkeys; i++)
  {
    if (!rik->spaceReadintkeys[i].found)
    {
      env_error_set(env,"file %s: missing line beginning with \"%s=\"",
                         str_get(prjfile),
                         rik->spaceReadintkeys[i].keystring);
      return -1;
    }
    printf("%s=",rik->spaceReadintkeys[i].keystring);
    if (rik->spaceReadintkeys[i].ptrdefined)
    {
      if (rik->spaceReadintkeys[i].smallvalueptr != NULL)
      {
        printf("%u\n",(unsigned int) *(rik->spaceReadintkeys[i].smallvalueptr));
      } else
      {
        if (rik->spaceReadintkeys[i].bigvalueptr != NULL)
        {
          printf(Formatuint64_t "\n",
                 PRINTuint64_tcast(*(rik->spaceReadintkeys[i].bigvalueptr)));
        } else
        {
          assert(false);
        }
      }
    } else
    {
      printf("0\n");
    }
  }
  return 0;
}

static int scanprjfileviafileptr(Suffixarray *suffixarray,
                                 Seqpos *totallength,
                                 const Str *prjfile,FILE *fpin,Env *env)
{
  ArrayUchar linebuffer;
  uint32_t integersize, littleendian, linenum, lengthofkey;
  unsigned long i, numofsequences, numoffiles = 0, numofallocatedfiles = 0;
  size_t dbfilelen = strlen(DBFILEKEY);
  int retval;
  bool found, haserr = false;
  Smallorbigint smallorbigint;
  ArrayReadintkeys rik;

  env_error_check(env);
  INITARRAY(&rik,Readintkeys);
  SETREADINTKEYS("totallength",totallength);
  SETREADINTKEYS("specialcharacters",
                 &suffixarray->specialcharinfo.specialcharacters);
  SETREADINTKEYS("specialranges",
                 &suffixarray->specialcharinfo.specialranges);
  SETREADINTKEYS("lengthofspecialprefix",
                 &suffixarray->specialcharinfo.lengthofspecialprefix);
  SETREADINTKEYS("lengthofspecialsuffix",
                 &suffixarray->specialcharinfo.lengthofspecialsuffix);
  SETREADINTKEYS("numofsequences",
                 &numofsequences);
  SETREADINTKEYS("numofdbsequences",
                 &suffixarray->numofdbsequences);
  setreadintkeys(&rik,"numofquerysequences",&suffixarray->numofdbsequences,
                 0,env);
  SETREADINTKEYS("prefixlength",
                 &suffixarray->prefixlength);
  SETREADINTKEYS("integersize",
                 &integersize);
  SETREADINTKEYS("littleendian",
                 &littleendian);
  assert(rik.spaceReadintkeys != NULL);
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
                         PairSeqpos,numofallocatedfiles);
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
        haserr = true;
        break;
      }
      if (readint1 < (int64_t) 1 || readint2 < (int64_t) 1)
      {
        env_error_set(env,"need positive integers in line %*.*s",
                          (int) linebuffer.nextfreeUchar,
                          (int) linebuffer.nextfreeUchar,
                          (const char *) linebuffer.spaceUchar);
        haserr = true;
      }
      if (!haserr)
      {
        strarray_add_cstr(suffixarray->filenametab,tmpfilename,env);
        FREESPACE(tmpfilename);
        suffixarray->filelengthtab[numoffiles].uint0 = (Seqpos) readint1;
        suffixarray->filelengthtab[numoffiles].uint1 = (Seqpos) readint2;
        printf("%s%s " FormatSeqpos " " FormatSeqpos "\n",
                DBFILEKEY,
                strarray_get(suffixarray->filenametab,numoffiles),
                PRINTSeqposcast(suffixarray->filelengthtab[numoffiles].uint0),
                PRINTSeqposcast(suffixarray->filelengthtab[numoffiles].uint1));
        numoffiles++;
      }
    } else
    {
      retval = scanuintintline(&lengthofkey,
                               &smallorbigint,
                               &linebuffer,
                               env);
      if (retval < 0)
      {
        haserr = true;
        break;
      }
      found = false;
      for (i=0; i<rik.nextfreeReadintkeys; i++)
      {
        if (memcmp(linebuffer.spaceUchar,
                   rik.spaceReadintkeys[i].keystring,
                   (size_t) lengthofkey) == 0)
        {
          rik.spaceReadintkeys[i].found = true;
          if (rik.spaceReadintkeys[i].ptrdefined)
          {
            if(rik.spaceReadintkeys[i].smallvalueptr != NULL)
            {
              if(retval == 1)
              {
                env_error_set(env,"bigvalue %llu does not fit into %s",
                              smallorbigint.bigvalue,
                              rik.spaceReadintkeys[i].keystring);
                haserr = true;
                break;
              }
              *(rik.spaceReadintkeys[i].smallvalueptr) 
                = smallorbigint.smallvalue;
            } else
            {
              if(retval == 1)
              {
                *(rik.spaceReadintkeys[i].bigvalueptr) = smallorbigint.bigvalue;
              } else
              {
                *(rik.spaceReadintkeys[i].bigvalueptr) 
                  = (uint64_t) smallorbigint.smallvalue;
              }
            }
          }
          found = true;
          break;
        }
      }
      if (!found)
      {
        env_error_set(env,"file %s, line %u: cannot find key for \"%*.*s\"",
                           str_get(prjfile),
                           linenum,
                           (Fieldwidthtype) lengthofkey,
                           (Fieldwidthtype) lengthofkey,
                           (char  *) linebuffer.spaceUchar);
        haserr = true;
        break;
      }
    }
  }
  if (!haserr && allkeysdefined(prjfile,&rik,env) != 0)
  {
    haserr = true;
  }
  if (!haserr && 
      integersize != (uint32_t) 32 && 
      integersize != (uint32_t) 64)
  {
    env_error_set(env,"%s contains illegal line defining the integer size",
                  str_get(prjfile));
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
  FREEARRAY(&rik,Readintkeys);
  return haserr ? -1 : 0;
}

int mapsuffixarray(Suffixarray *suffixarray,
                   bool withsuftab,
                   bool withencseq,
                   const Str *indexname,
                   Env *env)
{
  Str *tmpfilename;
  FILE *fp;
  bool haserr = false;
  Seqpos totallength;

  env_error_check(env);
  suffixarray->suftab = NULL;
  suffixarray->alpha = NULL;
  tmpfilename = str_clone(indexname,env);
  str_append_cstr(tmpfilename,".prj",env);
  fp = env_fa_fopen(env,str_get(tmpfilename),"rb");
  if (fp == NULL)
  {
    env_error_set(env,"cannot open file \"%s\": %s",str_get(tmpfilename),
                                                    strerror(errno));
    haserr = true;
  }
  if (!haserr && scanprjfileviafileptr(suffixarray,&totallength,
                                      tmpfilename,fp,env) != 0)
  {
    haserr = true;
  }
  env_fa_xfclose(fp,env);
  str_delete(tmpfilename,env);
  if (!haserr)
  {
    tmpfilename = str_clone(indexname,env);
    str_append_cstr(tmpfilename,".al1",env);
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
  }
  if (!haserr && withencseq)
  {
    suffixarray->encseq = initencodedseq(true,
					 NULL,
					 indexname,
					 totallength,
					 &suffixarray->specialcharinfo,
					 suffixarray->alpha,
                                         NULL,
					 env);
    if (suffixarray->encseq == NULL)
    {
      haserr = true;
    }
  }
  if (!haserr && withsuftab)
  {
    size_t numofbytes;

    tmpfilename = str_clone(indexname,env);
    str_append_cstr(tmpfilename,".suf",env);
    suffixarray->suftab = env_fa_mmap_read(env,str_get(tmpfilename),
                                           &numofbytes);
    if (suffixarray->suftab == NULL)
    {
      env_error_set(env,"cannot map file \"%s\": %s",str_get(tmpfilename),
                    strerror(errno));
      haserr = true;
    }
    if (!haserr && totallength + 1 != (Seqpos) (numofbytes/sizeof (Seqpos)))
    {
      env_error_set(env,"number of mapped integers = %lu != " FormatSeqpos
                        " = expected number of integers",
                         (unsigned long) (numofbytes/sizeof (Seqpos)),
                         PRINTSeqposcast(totallength + 1));
      env_fa_xmunmap((void *) suffixarray->suftab,env);
      haserr = true;
    }
    str_delete(tmpfilename,env);
  }
  return haserr ? -1 : 0;
}

void freesuffixarray(Suffixarray *suffixarray,Env *env)
{
  env_fa_xmunmap((void *) suffixarray->suftab,env);
  freeAlphabet(&suffixarray->alpha,env);
  freeEncodedsequence(&suffixarray->encseq,env);
  strarray_delete(suffixarray->filenametab,env);
  FREESPACE(suffixarray->filelengthtab);
}
