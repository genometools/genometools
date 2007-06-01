#include <stdio.h>
#include <errno.h>
#include <stdbool.h>
#include "libgtcore/env.h"
#include "types.h"
#include "arraydef.h"
#include "sarr-def.h"
#include "encseq-def.h"

#include "readnextline.pr"
#include "compfilenm.pr"
#include "endianess.pr"
#include "alphabet.pr"
#include "dstrdup.pr"

#define DBFILEKEY "dbfile="

typedef struct
{
  const char *keystring;
  Uint *valueptr;
  bool found;
} Readintkeys;

DECLAREARRAYSTRUCT(Readintkeys);

static void setreadintkeys(ArrayReadintkeys *riktab,
                           const char *keystring,
                           Uint *valueptr,
                           Env *env)
{
  Readintkeys *riktabptr;

  GETNEXTFREEINARRAY(riktabptr,riktab,Readintkeys,1);
  riktabptr->keystring = keystring;
  riktabptr->valueptr = valueptr;
  riktabptr->found = false;
}

static int scanuintintline(unsigned int *lengthofkey,
                           Uint *value,
                           const ArrayUchar *linebuffer,
                           Env *env)
{
  Scaninteger readint;
  unsigned int i;
  bool found = false;

  env_error_check(env);
  for (i=0; i<(unsigned int) linebuffer->nextfreeUchar; i++)
  {
    if (linebuffer->spaceUchar[i] == '=')
    {
      *lengthofkey = i;
      found = true;
      if (sscanf((char *) (linebuffer->spaceUchar + i + 1),"%ld",
                &readint) != 1 ||
         readint < (Scaninteger) 0)
      {
        env_error_set(env,"cannot find non-negative integer in \"%*.*s\"",
                           (Fieldwidthtype) (linebuffer->nextfreeUchar - (i+1)),
                           (Fieldwidthtype) (linebuffer->nextfreeUchar - (i+1)),
                           linebuffer->spaceUchar + i + 1);
        return -1;
      }
      *value = (Uint) readint;
      break;
    }
  }
  if (!found)
  {
    env_error_set(env,"missing equality symbol in \"%*.*s\"",
                       linebuffer->nextfreeUchar,
                       linebuffer->nextfreeUchar,
                       linebuffer->spaceUchar);
    return -2;
  }
  return 0;
}

static int allkeysdefined(const char *prjfile,const ArrayReadintkeys *rik,
                          Env *env)
{
  Uint i;

  env_error_check(env);
  for (i=0; i<rik->nextfreeReadintkeys; i++)
  {
    if (!rik->spaceReadintkeys[i].found)
    {
      env_error_set(env,"file %s: missing line beginning with \"%s=\"",
                         prjfile,
                         rik->spaceReadintkeys[i].keystring);
      return -1;
    }
    if (rik->spaceReadintkeys[i].valueptr == NULL)
    {
      printf("%s=0\n",rik->spaceReadintkeys[i].keystring);
    } else
    {
      printf("%s=%lu\n",rik->spaceReadintkeys[i].keystring,
                        (Showuint) *(rik->spaceReadintkeys[i].valueptr));
    }
  }
  return 0;
}

static int scanprjfileviafileptr(Suffixarray *suffixarray,
                                 Uint *totallength,
                                 const char *prjfile,FILE *fpin,Env *env)
{
  ArrayUchar linebuffer;
  Uint uintvalue, numofsequences, integersize, littleendian;
  unsigned int linenum, lengthofkey, i, numofallocatedfiles = 0;
  size_t dbfilelen = strlen(DBFILEKEY);
  int retval;
  bool found, haserr = false;
  ArrayReadintkeys rik;

  env_error_check(env);
  INITARRAY(&rik,Readintkeys);
  setreadintkeys(&rik,"totallength",totallength,env);
  setreadintkeys(&rik,"specialcharacters",
                         &suffixarray->specialcharinfo.specialcharacters,env);
  setreadintkeys(&rik,"specialranges",
                         &suffixarray->specialcharinfo.specialranges,env);
  setreadintkeys(&rik,"lengthofspecialprefix",
                         &suffixarray->specialcharinfo.lengthofspecialprefix,
                         env);
  setreadintkeys(&rik,"lengthofspecialsuffix",
                         &suffixarray->specialcharinfo.lengthofspecialsuffix,
                         env);
  setreadintkeys(&rik,"numofsequences",
                         &numofsequences,env);
  setreadintkeys(&rik,"numofdbsequences",
                         &suffixarray->numofdbsequences,env);
  setreadintkeys(&rik,"numofquerysequences",
                         NULL,env);
  setreadintkeys(&rik,"prefixlength",
                         &suffixarray->prefixlength,env);
  setreadintkeys(&rik,"integersize",
                         &integersize,env);
  setreadintkeys(&rik,"littleendian",
                      &littleendian,env);
  assert(rik.spaceReadintkeys != NULL);
  INITARRAY(&linebuffer,Uchar);
  suffixarray->filenametab = NULL;
  suffixarray->filelengthtab = NULL;
  suffixarray->numoffiles = 0;
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
      Scaninteger readint1, readint2;

      if (suffixarray->numoffiles >= numofallocatedfiles)
      {
        numofallocatedfiles += 2;
        ALLOCASSIGNSPACE(suffixarray->filenametab,suffixarray->filenametab,
                         char *,numofallocatedfiles);
        ALLOCASSIGNSPACE(suffixarray->filelengthtab,suffixarray->filelengthtab,
                         PairUint,numofallocatedfiles);
      }
      assert(suffixarray->filenametab != NULL);
      assert(suffixarray->filelengthtab != NULL);
      ALLOCASSIGNSPACE(tmpfilename,NULL,char,linebuffer.nextfreeUchar);
      if (sscanf((const char *) linebuffer.spaceUchar,"dbfile=%s %ld %ld\n",
                   tmpfilename,
                   &readint1,
                   &readint2) != 3)
      {
        env_error_set(env,"cannot parse line %*.*s",
                          (int) linebuffer.nextfreeUchar,
                          (int) linebuffer.nextfreeUchar,
                          (const char *) linebuffer.spaceUchar);
        FREESPACE(tmpfilename);
        haserr = true;
        break;
      }
      if (readint1 < (Scaninteger) 1 || readint2 < (Scaninteger) 1)
      {
        env_error_set(env,"need positive integers in line %*.*s",
                          (int) linebuffer.nextfreeUchar,
                          (int) linebuffer.nextfreeUchar,
                          (const char *) linebuffer.spaceUchar);
        haserr = true;
      }
      if (!haserr)
      {
        suffixarray->filenametab[suffixarray->numoffiles] = tmpfilename;
        suffixarray->filelengthtab[suffixarray->numoffiles].uint0
          = (Uint) readint1;
        suffixarray->filelengthtab[suffixarray->numoffiles].uint1
          = (Uint) readint2;
        printf("%s%s %lu %lu\n",
                DBFILEKEY,
                suffixarray->filenametab[suffixarray->numoffiles],
                (Showuint)
                suffixarray->filelengthtab[suffixarray->numoffiles].uint0,
                (Showuint)
                suffixarray->filelengthtab[suffixarray->numoffiles].uint1);
        suffixarray->numoffiles++;
      }
    } else
    {
      retval = scanuintintline(&lengthofkey,
                               &uintvalue,
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
          if (rik.spaceReadintkeys[i].valueptr != NULL)
          {
            *(rik.spaceReadintkeys[i].valueptr) = uintvalue;
          }
          found = true;
          break;
        }
      }
      if (!found)
      {
        env_error_set(env,"file %s, line %u: cannot find key for \"%*.*s\"",
                           prjfile,
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
  if (!haserr && integersize != UintConst(32) && integersize != UintConst(64))
  {
    env_error_set(env,"%s contains illegal line defining the integer size",
                  prjfile);
    haserr = true;
  }
  if (!haserr && integersize != (Uint) (sizeof (Uint) * CHAR_BIT))
  {
    env_error_set(env,"index was generated for %lu-bit integers while "
                      "this program uses %lu-bit integers",
                      (Showuint) integersize,
                      (Showuint) (sizeof (Uint) * CHAR_BIT));
    haserr = true;
  }
  if (!haserr)
  {
    if (islittleendian())
    {
      if (littleendian != UintConst(1))
      {
        env_error_set(env,"computer has little endian byte order, while index "
                          "was build on computer with big endian byte order");
        haserr = true;
      }
    } else
    {
      if (littleendian == UintConst(1))
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
                   const char *indexname,
                   Env *env)
{
  char *filename;
  FILE *fp;
  bool haserr = false;
  Uint totallength;

  env_error_check(env);
  suffixarray->suftab = NULL;
  suffixarray->alpha = NULL;
  filename = COMPOSEFILENAME(indexname,"prj");
  fp = env_fa_fopen(env,filename,"rb");
  if (fp == NULL)
  {
    env_error_set(env,"cannot open file \"%s\": %s",filename,
                                                    strerror(errno));
    haserr = true;
  }
  if (!haserr && scanprjfileviafileptr(suffixarray,&totallength,
                                      filename,fp,env) != 0)
  {
    haserr = true;
  }
  env_fa_xfclose(fp,env);
  FREESPACE(filename);
  if (!haserr)
  {
    filename = COMPOSEFILENAME(indexname,"al1");
    suffixarray->alpha = assigninputalphabet(false,
                                             false,
                                             filename,
                                             NULL,
                                             env);
    if (suffixarray->alpha == NULL)
    {
      haserr = true;
    }
    FREESPACE(filename);
  }
  if (!haserr && withencseq)
  {
    suffixarray->encseq = initencodedseq(true,
					 NULL,
					 0,
					 indexname,
					 (Uint64) totallength,
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

    filename = COMPOSEFILENAME(indexname,"suf");
    suffixarray->suftab = env_fa_mmap_read(env,filename,&numofbytes);
    if (suffixarray->suftab == NULL)
    {
      env_error_set(env,"cannot map file \"%s\": %s",filename,strerror(errno));
      haserr = true;
    }
    if (!haserr &&
       (Uint) (numofbytes/sizeof (Uint)) != totallength + 1)
    {
      env_error_set(env,"number of mapped integers = %lu != %lu = "
                        "expected number of integers",
                         (Showuint) (numofbytes/sizeof (Uint)),
                         (Showuint) (totallength + 1));
      env_fa_xmunmap((void *) suffixarray->suftab,env);
      haserr = true;
    }
    FREESPACE(filename);
  }
  return haserr ? -1 : 0;
}

void freesuffixarray(Suffixarray *suffixarray,Env *env)
{
  unsigned int i;

  env_fa_xmunmap((void *) suffixarray->suftab,env);
  freeAlphabet(&suffixarray->alpha,env);
  for (i=0; i<suffixarray->numoffiles; i++)
  {
    FREESPACE(suffixarray->filenametab[i]);
    suffixarray->filenametab[i] = NULL;
  }
  freeEncodedsequence(&suffixarray->encseq,env);
  FREESPACE(suffixarray->filenametab);
  FREESPACE(suffixarray->filelengthtab);
}
