/*
  Copyright (c) 2007 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
  Copyright (c) 2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#include <stdio.h>
#include <inttypes.h>
#include "libgtcore/env.h"
#include "libgtcore/str.h"
#include "libgtcore/array.h"
#include "types.h"
#include "sfx-ri-def.h"

struct Readintkeys
{
  const char *keystring;
  uint32_t *smallvalueptr;
  uint64_t *bigvalueptr;
  bool ptrdefined, 
       found,
       *readflag;
};

size_t sizeofReadintkeys(void)
{
  return sizeof(Readintkeys);
}

const char *keystringreadintkeys(const Readintkeys *rikptr)
{
  return rikptr->keystring;
}

void setfoundreadintkeys(Readintkeys *rikptr)
{
  rikptr->found = true;
}

bool ptrdefinedreadintkeys(const Readintkeys *rikptr)
{
  return rikptr->ptrdefined;
}

bool smallvalueptrisnullreadintkeys(const Readintkeys *rikptr)
{
   return (rikptr->smallvalueptr == NULL) ? true : false;
}

void setbigvalueptrreadintkeys(Readintkeys *rikptr,uint64_t bigvalue)
{
  *(rikptr->bigvalueptr) = bigvalue;
}

void setsmallvalueptrreadintkeys(Readintkeys *rikptr,uint32_t smallvalue)
{
  *(rikptr->smallvalueptr) = smallvalue;
}

void setreadintkeys(Array *riktab,
                    const char *keystring,
                    void *valueptr,
                    size_t sizeval,
                    bool *readflag,
                    Env *env)
{
  Readintkeys rikvalue;

  rikvalue.keystring = keystring;
  rikvalue.readflag = readflag;
  switch(sizeval)
  {
    case 0: rikvalue.smallvalueptr = NULL;
            rikvalue.bigvalueptr = NULL;
            rikvalue.ptrdefined = false;
            break;
    case 4: assert(sizeof(uint32_t) == (size_t) 4);
            rikvalue.smallvalueptr = valueptr;
            rikvalue.bigvalueptr = NULL;
            rikvalue.ptrdefined = true;
            break;
    case 8: assert(sizeof(uint64_t) == (size_t) 8);
            rikvalue.bigvalueptr = valueptr;
            rikvalue.smallvalueptr = NULL;
            rikvalue.ptrdefined = true;
            break;
    default: fprintf(stderr,"sizeval must be 0 or 4 or 8\n");
             exit(EXIT_FAILURE);
  }
  rikvalue.found = false;
  array_add_elem(riktab,&rikvalue,sizeof(Readintkeys),env);
}

int scanuintintline(uint32_t *lengthofkey,
                    Smallorbigint *smallorbigint,
                    const Uchar *linebuffer,
                    unsigned long linelength,
                    Env *env)
{
  int64_t readint;
  unsigned long i;
  bool found = false;
  int retval = 0;

  env_error_check(env);
  for (i=0; i<linelength; i++)
  {
    if (linebuffer[i] == '=')
    {
      *lengthofkey = (uint32_t) i;
      found = true;

      if (sscanf((const char *) (linebuffer + i + 1),
                 FormatScanint64_t,
                 Scanuint64_tcast(&readint)) != 1 ||
         readint < (int64_t) 0)
      {
        env_error_set(env,"cannot find non-negative integer in \"%*.*s\"",
                           (Fieldwidthtype) (linelength - (i+1)),
                           (Fieldwidthtype) (linelength - (i+1)),
                           (const char *) (linebuffer + i + 1));
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
                       (Fieldwidthtype) linelength,
                       (Fieldwidthtype) linelength,
                       linebuffer);
    return -2;
  }
  return retval;
}

int allkeysdefined(const Str *indexname,const Array *riktab,Env *env)
{
  unsigned long i;
  Readintkeys *rikptr;

  env_error_check(env);
  for (i=0; i<array_size(riktab); i++)
  {
    rikptr = (Readintkeys *) array_get(riktab,i);
    if(rikptr->found)
    {
      printf("%s=",rikptr->keystring);
      if (rikptr->ptrdefined)
      {
        if (rikptr->smallvalueptr != NULL)
        {
          printf("%u\n",(unsigned int) *(rikptr->smallvalueptr));
        } else
        {
          if (rikptr->bigvalueptr != NULL)
          {
            printf(Formatuint64_t "\n",
                  PRINTuint64_tcast(*(rikptr->bigvalueptr)));
          } else
          {
            assert(false);
          }
        }
      } else
      {
        printf("0\n");
      }
      if(rikptr->readflag != NULL)
      {
        *(rikptr->readflag) = true;
      }
    } else
    {
      if (rikptr->readflag == NULL)
      {
        env_error_set(env,"file %s%s: missing line beginning with \"%s=\"",
                           str_get(indexname),
                           ".prj",
                           rikptr->keystring);
        return -1;
      }
      *(rikptr->readflag) = false;
    }
  }
  return 0;
}
