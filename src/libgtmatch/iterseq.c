/*
  Copyright (c) 2007 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
  Copyright (c) 2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#include "libgtcore/strarray.h"
#include "symboldef.h"
#include "arraydef.h"
#include "chardef.h"
#include "fbs-def.h"

#include "fbsadv.pr"

#include "readnextUchar.gen"

int overallquerysequences(int(*processsequence)(void *,
                                                unsigned long,
                                                const Uchar *,
                                                const char *,
                                                unsigned long),
                          void *info,
                          ArrayUchar *sequencebuffer,
                          const StrArray *filenametab,
                          const Uchar *symbolmap,
                          Env *env)
{
  Fastabufferstate fbs;
  Uchar charcode;
  int retval;
  unsigned long unitnum = 0;

  initformatbufferstate(&fbs,
                        filenametab,
                        symbolmap,
                        false,
                        NULL,
                        env);
  sequencebuffer->nextfreeUchar = 0;
  while(true)
  {
    retval = readnextUchar(&charcode,&fbs,env);
    if (retval < 0)
    {
      return -1;
    }
    if (retval == 0)
    {
      break;
    }
    if (ISSPECIAL(charcode))
    {
      if (charcode == (Uchar) SEPARATOR)
      {
        if(sequencebuffer->nextfreeUchar > 0)
        {
          if(processsequence(info,unitnum,sequencebuffer->spaceUchar,
                             NULL,sequencebuffer->nextfreeUchar) != 0)
          {
            return -2;
          }
          sequencebuffer->nextfreeUchar = 0;
        }
        unitnum++;
      }
    } else
    {
      STOREINARRAY(sequencebuffer,Uchar,1024,charcode);
    }
  }
  if(sequencebuffer->nextfreeUchar > 0)
  {
    if(processsequence(info,unitnum,sequencebuffer->spaceUchar,
                       NULL,sequencebuffer->nextfreeUchar) != 0)
    {
      return -3;
    }
    sequencebuffer->nextfreeUchar = 0;
  }
  return 0;
}
