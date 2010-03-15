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
#include <string.h>
#include <inttypes.h>
#include "core/error.h"
#include "core/str.h"
#include "core/array.h"
#include "core/symboldef.h"
#include "core/logger.h"
#include "format64.h"
#include "sfx-ri-def.h"

typedef union
{
  uint32_t smallvalue;
  uint64_t bigvalue;
} Smallorbigint;

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
  return sizeof (Readintkeys);
}

void setreadintkeys(GtArray *riktab,
                    const char *keystring,
                    void *valueptr,
                    size_t sizeval,
                    bool *readflag)
{
  Readintkeys rikvalue;

  rikvalue.keystring = keystring;
  rikvalue.readflag = readflag;
  gt_assert(sizeval == 0 || sizeval == (size_t) 4 || sizeval == (size_t) 8);
  switch (sizeval)
  {
    case 0: rikvalue.smallvalueptr = NULL;
            rikvalue.bigvalueptr = NULL;
            rikvalue.ptrdefined = false;
            break;
    case 4: gt_assert(sizeof (uint32_t) == (size_t) 4);
            rikvalue.smallvalueptr = valueptr;
            rikvalue.bigvalueptr = NULL;
            rikvalue.ptrdefined = true;
            break;
    case 8: gt_assert(sizeof (uint64_t) == (size_t) 8);
            rikvalue.bigvalueptr = valueptr;
            rikvalue.smallvalueptr = NULL;
            rikvalue.ptrdefined = true;
            break;
  }
  rikvalue.found = false;
  gt_array_add_elem(riktab,&rikvalue,sizeof (Readintkeys));
}

static int scanuintintline(uint32_t *lengthofkey,
                           Smallorbigint *smallorbigint,
                           const char *linebuffer,
                           unsigned long linelength,
                           GtError *err)
{
  int64_t readint;
  unsigned long i;
  bool found = false;
  int retval = 0;

  gt_error_check(err);
  for (i=0; i<linelength; i++)
  {
    if (linebuffer[i] == '=')
    {
      *lengthofkey = (uint32_t) i;
      found = true;

      if (sscanf((const char *) (linebuffer + i + 1),
                 FormatScanint64_t,
                 SCANint64_tcast(&readint)) != 1 ||
         readint < (int64_t) 0)
      {
        gt_error_set(err,"cannot find non-negative integer in \"%*.*s\"",
                           (int) (linelength - (i+1)),
                           (int) (linelength - (i+1)),
                           linebuffer + i + 1);
        return -1;
      }
      if (readint <= (int64_t) UINT32_MAX)
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
    gt_error_set(err,"missing equality symbol in \"%*.*s\"",
                       (int) linelength,
                       (int) linelength,
                       linebuffer);
    return -2;
  }
  return retval;
}

int allkeysdefined(const GtStr *indexname,const char *suffix,
                   const GtArray *riktab,GtLogger *logger,
                   GtError *err)
{
  unsigned long i;
  Readintkeys *rikptr;

  gt_error_check(err);
  for (i=0; i<gt_array_size(riktab); i++)
  {
    rikptr = (Readintkeys *) gt_array_get(riktab,i);
    if (rikptr->found)
    {
      if (rikptr->ptrdefined)
      {
        if (rikptr->smallvalueptr != NULL)
        {
          gt_logger_log(logger,"%s=%u",
                      rikptr->keystring,
                      (unsigned int) *(rikptr->smallvalueptr));
        } else
        {
          if (rikptr->bigvalueptr != NULL)
          {
            gt_logger_log(logger,"%s=" Formatuint64_t,
                        rikptr->keystring,
                        PRINTuint64_tcast(*(rikptr->bigvalueptr)));
          } else
          {
            gt_assert(false);
          }
        }
      } else
      {
        gt_logger_log(logger,"%s=0",rikptr->keystring);
      }
      if (rikptr->readflag != NULL)
      {
        *(rikptr->readflag) = true;
      }
    } else
    {
      if (rikptr->readflag == NULL)
      {
        gt_error_set(err,"file %s%s: missing line beginning with \"%s=\"",
                           gt_str_get(indexname),
                           suffix,
                           rikptr->keystring);
        return -1;
      }
      *(rikptr->readflag) = false;
    }
  }
  return 0;
}

int analyzeuintline(const GtStr *indexname,
                    const char *suffix,
                    unsigned int linenum,
                    const char *linebuffer,
                    unsigned long linelength,
                    GtArray *riktab,
                    GtError *err)
{
  Readintkeys *rikptr;
  bool found = false, haserr = false;
  unsigned long i;
  int retval;
  Smallorbigint smallorbigint;
  uint32_t lengthofkey;

  gt_error_check(err);
  retval = scanuintintline(&lengthofkey,
                           &smallorbigint,
                           linebuffer,
                           linelength,
                           err);
  if (retval < 0)
  {
    haserr = true;
  } else
  {
    for (i=0; i<gt_array_size(riktab); i++)
    {
      rikptr = gt_array_get(riktab,i);
      if (memcmp(linebuffer,
                 rikptr->keystring,
                 (size_t) lengthofkey) == 0)
      {
        rikptr->found = true;
        if (rikptr->ptrdefined)
        {
          if (rikptr->smallvalueptr == NULL)
          {
            if (retval == 1)
            {
              *(rikptr->bigvalueptr) = smallorbigint.bigvalue;
            } else
            {
              *(rikptr->bigvalueptr) = (uint64_t) smallorbigint.smallvalue;
            }
          } else
          {
            if (retval == 1)
            {
              gt_error_set(err,
                           "bigvalue " Formatuint64_t " does not fit into %s",
                           PRINTuint64_tcast(smallorbigint.bigvalue),
                           rikptr->keystring);
              haserr = true;
              break;
            }
            *(rikptr->smallvalueptr) = smallorbigint.smallvalue;
          }
        }
        found = true;
        break;
      }
    }
    if (!found)
    {
      gt_error_set(err,"file %s%s, line %u: cannot find key for \"%*.*s\"",
                    gt_str_get(indexname),
                    suffix,
                    linenum,
                    (int) lengthofkey,
                    (int) lengthofkey,
                    linebuffer);
      haserr = true;
    }
  }
  return haserr ? -1 : 0;
}
