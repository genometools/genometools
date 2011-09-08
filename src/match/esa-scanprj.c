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
#include "core/ma.h"
#include "core/array.h"
#include "core/types_api.h"
#include "core/logger.h"
#include "core/format64.h"
#include "esa-scanprj.h"

typedef union
{
  uint32_t uint32value;
  uint64_t uint64value;
  double doublevalue;
} GtScannedprjvalue;

typedef struct
{
  const char *keystring;
  uint32_t *uint32valueptr;
  uint64_t *uint64valueptr;
  double *doublevalueptr;
  bool readdouble,
       ptrdefined,
       found,
       *readflag;
} GtScannedprjkey;

struct GtScannedprjkeytable
{
  GtArray *arr;
};

GtScannedprjkeytable *gt_scannedprjkeytable_new(void)
{
  GtScannedprjkeytable *scannedprjkeytable;

  scannedprjkeytable = gt_malloc(sizeof (*scannedprjkeytable));
  scannedprjkeytable->arr = gt_array_new(sizeof (GtScannedprjkey));
  return scannedprjkeytable;
}

void gt_scannedprjkeytable_delete(GtScannedprjkeytable *scannedprjkeytable)
{
  gt_array_delete(scannedprjkeytable->arr);
  gt_free(scannedprjkeytable);
}

void gt_scannedprjkey_add(GtScannedprjkeytable *scannedprjkeytable,
                          const char *keystring,
                          void *valueptr,
                          size_t sizeval,
                          bool readdouble,
                          bool *readflag)
{
  GtScannedprjkey rikvalue;

  rikvalue.keystring = keystring;
  rikvalue.readflag = readflag;
  rikvalue.readdouble = readdouble;
  if (readdouble)
  {
    gt_assert(sizeval == sizeof (double));
    rikvalue.uint32valueptr = NULL;
    rikvalue.uint64valueptr = NULL;
    rikvalue.doublevalueptr = valueptr;
    rikvalue.ptrdefined = true;
  } else
  {
    gt_assert(sizeval == 0 || sizeval == (size_t) 4 || sizeval == (size_t) 8);
    switch (sizeval)
    {
      case 0: rikvalue.uint32valueptr = NULL;
              rikvalue.uint64valueptr = NULL;
              rikvalue.ptrdefined = false;
              break;
      case 4: gt_assert(sizeof (uint32_t) == (size_t) 4);
              rikvalue.uint32valueptr = valueptr;
              rikvalue.uint64valueptr = NULL;
              rikvalue.ptrdefined = true;
              break;
      case 8: gt_assert(sizeof (uint64_t) == (size_t) 8);
              rikvalue.uint64valueptr = valueptr;
              rikvalue.uint32valueptr = NULL;
              rikvalue.ptrdefined = true;
              break;
    }
  }
  rikvalue.found = false;
  gt_array_add_elem(scannedprjkeytable->arr,&rikvalue,
                    sizeof (GtScannedprjkey));
}

static int gt_scannedprjkey_scanline(uint32_t *lengthofkey,
                                     GtScannedprjvalue *scannedprjvalue,
                                     const char *linebuffer,
                                     unsigned long linelength,
                                     GtError *err)
{
  unsigned long idx;
  bool haserr = false, found = false;
  int retval = 0;

  gt_error_check(err);
  for (idx=0; idx<linelength; idx++)
  {
    if (linebuffer[idx] == '=')
    {
      *lengthofkey = (uint32_t) idx;
      found = true;
      if (strchr(linebuffer + idx + 1,'.') != NULL)
      {
        double readdouble;

        if (sscanf((const char *) (linebuffer + idx + 1),"%lf",
                   &readdouble) != 1)
        {
          gt_error_set(err,"cannot find floating point number in \"%*.*s\"",
                             (int) (linelength - (idx+1)),
                             (int) (linelength - (idx+1)),
                             linebuffer + idx + 1);
          return -1;
        }
        scannedprjvalue->doublevalue = readdouble;
        retval = 2;
      } else
      {
        int64_t readint;

        if (sscanf((const char *) (linebuffer + idx + 1),
                   FormatScanint64_t,
                   SCANint64_tcast(&readint)) != 1 || readint < (int64_t) 0)
        {
          gt_error_set(err,"cannot find non-negative integer in \"%*.*s\"",
                             (int) (linelength - (idx+1)),
                             (int) (linelength - (idx+1)),
                             linebuffer + idx + 1);
          return -1;
        }
        if (readint <= (int64_t) UINT32_MAX)
        {
          scannedprjvalue->uint32value = (uint32_t) readint;
          retval = 0;
        } else
        {
          scannedprjvalue->uint64value = (uint64_t) readint;
          retval = 1;
        }
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
    haserr = true;
  }
  return haserr ? -1 : retval;
}

int gt_scannedprjkey_allkeysdefined(
                               const char *indexname,
                               const char *suffix,
                               const GtScannedprjkeytable *scannedprjkeytable,
                               GtLogger *logger,
                               GtError *err)
{
  unsigned long idx;
  GtScannedprjkey *rikptr;

  gt_error_check(err);
  for (idx=0; idx<gt_array_size(scannedprjkeytable->arr); idx++)
  {
    rikptr = (GtScannedprjkey *) gt_array_get(scannedprjkeytable->arr,idx);
    if (rikptr->found)
    {
      if (rikptr->ptrdefined)
      {
        if (rikptr->uint32valueptr != NULL)
        {
          gt_logger_log(logger,"%s=%u",
                      rikptr->keystring,
                      (unsigned int) *(rikptr->uint32valueptr));
        } else
        {
          if (rikptr->uint64valueptr != NULL)
          {
            gt_logger_log(logger,"%s=" Formatuint64_t,
                        rikptr->keystring,
                        PRINTuint64_tcast(*(rikptr->uint64valueptr)));
          } else
          {
            if (rikptr->doublevalueptr != NULL)
            {
              gt_logger_log(logger,"%s=%.2f", rikptr->keystring,
                            *(rikptr->doublevalueptr));
            } else
            {
              gt_assert(false);
            }
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
                     indexname, suffix, rikptr->keystring);
        return -1;
      }
      *(rikptr->readflag) = false;
    }
  }
  return 0;
}

int gt_scannedprjkey_analyze(const char *indexname,
                             const char *suffix,
                             unsigned int linenum,
                             const char *linebuffer,
                             unsigned long linelength,
                             GtScannedprjkeytable *scannedprjkeytable,
                             GtError *err)
{
  GtScannedprjkey *rikptr;
  bool found = false, haserr = false;
  unsigned long i;
  int retval;
  GtScannedprjvalue scannedprjvalue;
  uint32_t lengthofkey;

  gt_error_check(err);
  retval = gt_scannedprjkey_scanline(&lengthofkey,
                                     &scannedprjvalue,
                                     linebuffer,
                                     linelength,
                                     err);
  if (retval < 0)
  {
    haserr = true;
  } else
  {
    for (i=0; i<gt_array_size(scannedprjkeytable->arr); i++)
    {
      rikptr = gt_array_get(scannedprjkeytable->arr,i);
      if (memcmp(linebuffer, rikptr->keystring,(size_t) lengthofkey) == 0)
      {
        rikptr->found = true;
        if (rikptr->ptrdefined)
        {
          if (rikptr->uint32valueptr == NULL)
          {
            if (retval == 1)
            {
              *(rikptr->uint64valueptr) = scannedprjvalue.uint64value;
            } else
            {
              if (retval == 2)
              {
                gt_assert(rikptr->readdouble);
                *(rikptr->doublevalueptr) = scannedprjvalue.doublevalue;
              } else
              {
                gt_assert(retval == 0);
                *(rikptr->uint64valueptr)
                  = (uint64_t) scannedprjvalue.uint32value;
              }
            }
          } else
          {
            if (retval == 1)
            {
              gt_error_set(err,"uint64value " Formatuint64_t
                           " does not fit into %s",
                           PRINTuint64_tcast(scannedprjvalue.uint64value),
                           rikptr->keystring);
              haserr = true;
              break;
            }
            if (retval == 2)
            {
              gt_error_set(err,"double %.2f does not fit into %s",
                           scannedprjvalue.doublevalue,
                           rikptr->keystring);
              haserr = true;
              break;
            }
            *(rikptr->uint32valueptr) = scannedprjvalue.uint32value;
          }
        }
        found = true;
        break;
      }
    }
    if (!found)
    {
      gt_error_set(err,"file %s%s, line %u: cannot find key for \"%*.*s\"",
                    indexname,
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
