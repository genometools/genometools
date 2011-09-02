/*
  Copyright (c) 2009 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
  Copyright (c) 2009 Center for Bioinformatics, University of Hamburg

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

#ifndef S_SPLINT_S
#include <ctype.h>
#endif
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "core/ma_api.h"
#include "core/error_api.h"
#include "prsqualint.h"

#define ERRORMSG\
        "argument \"%s\" of option %s must be positive number "\
        "possibly followed by character %c or %c; if the "\
        "number is followed by character %c, then it must be "\
        " <= 100"

#define ERRORLPARAM\
        if (err != NULL)\
        {\
          gt_error_set(err,ERRORMSG,lparam,option,BESTCHARACTER,\
                       PERCENTAWAYCHARACTER,\
                       PERCENTAWAYCHARACTER);\
        } else\
        {\
          fprintf(stderr,ERRORMSG,lparam,option,BESTCHARACTER,\
                         PERCENTAWAYCHARACTER,\
                         PERCENTAWAYCHARACTER);\
        }

Qualifiedinteger *gt_parsequalifiedinteger(const char *option,
                                        const char *lparam,
                                        GtError *err)
{
  long readint;
  size_t i;
  char *lparamcopy;
  bool haserr = false;
  Qualifiedinteger *qualint;

  lparamcopy = gt_malloc(sizeof (char) * (strlen(lparam)+1));
  qualint = gt_malloc(sizeof (*qualint));
  strcpy(lparamcopy,lparam);
  for (i=0; lparamcopy[i] != '\0'; i++)
  {
    if (!isdigit((int) lparamcopy[i]) &&
        lparamcopy[i] != BESTCHARACTER &&
        lparamcopy[i] != PERCENTAWAYCHARACTER)
    {
      ERRORLPARAM;
      haserr = true;
      break;
    }
  }
  if (!haserr && i == 0)
  {
    ERRORLPARAM;
    haserr = true;
  }
  if (!haserr)
  {
    if (lparamcopy[i-1] == BESTCHARACTER)
    {
      lparamcopy[i-1] = '\0';
      qualint->qualtag = Qualbestof;
    } else
    {
      if (lparamcopy[i-1] == PERCENTAWAYCHARACTER)
      {
        lparamcopy[i-1] = '\0';
        qualint->qualtag = Qualpercentaway;
      } else
      {
        qualint->qualtag = Qualabsolute;
      }
    }
    if (sscanf(lparamcopy,"%ld",&readint) != 1 || readint <= 0)
    {
      ERRORLPARAM;
      haserr = true;
    }
  }
  if (!haserr &&
      (qualint->qualtag == Qualpercentaway || qualint->qualtag == Qualbestof))
  {
    if (readint > 100L)
    {
      ERRORLPARAM;
      haserr = true;
    }
  }
  qualint->integervalue = (unsigned long) readint;
  gt_free(lparamcopy);
  if (haserr)
  {
    gt_free (qualint);
    return NULL;
  }
  return qualint;
}
