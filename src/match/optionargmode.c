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

#include <string.h>
#include "optionargmode.h"
#include "core/unused_api.h"

int gt_optionargaddbitmask(const Optionargmodedesc *modedesc,
                           size_t numberofentries,
                           unsigned int *mode,
                           const char *optname,
                           const char *optionargument,
                           GtError *err)
{
  size_t modecount;

  gt_error_check(err);
  for (modecount=0; modecount < numberofentries; modecount++)
  {
    if (strcmp(optionargument,modedesc[modecount].name) == 0)
    {
      if (*mode & modedesc[modecount].bitmask)
      {
        gt_error_set(err,"argument \"%s\" to option %s already specified",
                      modedesc[modecount].name,optname);
        return -1;
      }
      *mode |= modedesc[modecount].bitmask;
      return 0;
    }
  }
  gt_error_set(err,"illegal argument \"%s\" to option %s",
                optionargument,optname);
  return -2;
}

int gt_optionargsetsingle(const Optionargmodedesc *modedesc,
                          size_t numberofentries,
                          const char *optname,
                          const char *optionargument,
                          GtError *err)
{
  size_t modecount;

  gt_error_check(err);
  for (modecount=0; modecount < numberofentries; modecount++)
  {
    if (strcmp(optionargument,modedesc[modecount].name) == 0)
    {
      return (int) modedesc[modecount].bitmask;
    }
  }
  gt_error_set(err,"illegal argument \"%s\" to option %s",
                optionargument,optname);
  return -1;
}

GtStr *gt_getargmodekeywords(const Optionargmodedesc *modedesc,
                             size_t numberofentries,
                             const char *what)
{
  GtStr *helpstring;
  size_t idx, modecount, len, maxlen = 0;
  GT_UNUSED size_t spacelen;
  const char *space = "    ";

  for (modecount=0; modecount < numberofentries; modecount++)
  {
    len = strlen(modedesc[modecount].name);
    if (maxlen < len)
    {
      maxlen = len;
    }
  }
  helpstring = gt_str_new_cstr("use combination of the following keywords "
                               "to specify ");
  gt_str_append_cstr(helpstring,what);
  gt_str_append_cstr(helpstring,"\n");
  spacelen = strlen(space);
  for (modecount=0; modecount < numberofentries; modecount++)
  {
    gt_str_append_cstr(helpstring,modedesc[modecount].name);
    gt_str_append_cstr(helpstring,space);
    len = strlen(modedesc[modecount].name);
    for (idx=0; idx<maxlen+-len; idx++)
    {
      gt_str_append_cstr(helpstring," ");
    }
    gt_str_append_cstr(helpstring,"show ");
    gt_str_append_cstr(helpstring,modedesc[modecount].desc);
    gt_str_append_cstr(helpstring,"\n");
  }
  return helpstring;
}

void gt_getsetargmodekeywords(const Optionargmodedesc *modedesc,
                              size_t numberofentries,
                              unsigned int bitfield)
{
  size_t modecount;

  for (modecount=0; modecount < numberofentries; modecount++)
  {
    if (bitfield & modedesc[modecount].bitmask)
    {
      printf("%s ",modedesc[modecount].name);
    }
  }
  printf("\n");
}
