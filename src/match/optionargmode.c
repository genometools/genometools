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
#include "core/option.h"
#include "optionargmode.h"

int optionaddbitmask(const Optionargmodedesc *modedesc,
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

GtStr *getargmodekeywords(const Optionargmodedesc *modedesc,
                          size_t numberofentries,
                          const char *final)
{
  GtStr *helpstring;
  size_t modecount;

  helpstring = gt_str_new_cstr("use combination of the following keywords:\n");
  for (modecount=0; modecount < numberofentries; modecount++)
  {
    gt_str_append_cstr(helpstring,modedesc[modecount].name);
    gt_str_append_cstr(helpstring,"\n");
  }
  gt_str_append_cstr(helpstring,"to specify ");
  gt_str_append_cstr(helpstring,final);
  return helpstring;
}
