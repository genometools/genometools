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

#ifndef OPTIONARGMODE_H
#define OPTIONARGMODE_H
#include "core/str_api.h"
#include "core/error_api.h"

typedef struct
{
  const char *name,
             *desc;
  unsigned int bitmask;
} Optionargmodedesc;

int gt_optionargaddbitmask(const Optionargmodedesc *modedesc,
                           size_t numberofentries,
                           unsigned int *mode,
                           const char *optname,
                           const char *optionargument,
                           GtError *err);

int gt_optionargsetsingle(const Optionargmodedesc *modedesc,
                          size_t numberofentries,
                          const char *optname,
                          const char *optionargument,
                          GtError *err);

GtStr *gt_getargmodekeywords(const Optionargmodedesc *modedesc,
                             size_t numberofentries,
                             const char *what);

void gt_getsetargmodekeywords(const Optionargmodedesc *modedesc,
                              size_t numberofentries,
                              unsigned int bitfield);

#endif
