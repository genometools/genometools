/*
  Copyright (c) 2008 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2008 Center for Bioinformatics, University of Hamburg

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

#include <stdlib.h>
#include "gt_config.h"
#include "core/version_api.h"

const char* gt_version_check(unsigned int required_major,
                             unsigned int required_minor,
                             unsigned int required_micro)
{
  if (required_major > GT_MAJOR_VERSION)
    return "GenomeTools library version too old (major mismatch)";
  if (required_major < GT_MAJOR_VERSION)
    return "GenomeTools library version too new (major mismatch)";
  if (required_minor > GT_MINOR_VERSION)
    return "GenomeTools library version too old (minor mismatch)";
  if (required_micro > GT_MICRO_VERSION)
    return "GenomeTools library version too old (micro mismatch)";
  return NULL;
}
