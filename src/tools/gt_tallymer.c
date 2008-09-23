/*
  Copyright (c) 2007 Thomas Jahns <Thomas.Jahns@gmx.net>
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

#include "core/cstr_array.h"
#include "core/error.h"
#include "core/option.h"
#include "core/str.h"
#include "core/unused_api.h"
#include "extended/toolbox.h"
#include "core/versionfunc.h"

static int gt_tallymer_mkindex(GT_UNUSED int argc, 
                               GT_UNUSED const char *argv[], 
                               GT_UNUSED GtError *err)
{
  return 0;
}

void *gt_tallymer_arguments_new(void)
{
  GtToolbox *tallymer_toolbox = gt_toolbox_new();
  gt_toolbox_add(tallymer_toolbox, "mkindex", gt_tallymer_mkindex);
  return tallymer_toolbox;
}
