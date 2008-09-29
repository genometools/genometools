/*
  Copyright (C) 2007 Thomas Jahns <Thomas.Jahns@gmx.net>

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
#include <stdlib.h>

#include <time.h>
#include <sys/time.h>

#include "core/assert_api.h"
#include "core/bitpackstring.h"
#include "core/error.h"
#include "core/ensure.h"

int
gt_bitPackString_unit_test(GtError *err)
{
  return gt_bitPackStringInt_unit_test(err)
    || gt_bitPackStringInt8_unit_test(err)
    || gt_bitPackStringInt16_unit_test(err)
    || gt_bitPackStringInt32_unit_test(err)
    || gt_bitPackStringInt64_unit_test(err);
}
