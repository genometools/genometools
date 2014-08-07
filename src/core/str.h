/*
  Copyright (c) 2006-2008 Gordon Gremme <gordon@gremme.org>
  Copyright (c) 2006-2008 Center for Bioinformatics, University of Hamburg

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

#ifndef STR_H
#define STR_H

#include <stdio.h>
#include "core/file.h"
#include "core/str_api.h"

/* never returns NULL, not always '\0' terminated */
void*         gt_str_get_mem(const GtStr*);
/* Remove end of <s> beginning with first occurrence of <c> */
void          gt_str_clip_suffix(GtStr *s, char c);
/* Read the next line from file pointer <fpin> and store the result in <str>
   (without the terminal newline). If the end of file <fpin> is reached, <EOF>
   is returned, otherwise 0. */
int           gt_str_read_next_line(GtStr *str, FILE *fpin);
int           gt_str_read_next_line_generic(GtStr*, GtFile*);
int           gt_str_unit_test(GtError*);

#endif
