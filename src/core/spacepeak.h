/*
  Copyright (c) 2011 Sascha Steinbiss <steinbiss@zbh.uni-hamburg.de>
  Copyright (c) 2011 Center for Bioinformatics, University of Hamburg

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

#ifndef SPACEPEAK_H
#define SPACEPEAK_H

/* The spacepeak module allows access to a combined space peak across multiple
   allocator objects. */

void          gt_spacepeak_init(void);
void          gt_spacepeak_add(unsigned long size);
void          gt_spacepeak_free(unsigned long size);
unsigned long gt_spacepeak_get_space_peak(void);
void          gt_spacepeak_show_space_peak(FILE *outfp);
void          gt_spacepeak_clean(void);

#endif
