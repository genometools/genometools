/*
  Copyright (c) 2008-2010 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2008      Center for Bioinformatics, University of Hamburg

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

#ifndef SYMBOL_H
#define SYMBOL_H

void        gt_symbol_init(void);

/* Return a symbol (a canonical representation) for <cstr>.  An advantage of
   symbols is that they can be compared for equality by a simple pointer
   comparison, rather than using <strcmp()> (as it is done in <gt_strcmp()>).
   Furthermore, a symbol is stored only once in memory for equal <cstr>s, but
   keep in mind that this memory can never be freed safely during the lifetime
   of the calling program. */
const char* gt_symbol(const char *cstr);

/* Free (and thereby invalidate) all created symbols! */
void        gt_symbol_clean(void);

int         gt_symbol_unit_test(GtError*);

#endif
