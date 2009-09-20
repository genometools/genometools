/*
  Copyright (c) 2004-2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2004-2007 Center for Bioinformatics, University of Hamburg

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

#ifndef GTHALPHATYPE_H
#define GTHALPHATYPE_H

/* the following enum type defines the valid alphabet types */
typedef enum {
  DNA_ALPHA,
  PROTEIN_ALPHA,
  MIXED_ALPHA,   /* not a real alphabet type, but used to denote that reference
                    files with DNA_ALPHA alphabets and with PROTEIN_ALPHA
                    alphabets are used at the same time */
  UNDEF_ALPHA,
} GthAlphatype;

#endif
