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

#ifndef CHARDEF_H
#define CHARDEF_H

#include <limits.h>

/*
  This file defines some character values used when storing
  multiple sequences.
*/

/*
  separator symbol in multiple seq
*/

#define SEPARATOR       UCHAR_MAX

/*
  wildcard symbol in multiple seq
*/

#define WILDCARD        (SEPARATOR-1)

/*
  undefined character, only to be used in conjunction with symbol maps
*/

#define UNDEFCHAR       (SEPARATOR-2)

/*
  either WILDCARD or SEPARATOR
*/

#define ISSPECIAL(C)    ((C) >= (Uchar) WILDCARD)

/*
  neither WILDCARD nor SEPARATOR
*/

#define ISNOTSPECIAL(C) ((C) < (Uchar) WILDCARD)

/*
  undefined character, only to be used in conjunction with the Burrows-Wheeler
  transform
*/

#define UNDEFBWTCHAR    WILDCARD

/*
  Either special character or UNDEFBWTCHAR
*/

#define ISBWTSPECIAL(C) ((C) >= (Uchar) UNDEFBWTCHAR)

/*
  the size of the DNA alphabet
*/

#define DNAALPHASIZE       ((unsigned int) 4)

#endif
