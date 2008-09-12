/*
  Copyright (c) 2006-2008 Gordon Gremme <gremme@zbh.uni-hamburg.de>
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

#ifndef STR_API_H
#define STR_API_H

#include <stdio.h>
#include "core/error_api.h"

/* Objects of the <GT_Str> class are strings which grow on demand. */
typedef struct GT_Str GT_Str;

/* Return an empty <GT_Str*> object. */
GT_Str*       gt_str_new(void);
/* Return a new <GT_Str*> object whose content is set to <cstr>. */
GT_Str*       gt_str_new_cstr(const char *cstr);
/* Return a clone of <str>. */
GT_Str*       gt_str_clone(const GT_Str *str);
/* Increase the reference count for <str> and return it.
   If <str> is <NULL>, <NULL> is returned without any side effects. */
GT_Str*       gt_str_ref(GT_Str *str);
/* Return the content of <str>.  Never returns NULL, and the content is always
   <\0>-terminated */
char*         gt_str_get(const GT_Str *str);
/* Set the content of <str> to <cstr>. */
void          gt_str_set(GT_Str *str, const char *cstr);
/* Append the string <src> to <dest>. */
void          gt_str_append_str(GT_Str *dest, const GT_Str *src);
/* Append the <\0>-terminated <cstr> to <str>. */
void          gt_str_append_cstr(GT_Str *str, const char *cstr);
/* Append the non <\0>-terminated <cstr> with given <length> to <str>. */
void          gt_str_append_cstr_nt(GT_Str *str,
                                    const char *cstr, unsigned long length);
/* Append character <c> to <str>. */
void          gt_str_append_char(GT_Str *str, char c);
/* Append double <d> to <str> with given <precision>. */
void          gt_str_append_double(GT_Str*, double d, int precision);
/* Append <ulong> to <str>. */
void          gt_str_append_ulong(GT_Str*, unsigned long ulong);
/* Set length of <str> to <length>. <length> must be smaller or equal than
   <gt_str_length(str)>. */
void          gt_str_set_length(GT_Str*, unsigned long length);
/* Reset <str> to length 0. */
void          gt_str_reset(GT_Str *str);
/* Compare <str1> and <str2> and return the result (similar to <strcmp(3)>). */
int           gt_str_cmp(const GT_Str *str1, const GT_Str *str2);
/* Return the length of <str>. If <str> is <NULL>, 0 is returned. */
unsigned long gt_str_length(const GT_Str *str);
/* Decrease the reference count for <str> or delete it, if this was the last
   reference. */
void          gt_str_delete(GT_Str *str);

#endif
