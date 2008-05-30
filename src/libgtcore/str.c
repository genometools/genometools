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

#include <math.h>
#include <assert.h>
#include "libgtcore/cstr.h"
#include "libgtcore/dynalloc.h"
#include "libgtcore/ensure.h"
#include "libgtcore/ma.h"
#include "libgtcore/str.h"
#include "libgtcore/xansi.h"

struct Str {
  char *cstr;           /* the actual string (always '\0' terminated) */
  unsigned long length; /* currently used length (without trailing '\0') */
  size_t allocated;     /* currently allocated memory */
  unsigned int reference_count;
};

Str* str_new(void)
{
  Str *s = ma_malloc(sizeof (Str));      /* create new string object */
  s->cstr = ma_calloc(1, sizeof (char)); /* init string with '\0' */
  s->length = 0;                         /* set the initial length */
  s->allocated = 1;                      /* set allocated space */
  s->reference_count = 0;                /* set reference count */
  return s;                              /* return new string object */
}

Str* str_new_cstr(const char *cstr)
{
  Str *s = str_new();
  if (cstr)
    str_append_cstr(s, cstr);
  return s;
}

void str_set(Str *s, const char *cstr)
{
  size_t cstrlen;
  char *sptr;
  assert(s);
  if (!cstr)
    s->length = 0;
  else {
    cstrlen = strlen(cstr);
    s->cstr = dynalloc(s->cstr, &s->allocated, (cstrlen + 1) * sizeof (char));
    sptr = s->cstr;
    while (*cstr != '\0') *sptr++ = *cstr++;
    s->length = cstrlen;
  }
}

void str_append_str(Str *dest, const Str* src)
{
  unsigned long i;
  assert(dest && src);
  dest->cstr = dynalloc(dest->cstr, &dest->allocated,
                        (dest->length + src->length + 1) * sizeof (char));
  for (i = 0; i < src->length; i++)
    dest->cstr[dest->length + i] = src->cstr[i];
  dest->length += src->length;
}

void str_append_cstr(Str *dest, const char *cstr)
{
  size_t cstrlen;
  char *destptr;
  assert(dest && cstr);
  cstrlen = strlen(cstr);
  dest->cstr = dynalloc(dest->cstr, &dest->allocated,
                        (dest->length + cstrlen + 1) * sizeof (char));
  destptr = dest->cstr + dest->length;
  while (*cstr != '\0')
    *destptr++ = *cstr++;
  dest->length += cstrlen;
}

void str_append_cstr_nt(Str *dest, const char *cstr, unsigned long length)
{
  unsigned long i;
  char *destptr;
  assert(dest && cstr);
  dest->cstr = dynalloc(dest->cstr, &dest->allocated,
                        (dest->length + length + 1) * sizeof (char));
  destptr = dest->cstr + dest->length;
  for (i = 0; i < length; i++)
    *destptr++ = *cstr++;
  dest->length += length;
}

/* inspired by D. J. Bernstein's fmt_ulong() */
void str_append_ulong(Str *dest, unsigned long u)
{
  unsigned int ulength = 1;
  unsigned long q = u;
  char *s;
  assert(dest);
  /* determine length of u */
  while (q > 9) {
    ulength++;
    q /= 10;
  }
  /* make sure the string is long enough */
  dest->cstr = dynalloc(dest->cstr, &dest->allocated,
                        (dest->length + ulength + 1) * sizeof (char));
 /* format */
  s = dest->cstr + dest->length + ulength;
  do {
    *--s = '0' + (u % 10);
    u /= 10;
  }
  while (u); /* handles u == 0 */
  dest->length += ulength;
}

void str_append_char(Str *dest, char c)
{
  assert(dest);
  if (dest->length + 2 > dest->allocated) {
    dest->cstr = dynalloc(dest->cstr, &dest->allocated,
                          (dest->length + 2) * sizeof (char));
  }
  dest->cstr[dest->length++] = c;
}

void str_append_double(Str *dest, double d, int precision)
{
  char buf[BUFSIZ];
  int rval;
  assert(dest);
  rval = snprintf(buf, BUFSIZ, "%.*f", precision, d);
  assert(rval < BUFSIZ);
  str_append_cstr(dest, buf);
}

char* str_get(const Str *s)
{
  assert(s);
  s->cstr[s->length] = '\0';
  return s->cstr;
}

void* str_get_mem(const Str *s)
{
  assert(s);
  return s->cstr;
}

unsigned long str_length(const Str *s)
{
  return s ? s->length : 0;
}

void str_set_length(Str *s, unsigned long length)
{
  assert(s && length <= s->length);
  s->length = length;
}

void str_reset(Str *s)
{
  assert(s);
  s->length = 0;
}

/* does not handle embedded \0's */
int str_cmp(const Str *s1, const Str *s2)
{
  assert(s1 && s2);
  if (s1 == s2)
    return 0; /* a string is equal to itself */
  assert(s1->cstr && s2->cstr);
  s1->cstr[s1->length] = '\0';
  s2->cstr[s2->length] = '\0';
  return strcmp(s1->cstr, s2->cstr);
}

Str* str_clone(const Str *s)
{
  Str *s_copy;
  assert(s);
  s_copy = ma_malloc(sizeof (Str));
  s->cstr[s->length] = '\0';
  s_copy->cstr = cstr_dup(s->cstr);
  s_copy->length = s_copy->allocated = s->length;
  s_copy->reference_count = 0;
  return s_copy;
}

Str* str_ref(Str *s)
{
  if (!s) return NULL;
  s->reference_count++; /* increase the reference counter */
  return s;
}

int str_read_next_line(Str *s, FILE *fpin)
{
  int cc;
  char c;
  assert(s && fpin);
  for (;;) {
    cc = xfgetc(fpin);
    if (cc == EOF)
      return EOF;
    if (cc == '\n') {
      if ((s->length+1) * sizeof (char) > s->allocated)
        s->cstr = dynalloc(s->cstr, &s->allocated, (s->length+1)*sizeof (char));
      s->cstr[s->length] = '\0';
      return 0;
    }
    c = cc;
    if ((s->length+2) * sizeof (char) > s->allocated)
      s->cstr = dynalloc(s->cstr, &s->allocated, (s->length+2)*sizeof (char));
    s->cstr[s->length++] = c;
  }
}

int str_read_next_line_generic(Str *s, GenFile *fpin)
{
  int cc;
  char c;
  assert(s);
  for (;;) {
    cc = genfile_xfgetc(fpin);
    if (cc == EOF)
      return EOF;
    if (cc == '\n') {
      if ((s->length+1) * sizeof (char) > s->allocated)
        s->cstr = dynalloc(s->cstr, &s->allocated, (s->length+1)*sizeof (char));
      s->cstr[s->length] = '\0';
      return 0;
    }
    c = cc;
    if ((s->length+2) * sizeof (char) > s->allocated)
      s->cstr = dynalloc(s->cstr, &s->allocated, (s->length+2)*sizeof (char));
    s->cstr[s->length++] = c;
  }
}

int str_unit_test(Error *err)
{
  Str *s, *s1, *s2;
  static char cstring_1[] = "test_string"; /* l=11 */
  int had_err = 0;
  error_check(err);

  /* the empty string */
  s1 = str_new();
  ensure(had_err, str_length(s1) == 0);
  str_delete(s1);

  /* string testing */
  s1 = str_new();
  str_set(s1, cstring_1);
  ensure(had_err, str_length(s1) == 11);
  ensure(had_err, strcmp(str_get(s1), cstring_1) == 0);
  str_delete(s1);

  s1 = str_new_cstr(cstring_1);
  ensure(had_err, str_length(s1) == 11);
  str_delete(s1);

  s1 = str_new();
  s2 = str_new_cstr("foo");
  ensure(had_err, str_length(s2) == 3);
  str_append_str(s1, s2);
  ensure(had_err, str_length(s1) == 3);
  str_append_cstr(s1, "bar");
  ensure(had_err, str_length(s1) == 6);
  str_append_char(s1, 'b');
  str_append_char(s1, 'a');
  str_append_char(s1, 'z');
  ensure(had_err, str_length(s1) == 9);
  ensure(had_err, strcmp("foobarbaz", str_get(s1)) == 0);
  ensure(had_err, strcmp("foo", str_get(s2)) == 0);
  str_append_ulong(s1, 1984);
  ensure(had_err, strcmp("foobarbaz1984", str_get(s1)) == 0);
  str_delete(s1);
  str_delete(s2);

  /* testing str_append_ulong() and str_set_length() */
  s = str_new();
  str_append_ulong(s, 0);
  ensure(had_err, strcmp("0", str_get(s)) == 0);
  str_reset(s);
  ensure(had_err, strcmp("", str_get(s)) == 0);
  str_append_ulong(s, 6);
  ensure(had_err, strcmp("6", str_get(s)) == 0);
  str_append_ulong(s, 16);
  ensure(had_err, strcmp("616", str_get(s)) == 0);
  str_delete(s);

  /* make sure that str_get never returns NULL */
  s = str_new();
  ensure(had_err, str_get(s));
  ensure(had_err, str_length(s) == 0);
  ensure(had_err, strlen(str_get(s)) == 0);
  str_delete(s);

  s = str_new_cstr(NULL);
  ensure(had_err, str_get(s));
  ensure(had_err, str_length(s) == 0);
  ensure(had_err, strlen(str_get(s)) == 0);
  str_delete(s);

  /* test str_new() followed by str_append_cstr_nt() */
  s = str_new();
  str_append_cstr_nt(s, "foo", 3);
  ensure(had_err, strcmp("foo", str_get(s)) == 0);
  ensure(had_err, str_length(s) == 3);
  str_delete(s);

  return had_err;
}

void str_delete(Str *s)
{
  if (!s) return;           /* return without action if 's' is NULL */
  if (s->reference_count) { /* there are multiple references to this string */
    s->reference_count--;   /* decrement the reference counter */
    return;                 /* return without freeing the object */
  }
  ma_free(s->cstr);         /* free the stored the C string */
  ma_free(s);               /* free the actual string object */
}
