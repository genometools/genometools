/*
  Copyright (c) 2012 Dirk Willrodt <willrodt@zbh.uni-hamburg.de>
  Copyright (c) 2012 Center for Bioinformatics, University of Hamburg

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

#ifndef STRING_ITER_H
#define STRING_ITER_H

/* <GtStringIter> is an abstract iterator class for data-structures containing
   strings */
typedef struct GtStringIter GtStringIter;

typedef struct GtStringIterClass GtStringIterClass;
typedef struct GtStringIterMembers GtStringIterMembers;

typedef int  (*GtStringIterNextFunc)(GtStringIter*,
                                     const char**,
                                     GtError*);
typedef int  (*GtStringIterResetFunc)(GtStringIter*, GtError*);
typedef void (*GtStringIterDeleteFunc)(GtStringIter*);

/* Sets <string> to the next string, retains ownership, will be overwritten by
   next call. Returns negative (<0)  on error and sets err acordingly, returns 0
   if no more strings are available and >0 on success. */
int                gt_string_iter_next(GtStringIter *string_iter,
                                       const char **string,
                                       GtError *err);

/* resets the iterator, a call to the next functien will return the first string
   again */
int                gt_string_iter_reset(GtStringIter *string_iter,
                                        GtError *err);

void               gt_string_iter_delete(GtStringIter *string_iter);

GtStringIterClass* gt_string_iter_class_new(size_t size,
                                            GtStringIterNextFunc,
                                            GtStringIterResetFunc,
                                            GtStringIterDeleteFunc);

GtStringIter*      gt_string_iter_create(const GtStringIterClass*);

void*              gt_string_iter_cast(const GtStringIterClass*, GtStringIter*);
#endif
