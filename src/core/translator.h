/*
  Copyright (c) 2006-2010 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c)      2009 Sascha Steinbiss <steinbiss@zbh.uni-hamburg.de>
  Copyright (c) 2006-2009 Center for Bioinformatics, University of Hamburg

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

#ifndef TRANSLATOR_H
#define TRANSLATOR_H

#include "core/translator_api.h"
#include "core/error_api.h"

/* Returns the translation of the next codon. The currently translated
   character is put in <translated> while the current reading frame is put in
   <frame>. <start> is set to <true> if the current codon is a start codon and
   to <false> otherwise.
   Returns GT_TRANSLATOR_ERROR if an error occurred, see <err> for details.
   If the end of the sequence region to translate has been reached,
   GT_TRANSLATOR_END is returned.
   Otherwise, GT_TRANSLATOR_OK (equal to 0) is returned. */
GtTranslatorStatus gt_translator_next_with_start(GtTranslator *translator,
                                                 char *translated,
                                                 unsigned int *frame,
                                                 bool *start, GtError *err);

int gt_translator_unit_test(GtError *err);

#endif
