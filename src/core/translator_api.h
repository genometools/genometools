/*
  Copyright (c) 2006-2008 Gordon Gremme <gremme@zbh.uni-hamburg.de>
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

#ifndef TRANSLATOR_API_H
#define TRANSLATOR_API_H

#include "core/str.h"
#include "core/str_array_api.h"

/* Implements a translator class supporting all translation schemes from
   __http://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi__.
   The <GtTranslator> can be used to generate single-frame translations of
   a DNA string to a <GtStr> or to produce 3-frame translations via an interator
   interface. */
typedef struct GtTranslator GtTranslator;

/* Creates a new <GtTranslator> with the standard genetic code translation
   scheme preselected. */
GtTranslator* gt_translator_new();

/* Selects the translation scheme in <tr> to the one identified by <transnum>.
   <transnum> refers to the numbers as reported by
   <gt_translator_get_translation_table_descriptions()> or the list given
   at the NCBI web site.
   Returns a negative value if an error occurred, see <e> for details. */
int           gt_translator_set_translation_scheme(GtTranslator *tr,
                                                   unsigned int transnum,
                                                   GtError *e);

/* Returns a <GtStrArray> of translation scheme descriptions, each of the
   format "%d: %s" where the number is the translation scheme number (usable in
   <gt_translator_set_translation_scheme()> and the string is the scheme
   name. */
GtStrArray*   gt_translator_get_translation_table_descriptions();

/* Starts a translation iteration process at the string position pointed to by
   <dnaseq>, intended to go on for <dnalen> bases. The current translation
   scheme set in <tr> will be used. The currently translated character is put in
   <translated> while the current reading frame (in reference to <dnaseq>
   position) is put in <frame>.
   Returns a negative value if an error occurred, see <e> for details. */
int           gt_translator_start(GtTranslator *tr,
                                  const char *dnaseq,
                                  unsigned long dnalen,
                                  char *translated,
                                  unsigned int *frame,
                                  GtError *e);

/* Continues a translation iteration process started by gt_translator_start().
   The currently translated character is put in <translated> while the current
   reading frame is put in <frame>.
   Returns a negative value if an error occurred, see <e> for details. */
int           gt_translator_next(GtTranslator *tr,
                                 char *translated,
                                 unsigned int *frame,
                                 GtError *e);

/* Translates <dnaseq> of length <dnalen> in reading frame <frame> using the
   settings currently active in <tr>. The resulting amino acid sequence is
   appended to <protein>.
   Returns a negative value if an error occurred, see <e> for details. */
int           gt_translator_translate_string(GtTranslator *tr,
                                             GtStr *protein,
                                             const char *dnaseq,
                                             unsigned long dnalen,
                                             unsigned int frame,
                                             GtError *e);

/* Determines the offset of the beginning of the first codon in <dnaseq> (of
   length <dnalen>) which is a start codon according to the current translation
   scheme in <tr>. The offset is written to the location pointed to by <pos>.
   Returns a negative value if an error occurred, see <e> for details. */
int           gt_translator_find_startcodon(GtTranslator *tr,
                                            const char *dnaseq,
                                            unsigned long dnalen,
                                            unsigned long *pos,
                                            GtError *e);

/* Writes the translation for the codon <c1>,<c2>,<c3> to the position pointed
   to by <amino>. The current translation scheme set in <tr> is used.
   Returns a negative value if an error occurred, see <e> for details. */
int           gt_translator_codon2amino(GtTranslator *tr,
                                        char c1, char c2, char c3,
                                        char *amino,
                                        GtError *e);

/* Deletes <tr> and frees all associated memory. */
void          gt_translator_delete(GtTranslator *tr);

#endif
