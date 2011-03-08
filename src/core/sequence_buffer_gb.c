/*
  Copyright (c) 2009 Sascha Steinbiss <steinbiss@zbh.uni-hamburg.de>
  Copyright (c) 2009 Center for Bioinformatics, University of Hamburg

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

#include <ctype.h>
#include <string.h>
#include "core/cstr_api.h"
#include "core/error.h"
#include "core/minmax.h"
#include "core/sequence_buffer_gb.h"
#include "core/sequence_buffer_rep.h"
#include "core/sequence_buffer_inline.h"
#include "core/warning_api.h"

#define GB_LOCUS_STRING      "LOCUS"
#define GB_DEFINITION_STRING "DEFINITION"
#define GB_ORIGIN_STRING     "ORIGIN"
#define GB_ENTRY_TERMINATOR  "//"
#define NEWLINESYMBOL        '\n'

typedef enum {
  GB_OUT_OF_ENTRY,
  GB_AWAITING_DESCRIPTION,
  GB_IN_DESCRIPTION,
  GB_AWAITING_SEQUENCE,
  GB_IN_SEQUENCE
} GtGBParserState;

struct GtSequenceBufferGB {
  const GtSequenceBuffer parent_instance;
  GtStr *headerbuffer,
        *overflowbuffer,
        *keywordbuffer;
  bool eof_was_set,
       description_set,
       filenum_has_changed;
  GtGBParserState state;
};

#define gt_sequence_buffer_gb_cast(SB)\
        gt_sequence_buffer_cast(gt_sequence_buffer_gb_class(), SB)

static inline int eat_newline(GtSequenceBuffer *sb,
                              unsigned long *currentfileread)
{
  int currentchar;
  currentchar = inlinebuf_getchar(sb, sb->pvt->inputstream);
  if (currentchar == EOF)
    return EOF;
  (*currentfileread)++;
  gt_assert(currentchar == NEWLINESYMBOL);
  return 0;
}

static inline int eat_whitespace(GtSequenceBuffer *sb,
                                 unsigned long *currentfileread)
{
  int currentchar;
  do {
    currentchar = inlinebuf_getchar(sb, sb->pvt->inputstream);
    if (currentchar == EOF)
      return EOF;
    (*currentfileread)++;
  } while (currentchar != NEWLINESYMBOL && isspace(currentchar));
  inlinebuf_ungetchar(sb);
  (*currentfileread)--;
  return 0;
}

static inline int eat_digits(GtSequenceBuffer *sb,
                             unsigned long *currentfileread)
{
  int currentchar;
  unsigned long i = 0;
  do {
    currentchar = inlinebuf_getchar(sb, sb->pvt->inputstream);
    if (currentchar == EOF)
      return EOF;
    (*currentfileread)++;
    i++;
  } while (isdigit(currentchar));
  inlinebuf_ungetchar(sb);
  (*currentfileread)--;
  return i-1;
}

static inline int eat_line(GtSequenceBuffer *sb,
                           unsigned long *currentfileread)
{
  int currentchar;
  do {
    currentchar = inlinebuf_getchar(sb, sb->pvt->inputstream);
    if (currentchar == EOF)
      return EOF;
    (*currentfileread)++;
  } while (currentchar != NEWLINESYMBOL);
  inlinebuf_ungetchar(sb);
  (*currentfileread)--;
  return 0;
}

static inline int get_keyword(GtSequenceBuffer *sb, GtStr *kwbuf,
                              unsigned long *currentfileread)
{
  int currentchar;
  currentchar = inlinebuf_getchar(sb, sb->pvt->inputstream);
  if (currentchar == EOF)
      return EOF;
  (*currentfileread)++;
  while (currentchar != NEWLINESYMBOL && !isspace(currentchar)) {
    gt_str_append_char(kwbuf, currentchar);
    currentchar = inlinebuf_getchar(sb, sb->pvt->inputstream);
    if (currentchar == EOF)
      return EOF;
    (*currentfileread)++;
  }
  inlinebuf_ungetchar(sb);
  (*currentfileread)--;
  return 0;
}

static inline int get_description(GtSequenceBuffer *sb,
                                  unsigned long *currentfileread)
{
  int currentchar;
  currentchar = inlinebuf_getchar(sb, sb->pvt->inputstream);
  if (currentchar == EOF)
      return EOF;
  (*currentfileread)++;
  while (currentchar != NEWLINESYMBOL) {
    if (sb->pvt->descptr)
      gt_desc_buffer_append_char(sb->pvt->descptr, currentchar);
    currentchar = inlinebuf_getchar(sb, sb->pvt->inputstream);
    if (currentchar == EOF)
      return EOF;
    (*currentfileread)++;
  }
  inlinebuf_ungetchar(sb);
  (*currentfileread)--;
  return 0;
}

static inline int get_sequence(GtSequenceBuffer *sb,
                               GtStr *linebuf,
                               unsigned long *currentoutpos,
                               unsigned long *currentfileread,
                               unsigned long *currentfileadd,
                               GtError *err)
{
  int currentchar, ret = 0, num_chars;
  /* expect digits */
  if ((num_chars = eat_digits(sb, currentfileread)) == EOF)
    return EOF;
  if (num_chars == 0) {
    gt_error_set(err, "sequence offset numbers missing in line %lu of file %s",
                      (unsigned long) sb->pvt->linenum,
                      gt_str_array_get(sb->pvt->filenametab,
                                       (unsigned long) sb->pvt->filenum));
    return -2;
  }
  currentchar = inlinebuf_getchar(sb, sb->pvt->inputstream);
  if (currentchar == EOF)
    return EOF;
  (*currentfileread)++;
  /* expect one blank */
  if (currentchar != ' ') {
    gt_error_set(err, "blank expected between offset and sequence in line %lu "
                      "of file %s",
                      (unsigned long) sb->pvt->linenum,
                      gt_str_array_get(sb->pvt->filenametab,
                                       (unsigned long) sb->pvt->filenum));
    return -2;
  }
  /* read sequence */
  do {
    currentchar = inlinebuf_getchar(sb, sb->pvt->inputstream);
    if (currentchar == EOF)
      return EOF;
    (*currentfileread)++;
    if (!isspace(currentchar)) {
      if (*currentoutpos >= (unsigned long) OUTBUFSIZE) {
        /* if outbuffer is full, keep rest of sequence in overflow buffer */
        gt_str_append_char(linebuf, currentchar);
      } else {
         /* fill output buffer w/ current char */
         if ((ret = process_char(sb, *currentoutpos, currentchar, err)))
             return ret;
        (*currentoutpos)++;
        (*currentfileadd)++;
      }
    }
  } while (currentchar != NEWLINESYMBOL);
  inlinebuf_ungetchar(sb);
  (*currentfileread)--;
  return 0;
}

static int gt_sequence_buffer_gb_advance(GtSequenceBuffer *sb, GtError *err)
{
  unsigned long currentoutpos = 0, currentfileadd = 0, currentfileread = 0;
  GtSequenceBufferMembers *pvt;
  GtSequenceBufferGB *sbe;
  int had_err = 0;
  gt_error_check(err);

  sbe = gt_sequence_buffer_gb_cast(sb);
  pvt = sb->pvt;

  /* open first stream */
  if (!pvt->inputstream) {
    sbe->state = GB_OUT_OF_ENTRY;
    pvt->linenum = (uint64_t) 1;
    pvt->inputstream = gt_file_xopen(gt_str_array_get(pvt->filenametab,
                                                  (unsigned long) pvt->filenum),
                                        "rb");
    pvt->currentinpos = 0;
    pvt->currentfillpos = 0;
  }

  if (gt_str_length(sbe->overflowbuffer) > 0) {
    /* we still have surplus sequence from the last line, process that first */
    unsigned long i;
    const char *overflowedstring;
    overflowedstring = gt_str_get_mem(sbe->overflowbuffer);
    for (i=0;i<gt_str_length(sbe->overflowbuffer);i++) {
      process_char(sb, currentoutpos, overflowedstring[i], err);
      currentoutpos++;
      currentfileadd++;
    }
    gt_str_reset(sbe->overflowbuffer);
    gt_assert(gt_str_length(sbe->overflowbuffer) == 0);
  }

  while (true) {
    if (had_err && had_err != EOF)
      break;

    if (currentoutpos >= (unsigned long) OUTBUFSIZE) {
      if (pvt->filelengthtab) {
        pvt->filelengthtab[pvt->filenum].length
          += (uint64_t) currentfileread;
        pvt->filelengthtab[pvt->filenum].effectivelength
          += (uint64_t) currentfileadd;
      }
      break;
    }

    had_err = get_keyword(sb, sbe->keywordbuffer, &currentfileread);

    if (had_err && had_err != EOF)
      break;

    /* terminators may occur in any line */
    if (!had_err && strcmp(gt_str_get(sbe->keywordbuffer),
                           GB_ENTRY_TERMINATOR) == 0) {
      pvt->outbuf[currentoutpos++] = (GtUchar) SEPARATOR;
      currentfileadd++;
      pvt->lastspeciallength++;
      if (!sbe->description_set && pvt->descptr)
          gt_desc_buffer_finish(pvt->descptr);
      sbe->description_set = false;
      if ((eat_line(sb, &currentfileread)) == EOF)
        break;
      sbe->state = GB_OUT_OF_ENTRY;
    }

    switch (sbe->state) {
      case GB_OUT_OF_ENTRY:
        /* eat everything */
        had_err = eat_line(sb, &currentfileread);
        if (!had_err && strcmp(gt_str_get(sbe->keywordbuffer),
                               GB_LOCUS_STRING) == 0) {
          sbe->state = GB_AWAITING_DESCRIPTION;
        }
        break;
      case GB_AWAITING_DESCRIPTION:
        if (strcmp(gt_str_get(sbe->keywordbuffer),
                   GB_DEFINITION_STRING) == 0) {
          had_err = eat_whitespace(sb, &currentfileread);
          if (!had_err)
            had_err = get_description(sb, &currentfileread);
          sbe->state = GB_IN_DESCRIPTION;
        } else if (strcmp(gt_str_get(sbe->keywordbuffer),
                   GB_ORIGIN_STRING) == 0) {
          gt_warning("sequence started without prior DEFINITION line in entry "
                     "in line %lu of file %s",
                     (unsigned long) pvt->linenum-1,
                     gt_str_array_get(pvt->filenametab,
                                     (unsigned long) pvt->filenum));
          had_err = eat_line(sb, &currentfileread);
          sbe->state = GB_IN_SEQUENCE;
        } else {
          had_err = eat_line(sb, &currentfileread);
        }
        break;
      case GB_IN_DESCRIPTION:
        if (gt_str_length(sbe->keywordbuffer) == 0)
        {
          had_err = eat_whitespace(sb, &currentfileread);
          if (pvt->descptr)
            gt_desc_buffer_append_char(pvt->descptr, ' ');
          if (!had_err)
            had_err = get_description(sb, &currentfileread);
        } else {
          if (strcmp(gt_str_get(sbe->keywordbuffer),
                     GB_DEFINITION_STRING) == 0) {
            gt_error_set(err, "encountered another DEFINITION line within one "
                              "entry in line %lu of file %s",
                              (unsigned long) pvt->linenum-1,
                              gt_str_array_get(pvt->filenametab,
                                                (unsigned long) pvt->filenum));
            return -1;
          } else {
            if (pvt->descptr) {
              gt_desc_buffer_finish(pvt->descptr);
            }
            sbe->description_set = true;
            if (strcmp(gt_str_get(sbe->keywordbuffer),
                       GB_ORIGIN_STRING) == 0) {
              had_err = eat_line(sb, &currentfileread);
              sbe->state = GB_IN_SEQUENCE;
            } else {
              had_err = eat_line(sb, &currentfileread);
              sbe->state = GB_AWAITING_SEQUENCE;
            }
          }
        }
        break;
      case GB_AWAITING_SEQUENCE:
        if (strcmp(gt_str_get(sbe->keywordbuffer),
                   GB_ORIGIN_STRING) == 0) {
          had_err = eat_line(sb, &currentfileread);
          sbe->state = GB_IN_SEQUENCE;
        } else {
          had_err = eat_line(sb, &currentfileread);
        }
        break;
      case GB_IN_SEQUENCE:
        if (gt_str_length(sbe->keywordbuffer) != 0) {
          gt_error_set(err, "only terminators allowed after a "
                            "sequence section, but found '%s' instead "
                            "in line %lu of file %s",
                            gt_str_get(sbe->keywordbuffer),
                            (unsigned long) pvt->linenum-1,
                            gt_str_array_get(pvt->filenametab,
                                                (unsigned long) pvt->filenum));
          return -1;
        } else {
          had_err = eat_whitespace(sb, &currentfileread);
          if (!had_err)
            had_err = get_sequence(sb, sbe->overflowbuffer, &currentoutpos,
                                   &currentfileread, &currentfileadd, err);
        }
        break;
    }
    /* ensure that line has been read completely */
    if (!had_err)
      had_err = eat_newline(sb, &currentfileread);
    gt_str_reset(sbe->keywordbuffer);
    /* handle end of file */
    if (had_err == EOF) {
      if (pvt->filelengthtab) {
        pvt->filelengthtab[pvt->filenum].length
          += (uint64_t) currentfileread;
        pvt->filelengthtab[pvt->filenum].effectivelength
          += (uint64_t) currentfileadd;
      }
      if (pvt->filenum+1 < gt_str_array_size(pvt->filenametab)) {
        pvt->filenum++;
        /* still files left, open next one */
        gt_file_delete(pvt->inputstream);
        sbe->state = GB_OUT_OF_ENTRY;
        pvt->linenum = (uint64_t) 1;
        pvt->inputstream = gt_file_xopen(gt_str_array_get(pvt->filenametab,
                                                  (unsigned long) pvt->filenum),
                                            "rb");
        if (pvt->filelengthtab) {
          pvt->filelengthtab[pvt->filenum].length = 0;
          pvt->filelengthtab[pvt->filenum].effectivelength = 0;
        }
        pvt->currentinpos = 0;
        pvt->currentfillpos = 0;
        currentfileread = currentfileadd = 0;
      } else {
        /* all files exhausted */
        pvt->complete = true;
        /* remove last separator */
        pvt->outbuf[--currentoutpos] = (GtUchar) '\0';
        if (pvt->filelengthtab) {
          pvt->filelengthtab[pvt->filenum].effectivelength--;
        }
        had_err = 0;
        break;
      }
    }
    pvt->linenum++;
  }
  pvt->nextfree = MIN(currentoutpos, OUTBUFSIZE);
  return had_err;
}

static void gt_sequence_buffer_gb_free(GtSequenceBuffer *sb)
{
  GtSequenceBufferGB *sbe = gt_sequence_buffer_gb_cast(sb);
  if (sb->pvt->inputstream)
    gt_file_delete(sb->pvt->inputstream);
  gt_str_delete(sbe->headerbuffer);
  gt_str_delete(sbe->overflowbuffer);
  gt_str_delete(sbe->keywordbuffer);
}

static unsigned long
gt_sequence_buffer_gb_get_file_index(GtSequenceBuffer *sb)
{
  gt_assert(sb);
  return sb->pvt->filenum;
}

const GtSequenceBufferClass* gt_sequence_buffer_gb_class(void)
{
  static const GtSequenceBufferClass sbc = { sizeof (GtSequenceBufferGB),
                                         gt_sequence_buffer_gb_advance,
                                         gt_sequence_buffer_gb_get_file_index,
                                         gt_sequence_buffer_gb_free };
  return &sbc;
}

bool gt_sequence_buffer_gb_guess(const char* txt)
{
  char *hit = NULL;
  if (!(hit = strstr(txt, "LOCUS ")))
    return false;
  /* LOCUS keyword must be at beginning of line */
  return (hit == txt || *(hit-1) == '\n');
}

GtSequenceBuffer* gt_sequence_buffer_gb_new(const GtStrArray *sequences)
{
  GtSequenceBuffer *sb;
  GtSequenceBufferGB *sbe;
  sb = gt_sequence_buffer_create(gt_sequence_buffer_gb_class());
  sb->pvt->filenametab = sequences;
  sb->pvt->filenum = 0;
  sb->pvt->nextread = sb->pvt->nextfree = 0;
  sb->pvt->lastspeciallength = 0;
  sbe = gt_sequence_buffer_gb_cast(sb);
  sbe->headerbuffer = gt_str_new();
  sbe->overflowbuffer = gt_str_new();
  sbe->keywordbuffer = gt_str_new();
  return sb;
}
