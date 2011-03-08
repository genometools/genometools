/*
  Copyright (c) 2009-2010 Sascha Steinbiss <steinbiss@zbh.uni-hamburg.de>
  Copyright (c) 2009-2010 Center for Bioinformatics, University of Hamburg

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
#include "core/sequence_buffer_embl.h"
#include "core/sequence_buffer_rep.h"
#include "core/sequence_buffer_inline.h"

#define EMBL_ID_LINE_STRING     "ID"
#define EMBL_DESCR_LINE_STRING  "DE"
#define EMBL_SEQ_LINE_STRING    "  "
#define EMBL_SPACER             "XX"
#define EMBL_ENTRY_TERMINATOR   "//"
#define NEWLINESYMBOL           '\n'

typedef enum {
  EMBL_UNDEFINED,
  EMBL_IN_SEQUENCE,
  EMBL_IN_DESCRIPTION
} GtEMBLParserState;

typedef enum {
  DESCRIPTION,
  ID,
  SEQUENCE,
  SPACER,
  TERMINATOR,
  OTHER
} GtEMBLParserLineCode;

struct GtSequenceBufferEMBL {
  const GtSequenceBuffer parent_instance;
  GtStr *headerbuffer,
        *overflowbuffer;
  bool firstentryinfile,
       firstoverallentry,
       nextfile,
       eof_was_set,
       description_set;
  GtEMBLParserState state;
};

#define gt_sequence_buffer_embl_cast(SB)\
        gt_sequence_buffer_cast(gt_sequence_buffer_embl_class(), SB)

static inline int
parse_next_line(GtSequenceBuffer *sb, GtEMBLParserLineCode *lc,
                GtStr *linebuf, const char **content,
                unsigned long *linelen, unsigned long *currentfileread,
                unsigned long *currentfileadd, unsigned long *currentoutpos,
                GtError *err)
{
  int currentchar, i, ret = 0;
  char linecode[2];
  GtSequenceBufferMembers *pvt;
  *lc = OTHER;

  pvt = sb->pvt;
  pvt->linenum++;
  *linelen = 0;

  /* determine line code */
  currentchar = inlinebuf_getchar(sb, pvt->inputstream);
  if (currentchar == EOF)
    return EOF;
  (*currentfileread)++;
  if (currentchar == NEWLINESYMBOL)
    return 0;
  linecode[0] = currentchar;
  currentchar = inlinebuf_getchar(sb, pvt->inputstream);
  if (currentchar == EOF)
    return EOF;
  (*currentfileread)++;
  if (currentchar == NEWLINESYMBOL) {
    gt_error_set(err, "2-character line code not found in line %lu",
                 (unsigned long) pvt->linenum - 1);
    return -2;
  }
  linecode[1] = currentchar;

  /* determine current line type */
  if (memcmp(linecode,      EMBL_DESCR_LINE_STRING, 2*sizeof (char)) == 0)
    *lc = DESCRIPTION;
  else if (memcmp(linecode, EMBL_SEQ_LINE_STRING ,  2*sizeof (char)) == 0)
    *lc = SEQUENCE;
  else if (memcmp(linecode, EMBL_SPACER,            2*sizeof (char)) == 0)
    *lc = SPACER;
  else if (memcmp(linecode, EMBL_ENTRY_TERMINATOR,  2*sizeof (char)) == 0)
    *lc = TERMINATOR;

  /* expect 3 blanks, except in spacer lines */
  if (*lc != SPACER) {
    for (i = 0; i < 3; i++) {
      currentchar = inlinebuf_getchar(sb, pvt->inputstream);
      if (currentchar == EOF)
        return EOF;
      (*currentfileread)++;
      if (currentchar == NEWLINESYMBOL)
        return 0;
      if (!isspace(currentchar)) {
        gt_error_set(err, "3 blanks expected between line code and content "
                          "in line %lu",
                     (unsigned long) pvt->linenum - 1);
        return -2;
      }
    }
  }

  /* read line body */
  i = 0;
  currentchar = inlinebuf_getchar(sb, pvt->inputstream);
  if (currentchar == EOF)
    return EOF;
  (*currentfileread)++;
  while (currentchar != NEWLINESYMBOL) {
    switch (*lc) {
      case SEQUENCE:
      /* sequences are 60 characters per line with 5 spaces in between */
        if (i++ < 65 && !isspace(currentchar)) {
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
          (*linelen)++;
        }
        break;
      case DESCRIPTION:
        /* write description into buffer */
        if (pvt->descptr)
          gt_desc_buffer_append_char(pvt->descptr, currentchar);
        break;
      default:
        break;
    }
    currentchar = inlinebuf_getchar(sb, pvt->inputstream);
    if (currentchar == EOF)
      return EOF;
    (*currentfileread)++;
  }
  *content = gt_str_get(linebuf);
  return 0;
}

static int gt_sequence_buffer_embl_advance(GtSequenceBuffer *sb, GtError *err)
{
  unsigned long currentoutpos = 0, currentfileadd = 0, currentfileread = 0,
                linelen = 0;
  GtSequenceBufferMembers *pvt;
  GtSequenceBufferEMBL *sbe;
  GtEMBLParserLineCode lc = OTHER;
  const char *line = "";
  int had_err = 0;

  gt_error_check(err);

  sbe = gt_sequence_buffer_embl_cast(sb);
  pvt = sb->pvt;
  /* open first stream */
  if (!pvt->inputstream) {
    sbe->firstentryinfile = true;
    sbe->state = EMBL_UNDEFINED;
    pvt->linenum = (uint64_t) 1;
    pvt->inputstream = gt_file_xopen(gt_str_array_get(pvt->filenametab,
                                                  (unsigned long) pvt->filenum),
                                         "rb");
    pvt->currentinpos = 0;
    pvt->currentfillpos = 0;
  }

  if (gt_str_length(sbe->overflowbuffer) > 0) {
    /* we still have surplus sequence from the last line, process that first */
    unsigned long i, len;
    const char *overflowedstring;
    overflowedstring = gt_str_get_mem(sbe->overflowbuffer);
    len = gt_str_length(sbe->overflowbuffer);
    for (i=0;i<len;i++)
    {
      process_char(sb, currentoutpos, overflowedstring[i], err);
      currentoutpos++;
      currentfileadd++;
    }
    gt_str_reset(sbe->overflowbuffer);
    gt_assert(gt_str_length(sbe->overflowbuffer) == 0);
  }

  /* parse file line-by-line */
  while (true) {
    had_err = parse_next_line(sb, &lc, sbe->overflowbuffer, &line, &linelen,
                              &currentfileread, &currentfileadd, &currentoutpos,
                              err);

    if (had_err && had_err != EOF) {
      break;
    }
    if (currentoutpos >= (unsigned long) OUTBUFSIZE) {
      if (pvt->filelengthtab) {
        pvt->filelengthtab[pvt->filenum].length
          += (uint64_t) currentfileread;
        pvt->filelengthtab[pvt->filenum].effectivelength
          += (uint64_t) currentfileadd;
      }
      pvt->nextfree = MIN(currentoutpos, OUTBUFSIZE);
      return 0; /* buffer full, finished */
    }
    if (lc == TERMINATOR) {
      pvt->outbuf[currentoutpos++] = (GtUchar) SEPARATOR;
      currentfileadd++;
      pvt->lastspeciallength++;
      sbe->state = EMBL_UNDEFINED;
      if (!sbe->description_set && pvt->descptr)
          gt_desc_buffer_finish(pvt->descptr);
      sbe->description_set = false;
    }
    /* FSM transitions begin here */
    switch (sbe->state) {
      case EMBL_IN_DESCRIPTION:
        if (lc != DESCRIPTION) {
          /* save description */
          if (pvt->descptr) {
            gt_desc_buffer_finish(pvt->descptr);
          }
          sbe->description_set = true;
          /* gt_str_reset(sbe->headerbuffer); */
          sbe->state = EMBL_UNDEFINED;
        }
        break;
      case EMBL_IN_SEQUENCE:
        /* only a terminator may come after a sequence */
        if (lc != SEQUENCE) {
          gt_error_set(err, "unterminated sequence in line %lu of file %s",
                            (unsigned long) pvt->linenum,
                            gt_str_array_get(pvt->filenametab,
                                                (unsigned long) pvt->filenum));
          return -1;
        }
        break;
      case EMBL_UNDEFINED:
        switch (lc) {
          case DESCRIPTION:
            sbe->state = EMBL_IN_DESCRIPTION;
            if (pvt->descptr) {
              gt_desc_buffer_append_char(pvt->descptr, ' ');
            }
            break;
          case SEQUENCE:
            sbe->state = EMBL_IN_SEQUENCE;
            break;
          default:
            break;
        }
        break;
    }
    if (had_err == EOF) {
      /* save length table entries if needed */
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
        sbe->state = EMBL_UNDEFINED;
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
        gt_assert(pvt->outbuf[currentoutpos-1] == SEPARATOR);
        pvt->outbuf[--currentoutpos] = (GtUchar) '\0';
        if (pvt->filelengthtab) {
          pvt->filelengthtab[pvt->filenum].effectivelength--;
        }
        had_err = 0;
        break;
      }
    }
  }
  pvt->nextfree = MIN(currentoutpos, OUTBUFSIZE);
  return had_err;
}

static void gt_sequence_buffer_embl_free(GtSequenceBuffer *sb)
{
  GtSequenceBufferEMBL *sbe = gt_sequence_buffer_embl_cast(sb);
  if (sb->pvt->inputstream)
    gt_file_delete(sb->pvt->inputstream);
  gt_str_delete(sbe->headerbuffer);
  gt_str_delete(sbe->overflowbuffer);
}

static unsigned long
gt_sequence_buffer_embl_get_file_index(GtSequenceBuffer *sb)
{
  gt_assert(sb);
  return sb->pvt->filenum;
}

const GtSequenceBufferClass* gt_sequence_buffer_embl_class(void)
{
  static const GtSequenceBufferClass sbc = { sizeof (GtSequenceBufferEMBL),
                                         gt_sequence_buffer_embl_advance,
                                         gt_sequence_buffer_embl_get_file_index,
                                         gt_sequence_buffer_embl_free };
  return &sbc;
}

bool gt_sequence_buffer_embl_guess(const char* txt)
{
  return !strncmp(txt, EMBL_ID_LINE_STRING, strlen(EMBL_ID_LINE_STRING));
}

GtSequenceBuffer* gt_sequence_buffer_embl_new(const GtStrArray *sequences)
{
  GtSequenceBuffer *sb;
  GtSequenceBufferEMBL *sbe;
  sb = gt_sequence_buffer_create(gt_sequence_buffer_embl_class());
  sbe = gt_sequence_buffer_embl_cast(sb);
  sb->pvt->filenametab = sequences;
  sbe->headerbuffer = gt_str_new();
  sbe->overflowbuffer = gt_str_new();
  sb->pvt->filenum = 0;
  sbe->firstoverallentry = true;
  sbe->firstentryinfile = true;
  sbe->nextfile = true;
  sb->pvt->nextread = sb->pvt->nextfree = 0;
  sb->pvt->lastspeciallength = 0;
  return sb;
}
