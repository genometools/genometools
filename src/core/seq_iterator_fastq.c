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

#include "core/chardef.h"
#include "core/class_alloc_lock.h"
#include "core/colorspace.h"
#include "core/cstr_api.h"
#include "core/file.h"
#include "core/filelengthvalues.h"
#include "core/seq_iterator_fastq_api.h"
#include "core/seq_iterator_rep.h"
#include "core/str_array.h"
#include "core/unused_api.h"

#define GT_SEQIT_QUAL_INBUFSIZE  8192

struct GtSeqIteratorFastQ
{
  const GtSeqIterator parent_instance;
  unsigned int filenum;
  uint64_t linenum;
  GtFilelengthvalues *filelengthtab;
  bool complete,
       use_ungetchar,
       is_color_space,
       relax_qualdesc_check;
  GtStr *sequencebuffer,
        *descbuffer,
        *qualsbuffer;
  GtStr *qdescbuffer;
  GtFile *curfile;
  unsigned long *chardisttab,
                currentfillpos,
                currentinpos,
                curline,
                reference_count;
  uint64_t lastspeciallength;
  unsigned long long maxread,
                     currentread;
  const GtStrArray *filenametab;
  unsigned char ungetchar,
                inbuf[GT_SEQIT_QUAL_INBUFSIZE];
  const GtUchar *symbolmap, **qualities;
};

#define gt_seq_iterator_fastq_cast(SI)\
        gt_seq_iterator_cast(gt_seq_iterator_fastq_class(), SI);

const GtSeqIteratorClass* gt_seq_iterator_fastq_class(void);

#define GT_FASTQ_BLOCK_START_CHAR      '@'
#define GT_FASTQ_QUAL_SEPARATOR_CHAR   '+'
#define GT_FASTQ_NEWLINESYMBOL         '\n'

static inline int fastq_buf_getchar(GtSeqIteratorFastQ *seqit)
{
  if (seqit->use_ungetchar) {
    seqit->use_ungetchar = false;
    return seqit->ungetchar;
  } else {
    if (seqit->currentinpos >= seqit->currentfillpos) {
      seqit->currentfillpos = gt_file_xread(seqit->curfile, seqit->inbuf,
                                             GT_SEQIT_QUAL_INBUFSIZE);
      if (seqit->currentfillpos == 0)
         return EOF;
      seqit->currentinpos = 0;
    }
    seqit->ungetchar = seqit->inbuf[seqit->currentinpos++];
    return seqit->ungetchar;
  }
}

static inline void fastq_buf_ungetchar(GtSeqIteratorFastQ *seqit)
{
  gt_assert(!seqit->use_ungetchar);
  seqit->use_ungetchar = true;
}

static inline int parse_fastq_seqname(GtSeqIteratorFastQ *seqit,
                                      GtStr *buffer,
                                      char startchar,
                                      GtError *err)
{
  char currentchar;
  bool firstsymbol = true;
  gt_error_check(err);
  gt_assert(seqit && buffer);
  gt_assert(gt_str_length(buffer) == 0);
  if ((currentchar = fastq_buf_getchar(seqit)) == EOF)
    return EOF;
  seqit->currentread++;
  if (currentchar != startchar) {
    gt_error_set(err, "'%c' expected, '%c' encountered instead in line %lu",
                      startchar,
                      currentchar,
                      seqit->curline);
    return -2;
  }
  while (currentchar != GT_FASTQ_NEWLINESYMBOL) {
    if (!firstsymbol)
      gt_str_append_char(buffer, currentchar);
    else
      firstsymbol = false;
    if ((currentchar = fastq_buf_getchar(seqit)) == EOF)
      return EOF;
    seqit->currentread++;
  }
  seqit->curline++;
  return 0;
}

static int parse_fastq_sequence(GtSeqIteratorFastQ *seqit,
                                GtError *err)
{
  int had_err = 0;
  char currentchar;
  GtStr *tmp_str = gt_str_new();

  gt_error_check(err);
  gt_assert(seqit);
  gt_assert(gt_str_length(seqit->sequencebuffer) == 0);
  /* read sequence */
  if ((currentchar = fastq_buf_getchar(seqit)) == EOF) {
    gt_str_delete(tmp_str);
    return EOF;
  }
  while (currentchar != GT_FASTQ_QUAL_SEPARATOR_CHAR) {
    if (currentchar != '\n' && currentchar != ' ') {
      gt_str_append_char(tmp_str, currentchar);
    } else if (currentchar == '\n') {
      seqit->curline++;
    }
    if ((currentchar = fastq_buf_getchar(seqit)) == EOF) {
      gt_str_delete(tmp_str);
      return EOF;
    }
    seqit->currentread++;
  }
  if (!gt_str_length(tmp_str)) {
    gt_error_set(err, "empty sequence given in file '%s', line %lu",
                      gt_str_array_get(seqit->filenametab,
                                       seqit->filenum),
                      seqit->curline-1);
    had_err = -2;
  }
  if (!had_err && seqit->is_color_space)
  {
    GtStr *translated = gt_str_new();
    had_err = gt_colorspace_decode_string(tmp_str,
                                          translated,
                                          err);
    gt_str_delete(tmp_str);
    tmp_str = translated;
  }
  if (!had_err)
  {
    if (seqit->symbolmap)
    {
      int charcode;
      char *input_str = gt_str_get(tmp_str);
      unsigned long str_len = gt_str_length(tmp_str),
                    idx;
      for (idx = 0; !had_err && idx < str_len; idx++)
      {
        charcode = seqit->symbolmap[(unsigned int) input_str[idx]];
        if (charcode == UNDEFCHAR) {
          gt_error_set(err, "illegal character '%c': file \"%s\", line %lu",
                            input_str[idx],
                            gt_str_array_get(seqit->filenametab,
                                             seqit->filenum),
                            (unsigned long) seqit->curline);
          had_err = -2;
        }
        if (ISSPECIAL(charcode)) {
          seqit->lastspeciallength++;
        } else {
          if (seqit->lastspeciallength > 0)
            seqit->lastspeciallength = 0;
          if (seqit->chardisttab)
            seqit->chardisttab[(int) charcode]++;
        }
        gt_str_append_char(seqit->sequencebuffer, charcode);
      }
    } else {
      gt_str_set(seqit->sequencebuffer, gt_str_get(tmp_str));
    }
  }
  fastq_buf_ungetchar(seqit);
  gt_str_delete(tmp_str);
  return had_err;
}

static inline int parse_fastq_qualities(GtSeqIteratorFastQ *seqit,
                                        GT_UNUSED GtError *err)
{
  char currentchar;
  unsigned long i = 0;
  gt_assert(gt_str_length(seqit->sequencebuffer) > 0);
  if ((currentchar = fastq_buf_getchar(seqit)) == EOF)
    return EOF;
  seqit->currentread++;

  for (i=0;i<gt_str_length(seqit->sequencebuffer);i++) {
    if (currentchar != '\n' && currentchar != ' ') {
      gt_str_append_char(seqit->qualsbuffer, currentchar);
    } else if (currentchar == '\n') {
      seqit->curline++;
      i--;
    } else {
      i--;
    }
    if (i+1 == gt_str_length(seqit->sequencebuffer)) {
      seqit->curline++;
    }
    if ((currentchar = fastq_buf_getchar(seqit)) == EOF)
      return EOF;
    seqit->currentread++;
  }
  /* expect newline at end of qualities */
  if (currentchar != GT_FASTQ_NEWLINESYMBOL) {
    gt_error_set(err, "qualities string of sequence length %lu is not ended "
                      "by newline in file '%s', line %lu -- this may be a "
                      "sign for sequence and qualities strings of different "
                      "length",
                      gt_str_length(seqit->sequencebuffer),
                      gt_str_array_get(seqit->filenametab,
                                       seqit->filenum),
                      seqit->curline-1);
    return -2;
  }
  return 0;
}

#define gt_fastq_premature_end_check(had_err, seqit) \
  if (had_err == EOF) { \
    gt_error_set(err, "premature end of file '%s' in line %lu: " \
                      "file ended before end of block", \
                      gt_str_array_get((seqit)->filenametab, \
                                       (seqit)->filenum), \
                      (unsigned long) (seqit)->curline-1); \
    return -2; \
  }

static inline int parse_fastq_block(GtSeqIteratorFastQ *seqit, GtError *err)
{
  int had_err = 0;
  gt_assert(seqit);
  gt_error_check(err);

  /* parse @<seqname> */
  had_err = parse_fastq_seqname(seqit,
                                seqit->descbuffer,
                                GT_FASTQ_BLOCK_START_CHAR,
                                err);
  if (!had_err) {
    /* parse sequence */
    had_err = parse_fastq_sequence(seqit, err);
    gt_fastq_premature_end_check(had_err, seqit);
  }
  if (!had_err) {
    /* parse +[seqname] */
    had_err = parse_fastq_seqname(seqit,
                                  seqit->qdescbuffer,
                                  GT_FASTQ_QUAL_SEPARATOR_CHAR,
                                  err);
    gt_fastq_premature_end_check(had_err, seqit);
  }
  if (!had_err
      && !seqit->relax_qualdesc_check
      && gt_str_length(seqit->qdescbuffer)
      && gt_str_cmp(seqit->descbuffer, seqit->qdescbuffer) != 0)
  {
      gt_error_set(err, "sequence description '%s' is not equal to "
                        "qualities description '%s' in line %lu",
                        gt_str_get(seqit->descbuffer),
                        gt_str_get(seqit->qdescbuffer),
                        seqit->curline-1);
      return -2;
  }
  if (!had_err) {
    /* parse qualities */
    had_err = parse_fastq_qualities(seqit, err);
    if (gt_str_length(seqit->qualsbuffer)
          != gt_str_length(seqit->sequencebuffer))
    {
      gt_error_set(err, "lengths of character sequence and qualities "
                        "sequence differ (%lu <-> %lu)",
                        gt_str_length(seqit->qualsbuffer),
                        gt_str_length(seqit->sequencebuffer));
      return -2;
    }
  }
  return had_err;
}

unsigned long gt_seq_iterator_fastq_get_file_index(GtSeqIteratorFastQ *seqit)
{
  gt_assert(seqit);
  return seqit->filenum;
}

void gt_seq_iterator_fastq_set_quality_buffer(GtSeqIterator *si,
                                             const GtUchar **qualities)
{
  GtSeqIteratorFastQ *seqit;
  gt_assert(si);
  seqit = gt_seq_iterator_fastq_cast(si);
  seqit->qualities = qualities;
}

void gt_seq_iterator_fastq_set_symbolmap(GtSeqIterator *si,
                                        const GtUchar *symbolmap)
{
  GtSeqIteratorFastQ *seqit;
  gt_assert(si);
  seqit = gt_seq_iterator_fastq_cast(si);
  seqit->symbolmap = symbolmap;
}

void gt_seq_iterator_fastq_set_chardisttab(GtSeqIterator *si,
                                          unsigned long *chardist)
{
  GtSeqIteratorFastQ *seqit;
  gt_assert(si && chardist);
  seqit = gt_seq_iterator_fastq_cast(si);
  seqit->chardisttab = chardist;
}

uint64_t gt_seq_iterator_fastq_get_lastspeciallength(const GtSeqIterator *si)
{
  GtSeqIteratorFastQ *seqit;
  gt_assert(si);
  seqit = gt_seq_iterator_fastq_cast((GtSeqIterator*) si);
  return seqit->lastspeciallength;
}

int gt_seq_iterator_fastq_next(GtSeqIterator *seqit,
                              const GtUchar **sequence,
                              unsigned long *len,
                              char **desc,
                              GtError *err)
{
  int errstatus = 0;
  GtSeqIteratorFastQ *seqitf;
  gt_assert(seqit);
  seqitf = gt_seq_iterator_fastq_cast((GtSeqIterator*) seqit);
  gt_assert(seqit && len && desc);

  seqitf = gt_seq_iterator_fastq_cast(seqit);
  gt_str_reset(seqitf->qualsbuffer);
  gt_str_reset(seqitf->qdescbuffer);
  gt_str_reset(seqitf->sequencebuffer);
  gt_str_reset(seqitf->descbuffer);

  /* parse file */
  errstatus = parse_fastq_block(seqitf, err);

  if (!errstatus) {
    *sequence = (GtUchar*) gt_str_get(seqitf->sequencebuffer);
    *len = gt_str_length(seqitf->sequencebuffer);
    *desc = gt_str_get(seqitf->descbuffer);
    if (seqitf->qualities)
      *seqitf->qualities = (GtUchar*) gt_str_get(seqitf->qualsbuffer);
    errstatus = 1;
  } else {
    if (errstatus == EOF) {
      /* we could not get a next entry from this file */
      /* can we open another? */
      if (seqitf->filenum+1 < gt_str_array_size(seqitf->filenametab)) {
        const char *filename;
        filename = gt_str_array_get(seqitf->filenametab, ++seqitf->filenum);
        gt_file_delete(seqitf->curfile);
        seqitf->curfile = gt_file_xopen(filename, "r");
        seqitf->curline = 1;
        /* get first entry from next file*/
        errstatus = parse_fastq_block(seqitf, err);
        if (!errstatus) {
          *sequence = (GtUchar*) gt_str_get(seqitf->sequencebuffer);
          *len = gt_str_length(seqitf->sequencebuffer);
          *desc = gt_str_get(seqitf->descbuffer);
          if (seqitf->qualities)
            *seqitf->qualities = (GtUchar*) gt_str_get(seqitf->qualsbuffer);
          errstatus = 1;
        } else {
          errstatus = -1;
        }
      } else {
        /* all entries read from all files */
        errstatus = 0;
      }
    } else {
      errstatus = -1;
    }
  }
  return errstatus;
}

const unsigned long long*
gt_seq_iterator_fastq_getcurrentcounter(GtSeqIterator *si,
                                       unsigned long long maxread)
{
  GtSeqIteratorFastQ *seqit;
  gt_assert(si);
  seqit = gt_seq_iterator_fastq_cast(si);
  seqit->maxread = maxread;
  return &seqit->currentread;
}

void gt_seq_iterator_fastq_delete(GtSeqIterator *si)
{
  GtSeqIteratorFastQ *seqit;
  if (!si) return;
  seqit = gt_seq_iterator_fastq_cast(si);
  gt_str_delete(seqit->qdescbuffer);
  gt_str_delete(seqit->sequencebuffer);
  gt_str_delete(seqit->qualsbuffer);
  gt_str_delete(seqit->descbuffer);
  if (seqit->curfile)
    gt_file_delete(seqit->curfile);
  seqit->currentread = seqit->maxread;
}

const GtSeqIteratorClass* gt_seq_iterator_fastq_class(void)
{
  static const GtSeqIteratorClass *sic = NULL;
  gt_class_alloc_lock_enter();
  if (!sic) {
    sic = gt_seq_iterator_class_new(sizeof (GtSeqIteratorFastQ),
                             gt_seq_iterator_fastq_set_symbolmap,
                             NULL,
                             gt_seq_iterator_fastq_next,
                             gt_seq_iterator_fastq_getcurrentcounter,
                             gt_seq_iterator_fastq_set_quality_buffer,
                             gt_seq_iterator_fastq_delete);
  }
  gt_class_alloc_lock_leave();
  return sic;
}

static GtSeqIterator* seqiterator_fastq_new_gen(const GtStrArray *filenametab,
                                                bool is_color_space,
                                                GT_UNUSED GtError *err)
{
  GtSeqIterator *seqit;
  GtSeqIteratorFastQ *seqitf;
  gt_assert(filenametab);
  seqit = gt_seq_iterator_create(gt_seq_iterator_fastq_class());
  seqitf = gt_seq_iterator_fastq_cast(seqit);
  seqitf->qdescbuffer = gt_str_new();
  seqitf->curfile = gt_file_xopen(gt_str_array_get(filenametab, 0), "r");
  seqitf->filenametab = filenametab;
  seqitf->curline = 1;
  seqitf->sequencebuffer = gt_str_new();
  seqitf->qualsbuffer = gt_str_new();
  seqitf->descbuffer = gt_str_new();
  seqitf->is_color_space = is_color_space;
  seqitf->relax_qualdesc_check = false;
  return seqit;
}

void gt_seq_iterator_fastq_relax_check_of_quality_description(
    GtSeqIteratorFastQ *seqitf)
{
  gt_assert(seqitf != NULL);
  seqitf->relax_qualdesc_check = true;
}

GtSeqIterator* gt_seq_iterator_fastq_new(const GtStrArray *filenametab,
                                     GtError *err)
{
  return seqiterator_fastq_new_gen(filenametab,
                                   false,
                                   err);
}

GtSeqIterator* gt_seq_iterator_colorspace_fastq_new(
                                                const GtStrArray *filenametab,
                                                GtError *err)
{
  return seqiterator_fastq_new_gen(filenametab,
                                   true,
                                   err);
}
