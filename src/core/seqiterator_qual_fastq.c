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

#include "core/cstr.h"
#include "core/genfile.h"
#include "core/str_array.h"
#include "core/seqiterator_qual_fastq.h"
#include "core/seqiterator_qual_rep.h"
#include "core/unused_api.h"

#define GT_FASTQ_BLOCK_START_CHAR      '@'
#define GT_FASTQ_QUAL_SEPARATOR_CHAR   '+'
#define GT_FASTQ_NEWLINESYMBOL         '\n'

#define gt_seqiterator_qual_fastq_cast(SI)\
        gt_seqiterator_qual_cast(gt_seqiterator_qual_fastq_class(), SI)

struct GtSeqIteratorQualFASTQ {
  const GtSeqIteratorQual parent_instance;
  GtStr *qdescbuffer;
};

static inline int fastq_buf_getchar(GtSeqIteratorQual *seqit)
{
  GtSeqIteratorQualMembers *pvt = seqit->pvt;
  if (pvt->use_ungetchar) {
    pvt->use_ungetchar = false;
    return pvt->ungetchar;
  } else {
    if (pvt->currentinpos >= pvt->currentfillpos) {
      pvt->currentfillpos = gt_genfile_xread(pvt->curfile, pvt->inbuf,
                                             GT_SEQIT_QUAL_INBUFSIZE);
      if (pvt->currentfillpos == 0)
         return EOF;
      pvt->currentinpos = 0;
    }
    pvt->ungetchar = pvt->inbuf[pvt->currentinpos++];
    return pvt->ungetchar;
  }
}

static inline void fastq_buf_ungetchar(GtSeqIteratorQual *seqit)
{
  gt_assert(!seqit->pvt->use_ungetchar);
  seqit->pvt->use_ungetchar = true;
}

static inline int parse_fastq_seqname(GtSeqIteratorQual *seqit,
                                      GtStr *buffer,
                                      char startchar,
                                      GtError *err)
{
  char currentchar;
  bool firstsymbol = true;
  GtSeqIteratorQualMembers *pvt = seqit->pvt;
  gt_error_check(err);
  gt_assert(seqit && buffer);
  gt_assert(gt_str_length(buffer) == 0);
  if ((currentchar = fastq_buf_getchar(seqit)) == EOF)
    return EOF;
  pvt->currentread++;
  if (currentchar != startchar) {
    gt_error_set(err, "'%c' expected, '%c' encountered instead in line %lu",
                      startchar,
                      currentchar,
                      pvt->curline);
    return -2;
  }
  while (currentchar != GT_FASTQ_NEWLINESYMBOL) {
    if (!firstsymbol)
      gt_str_append_char(buffer, currentchar);
    else
      firstsymbol = false;
    if ((currentchar = fastq_buf_getchar(seqit)) == EOF)
      return EOF;
    pvt->currentread++;
  }
  pvt->curline++;
  return 0;
}

static int parse_fastq_sequence(GtSeqIteratorQual *seqit,
                                GtError *err)
{
  char currentchar;
  GtSeqIteratorQualMembers *pvt = seqit->pvt;
  gt_error_check(err);
  gt_assert(seqit);
  gt_assert(gt_str_length(pvt->sequencebuffer) == 0);
  /* read sequence */
  if ((currentchar = fastq_buf_getchar(seqit)) == EOF)
    return EOF;
  while (currentchar != GT_FASTQ_QUAL_SEPARATOR_CHAR) {
    if (currentchar != '\n' && currentchar != ' ') {
      if (pvt->symbolmap) {
        int charcode;
        charcode = pvt->symbolmap[(unsigned int) currentchar];
        if (charcode == UNDEFCHAR) {
          gt_error_set(err, "illegal character '%c': file \"%s\", line %lu",
                            currentchar,
                            gt_str_array_get(pvt->filenametab,
                                             pvt->filenum),
                            (unsigned long) pvt->curline);
          return -2;
        }
        if (ISSPECIAL(charcode)) {
          pvt->lastspeciallength++;
        } else {
          if (pvt->lastspeciallength > 0)
            pvt->lastspeciallength = 0;
          if (pvt->chardisttab)
            pvt->chardisttab[(int) charcode]++;
        }
        gt_str_append_char(pvt->sequencebuffer, charcode);
      } else {
        gt_str_append_char(pvt->sequencebuffer, currentchar);
      }
    } else if (currentchar == '\n') {
      pvt->curline++;
    }
    if ((currentchar = fastq_buf_getchar(seqit)) == EOF)
      return EOF;
    pvt->currentread++;
  }
  if (!gt_str_length(pvt->sequencebuffer)) {
    gt_error_set(err, "empty sequence given in file '%s', line %lu",
                      gt_str_array_get(pvt->filenametab,
                                       pvt->filenum),
                      pvt->curline-1);
    return -2;
  }
  fastq_buf_ungetchar(seqit);
  return 0;
}

static inline int parse_fastq_qualities(GtSeqIteratorQual *seqit,
                                        GT_UNUSED GtError *err)
{
  char currentchar;
  unsigned long i = 0;
  GtSeqIteratorQualMembers *pvt = seqit->pvt;
  gt_assert(gt_str_length(pvt->sequencebuffer) > 0);
  if ((currentchar = fastq_buf_getchar(seqit)) == EOF)
    return EOF;
  pvt->currentread++;

  for (i=0;i<gt_str_length(pvt->sequencebuffer);i++) {
    if (currentchar != '\n' && currentchar != ' ') {
      gt_str_append_char(pvt->qualsbuffer, currentchar);
    } else if (currentchar == '\n') {
      pvt->curline++;
      i--;
    } else {
      i--;
    }
    if (i+1 == gt_str_length(pvt->sequencebuffer)) {
      pvt->curline++;
    }
    if ((currentchar = fastq_buf_getchar(seqit)) == EOF)
      return EOF;
    pvt->currentread++;
  }
  /* expect newline at end of qualities */
  if (currentchar != GT_FASTQ_NEWLINESYMBOL) {
    gt_error_set(err, "qualities string of sequence length %lu is not ended "
                      "by newline in file '%s', line %lu -- this may be a "
                      "sign for sequence and qualities strings of different "
                      "length",
                      gt_str_length(pvt->sequencebuffer),
                      gt_str_array_get(pvt->filenametab,
                                       pvt->filenum),
                      pvt->curline-1);
    return -2;
  }
  return 0;
}

#define gt_fastq_premature_end_check(had_err, pvt) \
  if (had_err == EOF) { \
    gt_error_set(err, "premature end of file '%s' in line %lu: " \
                      "file ended before end of block", \
                      gt_str_array_get((pvt)->filenametab, \
                                       (pvt)->filenum), \
                      (unsigned long) (pvt)->curline-1); \
    return -2; \
  }

static inline int parse_fastq_block(GtSeqIteratorQual *seqit, GtError *err)
{
  int had_err = 0;
  GtSeqIteratorQualMembers *pvt = seqit->pvt;
  GtSeqIteratorQualFASTQ *seqitf;
  gt_assert(seqit);
  gt_error_check(err);

  seqitf = gt_seqiterator_qual_fastq_cast(seqit);
  /* parse @<seqname> */
  had_err = parse_fastq_seqname(seqit,
                                pvt->descbuffer,
                                GT_FASTQ_BLOCK_START_CHAR,
                                err);
  if (!had_err) {
    /* parse sequence */
    had_err = parse_fastq_sequence(seqit, err);
    gt_fastq_premature_end_check(had_err, pvt);
  }
  if (!had_err) {
    /* parse +[seqname] */
    had_err = parse_fastq_seqname(seqit,
                                  seqitf->qdescbuffer,
                                  GT_FASTQ_QUAL_SEPARATOR_CHAR,
                                  err);
    gt_fastq_premature_end_check(had_err, pvt);
  }
  if (!had_err
      && gt_str_length(seqitf->qdescbuffer)
      && gt_str_cmp(pvt->descbuffer, seqitf->qdescbuffer) != 0)
  {
      gt_error_set(err, "sequence description '%s' is not equal to "
                        "qualities description '%s' in line %lu",
                        gt_str_get(pvt->descbuffer),
                        gt_str_get(seqitf->qdescbuffer),
                        pvt->curline-1);
      return -2;
  }
  if (!had_err) {
    /* parse qualities */
    had_err = parse_fastq_qualities(seqit, err);
    if (gt_str_length(pvt->qualsbuffer)
          != gt_str_length(pvt->sequencebuffer))
    {
      gt_error_set(err, "lengths of character sequence and qualities "
                        "sequence differ (%lu <-> %lu)",
                        gt_str_length(pvt->qualsbuffer),
                        gt_str_length(pvt->sequencebuffer));
      return -2;
    }
  }
  return had_err;
}

void gt_seqiterator_qual_fastq_set_symbolmap(GtSeqIteratorQual *seqit,
                                             const GtUchar *symbolmap)
{
  gt_assert(seqit);
  seqit->pvt->symbolmap = symbolmap;
}

void gt_seqiterator_qual_fastq_set_chardisttab(GtSeqIteratorQual *seqit,
                                               unsigned long *chardist)
{
  gt_assert(seqit && chardist);
  seqit->pvt->chardisttab = chardist;
}

uint64_t gt_seqiterator_qual_fastq_get_lastspeciallength(const
                                                         GtSeqIteratorQual *si)
{
  gt_assert(si);
  return si->pvt->lastspeciallength;
}

int gt_seqiterator_qual_fastq_next(GtSeqIteratorQual *seqit,
                                   const GtUchar **sequence,
                                   const GtUchar **qualities,
                                   unsigned long *len,
                                   char **desc,
                                   GtError *err)
{
  int errstatus = 0;
  GtSeqIteratorQualMembers *pvt = seqit->pvt;
  GtSeqIteratorQualFASTQ *seqitf;
  gt_assert(seqit && len && desc);

  seqitf = gt_seqiterator_qual_fastq_cast(seqit);
  gt_str_reset(seqitf->qdescbuffer);

  /* parse file */
  errstatus = parse_fastq_block(seqit, err);

  if (!errstatus) {
    *sequence = (GtUchar*) gt_str_get(pvt->sequencebuffer);
    *len = gt_str_length(pvt->sequencebuffer);
    *desc = gt_cstr_dup(gt_str_get(pvt->descbuffer));
    *qualities = (GtUchar*) gt_str_get(pvt->qualsbuffer);
    errstatus = 1;
  } else {
    if (errstatus == EOF) {
      /* we could not get a next entry from this file */
      /* can we open another? */
      if (pvt->filenum+1 < gt_str_array_size(pvt->filenametab)) {
        const char *filename;
        filename = gt_str_array_get(pvt->filenametab, ++pvt->filenum);
        gt_genfile_close(pvt->curfile);
        pvt->curfile = gt_genfile_xopen(filename, "r");
        pvt->curline = 1;
        /* get first entry from next file*/
        errstatus = parse_fastq_block(seqit, err);
        if (!errstatus) {
          *sequence = (GtUchar*) gt_str_get(pvt->sequencebuffer);
          *len = gt_str_length(pvt->sequencebuffer);
          *desc = gt_cstr_dup(gt_str_get(pvt->descbuffer));
          *qualities = (GtUchar*) gt_str_get(pvt->qualsbuffer);
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
gt_seqiterator_qual_fastq_getcurrentcounter(GtSeqIteratorQual *seqit,
                                            unsigned long long maxread)
{
  seqit->pvt->maxread = maxread;
  return &seqit->pvt->currentread;
}

void gt_seqiterator_qual_fastq_delete(GtSeqIteratorQual *seqit)
{
  GtSeqIteratorQualFASTQ *seqitf;
  if (!seqit) return;
  seqitf = gt_seqiterator_qual_fastq_cast(seqit);
  gt_str_delete(seqitf->qdescbuffer);
}

const GtSeqIteratorQualClass* gt_seqiterator_qual_fastq_class(void)
{
  static const GtSeqIteratorQualClass sbc = { sizeof (GtSeqIteratorQualFASTQ),
                                gt_seqiterator_qual_fastq_next,
                                gt_seqiterator_qual_fastq_set_symbolmap,
                                gt_seqiterator_qual_fastq_set_chardisttab,
                                gt_seqiterator_qual_fastq_get_lastspeciallength,
                                gt_seqiterator_qual_fastq_getcurrentcounter,
                                gt_seqiterator_qual_fastq_delete };
  return &sbc;
}

GtSeqIteratorQual* gt_seqiterator_qual_fastq_new(const GtStrArray *filenametab,
                                                 GT_UNUSED GtError *err)
{
  GtSeqIteratorQual *seqit;
  GtSeqIteratorQualFASTQ *seqitf;
  gt_assert(filenametab);
  seqit = gt_seqiterator_qual_create(gt_seqiterator_qual_fastq_class());
  seqitf = gt_seqiterator_qual_fastq_cast(seqit);
  seqitf->qdescbuffer = gt_str_new();
  seqit->pvt->curfile = gt_genfile_xopen(gt_str_array_get(filenametab, 0), "r");
  seqit->pvt->filenametab = filenametab;
  seqit->pvt->curline = 1;
  return seqit;
}
