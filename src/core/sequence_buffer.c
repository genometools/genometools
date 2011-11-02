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

#include <string.h>
#include "core/ensure.h"
#include "core/fa.h"
#include "core/fileutils_api.h"
#include "core/ma.h"
#include "core/sequence_buffer_rep.h"
#include "core/sequence_buffer_embl.h"
#include "core/sequence_buffer_fasta.h"
#include "core/sequence_buffer_fastq.h"
#include "core/sequence_buffer_gb.h"
#include "core/sequence_buffer_inline.h"
#include "core/unused_api.h"
#include "core/xansi_api.h"

GtSequenceBuffer*
gt_sequence_buffer_create(const GtSequenceBufferClass *sic)
{
  GtSequenceBuffer *si;
  gt_assert(sic && sic->size);
  si = gt_calloc(1, sic->size);
  si->c_class = sic;
  si->pvt = gt_calloc(1, sizeof (GtSequenceBufferMembers));
  return si;
}

void gt_sequence_buffer_delete(GtSequenceBuffer *si)
{
  if (!si) return;
  if (si->pvt->reference_count) {
    si->pvt->reference_count--;
    return;
  }
  gt_assert(si->c_class && si->c_class->free);
  si->c_class->free(si);
  gt_free(si->pvt);
  gt_free(si);
}

GtSequenceBuffer* gt_sequence_buffer_ref(GtSequenceBuffer *sb)
{
  gt_assert(sb);
  sb->pvt->reference_count++;
  return sb;
}

GtSequenceBuffer* gt_sequence_buffer_new_guess_type(const GtStrArray *seqs,
                                                    GtError *err)
{
  GtFile *file;
  GtSequenceBuffer *sb;
  char firstcontents[BUFSIZ];
  gt_assert(seqs);
  gt_error_check(err);

  if (gt_str_array_size(seqs) == 0) {
    gt_error_set(err, "cannot guess file type of file %s -- no sequence files "
                      "given",
                      gt_str_array_get(seqs, 0));
    return NULL;
  }

  memset(firstcontents, 0, BUFSIZ);
  file = gt_file_open(gt_file_mode_determine(gt_str_array_get(seqs, 0)),
                      gt_str_array_get(seqs, 0),
                      "rb", err);
  if (!file)
    return NULL;
  gt_file_xread(file, &firstcontents, BUFSIZ-1);
  gt_file_delete(file);

  if (gt_sequence_buffer_embl_guess(firstcontents)) {
    sb = gt_sequence_buffer_embl_new(seqs);
  } else if (gt_sequence_buffer_fasta_guess(firstcontents)) {
    sb = gt_sequence_buffer_fasta_new(seqs);
  } else if (gt_sequence_buffer_gb_guess(firstcontents)) {
    sb = gt_sequence_buffer_gb_new(seqs);
    } else if (gt_sequence_buffer_fastq_guess(firstcontents)) {
    sb = gt_sequence_buffer_fastq_new(seqs);
  } else {
    gt_error_set(err, "cannot guess file type of file %s -- unknown file "
                      "contents",
                      gt_str_array_get(seqs, 0));
    return NULL;
  }
  return sb;
}

unsigned long gt_sequence_buffer_get_file_index(GtSequenceBuffer *si)
{
  gt_assert(si && si->c_class && si->c_class->get_file_index);
  return si->c_class->get_file_index(si);
}

void* gt_sequence_buffer_cast(GT_UNUSED const GtSequenceBufferClass *sic,
                              GtSequenceBuffer *si)
{
  gt_assert(sic && si && si->c_class == sic);
  return si;
}

void gt_sequence_buffer_set_symbolmap(GtSequenceBuffer *si, const GtUchar *m)
{
  gt_assert(si && si->pvt);
  si->pvt->symbolmap = m;
}

void gt_sequence_buffer_set_desc_buffer(GtSequenceBuffer *si, GtDescBuffer *db)
{
  gt_assert(si && si->pvt && db);
  si->pvt->descptr = db;
}

void gt_sequence_buffer_set_filelengthtab(GtSequenceBuffer *si,
                                          GtFilelengthvalues *flv)
{
  gt_assert(si && si->pvt);
  si->pvt->filelengthtab = flv;
}

void gt_sequence_buffer_set_chardisttab(GtSequenceBuffer *si,
                                        unsigned long *chardisttab)
{
  gt_assert(si && si->pvt);
  si->pvt->chardisttab = chardisttab;
}

const unsigned long long*
gt_sequence_buffer_get_counter(const GtSequenceBuffer *si)
{
  gt_assert(si && si->pvt);
  return &si->pvt->counter;
}

int gt_sequence_buffer_advance(GtSequenceBuffer *sb, GtError *err)
{
  gt_assert(sb && sb->c_class && sb->c_class->advance);
  return sb->c_class->advance(sb, err);
}

int gt_sequence_buffer_next(GtSequenceBuffer *sb, GtUchar *val,
                            GtError *err)
{
  GtSequenceBufferMembers *pvt;
  pvt = sb->pvt;
  if (pvt->nextread >= pvt->nextfree)
  {
    if (pvt->complete)
    {
      return 0;
    }
    if (pvt->descptr && pvt->nextread > 0)
      gt_desc_buffer_reset(pvt->descptr);
    if (gt_sequence_buffer_advance(sb, err) != 0)
    {
      return -1;
    }
    pvt->nextread = 0;
    if (pvt->nextfree == 0)
    {
      return 0;
    }
  }
  *val = pvt->outbuf[pvt->nextread++];
  return 1;
}

int gt_sequence_buffer_next_with_original(GtSequenceBuffer *sb,
                                          GtUchar *val, char *orig,
                                          GtError *err)
{
  GtSequenceBufferMembers *pvt;
  pvt = sb->pvt;
  if (pvt->nextread >= pvt->nextfree)
  {
    if (pvt->complete)
    {
      return 0;
    }
    if (pvt->descptr && pvt->nextread > 0)
      gt_desc_buffer_reset(pvt->descptr);
    if (gt_sequence_buffer_advance(sb, err) != 0)
    {
      return -1;
    }
    pvt->nextread = 0;
    if (pvt->nextfree == 0)
    {
      return 0;
    }
  }
  *val = pvt->outbuf[pvt->nextread];
  *orig = pvt->outbuforig[pvt->nextread];
  pvt->nextread++;
  return 1;
}

int gt_sequence_buffer_unit_test(GtError *err)
{
  int had_err = 0;
  GtSequenceBuffer *sb;
  GtStrArray *testfiles;
  GtStr *tmpfilename;
  unsigned long i;
  FILE *tmpfp;

  gt_error_check(err);

  /* create test seqs */
  testfiles = gt_str_array_new();
  tmpfilename = gt_str_new();
  tmpfp = gt_xtmpfp(tmpfilename);
  fprintf(tmpfp, "ID   Foo\n"
                 "XX\n"
                 "DE   Bar\n"
                 "SQ\n"
                 "     gatcgcgcta\n"
                 "//\n");
  gt_fa_xfclose(tmpfp);
  gt_ensure(had_err, gt_file_exists(gt_str_get(tmpfilename)));
  gt_str_array_add(testfiles, tmpfilename);
  gt_str_reset(tmpfilename);
  tmpfp = gt_xtmpfp(tmpfilename);
  fprintf(tmpfp, ">test3\nseq3xx\n"
                 ">test4\nseq4xxx\n");
  gt_fa_xfclose(tmpfp);
  gt_ensure(had_err, gt_file_exists(gt_str_get(tmpfilename)));
  gt_str_array_add(testfiles, tmpfilename);
  gt_str_delete(tmpfilename);

  sb = gt_sequence_buffer_new_guess_type(testfiles, err);
  gt_ensure(had_err, sb != NULL);

  gt_sequence_buffer_delete(sb);
  for (i=0;i<gt_str_array_size(testfiles);i++) {
    const char *fn;
    fn = gt_str_array_get(testfiles, i);
    gt_xremove(fn);
    gt_ensure(had_err, !gt_file_exists(fn));
  }
  gt_str_array_delete(testfiles);

  return had_err;
}
