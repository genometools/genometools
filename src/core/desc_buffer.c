/*
  Copyright (c) 2011 Sascha Steinbiss <steinbiss@zbh.uni-hamburg.de>
  Copyright (c) 2011 Center for Bioinformatics, University of Hamburg

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
#include <string.h>
#include "core/assert_api.h"
#include "core/desc_buffer.h"
#include "core/dynalloc.h"
#include "core/ensure.h"
#include "core/log_api.h"
#include "core/ma.h"
#include "core/queue.h"
#include "core/unused_api.h"
#include "core/xansi_api.h"

#define GT_DESC_BUFFER_INIT_SIZE 8192

struct GtDescBuffer {
  char *buf;
  unsigned long length;
  size_t allocated;
  bool finished,
       dirty;
  GtQueue *startqueue;
  unsigned int reference_count;
};

GtDescBuffer* gt_desc_buffer_new(void)
{
  GtDescBuffer *db = gt_malloc(sizeof *db);
  db->buf = gt_calloc(GT_DESC_BUFFER_INIT_SIZE, sizeof (char));
  db->length = 0;
  db->allocated = GT_DESC_BUFFER_INIT_SIZE;
  db->finished = false;
  db->dirty = true;
  db->reference_count = 0;
  db->startqueue = gt_queue_new();
  gt_queue_add(db->startqueue, (void*) 0);
  return db;
}

void gt_desc_buffer_append_char(GtDescBuffer *db, char c)
{
  gt_assert(db);
  if (db->finished) {
    gt_queue_add(db->startqueue, (void*) (db->length));
    db->finished = false;
  }
  if (db->length + 2 > db->allocated) {
    db->buf = gt_dynalloc(db->buf, &db->allocated,
                          (db->length + 2) * sizeof (char));
  }
  db->buf[db->length++] = c;
}

const char* gt_desc_buffer_get_next(GtDescBuffer *db)
{
  gt_assert(db);
  unsigned long startpos = (unsigned long) gt_queue_get(db->startqueue);
  db->dirty = true;
  return db->buf + startpos;
}

void gt_desc_buffer_finish(GtDescBuffer *db)
{
  gt_assert(db);
  /* XXX: maybe do a gt_cstr_rtrim(..., ' ') equivalent? */
  gt_desc_buffer_append_char(db, '\0');
  db->finished = true;
}

unsigned long gt_desc_buffer_length(const GtDescBuffer *db)
{
  return db ? db->length : 0;
}

void gt_desc_buffer_reset(GtDescBuffer *db)
{
  unsigned long laststartpos;
  gt_assert(db);

  if (!db->dirty) return;
  if (gt_queue_size(db->startqueue) == 0) {
    db->length = 0;
    db->dirty = false;
    return;
  }
  laststartpos = (unsigned long) gt_queue_head(db->startqueue);
  if (laststartpos != 0) {
    laststartpos = (unsigned long) gt_queue_get(db->startqueue);
    db->length = db->length - laststartpos;
    if (db->length >= laststartpos) {
      /* strings overlap */
      memmove(db->buf, db->buf + laststartpos, db->length * sizeof (char));
    } else {
      /* no overlap */
      memcpy(db->buf, db->buf + laststartpos, db->length * sizeof (char));
    }
    gt_queue_add(db->startqueue, (void*) 0);
  }
  db->dirty = false;
}

GtDescBuffer* gt_desc_buffer_ref(GtDescBuffer *db)
{
  if (!db) return NULL;
  db->reference_count++;
  return db;
}

void gt_desc_buffer_delete(GtDescBuffer *db)
{
  if (!db) return;
  if (db->reference_count) {
    db->reference_count--;
    return;
  }
  gt_free(db->buf);
  gt_queue_delete(db->startqueue);
  gt_free(db);
}

static int queueprinter(void **elem, GT_UNUSED void *info,
                        GT_UNUSED GtError *err)
{
  printf("%lu ", *(unsigned long*) elem);
  return 0;
}

void gt_desc_buffer_show(GtDescBuffer *db)
{
  printf("str: %s\n", db->buf);
  printf("queue: ");
  gt_queue_iterate(db->startqueue, queueprinter, NULL, NULL);
  printf("\n");
}

int gt_desc_buffer_unit_test(GtError *err)
{
  GtDescBuffer *s;
  static char *strs[] = { "foo", "bar", "baz"};
  const char *ret;
  int had_err = 0;
  unsigned long i, j;
  gt_error_check(err);

  s = gt_desc_buffer_new();
  ret = gt_desc_buffer_get_next(s);
  gt_ensure(had_err, strcmp(ret, "") == 0);
  gt_ensure(had_err, ret == s->buf);
  gt_ensure(had_err, gt_desc_buffer_length(s) == 0);
  gt_desc_buffer_delete(s);

  s = gt_desc_buffer_new();
  for (i = 0; i < strlen(strs[0]); i++) {
    gt_desc_buffer_append_char(s, strs[0][i]);
  }
  gt_desc_buffer_finish(s);
  ret = gt_desc_buffer_get_next(s);
  gt_ensure(had_err, strcmp(ret, strs[0]) == 0);
  gt_ensure(had_err, ret == s->buf);
  gt_ensure(had_err, gt_desc_buffer_length(s) == 4);
  gt_desc_buffer_delete(s);

  s = gt_desc_buffer_new();
  for (j = 0; j < 2; j++) {
    for (i = 0; i < strlen(strs[j]); i++) {
      gt_desc_buffer_append_char(s, strs[j][i]);
    }
    gt_desc_buffer_finish(s);
  }
  ret = gt_desc_buffer_get_next(s);
  gt_ensure(had_err, strcmp(ret, strs[0]) == 0);
  gt_ensure(had_err, ret == s->buf);
  ret = gt_desc_buffer_get_next(s);
  gt_ensure(had_err, strcmp(ret, strs[1]) == 0);
  gt_ensure(had_err, ret == s->buf+4);
  gt_ensure(had_err, gt_desc_buffer_length(s) == 8);
  gt_desc_buffer_delete(s);

  s = gt_desc_buffer_new();
  for (j = 0; j < 3; j++) {
    for (i = 0; i < strlen(strs[j]); i++) {
      gt_desc_buffer_append_char(s, strs[j][i]);
    }
    gt_desc_buffer_finish(s);
  }
  ret = gt_desc_buffer_get_next(s);
  gt_ensure(had_err, strcmp(ret, strs[0]) == 0);
  gt_ensure(had_err, ret == s->buf);
  ret = gt_desc_buffer_get_next(s);
  gt_ensure(had_err, strcmp(ret, strs[1]) == 0);
  gt_ensure(had_err, ret == s->buf+4);
  ret = gt_desc_buffer_get_next(s);
  gt_ensure(had_err, strcmp(ret, strs[2]) == 0);
  gt_ensure(had_err, ret == s->buf+8);
  gt_ensure(had_err, gt_desc_buffer_length(s) == 12);
  gt_desc_buffer_delete(s);

  return had_err;
}
