/*
  Copyright (c) 2006-2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2006-2007 Center for Bioinformatics, University of Hamburg

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

#include <assert.h>
#include <stdlib.h>
#include "libgtcore/cstr.h"
#include "libgtcore/ma.h"
#include "libgtcore/unused.h"
#include "libgtext/comment.h"
#include "libgtext/genome_node_rep.h"

struct Comment
{
  const GenomeNode parent_instance;
  char *comment;
  Str *comment_str; /* used in comment_get_idstr() */
};

#define comment_cast(GN)\
        genome_node_cast(comment_class(), GN)

static void comment_free(GenomeNode *gn)
{
  Comment *c = comment_cast(gn);
  assert(c && c->comment);
  ma_free(c->comment);
  str_delete(c->comment_str);
}

static Str* comment_get_idstr(GenomeNode *gn)
{
  Comment *c;
  assert(gn);
  c = comment_cast(gn);
  return c->comment_str;
}

static Range comment_get_range(UNUSED GenomeNode *gn)
{
  Range range;
  range.start = 0;
  range.end = 0;
  return range;
}

static int comment_accept(GenomeNode *gn, GenomeVisitor *gv, Error *e)
{
  Comment *c;
  error_check(e);
  c = comment_cast(gn);
  return genome_visitor_visit_comment(gv, c, e);
}

const GenomeNodeClass* comment_class()
{
  static const GenomeNodeClass gnc = { sizeof (Comment),
                                       comment_free,
                                       NULL,
                                       comment_get_idstr,
                                       comment_get_range,
                                       NULL,
                                       NULL,
                                       comment_accept };
  return &gnc;
}

GenomeNode* comment_new(const char *comment, Str *filename,
                        unsigned long line_number)
{
  GenomeNode *gn = genome_node_create(comment_class(), filename, line_number);
  Comment *c = comment_cast(gn);
  assert(comment);
  c->comment = cstr_dup(comment);
  c->comment_str = str_new_cstr("");
  return gn;
}

const char* comment_get_comment(Comment *c)
{
  assert(c && c->comment);
  return c->comment;
}
