/*
  Copyright (c) 2006-2008 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2006-2008 Center for Bioinformatics, University of Hamburg

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
#include "core/cstr.h"
#include "core/ma.h"
#include "core/unused_api.h"
#include "extended/comment.h"
#include "extended/genome_node_rep.h"

struct GT_Comment
{
  const GT_GenomeNode parent_instance;
  char *comment;
  GT_Str *gt_comment_str; /* used in gt_comment_get_idstr() */
};

#define gt_comment_cast(GN)\
        gt_genome_node_cast(gt_comment_class(), GN)

static void gt_comment_free(GT_GenomeNode *gn)
{
  GT_Comment *c = gt_comment_cast(gn);
  assert(c && c->comment);
  gt_free(c->comment);
  gt_str_delete(c->gt_comment_str);
}

static GT_Str* gt_comment_get_idstr(GT_GenomeNode *gn)
{
  GT_Comment *c;
  assert(gn);
  c = gt_comment_cast(gn);
  return c->gt_comment_str;
}

static GT_Range gt_comment_get_range(GT_UNUSED GT_GenomeNode *gn)
{
  GT_Range range;
  range.start = 0;
  range.end = 0;
  return range;
}

static int gt_comment_accept(GT_GenomeNode *gn, GenomeVisitor *gv,
                             GT_Error *err)
{
  GT_Comment *c;
  gt_error_check(err);
  c = gt_comment_cast(gn);
  return genome_visitor_visit_comment(gv, c, err);
}

const GT_GenomeNodeClass* gt_comment_class()
{
  static const GT_GenomeNodeClass gnc = { sizeof (GT_Comment),
                                       gt_comment_free,
                                       NULL,
                                       gt_comment_get_idstr,
                                       gt_comment_get_range,
                                       NULL,
                                       NULL,
                                       gt_comment_accept };
  return &gnc;
}

GT_GenomeNode* gt_comment_new(const char *comment)
{
  GT_GenomeNode *gn = gt_genome_node_create(gt_comment_class());
  GT_Comment *c = gt_comment_cast(gn);
  assert(comment);
  c->comment = gt_cstr_dup(comment);
  c->gt_comment_str = gt_str_new_cstr("");
  return gn;
}

const char* gt_comment_get_comment(const GT_Comment *c)
{
  assert(c && c->comment);
  return c->comment;
}
