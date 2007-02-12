/*
  Copyright (c) 2006-2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2006-2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#include <assert.h>
#include <stdlib.h>
#include "comment.h"
#include "genome_node_rep.h"
#include "xansi.h"

struct Comment
{
  const GenomeNode parent_instance;
  char *comment;
};

#define comment_cast(GN)\
        genome_node_cast(comment_class(), GN)

static void comment_free(GenomeNode *gn)
{
  Comment *c = comment_cast(gn);
  assert(c && c->comment);
  free(c->comment);
}

static Str* comment_get_idstr(/*@unused@*/ GenomeNode *gn)
{
  static Str *comment_str;
  static unsigned int initialized = 0;

  if (!initialized) {
    comment_str = str_new_cstr("");
    initialized = 1;
  }
  return comment_str;
}

static Range comment_get_range(/*@unused@*/ GenomeNode *gn)
{
  Range range;
  range.start = 0;
  range.end = 0;
  return range;
}

static void comment_accept(GenomeNode *gn, Genome_visitor *gv, Log *l)
{
  Comment *c = comment_cast(gn);
  genome_visitor_visit_comment(gv, c, l);
}

const GenomeNodeClass* comment_class()
{
  static const GenomeNodeClass gnc = { sizeof(Comment), comment_free, NULL,
                                         comment_get_idstr, comment_get_range,
                                         NULL, NULL, NULL, NULL,
                                         comment_accept };
  return &gnc;
}

GenomeNode* comment_new(const char *comment,
                         const char *filename,
                         unsigned long line_number)
{
  GenomeNode *gn = genome_node_create(comment_class(), filename, line_number);
  Comment *c = comment_cast(gn);
  assert(comment);
  c->comment = xstrdup(comment);
  return gn;
}

const char* comment_get_comment(Comment *c)
{
  assert(c && c->comment);
  return c->comment;
}
