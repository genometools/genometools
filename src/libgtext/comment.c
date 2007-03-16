/*
  Copyright (c) 2006-2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2006-2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#include <assert.h>
#include <stdlib.h>
#include <libgtext/comment.h>
#include <libgtext/genome_node_rep.h>

struct Comment
{
  const GenomeNode parent_instance;
  char *comment;
  Str *comment_str; /* used in comment_get_idstr() */
};

#define comment_cast(GN)\
        genome_node_cast(comment_class(), GN)

static void comment_free(GenomeNode *gn, Env *env)
{
  Comment *c = comment_cast(gn);
  assert(c && c->comment);
  env_ma_free(c->comment, env);
  str_delete(c->comment_str, env);
}

static Str* comment_get_idstr(GenomeNode *gn)
{
  Comment *c;
  assert(gn);
  c = comment_cast(gn);
  return c->comment_str;
}

static Range comment_get_range(/*@unused@*/ GenomeNode *gn)
{
  Range range;
  range.start = 0;
  range.end = 0;
  return range;
}

static int comment_accept(GenomeNode *gn, GenomeVisitor *gv, Env *env)
{
  Comment *c;
  env_error_check(env);
  c = comment_cast(gn);
  return genome_visitor_visit_comment(gv, c, env);
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
                                       NULL,
                                       NULL,
                                       comment_accept };
  return &gnc;
}

GenomeNode* comment_new(const char *comment, const char *filename,
                        unsigned long line_number, Env *env)
{
  GenomeNode *gn = genome_node_create(comment_class(), filename, line_number,
                                      env);
  Comment *c = comment_cast(gn);
  assert(comment);
  c->comment = cstr_dup(comment, env);
  c->comment_str = str_new_cstr("", env);
  return gn;
}

const char* comment_get_comment(Comment *c)
{
  assert(c && c->comment);
  return c->comment;
}
