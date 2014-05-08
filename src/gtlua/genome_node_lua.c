/*
  Copyright (c) 2007-2008 Gordon Gremme <gordon@gremme.org>
  Copyright (c) 2007-2008 Center for Bioinformatics, University of Hamburg

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

#include "lauxlib.h"
#include "core/assert_api.h"
#include "core/symbol_api.h"
#include "extended/extract_feature_sequence.h"
#include "extended/feature_node.h"
#include "extended/genome_node.h"
#include "extended/region_node.h"
#include "extended/sequence_node_api.h"
#include "extended/meta_node_api.h"
#include "extended/gff3_output.h"
#include "extended/luahelper.h"
#include "gtlua/genome_node_lua.h"
#include "gtlua/genome_visitor_lua.h"
#include "gtlua/gt_lua.h"
#include "gtlua/range_lua.h"
#include "gtlua/region_mapping_lua.h"

static int feature_node_lua_new(lua_State *L)
{
  GtGenomeNode **gf;
  GtUword startpos, endpos;
  GtStrand strand;
  const char *seqid, *type, *strand_str;
  size_t length;
  GtStr *seqid_str;
  gt_assert(L);
  /* get/check parameters */
  seqid = luaL_checkstring(L, 1);
  type = luaL_checkstring(L, 2);
  startpos = luaL_checklong(L, 3);
  endpos   = luaL_checklong(L, 4);
  luaL_argcheck(L, startpos > 0, 3, "must be > 0");
  luaL_argcheck(L, endpos > 0, 4, "must be > 0");
  luaL_argcheck(L, startpos <= endpos, 3, "must be <= endpos");
  strand_str = luaL_checklstring(L, 5, &length);
  luaL_argcheck(L, length == 1, 5, "strand string must have length 1");
  luaL_argcheck(L, (strand = gt_strand_get(strand_str[0])) !=
                    GT_NUM_OF_STRAND_TYPES, 5, "invalid strand");
  /* construct object */
  gf = lua_newuserdata(L, sizeof (GtGenomeNode*));
  seqid_str = gt_str_new_cstr(seqid);
  *gf = gt_feature_node_new(seqid_str, type, startpos, endpos, strand);
  gt_str_delete(seqid_str);
  gt_assert(*gf);
  luaL_getmetatable(L, GENOME_NODE_METATABLE);
  lua_setmetatable(L, -2);
  return 1;
}

static int region_node_lua_new(lua_State *L)
{
  GtGenomeNode **rn;
  GtUword startpos, endpos;
  const char *seqid;
  GtStr *seqid_str;
  gt_assert(L);
  /* get_check parameters */
  seqid = luaL_checkstring(L, 1);
  startpos = luaL_checklong(L, 2);
  endpos   = luaL_checklong(L, 3);
  luaL_argcheck(L, startpos > 0, 2, "must be > 0");
  luaL_argcheck(L, endpos > 0, 3, "must be > 0");
  luaL_argcheck(L, startpos <= endpos, 2, "must be <= endpos");
  /* construct object */
  rn = lua_newuserdata(L, sizeof (GtGenomeNode*));
  seqid_str = gt_str_new_cstr(seqid);
  *rn = gt_region_node_new(seqid_str, startpos, endpos);
  gt_str_delete(seqid_str);
  gt_assert(*rn);
  luaL_getmetatable(L, GENOME_NODE_METATABLE);
  lua_setmetatable(L, -2);
  return 1;
}

static int meta_node_lua_new(lua_State *L)
{
  GtGenomeNode **mn;
  const char *directive, *data;
  gt_assert(L);
  /* get_check parameters */
  directive = luaL_checkstring(L, 1);
  data = luaL_checkstring(L, 2);
  /* construct object */
  mn = lua_newuserdata(L, sizeof (GtGenomeNode*));
  *mn = gt_meta_node_new(directive, data);
  gt_assert(*mn);
  luaL_getmetatable(L, GENOME_NODE_METATABLE);
  lua_setmetatable(L, -2);
  return 1;
}

static int comment_node_lua_new(lua_State *L)
{
  GtGenomeNode **cn;
  const char *comment;
  gt_assert(L);
  /* get_check parameters */
  comment = luaL_checkstring(L, 1);
  /* construct object */
  cn = lua_newuserdata(L, sizeof (GtGenomeNode*));
  *cn = gt_comment_node_new(comment);
  gt_assert(*cn);
  luaL_getmetatable(L, GENOME_NODE_METATABLE);
  lua_setmetatable(L, -2);
  return 1;
}

static int sequence_node_lua_new(lua_State *L)
{
  GtGenomeNode **sn;
  GtStr *seq_str;
  const char *desc, *seq;
  gt_assert(L);
  /* get_check parameters */
  desc = luaL_checkstring(L, 1);
  seq = luaL_checkstring(L, 2);
  /* construct object */
  sn = lua_newuserdata(L, sizeof (GtGenomeNode*));
  seq_str = gt_str_new_cstr(seq);
  *sn = gt_sequence_node_new(desc, seq_str);
  gt_assert(*sn);
  gt_str_delete(seq_str);
  luaL_getmetatable(L, GENOME_NODE_METATABLE);
  lua_setmetatable(L, -2);
  return 1;
}

static int genome_node_lua_get_filename(lua_State *L)
{
  GtGenomeNode **gn = check_genome_node(L, 1);
  lua_pushstring(L, gt_genome_node_get_filename(*gn));
  return 1;
}

static int genome_node_lua_get_line_number(lua_State *L)
{
  GtGenomeNode **gn = check_genome_node(L, 1);
  lua_pushnumber(L, gt_genome_node_get_line_number(*gn));
  return 1;
}

static int genome_node_lua_get_range(lua_State *L)
{
  GtGenomeNode **gn = check_genome_node(L, 1);
  return gt_lua_range_push(L, gt_genome_node_get_range(*gn));
}

static int genome_node_lua_get_seqid(lua_State *L)
{
  GtStr *seqid;
  GtGenomeNode **gn = check_genome_node(L, 1);
  if ((seqid = gt_genome_node_get_seqid(*gn)))
    lua_pushstring(L, gt_str_get(seqid));
  else
    lua_pushnil(L);
  return 1;
}

static int feature_node_lua_get_strand(lua_State *L)
{
  GtGenomeNode **gn = check_genome_node(L, 1);
  GtFeatureNode *fn;
  char strand_char[2];
  /* make sure we get a feature node */
  fn = gt_feature_node_try_cast(*gn);
  luaL_argcheck(L, fn, 1, "not a feature node");
  strand_char[0] = GT_STRAND_CHARS[gt_feature_node_get_strand(fn)];
  strand_char[1] = '\0';
  lua_pushstring(L, strand_char);
  return 1;
}

static int feature_node_lua_get_source(lua_State *L)
{
  GtGenomeNode **gn = check_genome_node(L, 1);
  GtFeatureNode *fn;
  /* make sure we get a feature node */
  fn = gt_feature_node_try_cast(*gn);
  luaL_argcheck(L, fn, 1, "not a feature node");
  lua_pushstring(L, gt_feature_node_get_source(fn));
  return 1;
}

static int feature_node_lua_get_score(lua_State *L)
{
  GtGenomeNode **gn = check_genome_node(L, 1);
  GtFeatureNode *fn;
  /* make sure we get a feature node */
  fn = gt_feature_node_try_cast(*gn);
  luaL_argcheck(L, fn, 1, "not a feature node");
  if (gt_feature_node_score_is_defined(fn))
    lua_pushnumber(L, gt_feature_node_get_score(fn));
  else
    lua_pushnil(L);
  return 1;
}

static int feature_node_lua_get_attribute(lua_State *L)
{
  GtGenomeNode **gn = check_genome_node(L, 1);
  const char *attr = NULL, *attrval = NULL;
  attr = luaL_checkstring(L, 2);
  GtFeatureNode *fn;
  /* make sure we get a feature node */
  fn = gt_feature_node_try_cast(*gn);
  luaL_argcheck(L, fn, 1, "not a feature node");
  attrval = gt_feature_node_get_attribute(fn, attr);
  if (attrval)
    lua_pushstring(L, attrval);
  else
    lua_pushnil(L);
  return 1;
}

static int feature_node_lua_get_exons(lua_State *L)
{
  GtGenomeNode **gn = check_genome_node(L, 1);
  GtArray *exons = gt_array_new(sizeof (GtGenomeNode*));
  GtUword i = 0;
  GtFeatureNode *fn;
  /* make sure we get a feature node */
  fn = gt_feature_node_try_cast(*gn);
  luaL_argcheck(L, fn, 1, "not a feature node");
  gt_feature_node_get_exons(fn, exons);
  lua_newtable(L);
  for (i = 0; i < gt_array_size(exons); i++) {
    lua_pushnumber(L, i+1);
    gt_lua_genome_node_push(L, (GtGenomeNode*)
                            gt_genome_node_ref(*(GtGenomeNode**)
                                               gt_array_get(exons, i)));
    lua_rawset(L, -3);
  }
  gt_array_delete(exons);
  return 1;
}

static int feature_node_lua_set_source(lua_State *L)
{
  const char *source;
  GtStr *source_str;
  GtGenomeNode **gn = check_genome_node(L, 1);
  GtFeatureNode *fn;
  /* make sure we get a feature node */
  fn = gt_feature_node_try_cast(*gn);
  luaL_argcheck(L, fn, 1, "not a feature node");
  source = luaL_checkstring(L, 2);
  source_str = gt_str_new_cstr(source);
  gt_feature_node_set_source(fn, source_str);
  gt_str_delete(source_str);
  return 0;
}

static int genome_node_lua_accept(lua_State *L)
{
  GtGenomeNode **gn;
  GtNodeVisitor **gv;
  GtError *err;
  gn = check_genome_node(L, 1);
  gv = check_genome_visitor(L, 2);
  err = gt_error_new();
  if (gt_genome_node_accept(*gn, *gv, err))
    return gt_lua_error(L, err);
  gt_error_delete(err);
  return 0;
}

static int genome_node_lua_add_child(lua_State *L)
{
  GtGenomeNode **parent, **child;
  GtFeatureNode *pf, *cf;
  parent = check_genome_node(L, 1);
  child  = check_genome_node(L, 2);
  pf = gt_feature_node_try_cast(*parent);
  luaL_argcheck(L, pf, 1, "not a feature node");
  cf = gt_feature_node_try_cast(*child);
  luaL_argcheck(L, cf, 2, "not a feature node");
  gt_feature_node_add_child(pf, (GtFeatureNode*)
                                gt_genome_node_ref((GtGenomeNode*) cf));
  return 0;
}

static int genome_node_lua_mark(lua_State *L)
{
  GtGenomeNode **gn = check_genome_node(L, 1);
  gt_feature_node_mark(gt_feature_node_cast(*gn));
  return 0;
}

static int genome_node_lua_is_marked(lua_State *L)
{
  GtGenomeNode **gn = check_genome_node(L, 1);
  lua_pushboolean(L, gt_feature_node_is_marked(gt_feature_node_cast(*gn)));
  return 1;
}

static int genome_node_lua_contains_marked(lua_State *L)
{
  GtGenomeNode **gn;
  gn = check_genome_node(L, 1);
  lua_pushboolean(L,
                  gt_feature_node_contains_marked(gt_feature_node_cast(*gn)));
  return 1;
}

static int feature_node_lua_output_leading(lua_State *L)
{
  GtGenomeNode **gn;
  GtFeatureNode *fn;
  gn = check_genome_node(L, 1);
  /* make sure we get a feature node */
  fn = gt_feature_node_try_cast(*gn);
  luaL_argcheck(L, fn, 1, "not a feature node");
  gt_gff3_output_leading(fn, NULL);
  return 0;
}

static int feature_node_lua_get_type(lua_State *L)
{
  GtGenomeNode **gn;
  GtFeatureNode *fn;
  gn = check_genome_node(L, 1);
  /* make sure we get a feature node */
  fn = gt_feature_node_try_cast(*gn);
  luaL_argcheck(L, fn, 1, "not a feature node");
  lua_pushstring(L, gt_feature_node_get_type(fn));
  return 1;
}

static int meta_node_lua_get_directive(lua_State *L)
{
  GtGenomeNode **gn;
  GtMetaNode *mn;
  gn = check_genome_node(L, 1);
  /* make sure we get a meta node */
  mn = gt_meta_node_try_cast(*gn);
  luaL_argcheck(L, mn, 1, "not a meta node");
  lua_pushstring(L, gt_meta_node_get_directive(mn));
  return 1;
}

static int meta_node_lua_get_data(lua_State *L)
{
  GtGenomeNode **gn;
  GtMetaNode *mn;
  gn = check_genome_node(L, 1);
  /* make sure we get a meta node */
  mn = gt_meta_node_try_cast(*gn);
  luaL_argcheck(L, mn, 1, "not a meta node");
  lua_pushstring(L, gt_meta_node_get_data(mn));
  return 1;
}

static int comment_node_lua_get_comment(lua_State *L)
{
  GtGenomeNode **gn;
  GtCommentNode *cn;
  gn = check_genome_node(L, 1);
  /* make sure we get a meta node */
  cn = gt_comment_node_try_cast(*gn);
  luaL_argcheck(L, cn, 1, "not a comment node");
  lua_pushstring(L, gt_comment_node_get_comment(cn));
  return 1;
}

static int sequence_node_lua_get_description(lua_State *L)
{
  GtGenomeNode **gn;
  GtSequenceNode *sn;
  gn = check_genome_node(L, 1);
  /* make sure we get a sequence node */
  sn = gt_sequence_node_try_cast(*gn);
  luaL_argcheck(L, sn, 1, "not a sequence node");
  lua_pushstring(L, gt_sequence_node_get_description(sn));
  return 1;
}

static int sequence_node_lua_get_sequence(lua_State *L)
{
  GtGenomeNode **gn;
  GtSequenceNode *sn;
  gn = check_genome_node(L, 1);
  /* make sure we get a sequence node */
  sn = gt_sequence_node_try_cast(*gn);
  luaL_argcheck(L, sn, 1, "not a sequence node");
  lua_pushstring(L, gt_sequence_node_get_sequence(sn));
  return 1;
}

static int sequence_node_lua_get_sequence_length(lua_State *L)
{
  GtGenomeNode **gn;
  GtSequenceNode *sn;
  gn = check_genome_node(L, 1);
  /* make sure we get a sequence node */
  sn = gt_sequence_node_try_cast(*gn);
  luaL_argcheck(L, sn, 1, "not a sequence node");
  lua_pushnumber(L, gt_sequence_node_get_sequence_length(sn));
  return 1;
}

static int feature_node_lua_extract_sequence(lua_State *L)
{
  GtGenomeNode **gn;
  GtFeatureNode *fn;
  const char *type;
  bool join;
  GtRegionMapping **region_mapping;
  GtStr *sequence;
  GtError *err;
  gn = check_genome_node(L, 1);
  /* make sure we get a feature node */
  fn = gt_feature_node_try_cast(*gn);
  luaL_argcheck(L, fn, 1, "not a feature node");
  type = luaL_checkstring(L, 2);
  join = lua_toboolean(L, 3);
  region_mapping = check_region_mapping(L, 4);
  err = gt_error_new();
  sequence = gt_str_new();
  if (gt_extract_feature_sequence(sequence, *gn, gt_symbol(type), join,
                                  NULL, NULL, *region_mapping, err)) {
    gt_str_delete(sequence);
    return gt_lua_error(L, err);
  }
  if (gt_str_length(sequence))
    lua_pushstring(L, gt_str_get(sequence));
  else
    lua_pushnil(L);
  gt_str_delete(sequence);
  gt_error_delete(err);
  return 1;
}

static int feature_node_lua_remove_leaf(lua_State *L)
{
  GtGenomeNode **parent, **leaf;
  GtFeatureNode *pf, *lf;
  parent = check_genome_node(L, 1);
  leaf  = check_genome_node(L, 2);
  pf = gt_feature_node_try_cast(*parent);
  luaL_argcheck(L, pf, 1, "not a feature node");
  lf = gt_feature_node_try_cast(*leaf);
  luaL_argcheck(L, lf, 2, "not a feature node");
  gt_feature_node_remove_leaf(pf, lf);
  return 0;
}

static int genome_node_lua_delete(lua_State *L)
{
  GtGenomeNode **gn;
  gn = check_genome_node(L, 1);
  gt_genome_node_delete(*gn);
  return 0;
}

static const struct luaL_Reg genome_node_lib_f [] = {
  { "feature_node_new", feature_node_lua_new },
  { "region_node_new", region_node_lua_new },
  { "meta_node_new", meta_node_lua_new },
  { "comment_node_new", comment_node_lua_new },
  { "sequence_node_new", sequence_node_lua_new },
  { NULL, NULL }
};

static const struct luaL_Reg genome_node_lib_m [] = {
  { "get_filename", genome_node_lua_get_filename },
  { "get_line_number", genome_node_lua_get_line_number },
  { "get_range", genome_node_lua_get_range },
  { "get_seqid", genome_node_lua_get_seqid },
  { "get_strand", feature_node_lua_get_strand },
  { "get_source", feature_node_lua_get_source },
  { "set_source", feature_node_lua_set_source },
  { "get_score", feature_node_lua_get_score },
  { "get_attribute", feature_node_lua_get_attribute },
  { "get_exons", feature_node_lua_get_exons },
  { "accept", genome_node_lua_accept },
  { "add_child", genome_node_lua_add_child },
  { "mark", genome_node_lua_mark },
  { "is_marked", genome_node_lua_is_marked },
  { "contains_marked", genome_node_lua_contains_marked },
  { "output_leading", feature_node_lua_output_leading },
  { "get_type", feature_node_lua_get_type },
  { "extract_sequence", feature_node_lua_extract_sequence },
  { "remove_leaf", feature_node_lua_remove_leaf },
  { "get_data", meta_node_lua_get_data },
  { "get_directive", meta_node_lua_get_directive },
  { "get_comment", comment_node_lua_get_comment },
  { "get_description", sequence_node_lua_get_description },
  { "get_sequence", sequence_node_lua_get_sequence },
  { "get_sequence_length", sequence_node_lua_get_sequence_length },
  { NULL, NULL }
};

int gt_lua_open_genome_node(lua_State *L)
{
#ifndef NDEBUG
  int stack_size;
#endif
  gt_assert(L);
#ifndef NDEBUG
  stack_size = lua_gettop(L);
#endif
  luaL_newmetatable(L, GENOME_NODE_METATABLE);
  /* metatable.__index = metatable */
  lua_pushvalue(L, -1); /* duplicate the metatable */
  lua_setfield(L, -2, "__index");
  /* set its _gc field */
  lua_pushstring(L, "__gc");
  lua_pushcfunction(L, genome_node_lua_delete);
  lua_settable(L, -3);
  /* register functions */
  luaL_register(L, NULL, genome_node_lib_m);
  gt_lua_export_metatable(L, GENOME_NODE_METATABLE);
  luaL_register(L, "gt", genome_node_lib_f);
  lua_pop(L, 1);
  gt_assert(lua_gettop(L) == stack_size);
  return 1;
}

void gt_lua_genome_node_push(lua_State *L, GtGenomeNode *gn)
{
  GtGenomeNode **gn_lua;
  gt_assert(L && gn);
  gn_lua = lua_newuserdata(L, sizeof (GtGenomeNode**));
  *gn_lua = gn;
  luaL_getmetatable(L, GENOME_NODE_METATABLE);
  lua_setmetatable(L, -2);
}
