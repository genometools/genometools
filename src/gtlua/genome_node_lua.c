/*
  Copyright (c) 2007-2008 Gordon Gremme <gordon@gremme.org>
  Copyright (c) 2014 Sascha Steinbiss <ss34@sanger.ac.uk>
  Copyright (c) 2016 Daniel Standage <daniel.standage@gmail.com>
  Copyright (c) 2007-2008 Center for Bioinformatics, University of Hamburg
  Copyright (c) 2014 Genome Research Ltd.

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
#include "lauxlib.h"
#include "core/assert_api.h"
#include "core/ma.h"
#include "core/phase_api.h"
#include "core/str_array_api.h"
#include "core/symbol_api.h"
#include "core/unused_api.h"
#include "extended/extract_feature_sequence.h"
#include "extended/feature_node.h"
#include "extended/feature_node_iterator_api.h"
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
  const char *directive = NULL, *data = NULL;
  gt_assert(L);
  /* get_check parameters */
  directive = luaL_checkstring(L, 1);
  if (!lua_isnil(L, 2))
    data = luaL_checkstring(L, 2);
  /* construct object */
  mn = lua_newuserdata(L, sizeof (GtGenomeNode*));
  gt_assert(directive);
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

static int feature_node_lua_get_phase(lua_State *L)
{
  GtGenomeNode **gn = check_genome_node(L, 1);
  GtFeatureNode *fn;
  GtPhase phase;
  char phasebuf[2];
  /* make sure we get a feature node */
  fn = gt_feature_node_try_cast(*gn);
  luaL_argcheck(L, fn, 1, "not a feature node");
  phase = gt_feature_node_get_phase(fn);
  phasebuf[0] = GT_PHASE_CHARS[(int) phase];
  phasebuf[1] = '\0';
  lua_pushstring(L, phasebuf);
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

static int genome_node_lua_change_seqid(lua_State *L)
{
  const char *seqid;
  GtStr *seqid_str;
  GtGenomeNode **gn = check_genome_node(L, 1);
  seqid = luaL_checkstring(L, 2);
  seqid_str = gt_str_new_cstr(seqid);
  gt_genome_node_change_seqid(*gn, seqid_str);
  gt_str_delete(seqid_str);
  return 0;
}

static int genome_node_lua_set_range(lua_State *L)
{
  GtRange *rng;
  GtGenomeNode **gn = check_genome_node(L, 1);
  rng = check_range(L, 2);
  gt_genome_node_set_range(*gn, rng);
  return 0;
}

static int feature_node_lua_set_strand(lua_State *L)
{
  const char *str;
  GtGenomeNode **gn = check_genome_node(L, 1);
  GtFeatureNode *fn;
  /* make sure we get a feature node */
  fn = gt_feature_node_try_cast(*gn);
  luaL_argcheck(L, fn, 1, "not a feature node");
  str = luaL_checkstring(L, 2);
  luaL_argcheck(L, strlen(str) == 1 && strchr(GT_STRAND_CHARS, str[0]),
                2, "must be one of '" GT_STRAND_CHARS "'");
  gt_feature_node_set_strand(fn, gt_strand_get(str[0]));
  return 0;
}

static int feature_node_lua_set_score(lua_State *L)
{
  float sc;
  GtGenomeNode **gn = check_genome_node(L, 1);
  GtFeatureNode *fn;
  /* make sure we get a feature node */
  fn = gt_feature_node_try_cast(*gn);
  luaL_argcheck(L, fn, 1, "not a feature node");
  sc = luaL_checknumber(L, 2);
  gt_feature_node_set_score(fn, sc);
  return 0;
}

static int feature_node_lua_set_phase(lua_State *L)
{
  const char *p;
  GtGenomeNode **gn = check_genome_node(L, 1);
  GtFeatureNode *fn;
  /* make sure we get a feature node */
  fn = gt_feature_node_try_cast(*gn);
  luaL_argcheck(L, fn, 1, "not a feature node");
  p = luaL_checkstring(L, 2);
  luaL_argcheck(L, strlen(p) == 1 && strchr(GT_PHASE_CHARS, p[0]),
                2, "must be one of '" GT_PHASE_CHARS "'");
  gt_feature_node_set_phase(fn, gt_phase_get(p[0]));
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
  luaL_argcheck(L,
               !strcmp(gt_str_get(gt_genome_node_get_seqid((GtGenomeNode*) pf)),
                      gt_str_get(gt_genome_node_get_seqid((GtGenomeNode*) cf))),
                2, "has a sequence ID different to the one of its parent");
  gt_feature_node_add_child(pf, (GtFeatureNode*)
                                gt_genome_node_ref((GtGenomeNode*) cf));
  return 0;
}

static int genome_node_lua_add_attribute(lua_State *L)
{
  GtGenomeNode **node;
  const char *key, *val;
  GtFeatureNode *f;
  node = check_genome_node(L, 1);
  f = gt_feature_node_try_cast(*node);
  luaL_argcheck(L, f, 1, "not a feature node");
  key = luaL_checkstring(L, 2);
  luaL_argcheck(L, !gt_feature_node_get_attribute(f, key), 2,
                "attribute already present");
  val = luaL_checkstring(L, 3);
  gt_feature_node_add_attribute(f, key, val);
  return 0;
}

static int genome_node_lua_remove_attribute(lua_State *L)
{
  GtGenomeNode **node;
  const char *key;
  GtFeatureNode *f;
  node = check_genome_node(L, 1);
  f = gt_feature_node_try_cast(*node);
  luaL_argcheck(L, f, 1, "not a feature node");
  key = luaL_checkstring(L, 2);
  luaL_argcheck(L, gt_feature_node_get_attribute(f, key), 2,
                "attribute not present in node");
  gt_feature_node_remove_attribute(f, key);
  return 0;
}

static int genome_node_lua_set_attribute(lua_State *L)
{
  GtGenomeNode **node;
  const char *key, *val;
  GtFeatureNode *f;
  node = check_genome_node(L, 1);
  f = gt_feature_node_try_cast(*node);
  luaL_argcheck(L, f, 1, "not a feature node");
  key = luaL_checkstring(L, 2);
  luaL_argcheck(L, strlen(key) > 0, 2,
                "key must have length > 0");
  val = luaL_checkstring(L, 3);
  luaL_argcheck(L, strlen(key) > 0, 3,
                "value must have length > 0");
  gt_feature_node_set_attribute(f, key, val);
  return 0;
}

typedef struct {
  GtUword cur_attr;
  GtStrArray *attribs;
  GtFeatureNode *fn;
} GtFeatureNodeLuaEachAttribInfo;

static int feature_node_lua_each_attribute_iter(lua_State *L) {
  GtFeatureNodeLuaEachAttribInfo *info;
  info = *(GtFeatureNodeLuaEachAttribInfo**) lua_touserdata(L,
                                                           lua_upvalueindex(1));
  if (!info->attribs) return 0;
  if (info->cur_attr < gt_str_array_size(info->attribs)) {
    const char *attr, *val;
    attr = gt_str_array_get(info->attribs, info->cur_attr);
    gt_assert(attr);
    val = gt_feature_node_get_attribute(info->fn, attr);
    gt_assert(val);
    lua_pushstring(L, attr);
    lua_pushstring(L, val);
    info->cur_attr++;
    return 2;
  } else
    return 0;
}

static int feature_node_lua_each_attribute_gc(lua_State *L) {
  GtFeatureNodeLuaEachAttribInfo *info;
  info = *(GtFeatureNodeLuaEachAttribInfo**) lua_touserdata(L, 1);
  if (info) {
    gt_str_array_delete(info->attribs);
    gt_genome_node_delete((GtGenomeNode*) info->fn);
    gt_free(info);
  }
  return 0;
}

static int feature_node_lua_each_attribute(lua_State *L) {
  GtGenomeNode **gn;
  GtFeatureNode *fn;
  GtFeatureNodeLuaEachAttribInfo **info;
  gn = check_genome_node(L, 1);
  /* make sure we get a feature node */
  fn = gt_feature_node_try_cast(*gn);
  luaL_argcheck(L, fn, 1, "not a feature node");
  info = (GtFeatureNodeLuaEachAttribInfo**) lua_newuserdata(L,
                                     sizeof (GtFeatureNodeLuaEachAttribInfo *));
  luaL_getmetatable(L, "GenomeTools.each_attrib");
  lua_setmetatable(L, -2);
  *info = gt_calloc(1, sizeof (GtFeatureNodeLuaEachAttribInfo));
  gt_assert(*info);
  (*info)->fn = (GtFeatureNode*) gt_genome_node_ref(*gn);
  (*info)->cur_attr = 0;
  (*info)->attribs = gt_feature_node_get_attribute_list(fn);
  lua_pushcclosure(L, feature_node_lua_each_attribute_iter, 1);
  return 1;
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

static int feature_node_lua_set_type(lua_State *L)
{
  GtGenomeNode **gn;
  const char *type = NULL;
  GtFeatureNode *fn;
  gn = check_genome_node(L, 1);
  /* make sure we get a feature node */
  fn = gt_feature_node_try_cast(*gn);
  luaL_argcheck(L, fn, 1, "not a feature node");
  type = luaL_checkstring(L, 2);
  gt_assert(type);
  luaL_argcheck(L, strlen(type) > 0, 2, "must not be empty");
  gt_feature_node_set_type(fn, type);
  return 0;
}

static int meta_node_lua_get_directive(lua_State *L)
{
  GtGenomeNode **gn;
  GtMetaNode *mn;
  const char *meta_directive = NULL;
  gn = check_genome_node(L, 1);
  /* make sure we get a meta node */
  mn = gt_meta_node_try_cast(*gn);
  luaL_argcheck(L, mn, 1, "not a meta node");
  meta_directive = gt_meta_node_get_directive(mn);
  gt_assert(meta_directive);
  lua_pushstring(L, meta_directive);
  return 1;
}

static int meta_node_lua_get_data(lua_State *L)
{
  GtGenomeNode **gn;
  GtMetaNode *mn;
  const char *meta_data = NULL;
  gn = check_genome_node(L, 1);
  /* make sure we get a meta node */
  mn = gt_meta_node_try_cast(*gn);
  luaL_argcheck(L, mn, 1, "not a meta node");
  meta_data = gt_meta_node_get_data(mn);
  if (meta_data) {
    lua_pushstring(L, meta_data);
  } else {
    lua_pushnil(L);
  }
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
  lua_pushstring(L, gt_str_get(sequence));
  gt_str_delete(sequence);
  gt_error_delete(err);
  return 1;
}

static int feature_node_lua_extract_and_translate_sequence(lua_State *L)
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
  if (gt_extract_and_translate_feature_sequence(fn, gt_symbol(type), join,
                                                NULL, NULL, *region_mapping,
                                                NULL,
                                                sequence, NULL, NULL, err)) {
    gt_str_delete(sequence);
    return gt_lua_error(L, err);
  }
  lua_pushstring(L, gt_str_get(sequence));
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

static int feature_node_lua_get_children_iter(lua_State *L)
{
  GtFeatureNodeIterator *it;
  GtFeatureNode *fn;
  it = *(GtFeatureNodeIterator**) lua_touserdata(L, lua_upvalueindex(1));
  if ((fn = gt_feature_node_iterator_next(it))) {
    gt_lua_genome_node_push(L, gt_genome_node_ref((GtGenomeNode*) fn));
    return 1;
  } else
    return 0;
}

static int feature_node_lua_get_children_gc(lua_State *L)
{
  GtFeatureNodeIterator *it;
  it = *(GtFeatureNodeIterator**) lua_touserdata(L, 1);
  if (it)
    gt_feature_node_iterator_delete(it);
  return 0;
}

static int feature_node_lua_get_children_generic(lua_State *L, bool direct)
{
  GtGenomeNode **gn;
  GtFeatureNode *fn;
  GtFeatureNodeIterator **it;
  gn = check_genome_node(L, 1);
  /* make sure we get a feature node */
  fn = gt_feature_node_try_cast(*gn);
  luaL_argcheck(L, fn, 1, "not a feature node");
  it = (GtFeatureNodeIterator**) lua_newuserdata(L,
                                              sizeof (GtFeatureNodeIterator *));
  luaL_getmetatable(L, "GenomeTools.get_children");
  lua_setmetatable(L, -2);
  if (direct)
    *it = gt_feature_node_iterator_new_direct(fn);
  else
    *it = gt_feature_node_iterator_new(fn);
  gt_assert(*it);
  lua_pushcclosure(L, feature_node_lua_get_children_iter, 1);
  return 1;
}

static int feature_node_lua_get_children(lua_State *L)
{
  return feature_node_lua_get_children_generic(L, false);
}

static int feature_node_lua_get_direct_children(lua_State *L)
{
  return feature_node_lua_get_children_generic(L, true);
}

static int feature_node_lua_has_child_of_type(lua_State *L)
{
  GtGenomeNode **gn;
  GtFeatureNode *fn, GT_UNUSED *fn2 = NULL;
  GtFeatureNodeIterator *it;
  bool found = false;
  const char *type;

  gn = check_genome_node(L, 1);
  /* make sure we get a feature node */
  fn = gt_feature_node_try_cast(*gn);
  luaL_argcheck(L, fn, 1, "not a feature node");
  type = gt_symbol(luaL_checkstring(L, 2));
  it = gt_feature_node_iterator_new(fn);
  /* skip parent node itself */
  fn2 = gt_feature_node_iterator_next(it);
  gt_assert(fn2);
  while (!found && (fn2 = gt_feature_node_iterator_next(it))) {
    found = (gt_feature_node_get_type(fn2) == type);
  }
  gt_feature_node_iterator_delete(it);
  lua_pushboolean(L, found);
  return 1;
}

static int feature_node_lua_count_children_of_type(lua_State *L)
{
  GtGenomeNode **gn;
  GtFeatureNode *fn, GT_UNUSED *fn2 = NULL;
  GtFeatureNodeIterator *it;
  bool count = 0;
  const char *type;

  gn = check_genome_node(L, 1);
  /* make sure we get a feature node */
  fn = gt_feature_node_try_cast(*gn);
  luaL_argcheck(L, fn, 1, "not a feature node");
  type = gt_symbol(luaL_checkstring(L, 2));
  it = gt_feature_node_iterator_new(fn);
  /* skip parent node itself */
  fn2 = gt_feature_node_iterator_next(it);
  gt_assert(fn2);
  while ((fn2 = gt_feature_node_iterator_next(it))) {
    if (gt_feature_node_get_type(fn2) == type) {
      count++;
    }
  }
  gt_feature_node_iterator_delete(it);
  lua_pushnumber(L, count);
  return 1;
}

static int genome_node_lua_tostring (lua_State *L)
{
  GtGenomeNode **gn;
  char buf[BUFSIZ];
  gn = check_genome_node(L, 1);
  if (gt_feature_node_try_cast(*gn)) {
    GtFeatureNode *fn = gt_feature_node_try_cast(*gn);
    (void) snprintf(buf, BUFSIZ, "feature: %s "GT_WU"-"GT_WU" %c",
                    gt_feature_node_get_type(fn),
                    gt_genome_node_get_start(*gn),
                    gt_genome_node_get_end(*gn),
                    GT_STRAND_CHARS[gt_feature_node_get_strand(fn)]);
  } else if (gt_region_node_try_cast(*gn)) {
    GtRange rng = gt_genome_node_get_range(*gn);
    (void) snprintf(buf, BUFSIZ, "region: %s "GT_WU"-"GT_WU,
                    gt_str_get(gt_genome_node_get_seqid(*gn)),
                    rng.start, rng.end);
  }
  lua_pushfstring(L, "%s", buf);
  return 1;
}

static int genome_node_lua_eq(lua_State *L)
{
  GtGenomeNode **gn1, **gn2;
  bool ret = false;
  gt_assert(L);
  gn1 = check_genome_node(L, 1);
  gn2 = check_genome_node(L, 2);
  if (*gn1 == *gn2)
    ret = true;
  lua_pushboolean(L, ret);
  return 1;
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
  { "set_range", genome_node_lua_set_range },
  { "get_seqid", genome_node_lua_get_seqid },
  { "change_seqid", genome_node_lua_change_seqid },
  { "get_strand", feature_node_lua_get_strand },
  { "set_strand", feature_node_lua_set_strand },
  { "get_source", feature_node_lua_get_source },
  { "set_source", feature_node_lua_set_source },
  { "get_score", feature_node_lua_get_score },
  { "set_score", feature_node_lua_set_score },
  { "get_phase", feature_node_lua_get_phase },
  { "set_phase", feature_node_lua_set_phase },
  { "get_attribute", feature_node_lua_get_attribute },
  { "get_exons", feature_node_lua_get_exons },
  { "accept", genome_node_lua_accept },
  { "add_attribute", genome_node_lua_add_attribute },
  { "remove_attribute", genome_node_lua_remove_attribute },
  { "set_attribute", genome_node_lua_set_attribute },
  { "add_child", genome_node_lua_add_child },
  { "mark", genome_node_lua_mark },
  { "is_marked", genome_node_lua_is_marked },
  { "contains_marked", genome_node_lua_contains_marked },
  { "output_leading", feature_node_lua_output_leading },
  { "get_type", feature_node_lua_get_type },
  { "set_type", feature_node_lua_set_type },
  { "get_children", feature_node_lua_get_children },
  { "get_direct_children", feature_node_lua_get_direct_children },
  { "children", feature_node_lua_get_children },
  { "direct_children", feature_node_lua_get_direct_children },
  { "attribute_pairs", feature_node_lua_each_attribute },
  { "has_child_of_type", feature_node_lua_has_child_of_type },
  { "count_children_of_type", feature_node_lua_count_children_of_type },
  { "extract_sequence", feature_node_lua_extract_sequence },
  { "extract_and_translate_sequence",
                              feature_node_lua_extract_and_translate_sequence },
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
  /* set its _tostring field */
  lua_pushstring(L, "__tostring");
  lua_pushcfunction(L, genome_node_lua_tostring);
  lua_settable(L, -3);
  /* set its _eq field */
  lua_pushstring(L, "__eq");
  lua_pushcfunction(L, genome_node_lua_eq);
  lua_settable(L, -3);
  /* register functions */
  luaL_register(L, NULL, genome_node_lib_m);
  gt_lua_export_metatable(L, GENOME_NODE_METATABLE);
  luaL_register(L, "gt", genome_node_lib_f);
  /* child node iterator upvalue needs custom gc callback */
  luaL_newmetatable(L, "GenomeTools.get_children");
  lua_pushstring(L, "__gc");
  lua_pushcfunction(L, feature_node_lua_get_children_gc);
  lua_settable(L, -3);
  lua_pop(L, 1);
  luaL_newmetatable(L, "GenomeTools.each_attrib");
  lua_pushstring(L, "__gc");
  lua_pushcfunction(L, feature_node_lua_each_attribute_gc);
  lua_settable(L, -3);
  lua_pop(L, 2);
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
