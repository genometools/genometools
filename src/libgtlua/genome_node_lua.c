/*
  Copyright (c) 2007-2008 Gordon Gremme <gremme@zbh.uni-hamburg.de>
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
#include "libgtext/extract_feat_sequence.h"
#include "libgtext/genome_node.h"
#include "libgtext/gff3_output.h"
#include "libgtext/luahelper.h"
#include "libgtlua/genome_node_lua.h"
#include "libgtlua/genome_visitor_lua.h"
#include "libgtlua/range_lua.h"
#include "libgtlua/region_mapping_lua.h"

static int genome_feature_lua_new(lua_State *L)
{
  GenomeNode **gf;
  GenomeFeatureType type;
  Range *range;
  Strand strand;
  const char *type_str, *strand_str;
  size_t length;
  Str *filename;
  assert(L);
  /* get/check parameters */
  type_str = luaL_checkstring(L, 1);
  luaL_argcheck(L, !genome_feature_type_get(&type, type_str), 1,
                "invalid feature type");
  range = check_range(L, 2);
  strand_str = luaL_checklstring(L, 3, &length);
  luaL_argcheck(L, length == 1, 3, "strand string must have length 1");
  luaL_argcheck(L, (strand = strand_get(strand_str[0])) != NUM_OF_STRAND_TYPES,
                3, "invalid strand");
  /* construct object */
  gf = lua_newuserdata(L, sizeof (GenomeNode*));
  filename = str_new_cstr("Lua");
  *gf = genome_feature_new(type, *range, strand, filename, 0);
  str_delete(filename);
  assert(*gf);
  luaL_getmetatable(L, GENOME_NODE_METATABLE);
  lua_setmetatable(L, -2);
  return 1;
}

static int sequence_region_lua_new(lua_State *L)
{
  GenomeNode **sr;
  const char *seqid;
  Str *seqid_str, *filename;
  Range *range;
  assert(L);
  /* get_check parameters */
  seqid = luaL_checkstring(L, 1);
  range = check_range(L, 2);
  /* construct object */
  sr = lua_newuserdata(L, sizeof (GenomeNode*));
  seqid_str = str_new_cstr(seqid);
  filename = str_new_cstr("Lua");
  *sr = sequence_region_new(seqid_str, *range, filename, 0);
  str_delete(filename);
  str_delete(seqid_str);
  assert(*sr);
  luaL_getmetatable(L, GENOME_NODE_METATABLE);
  lua_setmetatable(L, -2);
  return 1;
}

static int genome_node_lua_get_filename(lua_State *L)
{
  GenomeNode **gn = check_genome_node(L, 1);
  lua_pushstring(L, genome_node_get_filename(*gn));
  return 1;
}

static int genome_node_lua_get_range(lua_State *L)
{
  GenomeNode **gn = check_genome_node(L, 1);
  return range_lua_push(L, genome_node_get_range(*gn));
}

static int genome_node_lua_get_seqid(lua_State *L)
{
  Str *seqid;
  GenomeNode **gn = check_genome_node(L, 1);
  if ((seqid = genome_node_get_seqid(*gn)))
    lua_pushstring(L, str_get(seqid));
  else
    lua_pushnil(L);
  return 1;
}

static int genome_node_lua_set_seqid(lua_State *L)
{
  const char *seqid;
  Str *seqid_str;
  GenomeNode **gn = check_genome_node(L, 1);
  seqid = luaL_checkstring(L, 2);
  seqid_str = str_new_cstr(seqid);
  genome_node_set_seqid(*gn, seqid_str);
  str_delete(seqid_str);
  return 0;
}

static int genome_feature_lua_get_strand(lua_State *L)
{
  GenomeNode **gn = check_genome_node(L, 1);
  GenomeFeature *gf;
  char strand_char[2];
  /* make sure we get a genome feature */
  gf = genome_node_cast(genome_feature_class(), *gn);
  luaL_argcheck(L, gf, 1, "not a genome feature");
  strand_char[0] = STRANDCHARS[genome_feature_get_strand(gf)];
  strand_char[1] = '\0';
  lua_pushstring(L, strand_char);
  return 1;
}

static int genome_feature_lua_get_source(lua_State *L)
{
  GenomeNode **gn = check_genome_node(L, 1);
  GenomeFeature *gf;
  /* make sure we get a genome feature */
  gf = genome_node_cast(genome_feature_class(), *gn);
  luaL_argcheck(L, gf, 1, "not a genome feature");
  lua_pushstring(L, genome_feature_get_source(gf));
  return 1;
}

static int genome_feature_lua_set_source(lua_State *L)
{
  const char *source;
  Str *source_str;
  GenomeNode **gn = check_genome_node(L, 1);
  GenomeFeature *gf;
  /* make sure we get a genome feature */
  gf = genome_node_cast(genome_feature_class(), *gn);
  luaL_argcheck(L, gf, 1, "not a genome feature");
  source = luaL_checkstring(L, 2);
  source_str = str_new_cstr(source);
  genome_feature_set_source(*gn, source_str);
  str_delete(source_str);
  return 0;
}

static int genome_node_lua_accept(lua_State *L)
{
  GenomeNode **gn;
  GenomeVisitor **gv;
  Error *e;
  gn = check_genome_node(L, 1);
  gv = check_genome_visitor(L, 2);
  e = error_new();
  if (genome_node_accept(*gn, *gv, e))
    return lua_gt_error(L, e);
  error_delete(e);
  return 0;
}

static int genome_node_lua_is_part_of_genome_node(lua_State *L)
{
  GenomeNode **parent, **child;
  parent = check_genome_node(L, 1);
  child  = check_genome_node(L, 2);
  genome_node_is_part_of_genome_node(*parent, genome_node_rec_ref(*child));
  return 0;
}

static int genome_node_lua_mark(lua_State *L)
{
  GenomeNode **gn = check_genome_node(L, 1);
  genome_node_mark(*gn);
  return 0;
}

static int genome_node_lua_is_marked(lua_State *L)
{
  GenomeNode **gn = check_genome_node(L, 1);
  lua_pushboolean(L, genome_node_is_marked(*gn));
  return 1;
}

static int genome_node_lua_contains_marked(lua_State *L)
{
  GenomeNode **gn;
  gn = check_genome_node(L, 1);
  lua_pushboolean(L, genome_node_contains_marked(*gn));
  return 1;
}

static int genome_feature_lua_output_leading(lua_State *L)
{
  GenomeNode **gn;
  GenomeFeature *gf;
  gn = check_genome_node(L, 1);
  /* make sure we get a genome feature */
  gf = genome_node_cast(genome_feature_class(), *gn);
  luaL_argcheck(L, gf, 1, "not a genome feature");
  gff3_output_leading(gf, NULL);
  return 0;
}

static int genome_feature_lua_get_type(lua_State *L)
{
  GenomeNode **gn;
  GenomeFeature *gf;
  gn = check_genome_node(L, 1);
  /* make sure we get a genome feature */
  gf = genome_node_cast(genome_feature_class(), *gn);
  luaL_argcheck(L, gf, 1, "not a genome feature");
  lua_pushstring(L, genome_feature_type_get_cstr(genome_feature_get_type(gf)));
  return 1;
}

static int genome_feature_lua_extract_sequence(lua_State *L)
{
  GenomeNode **gn;
  GenomeFeature *gf;
  const char *typestr;
  GenomeFeatureType type;
  bool join;
  RegionMapping **region_mapping;
  Str *sequence;
  Error *err;
  gn = check_genome_node(L, 1);
  /* make sure we get a genome feature */
  gf = genome_node_cast(genome_feature_class(), *gn);
  luaL_argcheck(L, gf, 1, "not a genome feature");
  typestr = lua_tostring(L, 2);
  luaL_argcheck(L, !genome_feature_type_get(&type, typestr), 2,
                "not a valid type");
  join = lua_toboolean(L, 3);
  region_mapping = check_region_mapping(L, 4);
  err = error_new();
  sequence = str_new();
  if (extract_feat_sequence(sequence, *gn, type, join, *region_mapping, err)) {
    str_delete(sequence);
    return lua_gt_error(L, err);
  }
  if (str_length(sequence))
    lua_pushstring(L, str_get(sequence));
  else
    lua_pushnil(L);
  str_delete(sequence);
  error_delete(err);
  return 1;
}

static int genome_node_lua_delete(lua_State *L)
{
  GenomeNode **gn;
  gn = check_genome_node(L, 1);
  genome_node_rec_delete(*gn);
  return 0;
}

static const struct luaL_Reg genome_node_lib_f [] = {
  { "genome_feature_new", genome_feature_lua_new },
  { "sequence_region_new", sequence_region_lua_new },
  { NULL, NULL }
};

static const struct luaL_Reg genome_node_lib_m [] = {
  { "get_filename", genome_node_lua_get_filename },
  { "get_range", genome_node_lua_get_range },
  { "get_seqid", genome_node_lua_get_seqid },
  { "set_seqid", genome_node_lua_set_seqid },
  { "get_strand", genome_feature_lua_get_strand },
  { "get_source", genome_feature_lua_get_source },
  { "set_source", genome_feature_lua_set_source },
  { "accept", genome_node_lua_accept },
  { "is_part_of_genome_node", genome_node_lua_is_part_of_genome_node },
  { "mark", genome_node_lua_mark },
  { "is_marked", genome_node_lua_is_marked },
  { "contains_marked", genome_node_lua_contains_marked },
  { "output_leading", genome_feature_lua_output_leading },
  { "get_type", genome_feature_lua_get_type },
  { "extract_sequence", genome_feature_lua_extract_sequence },
  { NULL, NULL }
};

int luaopen_genome_node(lua_State *L)
{
  assert(L);
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
  lua_export_metatable(L, GENOME_NODE_METATABLE);
  luaL_register(L, "gt", genome_node_lib_f);
  return 1;
}

void genome_node_lua_push(lua_State *L, GenomeNode *gn)
{
  GenomeNode **gn_lua;
  assert(L && gn);
  gn_lua = lua_newuserdata(L, sizeof (GenomeNode**));
  *gn_lua = gn;
  luaL_getmetatable(L, GENOME_NODE_METATABLE);
  lua_setmetatable(L, -2);
}
