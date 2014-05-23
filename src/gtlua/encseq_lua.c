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

#include <string.h>
#include "lauxlib.h"
#include "core/error.h"
#include "core/ma.h"
#include "core/str_array.h"
#include "core/unused_api.h"
#include "extended/luahelper.h"
#include "gtlua/alphabet_lua.h"
#include "gtlua/encseq_lua.h"
#include "gtlua/gtcore_lua.h"

#define ENCSEQ_BUFFER_METATABLE  "GenomeTools.encseq_buffer"

/* GtEncseqReader */

typedef struct {
  unsigned char *buf;
  GtUword length;
} GtEncseqExtractedBuffer;

static void gt_lua_encseq_reader_push(lua_State *L, GtEncseqReader *reader)
{
  GtEncseqReader **readerptr;
  readerptr = lua_newuserdata(L, sizeof (GtEncseqReader*));
  *readerptr = reader;
  luaL_getmetatable(L, ENCSEQ_READER_METATABLE);
  lua_setmetatable(L, -2);
}

static int encseq_reader_lua_next_encoded_char(lua_State *L)
{
  GtEncseqReader **reader;
  unsigned char cc;
  reader = check_encseq_reader(L, 1);
  cc = gt_encseq_reader_next_encoded_char(*reader);
  lua_pushnumber(L, cc);
  return 1;
}

static int encseq_reader_lua_next_decoded_char(lua_State *L)
{
  GtEncseqReader **reader;
  char cc;
  reader = check_encseq_reader(L, 1);
  cc = gt_encseq_reader_next_decoded_char(*reader);
  lua_pushlstring(L, &cc, sizeof (char));
  return 1;
}

static int encseq_reader_lua_delete(lua_State *L)
{
  GtEncseqReader **reader;
  reader = check_encseq_reader(L, 1);
  gt_encseq_reader_delete(*reader);
  return 0;
}

static int encseq_reader_lua_reinit_with_readmode(lua_State *L)
{
  GtEncseq **encseq;
  GtEncseqReader **reader;
  GtUword startpos;
  GtReadmode readmode;
  reader = check_encseq_reader(L, 1);
  encseq = check_encseq(L, 2);
  readmode = luaL_checklong(L, 3);
  startpos = luaL_checklong(L, 4);
  luaL_argcheck(L, readmode <= 3, 3,
                "invalid readmode value, must be <= 3");
  luaL_argcheck(L, startpos < gt_encseq_total_length(*encseq), 4,
                "cannot exceed total length of encoded sequence");
  gt_encseq_reader_reinit_with_readmode(*reader, *encseq, readmode, startpos);
  return 0;
}

/* GtEncseq */

static int encseq_lua_total_length(lua_State *L)
{
  GtEncseq **encseq;
  encseq = check_encseq(L, 1);
  lua_pushnumber(L, gt_encseq_total_length(*encseq));
  return 1;
}

static int encseq_lua_num_of_sequences(lua_State *L)
{
  GtEncseq **encseq;
  encseq = check_encseq(L, 1);
  lua_pushnumber(L, gt_encseq_num_of_sequences(*encseq));
  return 1;
}

static int encseq_lua_get_encoded_char(lua_State *L)
{
  GtEncseq **encseq;
  GtUword pos;
  int readmode;
  unsigned char cc;
  encseq = check_encseq(L, 1);
  pos = luaL_checklong(L, 2);
  readmode = luaL_checklong(L, 3);
  luaL_argcheck(L, pos < gt_encseq_total_length(*encseq), 2,
                "cannot exceed total length of encoded sequence");
  luaL_argcheck(L, readmode <= 3, 3,
                "invalid readmode value, must be <= 3");
  cc = gt_encseq_get_encoded_char(*encseq, pos, readmode);
  lua_pushnumber(L, cc);
  return 1;
}

static int encseq_lua_get_decoded_char(lua_State *L)
{
  GtEncseq **encseq;
  GtUword pos;
  int readmode;
  char cc;
  encseq = check_encseq(L, 1);
  pos = luaL_checklong(L, 2);
  readmode = luaL_checklong(L, 3);
  luaL_argcheck(L, pos < gt_encseq_total_length(*encseq), 2,
                "cannot exceed total length of encoded sequence");
  luaL_argcheck(L, readmode <= 3, 3,
                "invalid readmode value, must be <= 3");
  cc = gt_encseq_get_decoded_char(*encseq, pos, readmode);
  lua_pushlstring(L, &cc, sizeof (char));
  return 1;
}

static int encseq_lua_push_buffer(lua_State *L, unsigned char *arr,
                                  GtUword len)
{
  GtEncseqExtractedBuffer** buf;
  buf = lua_newuserdata(L, sizeof (GtEncseqExtractedBuffer*));
  gt_assert(buf);
  *buf = gt_malloc(sizeof (GtEncseqExtractedBuffer));
  gt_assert(*buf);
  (*buf)->buf = arr;
  (*buf)->length = len;
  luaL_getmetatable(L, ENCSEQ_BUFFER_METATABLE);
  lua_setmetatable(L, -2);
  return 1;
}

static int encseq_lua_index_buffer(lua_State *L)
{
  GtEncseqExtractedBuffer **buf = luaL_checkudata(L, 1,
                                                  ENCSEQ_BUFFER_METATABLE);
  GtUword num = luaL_checknumber(L, 2);
  luaL_argcheck(L, num <= (*buf)->length, 2,
                "must be inside extracted substring");
  lua_pushinteger(L, ((*buf)->buf)[num-1]);
  return 1;
}

static int encseq_lua_delete_buffer(lua_State *L)
{
  GtEncseqExtractedBuffer **buf = luaL_checkudata(L, 1,
                                                  ENCSEQ_BUFFER_METATABLE);
  gt_free((*buf)->buf);
  gt_free(*buf);
  return 0;
}

static int encseq_lua_extract_encoded(lua_State *L)
{
  GtEncseq **encseq;
  GtUword from, to;
  unsigned char *string;
  encseq = check_encseq(L, 1);
  from = luaL_checklong(L, 2);
  to = luaL_checklong(L, 3);
  luaL_argcheck(L, from <= to, 2, "must be <= range endposition");
  luaL_argcheck(L, to < gt_encseq_total_length(*encseq), 3,
                "cannot exceed total length of encoded sequence");
  string = gt_malloc((to - from + 1) * sizeof (unsigned char));
  gt_encseq_extract_encoded(*encseq, string, from, to);
  encseq_lua_push_buffer(L, string, (to - from + 1));
  return 1;
}

static int encseq_lua_extract_decoded(lua_State *L)
{
  GtEncseq **encseq;
  GtUword from, to;
  char *string;
  encseq = check_encseq(L, 1);
  from = luaL_checklong(L, 2);
  to = luaL_checklong(L, 3);
  luaL_argcheck(L, from <= to, 2, "must be <= range endposition");
  luaL_argcheck(L, to < gt_encseq_total_length(*encseq), 3,
                "cannot exceed total length of encoded sequence");
  string = gt_malloc((to - from + 1) * sizeof (char));
  gt_encseq_extract_decoded(*encseq, string, from, to);
  lua_pushlstring(L, string, (to - from + 1));
  gt_free(string);
  return 1;
}

static int encseq_lua_seqlength(lua_State *L)
{
  GtEncseq **encseq;
  GtUword pos;
  encseq = check_encseq(L, 1);
  pos = luaL_checklong(L, 2);
  luaL_argcheck(L, pos < gt_encseq_num_of_sequences(*encseq), 2,
                "cannot exceed number of sequences");
  lua_pushnumber(L, gt_encseq_seqlength(*encseq, pos));
  return 1;
}

static int encseq_lua_seqstartpos(lua_State *L)
{
  GtEncseq **encseq;
  GtUword pos;
  encseq = check_encseq(L, 1);
  pos = luaL_checklong(L, 2);
  luaL_argcheck(L, pos < gt_encseq_num_of_sequences(*encseq), 2,
                "cannot exceed number of sequences");
  lua_pushnumber(L, gt_encseq_seqstartpos(*encseq, pos));
  return 1;
}

static int encseq_lua_seqnum(lua_State *L)
{
  GtEncseq **encseq;
  GtUword pos;
  encseq = check_encseq(L, 1);
  pos = luaL_checklong(L, 2);
  luaL_argcheck(L, pos < gt_encseq_total_length(*encseq), 2,
                "cannot exceed total length of encoded sequence");
  lua_pushnumber(L, gt_encseq_seqnum(*encseq, pos));
  return 1;
}

static int encseq_lua_has_multiseq_support(lua_State *L)
{
  GtEncseq **encseq;
  encseq = check_encseq(L, 1);
  lua_pushboolean(L, gt_encseq_has_multiseq_support(*encseq));
  return 1;
}

static int encseq_lua_has_description_support(lua_State *L)
{
  GtEncseq **encseq;
  encseq = check_encseq(L, 1);
  lua_pushboolean(L, gt_encseq_has_description_support(*encseq));
  return 1;
}

static int encseq_lua_description(lua_State *L)
{
  GtEncseq **encseq;
  GtUword seqno, desclen;
  const char *string;
  encseq = check_encseq(L, 1);
  seqno = luaL_checklong(L, 2);
  luaL_argcheck(L, seqno < gt_encseq_num_of_sequences(*encseq), 2,
                "cannot exceed number of sequences");
  string = gt_encseq_description(*encseq, &desclen, seqno);
  lua_pushlstring(L, string, desclen);
  return 1;
}

static int encseq_lua_num_of_files(lua_State *L)
{
  GtEncseq **encseq;
  encseq = check_encseq(L, 1);
  lua_pushnumber(L, gt_encseq_num_of_files(*encseq));
  return 1;
}

static int encseq_lua_effective_filelength(lua_State *L)
{
  GtEncseq **encseq;
  GtUword fileno;
  encseq = check_encseq(L, 1);
  fileno = luaL_checklong(L, 2);
  luaL_argcheck(L, fileno < gt_encseq_num_of_files(*encseq), 2,
                "cannot exceed number of files");
  lua_pushnumber(L, gt_encseq_effective_filelength(*encseq, fileno));
  return 1;
}

static int encseq_lua_filestartpos(lua_State *L)
{
  GtEncseq **encseq;
  GtUword fileno;
  encseq = check_encseq(L, 1);
  fileno = luaL_checklong(L, 2);
  luaL_argcheck(L, fileno < gt_encseq_num_of_files(*encseq), 2,
                "cannot exceed number of files");
  lua_pushnumber(L, gt_encseq_filestartpos(*encseq, fileno));
  return 1;
}

static int encseq_lua_filenum(lua_State *L)
{
  GtEncseq **encseq;
  GtUword pos;
  encseq = check_encseq(L, 1);
  pos = luaL_checklong(L, 2);
  luaL_argcheck(L, pos < gt_encseq_total_length(*encseq), 2,
                "cannot exceed total length of encoded sequence");
  lua_pushnumber(L, gt_encseq_filenum(*encseq, pos));
  return 1;
}

static int encseq_lua_filenames(lua_State *L)
{
  GtEncseq **encseq;
  const GtStrArray *filenames;
  GtUword i;
  encseq = check_encseq(L, 1);
  filenames = gt_encseq_filenames(*encseq);
  lua_newtable(L);
  for (i = 0; i < gt_str_array_size(filenames); i++) {
    lua_pushinteger(L, i+1); /* in Lua we index from 1 on */
    lua_pushstring(L, gt_str_array_get(filenames, i));
    lua_rawset(L, -3);
  }
  return 1;
}

static int encseq_lua_alphabet(lua_State *L)
{
  GtEncseq **encseq;
  GtAlphabet *alpha;
  encseq = check_encseq(L, 1);
  gt_assert(*encseq);
  alpha = gt_alphabet_ref(gt_encseq_alphabet(*encseq));
  gt_lua_alphabet_push(L, alpha);
  return 1;
}

static int encseq_lua_mirror(lua_State *L)
{
  GtEncseq **encseq;
  GtError *err = gt_error_new();
  encseq = check_encseq(L, 1);
  gt_assert(*encseq);
  luaL_argcheck(L, !gt_encseq_is_mirrored(*encseq), 1, "is already mirrored");
  if (gt_encseq_mirror(*encseq, err) != 0)
    gt_lua_error(L, err);
  gt_error_delete(err);
  return 0;
}

static int encseq_lua_unmirror(lua_State *L)
{
  GtEncseq **encseq;
  encseq = check_encseq(L, 1);
  gt_assert(*encseq);
  luaL_argcheck(L, gt_encseq_is_mirrored(*encseq), 1, "is not mirrored");
  gt_encseq_unmirror(*encseq);
  return 0;
}

static int encseq_lua_is_mirrored(lua_State *L)
{
  GtEncseq **encseq;
  encseq = check_encseq(L, 1);
  gt_assert(*encseq);
  lua_pushboolean(L, gt_encseq_is_mirrored(*encseq));
  return 1;
}

static int encseq_lua_create_reader_with_readmode(lua_State *L)
{
  GtEncseq **encseq;
  GtEncseqReader *reader;
  GtUword startpos;
  GtReadmode readmode;
  encseq = check_encseq(L, 1);
  readmode = luaL_checklong(L, 2);
  startpos = luaL_checklong(L, 3);
  luaL_argcheck(L, readmode <= 3, 2,
                "invalid readmode value, must be <= 3");
  luaL_argcheck(L, startpos < gt_encseq_total_length(*encseq), 3,
                "cannot exceed total length of encoded sequence");
  reader = gt_encseq_create_reader_with_readmode(*encseq, readmode, startpos);
  gt_assert(reader);
  gt_lua_encseq_reader_push(L, reader);
  return 1;
}

static int encseq_lua_delete(lua_State *L)
{
  GtEncseq **encseq;
  encseq = check_encseq(L, 1);
  gt_encseq_delete(*encseq);
  return 0;
}

/* GtEncseqEncoder */

static int encseq_encoder_lua_new(lua_State *L)
{
  GtEncseqEncoder **encoder;
  encoder = lua_newuserdata(L, sizeof (GtEncseqEncoder*));
  gt_assert(encoder);
  *encoder = gt_encseq_encoder_new();
  gt_assert(*encoder);
  luaL_getmetatable(L, ENCSEQ_ENCODER_METATABLE);
  lua_setmetatable(L, -2);
  return 1;
}

static int encseq_encoder_lua_use_representation(lua_State *L)
{
  GtEncseqEncoder **encoder;
  const char *repr;
  GtError *err;
  encoder = check_encseq_encoder(L, 1);
  repr = luaL_checkstring(L, 2);
  gt_assert(*encoder);
  err = gt_error_new();
  if (gt_encseq_encoder_use_representation(*encoder, repr, err) != 0)
    gt_lua_error(L, err);
  gt_error_delete(err);
  return 0;
}

static int encseq_encoder_lua_representation(lua_State *L)
{
  GtEncseqEncoder **encoder;
  encoder = check_encseq_encoder(L, 1);
  gt_assert(*encoder);
  lua_pushstring(L, gt_str_get(gt_encseq_encoder_representation(*encoder)));
  return 1;
}

static int encseq_encoder_lua_use_symbolmap_file(lua_State *L)
{
  GtEncseqEncoder **encoder;
  const char *smapfile;
  GtError *err;
  encoder = check_encseq_encoder(L, 1);
  smapfile = luaL_checkstring(L, 2);
  gt_assert(*encoder);
  err = gt_error_new();
  if (gt_encseq_encoder_use_symbolmap_file(*encoder, smapfile, err) != 0)
    gt_lua_error(L, err);
  gt_error_delete(err);
  return 0;
}

static int encseq_encoder_lua_symbolmap_file(lua_State *L)
{
  GtEncseqEncoder **encoder;
  encoder = check_encseq_encoder(L, 1);
  gt_assert(*encoder);
  lua_pushstring(L, gt_encseq_encoder_symbolmap_file(*encoder));
  return 1;
}

static int encseq_encoder_lua_enable_description_support(lua_State *L)
{
  GtEncseqEncoder **encoder;
  encoder = check_encseq_encoder(L, 1);
  gt_assert(*encoder);
  gt_encseq_encoder_enable_description_support(*encoder);
  return 0;
}

static int encseq_encoder_lua_disable_description_support(lua_State *L)
{
  GtEncseqEncoder **encoder;
  encoder = check_encseq_encoder(L, 1);
  gt_assert(*encoder);
  gt_encseq_encoder_disable_description_support(*encoder);
  return 0;
}

static int encseq_encoder_lua_enable_multiseq_support(lua_State *L)
{
  GtEncseqEncoder **encoder;
  encoder = check_encseq_encoder(L, 1);
  gt_assert(*encoder);
  gt_encseq_encoder_enable_multiseq_support(*encoder);
  return 0;
}

static int encseq_encoder_lua_disable_multiseq_support(lua_State *L)
{
  GtEncseqEncoder **encoder;
  encoder = check_encseq_encoder(L, 1);
  gt_assert(*encoder);
  gt_encseq_encoder_disable_multiseq_support(*encoder);
  return 0;
}

static int encseq_encoder_lua_enable_lossless_support(lua_State *L)
{
  GtEncseqEncoder **encoder;
  encoder = check_encseq_encoder(L, 1);
  gt_assert(*encoder);
  gt_encseq_encoder_enable_lossless_support(*encoder);
  return 0;
}

static int encseq_encoder_lua_disable_lossless_support(lua_State *L)
{
  GtEncseqEncoder **encoder;
  encoder = check_encseq_encoder(L, 1);
  gt_assert(*encoder);
  gt_encseq_encoder_disable_lossless_support(*encoder);
  return 0;
}

static int encseq_encoder_lua_set_input_dna(lua_State *L)
{
  GtEncseqEncoder **encoder;
  encoder = check_encseq_encoder(L, 1);
  gt_assert(*encoder);
  gt_encseq_encoder_set_input_dna(*encoder);
  return 0;
}

static int encseq_encoder_lua_is_input_dna(lua_State *L)
{
  GtEncseqEncoder **encoder;
  encoder = check_encseq_encoder(L, 1);
  gt_assert(*encoder);
  lua_pushboolean(L, gt_encseq_encoder_is_input_dna(*encoder));
  return 1;
}

static int encseq_encoder_lua_set_input_protein(lua_State *L)
{
  GtEncseqEncoder **encoder;
  encoder = check_encseq_encoder(L, 1);
  gt_assert(*encoder);
  gt_encseq_encoder_set_input_protein(*encoder);
  return 0;
}

static int encseq_encoder_lua_is_input_protein(lua_State *L)
{
  GtEncseqEncoder **encoder;
  encoder = check_encseq_encoder(L, 1);
  gt_assert(*encoder);
  lua_pushboolean(L, gt_encseq_encoder_is_input_protein(*encoder));
  return 1;
}

static int encseq_encoder_lua_encode(lua_State *L)
{
  GtEncseqEncoder **encoder;
  GtStrArray *seqfiles = gt_str_array_new();
  GtError *err;
  const char *indexname;
  encoder = check_encseq_encoder(L, 1);
  err = gt_error_new();
  if (gt_lua_get_table_as_strarray(L, 2, seqfiles, err) != 0) {
    gt_str_array_delete(seqfiles);
    gt_lua_error(L, err);
  }
  indexname = luaL_checkstring(L, 3);
  gt_assert(*encoder);
  if (gt_encseq_encoder_encode(*encoder, seqfiles, indexname, err) != 0) {
    gt_str_array_delete(seqfiles);
    gt_lua_error(L, err);
  }
  gt_str_array_delete(seqfiles);
  gt_error_delete(err);
  return 0;
}

static int encseq_encoder_lua_delete(lua_State *L)
{
  GtEncseqEncoder **encoder;
  encoder = check_encseq_encoder(L, 1);
  gt_encseq_encoder_delete(*encoder);
  return 0;
}

/* GtEncseqLoader */

static int encseq_loader_lua_new(lua_State *L)
{
  GtEncseqLoader **loader;
  loader = lua_newuserdata(L, sizeof (GtEncseqLoader*));
  gt_assert(loader);
  *loader = gt_encseq_loader_new();
  gt_assert(*loader);
  luaL_getmetatable(L, ENCSEQ_LOADER_METATABLE);
  lua_setmetatable(L, -2);
  return 1;
}

static int encseq_loader_lua_enable_autosupport(lua_State *L)
{
  GtEncseqLoader **loader;
  loader = check_encseq_loader(L, 1);
  gt_assert(*loader);
  gt_encseq_loader_enable_autosupport(*loader);
  return 0;
}

static int encseq_loader_lua_disable_autosupport(lua_State *L)
{
  GtEncseqLoader **loader;
  loader = check_encseq_loader(L, 1);
  gt_assert(*loader);
  gt_encseq_loader_disable_autosupport(*loader);
  return 0;
}

static int encseq_loader_lua_require_description_support(lua_State *L)
{
  GtEncseqLoader **loader;
  loader = check_encseq_loader(L, 1);
  gt_assert(*loader);
  gt_encseq_loader_require_description_support(*loader);
  return 0;
}

static int encseq_loader_lua_drop_description_support(lua_State *L)
{
  GtEncseqLoader **loader;
  loader = check_encseq_loader(L, 1);
  gt_assert(*loader);
  gt_encseq_loader_drop_description_support(*loader);
  return 0;
}

static int encseq_loader_lua_require_multiseq_support(lua_State *L)
{
  GtEncseqLoader **loader;
  loader = check_encseq_loader(L, 1);
  gt_assert(*loader);
  gt_encseq_loader_require_multiseq_support(*loader);
  return 0;
}

static int encseq_loader_lua_drop_multiseq_support(lua_State *L)
{
  GtEncseqLoader **loader;
  loader = check_encseq_loader(L, 1);
  gt_assert(*loader);
  gt_encseq_loader_drop_multiseq_support(*loader);
  return 0;
}

static int encseq_loader_lua_require_lossless_support(lua_State *L)
{
  GtEncseqLoader **loader;
  loader = check_encseq_loader(L, 1);
  gt_assert(*loader);
  gt_encseq_loader_require_lossless_support(*loader);
  return 0;
}

static int encseq_loader_lua_drop_lossless_support(lua_State *L)
{
  GtEncseqLoader **loader;
  loader = check_encseq_loader(L, 1);
  gt_assert(*loader);
  gt_encseq_loader_drop_lossless_support(*loader);
  return 0;
}

static int encseq_loader_lua_mirror(lua_State *L)
{
  GtEncseqLoader **loader;
  loader = check_encseq_loader(L, 1);
  gt_assert(*loader);
  gt_encseq_loader_mirror(*loader);
  return 0;
}

static int encseq_loader_lua_do_not_mirror(lua_State *L)
{
  GtEncseqLoader **loader;
  loader = check_encseq_loader(L, 1);
  gt_assert(*loader);
  gt_encseq_loader_do_not_mirror(*loader);
  return 0;
}

static int encseq_loader_lua_load(lua_State *L)
{
  GtEncseqLoader **loader;
  const char *indexname;
  GtError *err = gt_error_new();
  GtEncseq *encseq;
  loader = check_encseq_loader(L, 1);
  indexname = luaL_checkstring(L, 2);
  gt_assert(*loader);
  encseq = gt_encseq_loader_load(*loader, indexname, err);
  if (encseq == NULL)
    gt_lua_error(L, err);
  else {
    gt_lua_encseq_push(L, encseq);
  }
  gt_error_delete(err);
  return 1;
}

static int encseq_loader_lua_delete(lua_State *L)
{
  GtEncseqLoader **loader;
  loader = check_encseq_loader(L, 1);
  gt_encseq_loader_delete(*loader);
  return 0;
}

/* GtEncseqBuilder */

static int encseq_builder_lua_new(lua_State *L)
{
  GtEncseqBuilder **builder;
  GtAlphabet **alpha;
  builder = lua_newuserdata(L, sizeof (GtEncseqBuilder*));
  alpha = check_alphabet(L, 1);
  gt_assert(builder && *alpha);
  *builder = gt_encseq_builder_new(*alpha);
  gt_assert(*builder);
  luaL_getmetatable(L, ENCSEQ_BUILDER_METATABLE);
  lua_setmetatable(L, -2);
  return 1;
}

static int encseq_builder_lua_enable_description_support(lua_State *L)
{
  GtEncseqBuilder **builder;
  builder = check_encseq_builder(L, 1);
  gt_assert(*builder);
  gt_encseq_builder_enable_description_support(*builder);
  return 0;
}

static int encseq_builder_lua_disable_description_support(lua_State *L)
{
  GtEncseqBuilder **builder;
  builder = check_encseq_builder(L, 1);
  gt_assert(*builder);
  gt_encseq_builder_disable_description_support(*builder);
  return 0;
}

static int encseq_builder_lua_enable_multiseq_support(lua_State *L)
{
  GtEncseqBuilder **builder;
  builder = check_encseq_builder(L, 1);
  gt_assert(*builder);
  gt_encseq_builder_enable_multiseq_support(*builder);
  return 0;
}

static int encseq_builder_lua_disable_multiseq_support(lua_State *L)
{
  GtEncseqBuilder **builder;
  builder = check_encseq_builder(L, 1);
  gt_assert(*builder);
  gt_encseq_builder_disable_multiseq_support(*builder);
  return 0;
}

static int encseq_builder_lua_add_str(lua_State *L)
{
  GtEncseqBuilder **builder;
  const char *str, *desc;
  builder = check_encseq_builder(L, 1);
  str = luaL_checkstring(L, 2);
  if (lua_isnil(L, 3))
    desc = "";
  else
    desc = luaL_checkstring(L, 3);
  gt_assert(*builder);
  gt_encseq_builder_add_cstr(*builder, str, strlen(str), desc);
  return 0;
}

static inline int gt_lua_get_table_as_uchararray(lua_State *L, int index,
                                                 unsigned char **outarray,
                                                 GtUword *arrlen,
                                                 GtError *err)
{
  int had_err = 0, val;
  GtUword i;
  unsigned char *arr;
  if (!lua_istable(L, index)) {
    gt_error_set(err, "argument is not a table");
    return -1;
  }
  *arrlen = lua_objlen(L, index);
  arr = gt_malloc(*arrlen * sizeof (unsigned char));
  for (i = 1; i <= *arrlen; i++)
  {
    lua_rawgeti(L, index, i);
    if (!lua_isnumber(L, -1)) {
      had_err = -1;
      gt_error_set(err, "input contains non-numeric value");
      break;
    }
    if (!had_err && (val = lua_tointeger(L, -1)) > UCHAR_MAX) {
      had_err = -1;
      gt_error_set(err, "input contains oversized encoded value");
    }
    if (!had_err)
      arr[i-1] = (unsigned char) val;
    lua_pop(L, 1);
  }
  *outarray = arr;
  return had_err;
}

static int encseq_builder_lua_add_encoded(lua_State *L)
{
  GtEncseqBuilder **builder;
  const char *desc;
  unsigned char *arr = NULL;
  GtUword arrlen = 0;
  GtError *err;
  builder = check_encseq_builder(L, 1);
  err = gt_error_new();
  if (gt_lua_get_table_as_uchararray(L, 2, &arr, &arrlen, err) != 0) {
    gt_lua_error(L, err);
  }
  if (lua_isnil(L, 3))
    desc = "";
  else
    desc = luaL_checkstring(L, 3);
  gt_assert(*builder);
  gt_encseq_builder_add_encoded_own(*builder, arr, arrlen, desc);
  gt_free(arr);
  gt_error_delete(err);
  return 0;
}

static int encseq_builder_lua_reset(lua_State *L)
{
  GtEncseqBuilder **builder;
  builder = check_encseq_builder(L, 1);
  gt_encseq_builder_reset(*builder);
  return 0;
}

static int encseq_builder_lua_build(lua_State *L)
{
  GtEncseqBuilder **builder;
  GtEncseq *encseq;
  GtError *err;
  builder = check_encseq_builder(L, 1);
  err = gt_error_new();
  encseq = gt_encseq_builder_build(*builder, err);
  if (encseq == NULL) {
    gt_lua_error(L, err);
  } else
    gt_lua_encseq_push(L, encseq);
  gt_error_delete(err);
  return 1;
}

static int encseq_builder_lua_delete(lua_State *L)
{
  GtEncseqBuilder **builder;
  builder = check_encseq_builder(L, 1);
  gt_encseq_builder_delete(*builder);
  return 0;
}

static const struct luaL_Reg encseq_lib_f [] = {
  { "encseq_encoder_new", encseq_encoder_lua_new },
  { "encseq_loader_new", encseq_loader_lua_new },
  { "encseq_builder_new", encseq_builder_lua_new },
  { NULL, NULL }
};

static const struct luaL_Reg encseq_reader_lib_m [] = {
  { "next_encoded_char", encseq_reader_lua_next_encoded_char },
  { "next_decoded_char", encseq_reader_lua_next_decoded_char },
  { "reinit_with_readmode", encseq_reader_lua_reinit_with_readmode },
  { NULL, NULL }
};

static const struct luaL_Reg encseq_encoder_lib_m [] = {
  { "use_representation", encseq_encoder_lua_use_representation },
  { "representation", encseq_encoder_lua_representation },
  { "use_symbolmap_file", encseq_encoder_lua_use_symbolmap_file },
  { "symbolmap_file", encseq_encoder_lua_symbolmap_file },
  { "enable_description_support",
                                encseq_encoder_lua_enable_description_support },
  { "disable_description_support",
                               encseq_encoder_lua_disable_description_support },
  { "enable_multiseq_support", encseq_encoder_lua_enable_multiseq_support },
  { "disable_multiseq_support", encseq_encoder_lua_disable_multiseq_support },
  { "enable_lossless_support", encseq_encoder_lua_enable_lossless_support },
  { "disable_lossless_support", encseq_encoder_lua_disable_lossless_support },
  { "set_input_dna", encseq_encoder_lua_set_input_dna },
  { "is_input_dna", encseq_encoder_lua_is_input_dna },
  { "set_input_protein", encseq_encoder_lua_set_input_protein },
  { "is_input_protein", encseq_encoder_lua_is_input_protein },
  { "encode", encseq_encoder_lua_encode },
  { NULL, NULL }
};

static const struct luaL_Reg encseq_loader_lib_m [] = {
  { "enable_autosupport", encseq_loader_lua_enable_autosupport },
  { "disable_autosupport", encseq_loader_lua_disable_autosupport },
  { "require_description_support",
                                encseq_loader_lua_require_description_support },
  { "drop_description_support", encseq_loader_lua_drop_description_support },
  { "require_multiseq_support", encseq_loader_lua_require_multiseq_support },
  { "drop_multiseq_support", encseq_loader_lua_drop_multiseq_support },
  { "require_lossless_support", encseq_loader_lua_require_lossless_support },
  { "drop_lossless_support", encseq_loader_lua_drop_lossless_support },
  { "mirror", encseq_loader_lua_mirror },
  { "do_not_mirror", encseq_loader_lua_do_not_mirror },
  { "load", encseq_loader_lua_load },
  { NULL, NULL }
};

static const struct luaL_Reg encseq_builder_lib_m [] = {
  { "enable_description_support",
                                encseq_builder_lua_enable_description_support },
  { "disable_description_support",
                               encseq_builder_lua_disable_description_support },
  { "enable_multiseq_support", encseq_builder_lua_enable_multiseq_support },
  { "disable_multiseq_support", encseq_builder_lua_disable_multiseq_support },
  { "add_string", encseq_builder_lua_add_str },
  { "add_encoded", encseq_builder_lua_add_encoded },
  { "build", encseq_builder_lua_build },
  { "reset", encseq_builder_lua_reset },
  { NULL, NULL }
};

static const struct luaL_Reg encseq_lib_m [] = {
  { "total_length", encseq_lua_total_length },
  { "num_of_sequences", encseq_lua_num_of_sequences },
  { "get_encoded_char", encseq_lua_get_encoded_char },
  { "get_decoded_char", encseq_lua_get_decoded_char },
  { "create_reader_with_readmode", encseq_lua_create_reader_with_readmode },
  { "extract_encoded", encseq_lua_extract_encoded },
  { "extract_decoded", encseq_lua_extract_decoded },
  { "seqlength", encseq_lua_seqlength },
  { "seqstartpos", encseq_lua_seqstartpos },
  { "seqnum", encseq_lua_seqnum },
  { "has_multiseq_support", encseq_lua_has_multiseq_support },
  { "has_description_support", encseq_lua_has_description_support },
  { "description", encseq_lua_description },
  { "filenames", encseq_lua_filenames },
  { "num_of_files", encseq_lua_num_of_files },
  { "effective_filelength", encseq_lua_effective_filelength },
  { "filestartpos", encseq_lua_filestartpos },
  { "filenum", encseq_lua_filenum },
  { "alphabet", encseq_lua_alphabet },
  { "mirror", encseq_lua_mirror },
  { "unmirror", encseq_lua_unmirror },
  { "is_mirrored", encseq_lua_is_mirrored },
  { NULL, NULL }
};

int gt_lua_open_encseq(lua_State *L)
{
#ifndef NDEBUG
  int stack_size;
#endif
  gt_assert(L);
#ifndef NDEBUG
  stack_size = lua_gettop(L);
#endif
  luaL_newmetatable(L, ENCSEQ_METATABLE);
  lua_pushvalue(L, -1); /* duplicate the metatable */
  lua_setfield(L, -2, "__index");
  lua_pushstring(L, "__gc");
  lua_pushcfunction(L, encseq_lua_delete);
  lua_settable(L, -3);
  luaL_register(L, NULL, encseq_lib_m);
  gt_lua_export_metatable(L, ENCSEQ_METATABLE);

  luaL_newmetatable(L, ENCSEQ_BUFFER_METATABLE);
  lua_pushcfunction(L, encseq_lua_index_buffer);
  lua_setfield(L, -2, "__index");
  lua_pushstring(L, "__gc");
  lua_pushcfunction(L, encseq_lua_delete_buffer);
  lua_settable(L, -3);
  lua_pop(L, 1);

  luaL_newmetatable(L, ENCSEQ_READER_METATABLE);
  lua_pushvalue(L, -1); /* duplicate the metatable */
  lua_setfield(L, -2, "__index");
  lua_pushstring(L, "__gc");
  lua_pushcfunction(L, encseq_reader_lua_delete);
  lua_settable(L, -3);
  luaL_register(L, NULL, encseq_reader_lib_m);
  gt_lua_export_metatable(L, ENCSEQ_READER_METATABLE);

  luaL_newmetatable(L, ENCSEQ_ENCODER_METATABLE);
  lua_pushvalue(L, -1); /* duplicate the metatable */
  lua_setfield(L, -2, "__index");
  lua_pushstring(L, "__gc");
  lua_pushcfunction(L, encseq_encoder_lua_delete);
  lua_settable(L, -3);
  luaL_register(L, NULL, encseq_encoder_lib_m);
  gt_lua_export_metatable(L, ENCSEQ_ENCODER_METATABLE);

  luaL_newmetatable(L, ENCSEQ_LOADER_METATABLE);
  lua_pushvalue(L, -1); /* duplicate the metatable */
  lua_setfield(L, -2, "__index");
  lua_pushstring(L, "__gc");
  lua_pushcfunction(L, encseq_loader_lua_delete);
  lua_settable(L, -3);
  luaL_register(L, NULL, encseq_loader_lib_m);
  gt_lua_export_metatable(L, ENCSEQ_LOADER_METATABLE);

  luaL_newmetatable(L, ENCSEQ_BUILDER_METATABLE);
  lua_pushvalue(L, -1); /* duplicate the metatable */
  lua_setfield(L, -2, "__index");
  lua_pushstring(L, "__gc");
  lua_pushcfunction(L, encseq_builder_lua_delete);
  lua_settable(L, -3);
  luaL_register(L, NULL, encseq_builder_lib_m);
  gt_lua_export_metatable(L, ENCSEQ_BUILDER_METATABLE);

  luaL_register(L, "gt", encseq_lib_f);
  lua_pop(L, 1);
  gt_assert(lua_gettop(L) == stack_size);
  return 1;
}

void gt_lua_encseq_push(lua_State *L, GtEncseq *encseq)
{
  GtEncseq **encseqptr;
  gt_assert(L);
  encseqptr = lua_newuserdata(L, sizeof (GtEncseq*));
  *encseqptr = encseq;
  luaL_getmetatable(L, ENCSEQ_METATABLE);
  lua_setmetatable(L, -2);
}
