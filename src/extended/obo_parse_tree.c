/*
  Copyright (c) 2008 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2008 Center for Bioinformatics, University of Hamburg

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
#include "core/array.h"
#include "core/cstr.h"
#include "core/fa.h"
#include "core/hashmap.h"
#include "core/io.h"
#include "core/ma.h"
#include "extended/obo_parse_tree.h"

#define BLANK_CHAR         ' '
#define COMMENT_CHAR       '!'
#define SEPARATOR_CHAR     ':'
#define STANZA_OPEN_CHAR   '['
#define STANZA_CLOSE_CHAR  ']'
#define CARRIAGE_RETURN    '\r'
#define END_OF_LINE        '\n'
#define END_OF_FILE        EOF

typedef struct {
  char *tag,
       *value;
} OBOHeaderEntry;

typedef struct {
  GT_Array *content;
} OBOHeader;

static OBOHeader* obo_header_new(void)
{
  OBOHeader *obo_header = ma_malloc(sizeof *obo_header);
  obo_header->content = gt_array_new(sizeof (OBOHeaderEntry*));
  return obo_header;
}

static void obo_header_delete(OBOHeader *obo_header)
{
  unsigned long i;
  if (!obo_header) return;
  for (i = 0; i < gt_array_size(obo_header->content); i++) {
    OBOHeaderEntry *entry = *(OBOHeaderEntry**)
                            gt_array_get(obo_header->content, i);
    ma_free(entry->value);
    ma_free(entry->tag);
    ma_free(entry);
  }
  gt_array_delete(obo_header->content);
  ma_free(obo_header);
}

static void obo_header_add(OBOHeader *obo_header,
                           const char *tag, const char *value)
{
  OBOHeaderEntry *entry;
  assert(obo_header && tag && value);
  entry = ma_malloc(sizeof *entry);
  entry->tag = cstr_dup(tag);
  entry->value = cstr_dup(value);
  gt_array_add(obo_header->content, entry);
}

static const char* obo_header_get(OBOHeader *obo_header, const char *tag)
{
  unsigned long i;
  assert(obo_header && tag);
  for (i = 0; i < gt_array_size(obo_header->content); i++) {
    OBOHeaderEntry *entry = *(OBOHeaderEntry**)
                            gt_array_get(obo_header->content, i);
    if (!strcmp(entry->tag, tag))
      return entry->value;
  }
  return NULL;
}

static int obo_header_validate(OBOHeader *obo_header, const char *obo_file_name,
                               Error *err)
{
  error_check(err);
  assert(obo_header && obo_file_name);
  /* make sure header contains the required tag */
  if (!obo_header_get(obo_header, "format-version")) {
    error_set(err, "the header of OBO-file \"%s\" does not contain "
              "\"format-version\" tag", obo_file_name);
    return -1;
  }
  return 0;
}

typedef struct {
  char *type;
  Hashmap *content;
  unsigned long line;
  Str *filename;
} OBOStanza;

static OBOStanza* obo_stanza_new(const char *type, unsigned long line,
                                 Str *filename)
{
  OBOStanza *obo_stanza = ma_malloc(sizeof *obo_stanza);
  obo_stanza->type = cstr_dup(type);
  obo_stanza->content = hashmap_new(HASH_STRING, ma_free_func, ma_free_func);
  obo_stanza->line = line;
  obo_stanza->filename = str_ref(filename);
  return obo_stanza;
}

static void obo_stanza_delete(OBOStanza *obo_stanza)
{
  if (!obo_stanza) return;
  str_delete(obo_stanza->filename);
  hashmap_delete(obo_stanza->content);
  ma_free(obo_stanza->type);
  ma_free(obo_stanza);
}

static void obo_stanza_add(OBOStanza *obo_stanza,
                          const char *tag, const char *value)
{
  assert(obo_stanza && tag && value);
  /* XXX: currently duplicate tags are silently skipped */
  if (!hashmap_get(obo_stanza->content, tag))
    hashmap_add(obo_stanza->content, cstr_dup(tag), cstr_dup(value));
}

static const char* obo_stanza_get_type(const OBOStanza *obo_stanza)
{
  assert(obo_stanza);
  return obo_stanza->type;
}

static const char* obo_stanza_get_value(const OBOStanza *obo_stanza,
                                        const char *stanza_key)
{
  assert(obo_stanza);
  return hashmap_get(obo_stanza->content, stanza_key);
}

struct OBOParseTree {
  OBOHeader *obo_header;
  GT_Array *stanzas;
};

static void obo_parse_tree_add_stanza(OBOParseTree *obo_parse_tree,
                                      OBOStanza *obo_stanza)
{
  assert(obo_parse_tree && obo_stanza);
  gt_array_add(obo_parse_tree->stanzas, obo_stanza);
}

static int validate_value(const OBOStanza *obo_stanza, const char *value,
                          Error *err)
{
  error_check(err);
  assert(obo_stanza && value);
  if (!obo_stanza_get_value(obo_stanza, value)) {
    error_set(err, "%s stanza starting on line %lu in file \"%s\" lacks "
              "required \"%s\" tag",
              obo_stanza_get_type(obo_stanza), obo_stanza->line,
              str_get(obo_stanza->filename), value);
    return -1;
  }
  return 0;
}

static int obo_parse_tree_validate_stanzas(const OBOParseTree *obo_parse_tree,
                                            Error *err)
{
  unsigned long i;
  int had_err = 0;
  error_check(err);
  assert(obo_parse_tree);
  for (i = 0; !had_err && i < obo_parse_tree_num_of_stanzas(obo_parse_tree);
       i++) {
    OBOStanza *stanza = *(OBOStanza**) gt_array_get(obo_parse_tree->stanzas, i);
    if (!strcmp(obo_parse_tree_get_stanza_type(obo_parse_tree, i), "Term")) {
      had_err = validate_value(stanza, "id", err);
      if (!had_err)
        had_err = validate_value(stanza, "name", err);
    }
    else if (!strcmp(obo_parse_tree_get_stanza_type(obo_parse_tree, i),
                     "Typedef")) {
      had_err = validate_value(stanza, "id", err);
      if (!had_err)
        had_err = validate_value(stanza, "name", err);
    }
    else if (!strcmp(obo_parse_tree_get_stanza_type(obo_parse_tree, i),
                                                    "Instance")) {
      had_err = validate_value(stanza, "id", err);
      if (!had_err)
        had_err = validate_value(stanza, "name", err);
      if (!had_err)
        had_err = validate_value(stanza, "instance_of", err);
    }
  }
  return had_err;
}

static bool any_char(IO *obo_file, bool be_permissive)
{
  switch (io_peek(obo_file)) {
    case BLANK_CHAR:
    case SEPARATOR_CHAR:
    case STANZA_OPEN_CHAR:
    case STANZA_CLOSE_CHAR:
      if (be_permissive)
        return true;
    case COMMENT_CHAR:
    case CARRIAGE_RETURN:
    case END_OF_LINE:
    case END_OF_FILE:
      return false;
  }
  return true;
}

static bool ignored_char(IO *obo_file)
{
  char cc = io_peek(obo_file);
  if ((cc == BLANK_CHAR) || (cc == COMMENT_CHAR) || (cc == CARRIAGE_RETURN) ||
      (cc == END_OF_LINE))
    return true;
  return false;
}

static int expect(IO *obo_file, char expected_char, Error *err)
{
  char cc;
  error_check(err);
  cc = io_next(obo_file);
  if (cc != expected_char) {
    if (cc == CARRIAGE_RETURN) {
      if (io_peek(obo_file) == END_OF_LINE)
        io_next(obo_file);
      return 0;
    }
    if (expected_char == END_OF_FILE) {
      error_set(err, "file \"%s\": line %lu: expected end-of-file, got '%c'",
                io_get_filename(obo_file), io_get_line_number(obo_file), cc);
    }
    else if ((cc == CARRIAGE_RETURN) || (cc == END_OF_LINE)) {
      error_set(err, "file \"%s\": line %lu: expected character '%c', got "
                "newline", io_get_filename(obo_file),
                io_get_line_number(obo_file), expected_char);
    }
    else {
      error_set(err, "file \"%s\": line %lu: expected character '%c', got '%c'",
                io_get_filename(obo_file), io_get_line_number(obo_file),
                expected_char, cc);
    }
    return -1;
  }
  return 0;
}

static int comment_line(IO *obo_file, Error *err)
{
  int had_err;
  error_check(err);
  had_err = expect(obo_file, COMMENT_CHAR, err);
  while (!had_err) {
    switch (io_peek(obo_file)) {
      case CARRIAGE_RETURN:
        io_next(obo_file);
        if (io_peek(obo_file) == END_OF_LINE)
          io_next(obo_file);
        return had_err;
      case END_OF_LINE:
        io_next(obo_file);
      case END_OF_FILE:
        return had_err;
      default:
        io_next(obo_file);
    }
  }
  return had_err;
}

static int blank_line(IO *obo_file, Error *err)
{
  int had_err;
  error_check(err);
  had_err = expect(obo_file, BLANK_CHAR, err);
  while (!had_err) {
    char cc = io_peek(obo_file);
    if (cc == COMMENT_CHAR)
      return comment_line(obo_file, err);
    else if (cc == CARRIAGE_RETURN) {
      io_next(obo_file);
      if (io_peek(obo_file) == END_OF_LINE)
        io_next(obo_file);
      break;
    }
    else if ((cc == END_OF_LINE) || (cc == END_OF_FILE)) {
      io_next(obo_file);
      break;
    }
    else
      had_err = expect(obo_file, BLANK_CHAR, err);
  }
  return had_err;
}

static bool ignored_line(IO *obo_file, Error *err)
{
  error_check(err);
  if (io_peek(obo_file) == BLANK_CHAR)
    return blank_line(obo_file, err);
  return comment_line(obo_file, err);
}

static int proc_any_char(IO *obo_file, Str *capture, bool be_permissive,
                         Error *err)
{
  error_check(err);
  assert(obo_file && capture);
  if (!any_char(obo_file, be_permissive)) {
    if (io_peek(obo_file) == END_OF_FILE) {
      error_set(err, "file \"%s\": line %lu: unexpected end-of-file",
                io_get_filename(obo_file), io_get_line_number(obo_file));
    }
    else if ((io_peek(obo_file) == CARRIAGE_RETURN) ||
             (io_peek(obo_file) == END_OF_LINE)) {
      error_set(err, "file \"%s\": line %lu: unexpected newline",
                io_get_filename(obo_file), io_get_line_number(obo_file));
    }
    else {
      error_set(err, "file \"%s\": line %lu: unexpected character '%c'",
                io_get_filename(obo_file), io_get_line_number(obo_file),
                io_peek(obo_file));
    }
    return -1;
  }
  str_append_char(capture, io_next(obo_file));
  return 0;
}

static int tag_line(IO *obo_file, Str *tag, Str *value, Error *err)
{
  int had_err;
  error_check(err);
  assert(obo_file && tag && value);
  do {
    had_err = proc_any_char(obo_file, tag, false, err);
  } while (!had_err && any_char(obo_file, false));
  if (!had_err)
    had_err = expect(obo_file, SEPARATOR_CHAR, err);
  while (!had_err && io_peek(obo_file) == BLANK_CHAR)
    io_next(obo_file);
  if (!had_err) {
    do {
      had_err = proc_any_char(obo_file, value, true, err);
    } while (!had_err && any_char(obo_file, true));
  }
  if (!had_err) {
    if (io_peek(obo_file) == COMMENT_CHAR)
      had_err = comment_line(obo_file, err);
    else
      had_err = expect(obo_file, END_OF_LINE, err);
  }
  return had_err;
}

static int header(OBOParseTree *obo_parse_tree, IO *obo_file, Error *err)
{
  Str *tag, *value;
  int had_err;
  error_check(err);
  assert(obo_parse_tree && obo_file);
  tag = str_new();
  value = str_new();
  do {
    str_reset(tag);
    str_reset(value);
    had_err = tag_line(obo_file, tag, value, err);
    if (!had_err)
      obo_header_add(obo_parse_tree->obo_header, str_get(tag), str_get(value));
  } while (!had_err && any_char(obo_file, false));
  if (!had_err) {
    had_err = obo_header_validate(obo_parse_tree->obo_header,
                                  io_get_filename(obo_file), err);
  }
  str_delete(value);
  str_delete(tag);
  return had_err;
}

static int stanza_line(IO *obo_file, Str *type, Error *err)
{
  int had_err;
  error_check(err);
  assert(obo_file && type);
  had_err = expect(obo_file, STANZA_OPEN_CHAR, err);
  if (!had_err) {
    do {
      had_err = proc_any_char(obo_file, type, false, err);
    } while (!had_err && any_char(obo_file, false));
  }
  if (!had_err)
    had_err = expect(obo_file, STANZA_CLOSE_CHAR, err);
  if (!had_err)
    had_err = expect(obo_file, END_OF_LINE, err);
  return had_err;
}

static int stanza(OBOParseTree *obo_parse_tree, IO *obo_file, Error *err)
{
  unsigned long stanza_line_number;
  int had_err;
  Str *type, *tag, *value;
  error_check(err);
  assert(obo_parse_tree && obo_file);
  type = str_new();
  tag = str_new();
  value = str_new();
  stanza_line_number = io_get_line_number(obo_file);
  had_err = stanza_line(obo_file, type, err);
  if (!had_err) {
    OBOStanza *obo_stanza = obo_stanza_new(str_get(type), stanza_line_number,
                                           io_get_filename_str(obo_file));
    obo_parse_tree_add_stanza(obo_parse_tree, obo_stanza);
    while (!had_err &&
           (any_char(obo_file, false) || io_peek(obo_file) == COMMENT_CHAR)) {
      str_reset(tag);
      str_reset(value);
      if (io_peek(obo_file) == COMMENT_CHAR)
        had_err = comment_line(obo_file, err);
      else {
        had_err = tag_line(obo_file, tag, value, err);
        obo_stanza_add(obo_stanza, str_get(tag), str_get(value));
      }
    }
  }
  str_delete(value);
  str_delete(tag);
  str_delete(type);
  return had_err;
}

static int parse_obo_file(OBOParseTree *obo_parse_tree,
                          IO *obo_file, Error *err)
{
  int had_err = 0;
  error_check(err);
  assert(obo_parse_tree && obo_file);
  while (!had_err && ignored_char(obo_file)) {
    had_err = ignored_line(obo_file, err);
  }
  if (!had_err)
    had_err = header(obo_parse_tree, obo_file, err);
  while (!had_err && io_has_char(obo_file)) {
    switch (io_peek(obo_file)) {
      case BLANK_CHAR:
        had_err = blank_line(obo_file, err);
        break;
      case COMMENT_CHAR:
        had_err = comment_line(obo_file, err);
        break;
      case CARRIAGE_RETURN:
        io_next(obo_file);
        if (io_peek(obo_file) == END_OF_LINE)
          io_next(obo_file);
        break;
      case END_OF_LINE:
        io_next(obo_file);
        break;
      default:
        had_err = stanza(obo_parse_tree, obo_file, err);
    }
  }
  if (!had_err)
    had_err = expect(obo_file, END_OF_FILE, err);
  if (!had_err)
    had_err = obo_parse_tree_validate_stanzas(obo_parse_tree, err);
  return had_err;
}

OBOParseTree* obo_parse_tree_new(const char *obo_file_path, Error *err)
{
  OBOParseTree *obo_parse_tree;
  IO *obo_file;
  error_check(err);
  assert(obo_file_path);
  obo_file = io_new(obo_file_path, "r");
  obo_parse_tree = ma_malloc(sizeof *obo_parse_tree);
  obo_parse_tree->obo_header = obo_header_new();
  obo_parse_tree->stanzas = gt_array_new(sizeof (OBOStanza*));
  if (parse_obo_file(obo_parse_tree, obo_file, err)) {
    obo_parse_tree_delete(obo_parse_tree);
    io_delete(obo_file);
    return NULL;
  }
  io_delete(obo_file);
  return obo_parse_tree;
}

void obo_parse_tree_delete(OBOParseTree *obo_parse_tree)
{
  unsigned long i;
  if (!obo_parse_tree) return;
  for (i = 0; i < gt_array_size(obo_parse_tree->stanzas); i++)
    obo_stanza_delete(*(OBOStanza**) gt_array_get(obo_parse_tree->stanzas, i));
  gt_array_delete(obo_parse_tree->stanzas);
  obo_header_delete(obo_parse_tree->obo_header);
  ma_free(obo_parse_tree);
}

const char* obo_parse_tree_get_stanza_type(const OBOParseTree *obo_parse_tree,
                                           unsigned long stanza_num)
{
  assert(obo_parse_tree);
  return obo_stanza_get_type(*(OBOStanza**)
                             gt_array_get(obo_parse_tree->stanzas, stanza_num));
}

const char* obo_parse_tree_get_stanza_value(const OBOParseTree *obo_parse_tree,
                                            unsigned long stanza_num,
                                            const char *stanza_key)
{
  assert(obo_parse_tree);
  return obo_stanza_get_value(*(OBOStanza**)
                              gt_array_get(obo_parse_tree->stanzas, stanza_num),
                              stanza_key);
}

unsigned long obo_parse_tree_num_of_stanzas(const OBOParseTree *obo_parse_tree)
{
  assert(obo_parse_tree);
  return gt_array_size(obo_parse_tree->stanzas);
}
