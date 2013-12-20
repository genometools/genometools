/*
  Copyright (c) 2013 Sascha Steinbiss <ss34@sanger.ac.uk>

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
#include "core/cstr_api.h"
#include "core/fa.h"
#include "core/grep_api.h"
#include "core/hashmap_api.h"
#include "core/io.h"
#include "core/log_api.h"
#include "core/ma.h"
#include "core/unused_api.h"
#include "core/warning_api.h"
#include "extended/xrf_abbr_entry.h"
#include "extended/xrf_abbr_parse_tree.h"

#define XRF_BLANK_CHAR           ' '
#define XRF_COMMENT_CHAR         '!'
#define XRF_SEPARATOR_CHAR       ':'

#define XRF_LABEL_ABBREVIATION    "abbreviation"
#define XRF_LABEL_SHORTHAND_NAME  "shorthand_name"
#define XRF_LABEL_DATABASE        "database"
#define XRF_LABEL_OBJECT          "object"
#define XRF_LABEL_SYNONYM         "synonym"
#define XRF_LABEL_EXAMPLE_ID      "example_id"
#define XRF_LABEL_LOCAL_ID_SYNTAX "local_id_syntax"
#define XRF_LABEL_GENERIC_URL     "generic_url"
#define XRF_LABEL_URL_SYNTAX      "url_syntax"
#define XRF_LABEL_URL_EXAMPLE     "url_example"
#define XRF_LABEL_IS_OBSOLETE     "is_obsolete"
#define XRF_LABEL_CONSIDER        "consider"
#define XRF_LABEL_REPLACED_BY     "replaced_by"

struct GtXRFAbbrParseTree {
  GtArray *entries;
};

static bool gt_xrf_abbr_parse_tree_valid_label(const char *label)
{
  if (strcmp(label, XRF_LABEL_ABBREVIATION) == 0
        || strcmp(label, XRF_LABEL_SHORTHAND_NAME) == 0
        || strcmp(label, XRF_LABEL_DATABASE) == 0
        || strcmp(label, XRF_LABEL_OBJECT) == 0
        || strcmp(label, XRF_LABEL_SYNONYM) == 0
        || strcmp(label, XRF_LABEL_EXAMPLE_ID) == 0
        || strcmp(label, XRF_LABEL_LOCAL_ID_SYNTAX) == 0
        || strcmp(label, XRF_LABEL_GENERIC_URL) == 0
        || strcmp(label, XRF_LABEL_URL_SYNTAX) == 0
        || strcmp(label, XRF_LABEL_URL_EXAMPLE) == 0
        || strcmp(label, XRF_LABEL_IS_OBSOLETE) == 0
        || strcmp(label, XRF_LABEL_CONSIDER) == 0
        || strcmp(label, XRF_LABEL_REPLACED_BY) == 0)
    return true;
  return false;
}

static void gt_xrf_abbr_parse_tree_add_entry(GtXRFAbbrParseTree
                                             *xrf_abbr_parse_tree,
                                             GtXRFAbbrEntry *xrf_abbr_entry)
{
  gt_assert(xrf_abbr_parse_tree && xrf_abbr_entry);
  gt_array_add(xrf_abbr_parse_tree->entries, xrf_abbr_entry);
}

static int gt_xrf_abbr_parse_tree_validate_entries(const GtXRFAbbrParseTree
                                                           *xrf_abbr_parse_tree,
                                                   GtError *err)
{
  GtUword i;
  GtHashmap *abbrvs;
  const char *value;
  int had_err = 0;
  gt_error_check(err);
  gt_assert(xrf_abbr_parse_tree);

  abbrvs = gt_hashmap_new(GT_HASH_STRING, NULL, NULL);
  for (i = 0; !had_err
         && i < gt_xrf_abbr_parse_tree_num_of_entries(xrf_abbr_parse_tree);
       i++) {
    GtXRFAbbrEntry *entry = *(GtXRFAbbrEntry**)
                                gt_array_get(xrf_abbr_parse_tree->entries, i);
    if (!(value = gt_xrf_abbr_entry_get_value(entry, XRF_LABEL_ABBREVIATION))) {
      gt_error_set(err, "file \"%s\": line "GT_WU": required "
                        "label \"" XRF_LABEL_ABBREVIATION "\" missing",
                   gt_xrf_abbr_entry_filename(entry),
                   gt_xrf_abbr_entry_line(entry));
      had_err = -1;
    }
    if (!had_err) {
      gt_assert(value);
      if (gt_hashmap_get(abbrvs, value)) {
        gt_error_set(err, "file \"%s\": line "GT_WU": duplicate abbreviation "
                          "\"%s\", must be unique",
                     gt_xrf_abbr_entry_filename(entry),
                     gt_xrf_abbr_entry_line(entry),
                     value);
        had_err = -1;
      } else {
        gt_hashmap_add(abbrvs, (void*) value, (void*) value);
      }
    }
    if (!had_err && (value = gt_xrf_abbr_entry_get_value(entry,
                                                XRF_LABEL_SHORTHAND_NAME))) {
      if (strlen(value) >= 10) {
        gt_error_set(err, "file \"%s\": line "GT_WU": length of "
                          "shorthand name \"%s\" "
                          "is not less than 10 characters",
                     gt_xrf_abbr_entry_filename(entry),
                     gt_xrf_abbr_entry_line(entry), value);
        had_err = -1;
      }
    }
    if (!had_err && (value = gt_xrf_abbr_entry_get_value(entry,
                                              XRF_LABEL_LOCAL_ID_SYNTAX))) {
      GtError *regex_error = gt_error_new();
      bool match;
      if (gt_grep(&match, value, "", regex_error)) {
        gt_error_set(err, "file \"%s\": line "GT_WU": invalid "
                          "regular expression \"%s\" (%s)",
                     gt_xrf_abbr_entry_filename(entry),
                     gt_xrf_abbr_entry_line(entry), value,
                     gt_error_get(regex_error));
        had_err = -1;
      }
      gt_error_delete(regex_error);
    }
  }
  gt_hashmap_delete(abbrvs);
  return had_err;
}

static bool gt_xrf_abbr_parse_tree_any_char(GtIO *xrf_abbr_file,
                                            bool be_permissive)
{
  switch (gt_io_peek(xrf_abbr_file)) {
    case XRF_BLANK_CHAR:
    case XRF_SEPARATOR_CHAR:
      if (be_permissive)
        return true;
    case XRF_COMMENT_CHAR:
    case GT_CARRIAGE_RETURN:
    case GT_END_OF_LINE:
    case GT_END_OF_FILE:
      return false;
  }
  return true;
}

static bool gt_xrf_abbr_parse_tree_ignored_char(GtIO *xrf_abbr_file)
{
  char cc = gt_io_peek(xrf_abbr_file);
  if ((cc == XRF_BLANK_CHAR) || (cc == XRF_COMMENT_CHAR) ||
      (cc == GT_CARRIAGE_RETURN) || (cc == GT_END_OF_LINE))
    return true;
  return false;
}

static int gt_xrf_abbr_parse_tree_comment_line(GtIO *xrf_abbr_file,
                                               GtError *err)
{
  int had_err;
  gt_error_check(err);
  gt_log_log("comment");
  had_err = gt_io_expect(xrf_abbr_file, XRF_COMMENT_CHAR, err);
  while (!had_err) {
    switch (gt_io_peek(xrf_abbr_file)) {
      case GT_CARRIAGE_RETURN:
        gt_io_next(xrf_abbr_file);
        if (gt_io_peek(xrf_abbr_file) == GT_END_OF_LINE)
          gt_io_next(xrf_abbr_file);
        return had_err;
      case GT_END_OF_LINE:
        gt_io_next(xrf_abbr_file);
        /*@fallthrough@*/
      case GT_END_OF_FILE:
        return had_err;
      default:
        gt_io_next(xrf_abbr_file);
    }
  }
  return had_err;
}

static int gt_xrf_abbr_parse_tree_blank_line(GtIO *xrf_abbr_file, GtError *err)
{
  int had_err;
  gt_error_check(err);
  gt_log_log("blank");
  had_err = gt_io_expect(xrf_abbr_file, XRF_BLANK_CHAR, err);
  while (!had_err) {
    char cc = gt_io_peek(xrf_abbr_file);
    if (cc == XRF_COMMENT_CHAR)
      return gt_xrf_abbr_parse_tree_comment_line(xrf_abbr_file, err);
    else if (cc == GT_CARRIAGE_RETURN) {
      gt_io_next(xrf_abbr_file);
      if (gt_io_peek(xrf_abbr_file) == GT_END_OF_LINE)
        gt_io_next(xrf_abbr_file);
      break;
    }
    else if ((cc == GT_END_OF_LINE) || (cc == GT_END_OF_FILE)) {
      gt_io_next(xrf_abbr_file);
      break;
    }
    else
      had_err = gt_io_expect(xrf_abbr_file, XRF_BLANK_CHAR, err);
  }
  return had_err;
}

static bool gt_xrf_abbr_parse_tree_ignored_line(GtIO *xrf_abbr_file,
                                                GtError *err)
{
  gt_error_check(err);
    gt_log_log("ignored");
  if (gt_io_peek(xrf_abbr_file) == XRF_BLANK_CHAR)
    return gt_xrf_abbr_parse_tree_blank_line(xrf_abbr_file, err);
  if (gt_io_peek(xrf_abbr_file) == XRF_COMMENT_CHAR)
    return gt_xrf_abbr_parse_tree_comment_line(xrf_abbr_file, err);
  gt_io_next(xrf_abbr_file);
  return false;
}

static int gt_xrf_abbr_parse_tree_proc_any_char(GtIO *xrf_abbr_file,
                                                GtStr *capture,
                                                bool be_permissive,
                                                GtError *err)
{
  gt_error_check(err);
  gt_assert(xrf_abbr_file && capture);
  if (!gt_xrf_abbr_parse_tree_any_char(xrf_abbr_file, be_permissive)) {
    if (gt_io_peek(xrf_abbr_file) == GT_END_OF_FILE) {
      gt_error_set(err, "file \"%s\": line "GT_WU": unexpected end-of-file",
                gt_io_get_filename(xrf_abbr_file),
                gt_io_get_line_number(xrf_abbr_file));
    }
    else if ((gt_io_peek(xrf_abbr_file) == GT_CARRIAGE_RETURN) ||
             (gt_io_peek(xrf_abbr_file) == GT_END_OF_LINE)) {
      gt_error_set(err, "file \"%s\": line "GT_WU": unexpected newline",
                gt_io_get_filename(xrf_abbr_file),
                gt_io_get_line_number(xrf_abbr_file));
    }
    else {
      gt_error_set(err, "file \"%s\": line "GT_WU": unexpected character '%c'",
                gt_io_get_filename(xrf_abbr_file),
                gt_io_get_line_number(xrf_abbr_file),
                gt_io_peek(xrf_abbr_file));
    }
    return -1;
  }
  gt_str_append_char(capture, gt_io_next(xrf_abbr_file));
  return 0;
}

static int gt_xrf_abbr_parse_tree_tag_line(GtIO *xrf_abbr_file, GtStr *tag,
                                           GtStr *value, GtError *err)
{
  int had_err = 0;
  gt_error_check(err);
  gt_log_log("tag");
  gt_assert(xrf_abbr_file && tag && value);
  do {
    had_err = gt_xrf_abbr_parse_tree_proc_any_char(xrf_abbr_file, tag,
                                                   false, err);
  } while (!had_err && gt_xrf_abbr_parse_tree_any_char(xrf_abbr_file, false));
  if (!had_err)
    had_err = gt_io_expect(xrf_abbr_file, XRF_SEPARATOR_CHAR, err);
  while (!had_err && gt_io_peek(xrf_abbr_file) == XRF_BLANK_CHAR)
    gt_io_next(xrf_abbr_file);
  if (!had_err) {
    do {
      had_err = gt_xrf_abbr_parse_tree_proc_any_char(xrf_abbr_file, value,
                                                     true, err);
    } while (!had_err && gt_xrf_abbr_parse_tree_any_char(xrf_abbr_file, true));
  }
  if (!had_err) {
    if (gt_io_peek(xrf_abbr_file) == XRF_COMMENT_CHAR)
      had_err = gt_xrf_abbr_parse_tree_comment_line(xrf_abbr_file, err);
    else
      had_err = gt_io_expect(xrf_abbr_file, GT_END_OF_LINE, err);
  }
  if (!had_err && !gt_xrf_abbr_parse_tree_valid_label(gt_str_get(tag))) {
    gt_warning("file \"%s\": line "GT_WU": unknown label \"%s\"",
                gt_io_get_filename(xrf_abbr_file),
                gt_io_get_line_number(xrf_abbr_file),
                gt_str_get(tag));
  }
  gt_log_log("parsed line %s/%s", gt_str_get(tag), gt_str_get(value));
  return had_err;
}

static int gt_xrf_abbr_parse_tree_entry(GtXRFAbbrParseTree *xrf_abbr_parse_tree,
                 GtIO *xrf_abbr_file, GtError *err)
{
  GtUword entry_line_number;
  int had_err = 0;
  GtStr *tag, *value;
  gt_error_check(err);
  gt_assert(xrf_abbr_parse_tree && xrf_abbr_file);
  tag = gt_str_new();
  value = gt_str_new();
  entry_line_number = gt_io_get_line_number(xrf_abbr_file);
  if (!had_err) {
    GtXRFAbbrEntry *xrf_abbr_entry =
      gt_xrf_abbr_entry_new(entry_line_number,
                            gt_io_get_filename_str(xrf_abbr_file));
    gt_xrf_abbr_parse_tree_add_entry(xrf_abbr_parse_tree, xrf_abbr_entry);
    while (!had_err &&
           (gt_xrf_abbr_parse_tree_any_char(xrf_abbr_file, false) ||
            gt_io_peek(xrf_abbr_file) == XRF_COMMENT_CHAR)) {
      gt_str_reset(tag);
      gt_str_reset(value);
      if (gt_io_peek(xrf_abbr_file) == XRF_COMMENT_CHAR)
        had_err = gt_xrf_abbr_parse_tree_comment_line(xrf_abbr_file, err);
      else {
        had_err = gt_xrf_abbr_parse_tree_tag_line(xrf_abbr_file, tag, value,
                                                  err);
        gt_xrf_abbr_entry_add(xrf_abbr_entry, gt_str_get(tag),
                              gt_str_get(value));
      }
    }
  }
  gt_str_delete(value);
  gt_str_delete(tag);
  return had_err;
}

static int parse_xrf_abbr_file(GtXRFAbbrParseTree *xrf_abbr_parse_tree,
                               GtIO *xrf_abbr_file, GtError *err)
{
  int had_err = 0;
  gt_error_check(err);
  gt_assert(xrf_abbr_parse_tree && xrf_abbr_file);
  while (!had_err && gt_xrf_abbr_parse_tree_ignored_char(xrf_abbr_file)) {
    had_err = gt_xrf_abbr_parse_tree_ignored_line(xrf_abbr_file, err);
  }
  while (!had_err && gt_io_has_char(xrf_abbr_file)) {
    switch (gt_io_peek(xrf_abbr_file)) {
      case XRF_BLANK_CHAR:
        had_err = gt_xrf_abbr_parse_tree_blank_line(xrf_abbr_file, err);
        break;
      case XRF_COMMENT_CHAR:
        had_err = gt_xrf_abbr_parse_tree_comment_line(xrf_abbr_file, err);
        break;
      case GT_CARRIAGE_RETURN:
        gt_io_next(xrf_abbr_file);
        if (gt_io_peek(xrf_abbr_file) == GT_END_OF_LINE)
          gt_io_next(xrf_abbr_file);
        break;
      case GT_END_OF_LINE:
        gt_io_next(xrf_abbr_file);
        break;
      default:
        had_err = gt_xrf_abbr_parse_tree_entry(xrf_abbr_parse_tree,
                                               xrf_abbr_file, err);
    }
  }
  if (!had_err)
    had_err = gt_io_expect(xrf_abbr_file, GT_END_OF_FILE, err);
  if (!had_err)
    had_err = gt_xrf_abbr_parse_tree_validate_entries(xrf_abbr_parse_tree,
                                                      err);
  return had_err;
}

GtXRFAbbrParseTree* gt_xrf_abbr_parse_tree_new(const char *xrf_abbr_file_path,
                                               GtError *err)
{
  GtXRFAbbrParseTree *xrf_abbr_parse_tree;
  GtIO *xrf_abbr_file;
  gt_error_check(err);
  gt_assert(xrf_abbr_file_path);
  xrf_abbr_file = gt_io_new(xrf_abbr_file_path, "r");
  xrf_abbr_parse_tree = gt_malloc(sizeof *xrf_abbr_parse_tree);
  xrf_abbr_parse_tree->entries = gt_array_new(sizeof (GtXRFAbbrEntry*));
  if (parse_xrf_abbr_file(xrf_abbr_parse_tree, xrf_abbr_file, err)) {
    gt_xrf_abbr_parse_tree_delete(xrf_abbr_parse_tree);
    gt_io_delete(xrf_abbr_file);
    return NULL;
  }
  gt_io_delete(xrf_abbr_file);
  return xrf_abbr_parse_tree;
}

void gt_xrf_abbr_parse_tree_delete(GtXRFAbbrParseTree *xrf_abbr_parse_tree)
{
  GtUword i;
  if (!xrf_abbr_parse_tree) return;
  for (i = 0; i < gt_array_size(xrf_abbr_parse_tree->entries); i++) {
    gt_xrf_abbr_entry_delete(*(GtXRFAbbrEntry**)
                         gt_array_get(xrf_abbr_parse_tree->entries, i));
  }
  gt_array_delete(xrf_abbr_parse_tree->entries);
  gt_free(xrf_abbr_parse_tree);
}

const char* gt_xrf_abbr_parse_tree_get_entry_value(const GtXRFAbbrParseTree
                                                           *xrf_abbr_parse_tree,
                                                   GtUword entry_num,
                                                   const char *entry_key)
{
  gt_assert(xrf_abbr_parse_tree);
  return gt_xrf_abbr_entry_get_value(*(GtXRFAbbrEntry**)
                                      gt_array_get(xrf_abbr_parse_tree->entries,
                                                   entry_num), entry_key);
}

const GtXRFAbbrEntry* gt_xrf_abbr_parse_tree_get_entry(const GtXRFAbbrParseTree
                                                           *xrf_abbr_parse_tree,
                                                       GtUword entry_num)
{
  gt_assert(xrf_abbr_parse_tree);
  return *(GtXRFAbbrEntry**) gt_array_get(xrf_abbr_parse_tree->entries,
                                          entry_num);
}

GtUword gt_xrf_abbr_parse_tree_num_of_entries(const GtXRFAbbrParseTree
                                               *xrf_abbr_parse_tree)
{
  gt_assert(xrf_abbr_parse_tree);
  return gt_array_size(xrf_abbr_parse_tree->entries);
}
