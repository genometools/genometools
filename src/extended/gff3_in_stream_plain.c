/*
  Copyright (c) 2006-2012 Gordon Gremme <gremme@zbh.uni-hamburg.de>
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

#include <string.h>
#include "core/assert_api.h"
#include "core/class_alloc_lock.h"
#include "core/cstr_table.h"
#include "core/fileutils_api.h"
#include "core/queue.h"
#include "core/progressbar.h"
#include "core/str_array.h"
#include "extended/genome_node.h"
#include "extended/gff3_in_stream_plain.h"
#include "extended/gff3_parser.h"
#include "extended/node_stream_api.h"

struct GtGFF3InStreamPlain {
  const GtNodeStream parent_instance;
  unsigned long next_file;
  GtStrArray *files;
  GtStr *stdinstr;
  bool ensure_sorting,
       stdin_argument,
       stdin_processed,
       file_is_open,
       progress_bar;
  GtFile *fpin;
  unsigned long long line_number;
  GtQueue *genome_node_buffer;
  GtGFF3Parser *gff3_parser;
  GtCstrTable *used_types;
};

#define gff3_in_stream_plain_cast(NS)\
        gt_node_stream_cast(gt_gff3_in_stream_plain_class(), NS)

static int buffer_is_sorted(void **elem, void *info, GtError *err)
{
  GtGenomeNode *current_node, **last_node;

  gt_error_check(err);
  gt_assert(elem && info);

  current_node = *(GtGenomeNode**) elem,
  last_node = info;

  if (*last_node && gt_genome_node_compare(last_node, &current_node) > 0) {
    gt_assert(*last_node);
    gt_error_set(err, "the file %s is not sorted (example: line %u and %u)",
              gt_genome_node_get_filename(*last_node),
              gt_genome_node_get_line_number(*last_node),
              gt_genome_node_get_line_number(current_node));
    return -1;
  }
  else
    *last_node = current_node;
  return 0;
}

static int gff3_in_stream_plain_next(GtNodeStream *ns, GtGenomeNode **gn,
                                     GtError *err)
{
  GtGFF3InStreamPlain *is = gff3_in_stream_plain_cast(ns);
  GtStr *filenamestr;
  int had_err = 0, status_code;

  gt_error_check(err);

  if (gt_queue_size(is->genome_node_buffer) > 1) {
    /* we still have at least two nodes in the buffer -> serve from there */
    *gn = gt_queue_get(is->genome_node_buffer);
    return 0;
  }

  /* the buffer is empty or has one element */
  gt_assert(gt_queue_size(is->genome_node_buffer) <= 1);

  for (;;) {
    /* open file if necessary */
    if (!is->file_is_open) {
      if (gt_str_array_size(is->files) &&
          is->next_file == gt_str_array_size(is->files)) {
        break;
      }
      if (gt_str_array_size(is->files)) {
        if (strcmp(gt_str_array_get(is->files, is->next_file), "-") == 0) {
          if (is->stdin_argument) {
            gt_error_set(err, "multiple specification of argument file \"-\"");
            had_err = -1;
            break;
          }
          is->fpin = gt_file_xopen(NULL, "r");
          is->file_is_open = true;
          is->stdin_argument = true;
        }
        else {
          is->fpin = gt_file_xopen(gt_str_array_get(is->files,
                                                       is->next_file), "r");
          is->file_is_open = true;
        }
        is->next_file++;
      }
      else {
        if (is->stdin_processed)
          break;
        is->fpin = NULL;
        is->file_is_open = true;
      }
      is->line_number = 0;

      if (!had_err && is->progress_bar) {
        printf("processing file \"%s\"\n", gt_str_array_size(is->files)
               ? gt_str_array_get(is->files, is->next_file-1) : "stdin");
      }
      if (!had_err && is->fpin && is->progress_bar) {
        gt_progressbar_start(&is->line_number,
                            gt_file_number_of_lines(gt_str_array_get(is->files,
                                                             is->next_file-1)));
      }
    }

    gt_assert(is->file_is_open);

    filenamestr = gt_str_array_size(is->files)
                  ? gt_str_array_get_str(is->files, is->next_file-1)
                  : is->stdinstr;
    /* read two nodes */
    had_err = gt_gff3_parser_parse_genome_nodes(is->gff3_parser, &status_code,
                                                is->genome_node_buffer,
                                                is->used_types, filenamestr,
                                                &is->line_number, is->fpin,
                                                err);
    if (had_err)
      break;
    if (status_code != EOF) {
      had_err = gt_gff3_parser_parse_genome_nodes(is->gff3_parser, &status_code,
                                                  is->genome_node_buffer,
                                                  is->used_types, filenamestr,
                                                  &is->line_number, is->fpin,
                                                  err);
      if (had_err)
        break;
    }

    if (status_code == EOF) {
      /* end of current file */
      if (is->progress_bar) gt_progressbar_stop();
      gt_file_delete(is->fpin);
      is->fpin = NULL;
      is->file_is_open = false;
      gt_gff3_parser_reset(is->gff3_parser);
      if (!gt_str_array_size(is->files)) {
        is->stdin_processed = true;
        break;
      }
      continue;
    }

    gt_assert(gt_queue_size(is->genome_node_buffer));

    /* make sure the parsed nodes are sorted */
    if (is->ensure_sorting && gt_queue_size(is->genome_node_buffer) > 1) {
      GtGenomeNode *last_node = NULL;
      /* a sorted stream can have at most one input file */
      gt_assert(gt_str_array_size(is->files) == 0 ||
                gt_str_array_size(is->files) == 1);
      had_err = gt_queue_iterate(is->genome_node_buffer, buffer_is_sorted,
                                 &last_node, err);
    }
    if (!had_err) {
      *gn = gt_queue_get(is->genome_node_buffer);
    }
    return had_err;
  }
  gt_assert(!gt_queue_size(is->genome_node_buffer));
  *gn = NULL;
  return had_err;
}

static void gff3_in_stream_plain_free(GtNodeStream *ns)
{
  GtGFF3InStreamPlain *gff3_in_stream_plain = gff3_in_stream_plain_cast(ns);
  gt_str_array_delete(gff3_in_stream_plain->files);
  gt_str_delete(gff3_in_stream_plain->stdinstr);
  while (gt_queue_size(gff3_in_stream_plain->genome_node_buffer)) {
    gt_genome_node_delete(gt_queue_get(gff3_in_stream_plain
                                       ->genome_node_buffer));
  }
  gt_queue_delete(gff3_in_stream_plain->genome_node_buffer);
  gt_gff3_parser_delete(gff3_in_stream_plain->gff3_parser);
  gt_cstr_table_delete(gff3_in_stream_plain->used_types);
  gt_file_delete(gff3_in_stream_plain->fpin);
}

const GtNodeStreamClass* gt_gff3_in_stream_plain_class(void)
{
  static const GtNodeStreamClass *nsc = NULL;
  gt_class_alloc_lock_enter();
  if (!nsc) {
    nsc = gt_node_stream_class_new(sizeof (GtGFF3InStreamPlain),
                                   gff3_in_stream_plain_free,
                                   gff3_in_stream_plain_next);
  }
  gt_class_alloc_lock_leave();
  return nsc;
}

/* takes ownership of <files> */
static GtNodeStream* gff3_in_stream_plain_new(GtStrArray *files,
                                              bool ensure_sorting)
{
  GtNodeStream *ns = gt_node_stream_create(gt_gff3_in_stream_plain_class(),
                                           ensure_sorting);
  GtGFF3InStreamPlain *gff3_in_stream_plain = gff3_in_stream_plain_cast(ns);
  gff3_in_stream_plain->files               = files;
  gff3_in_stream_plain->stdinstr            = gt_str_new_cstr("stdin");
  gff3_in_stream_plain->ensure_sorting      = ensure_sorting;
  gff3_in_stream_plain->genome_node_buffer  = gt_queue_new();
  gff3_in_stream_plain->gff3_parser         = gt_gff3_parser_new(NULL);
  gff3_in_stream_plain->used_types          = gt_cstr_table_new();
  return ns;
}

void gt_gff3_in_stream_plain_check_id_attributes(GtGFF3InStreamPlain *is)
{
  gt_assert(is);
  gt_gff3_parser_check_id_attributes(is->gff3_parser);
}

void gt_gff3_in_stream_plain_check_region_boundaries(GtGFF3InStreamPlain *is)
{
  gt_assert(is);
  gt_gff3_parser_check_region_boundaries(is->gff3_parser);
}

void gt_gff3_in_stream_plain_do_not_check_region_boundaries(GtGFF3InStreamPlain
                                                            *is)
{
  gt_assert(is);
  gt_gff3_parser_do_not_check_region_boundaries(is->gff3_parser);
}

void gt_gff3_in_stream_plain_show_progress_bar(GtGFF3InStreamPlain *is)
{
  gt_assert(is);
  is->progress_bar = true;
}

void gt_gff3_in_stream_plain_set_type_checker(GtNodeStream *ns,
                                              GtTypeChecker *type_checker)
{
  GtGFF3InStreamPlain *is = gff3_in_stream_plain_cast(ns);
  gt_assert(is);
  gt_gff3_parser_set_type_checker(is->gff3_parser, type_checker);
}

GtStrArray* gt_gff3_in_stream_plain_get_used_types(GtNodeStream *ns)
{
  GtGFF3InStreamPlain *is = gff3_in_stream_plain_cast(ns);
  gt_assert(is);
  return gt_cstr_table_get_all(is->used_types);
}

void gt_gff3_in_stream_plain_set_offset(GtNodeStream *ns, long offset)
{
  GtGFF3InStreamPlain *is = gff3_in_stream_plain_cast(ns);
  gt_assert(is);
  gt_gff3_parser_set_offset(is->gff3_parser, offset);
}

int gt_gff3_in_stream_plain_set_offsetfile(GtNodeStream *ns, GtStr *offsetfile,
                                           GtError *err)
{
  GtGFF3InStreamPlain *is = gff3_in_stream_plain_cast(ns);
  gt_assert(is);
  return gt_gff3_parser_set_offsetfile(is->gff3_parser, offsetfile, err);
}

void gt_gff3_in_stream_plain_enable_strict_mode(GtNodeStream *ns)
{
  GtGFF3InStreamPlain *is = gff3_in_stream_plain_cast(ns);
  gt_assert(is);
  gt_gff3_parser_enable_strict_mode(is->gff3_parser);
}

void gt_gff3_in_stream_plain_enable_tidy_mode(GtNodeStream *ns)
{
  GtGFF3InStreamPlain *is = gff3_in_stream_plain_cast(ns);
  gt_assert(is);
  gt_gff3_parser_enable_tidy_mode(is->gff3_parser);
}

GtNodeStream* gt_gff3_in_stream_plain_new_unsorted(int num_of_files,
                                                   const char **filenames)
{
  int i;
  GtStrArray *files = gt_str_array_new();
  for (i = 0; i < num_of_files; i++)
    gt_str_array_add_cstr(files, filenames[i]);
  return gff3_in_stream_plain_new(files, false);
}

GtNodeStream* gt_gff3_in_stream_plain_new_sorted(const char *filename)
{
  GtStrArray *files = gt_str_array_new();
  if (filename)
    gt_str_array_add_cstr(files, filename);
  return gff3_in_stream_plain_new(files, true);
}
