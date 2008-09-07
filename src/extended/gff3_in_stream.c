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
#include <string.h>
#include "core/fileutils.h"
#include "core/progressbar.h"
#include "core/strarray.h"
#include "extended/feature_type_factory_any.h"
#include "extended/genome_stream_rep.h"
#include "extended/gff3_in_stream.h"
#include "extended/gff3_parser.h"

struct GFF3InStream
{
  const GenomeStream parent_instance;
  unsigned long next_file;
  GT_StrArray *files;
  Str *stdinstr;
  bool ensure_sorting,
       stdin_argument,
       file_is_open,
       be_verbose,
       checkids,
       own_factory;
  GenFile *fpin;
  unsigned long long line_number;
  Queue *genome_node_buffer;
  FeatureTypeFactory *feature_type_factory;
  GFF3Parser *gff3_parser;
};

#define gff3_in_stream_cast(GS)\
        genome_stream_cast(gff3_in_stream_class(), GS)

static int buffer_is_sorted(void **elem, void *info, GT_Error *err)
{
  GT_GenomeNode *current_node, **last_node;

  gt_error_check(err);
  assert(elem && info);

  current_node = *(GT_GenomeNode**) elem,
  last_node = info;

  if (*last_node && genome_node_compare(last_node, &current_node) > 0) {
    assert(*last_node);
    gt_error_set(err, "the file %s is not sorted (example: line %u and %u)",
              genome_node_get_filename(*last_node),
              genome_node_get_line_number(*last_node),
              genome_node_get_line_number(current_node));
    return -1;
  }
  else
    *last_node = current_node;
  return 0;
}

static int gff3_in_stream_next_tree(GenomeStream *gs, GT_GenomeNode **gn,
                                    GT_Error *err)
{
  GFF3InStream *is = gff3_in_stream_cast(gs);
  Str *filenamestr;
  int had_err = 0, status_code;

  gt_error_check(err);

  if (queue_size(is->genome_node_buffer) > 1) {
    /* we still have at least two nodes in the buffer -> serve from there */
    *gn = queue_get(is->genome_node_buffer);
    return 0;
  }

  /* the buffer is empty or has one element */
  assert(queue_size(is->genome_node_buffer) <= 1);

  for (;;) {
    /* open file if necessary */
    if (!is->file_is_open) {
      if (gt_strarray_size(is->files) &&
          is->next_file == gt_strarray_size(is->files)) {
        break;
      }
      if (gt_strarray_size(is->files)) {
        if (strcmp(gt_strarray_get(is->files, is->next_file), "-") == 0) {
          if (is->stdin_argument) {
            gt_error_set(err, "multiple specification of argument file \"-\"");
            had_err = -1;
            break;
          }
          is->fpin = genfile_xopen(NULL, "r");
          is->file_is_open = true;
          is->stdin_argument = true;
        }
        else {
          is->fpin = genfile_xopen(gt_strarray_get(is->files, is->next_file),
                                   "r");
          is->file_is_open = true;
        }
        is->next_file++;
      }
      else {
        is->fpin = NULL;
        is->file_is_open = true;
      }
      is->line_number = 0;

      if (!had_err && is->be_verbose) {
        printf("processing file \"%s\"\n", gt_strarray_size(is->files)
               ? gt_strarray_get(is->files, is->next_file-1) : "stdin");
      }
      if (!had_err && is->fpin && is->be_verbose) {
        progressbar_start(&is->line_number,
                          file_number_of_lines(gt_strarray_get(is->files,
                                                            is->next_file-1)));
      }
    }

    assert(is->file_is_open);

    filenamestr = gt_strarray_size(is->files)
                  ? gt_strarray_get_str(is->files, is->next_file-1)
                  : is->stdinstr;
    /* read two nodes */
    had_err = gff3parser_parse_genome_nodes(&status_code, is->gff3_parser,
                                            is->genome_node_buffer, filenamestr,
                                            &is->line_number, is->fpin, err);
    if (had_err)
      break;
    if (status_code != EOF) {
      had_err = gff3parser_parse_genome_nodes(&status_code, is->gff3_parser,
                                              is->genome_node_buffer,
                                              filenamestr, &is->line_number,
                                              is->fpin, err);
      if (had_err)
        break;
    }

    if (status_code == EOF) {
      /* end of current file */
      if (is->be_verbose) progressbar_stop();
      genfile_close(is->fpin);
      is->fpin = NULL;
      is->file_is_open = false;
      gff3parser_reset(is->gff3_parser);
      if (!gt_strarray_size(is->files))
        break;
      continue;
    }

    assert(queue_size(is->genome_node_buffer));

    /* make sure the parsed nodes are sorted */
    if (is->ensure_sorting && queue_size(is->genome_node_buffer) > 1) {
      GT_GenomeNode *last_node = NULL;
      /* a sorted stream can have at most one input file */
      assert(gt_strarray_size(is->files) == 0 ||
             gt_strarray_size(is->files) == 1);
      had_err = queue_iterate(is->genome_node_buffer, buffer_is_sorted,
                              &last_node, err);
    }
    if (!had_err) {
      *gn = queue_get(is->genome_node_buffer);
    }
    return had_err;
  }
  assert(!queue_size(is->genome_node_buffer));
  *gn = NULL;
  return had_err;
}

static void gff3_in_stream_free(GenomeStream *gs)
{
  GFF3InStream *gff3_in_stream = gff3_in_stream_cast(gs);
  gt_strarray_delete(gff3_in_stream->files);
  str_delete(gff3_in_stream->stdinstr);
  while (queue_size(gff3_in_stream->genome_node_buffer))
    genome_node_rec_delete(queue_get(gff3_in_stream->genome_node_buffer));
  queue_delete(gff3_in_stream->genome_node_buffer);
  gff3parser_delete(gff3_in_stream->gff3_parser);
  if (gff3_in_stream->own_factory)
    feature_type_factory_delete(gff3_in_stream->feature_type_factory);
  genfile_close(gff3_in_stream->fpin);
}

const GenomeStreamClass* gff3_in_stream_class(void)
{
  static const GenomeStreamClass gsc = { sizeof (GFF3InStream),
                                         gff3_in_stream_next_tree,
                                         gff3_in_stream_free };
  return &gsc;
}

/* takes ownership of <files> */
static GenomeStream* gff3_in_stream_new(GT_StrArray *files,
                                        bool ensure_sorting, bool be_verbose,
                                        bool checkids)
{
  GenomeStream *gs = genome_stream_create(gff3_in_stream_class(),
                                          ensure_sorting);
  GFF3InStream *gff3_in_stream           = gff3_in_stream_cast(gs);
  gff3_in_stream->next_file              = 0;
  gff3_in_stream->files                  = files;
  gff3_in_stream->stdinstr               = str_new_cstr("stdin");
  gff3_in_stream->ensure_sorting         = ensure_sorting;
  gff3_in_stream->stdin_argument         = false;
  gff3_in_stream->file_is_open           = false;
  gff3_in_stream->fpin                   = NULL;
  gff3_in_stream->line_number            = 0;
  gff3_in_stream->genome_node_buffer     = queue_new();
  gff3_in_stream->feature_type_factory   = feature_type_factory_any_new();
  gff3_in_stream->own_factory            = true;
  gff3_in_stream->checkids               = checkids;
  gff3_in_stream->gff3_parser            = gff3parser_new(checkids,
                                                          gff3_in_stream
                                                        ->feature_type_factory);
  gff3_in_stream->be_verbose             = be_verbose;
  return gs;
}

void gff3_in_stream_set_feature_type_factory(GenomeStream *gs,
                                             FeatureTypeFactory
                                             *feature_type_factory)
{
  GFF3InStream *is = gff3_in_stream_cast(gs);
  assert(is);
  gff3parser_delete(is->gff3_parser);
  if (is->own_factory) {
    feature_type_factory_delete(is->feature_type_factory);
    is->own_factory = false;
  }
  is->gff3_parser = gff3parser_new(is->checkids, feature_type_factory);
  is->feature_type_factory = feature_type_factory;
}

GT_StrArray* gff3_in_stream_get_used_types(GenomeStream *gs)
{
  GFF3InStream *is = gff3_in_stream_cast(gs);
  assert(is);
  return feature_type_factory_get_used_types(is->feature_type_factory);
}

void gff3_in_stream_set_offset(GenomeStream *gs, long offset)
{
  GFF3InStream *is = gff3_in_stream_cast(gs);
  assert(is);
  gff3parser_set_offset(is->gff3_parser, offset);
}

int gff3_in_stream_set_offsetfile(GenomeStream *gs, Str *offsetfile, GT_Error *err)
{
  GFF3InStream *is = gff3_in_stream_cast(gs);
  assert(is);
  return gff3parser_set_offsetfile(is->gff3_parser, offsetfile, err);
}

void gff3_in_stream_enable_tidy_mode(GenomeStream *gs)
{
  GFF3InStream *is = gff3_in_stream_cast(gs);
  assert(is);
  gff3parser_enable_tidy_mode(is->gff3_parser);
}

GenomeStream* gff3_in_stream_new_unsorted(int num_of_files,
                                          const char **filenames,
                                          bool be_verbose, bool checkids)
{
  int i;
  GT_StrArray *files = gt_strarray_new();
  for (i = 0; i < num_of_files; i++)
    gt_strarray_add_cstr(files, filenames[i]);
  return gff3_in_stream_new(files, false, be_verbose, checkids);
}

GenomeStream* gff3_in_stream_new_sorted(const char *filename, bool be_verbose)
{
  GT_StrArray *files = gt_strarray_new();
  if (filename)
    gt_strarray_add_cstr(files, filename);
  return gff3_in_stream_new(files, true, be_verbose, false);
}
