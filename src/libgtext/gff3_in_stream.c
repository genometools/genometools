/*
  Copyright (c) 2006-2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2006-2007 Center for Bioinformatics, University of Hamburg

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
#include "libgtcore/fileutils.h"
#include "libgtcore/progressbar.h"
#include "libgtcore/strarray.h"
#include "libgtext/genome_stream_rep.h"
#include "libgtext/gff3_in_stream.h"
#include "libgtext/gff3_parser.h"

struct GFF3InStream
{
  const GenomeStream parent_instance;
  int next_file;
  StrArray *files; /* contains char* to filenames */
  bool ensure_sorting,
       stdin_argument,
       file_is_open,
       be_verbose;
  GenFile *fpin;
  unsigned long line_number;
  Queue *genome_node_buffer;
  GFF3Parser *gff3_parser;
  GenomeNode *last_node;
};

#define gff3_in_stream_cast(GS)\
        genome_stream_cast(gff3_in_stream_class(), GS)

static int gff3_in_stream_next_tree(GenomeStream *gs, GenomeNode **gn, Env *env)
{
  GFF3InStream *is = gff3_in_stream_cast(gs);
  unsigned long i;
  Str *filenamestr;
  int had_err = 0, status_code;

  env_error_check(env);

  if (queue_size(is->genome_node_buffer)) {
    /* we still have a node in the buffer -> serve it from there */
    *gn = *(GenomeNode**) queue_get(is->genome_node_buffer);
    return 0;
  }

  /* the buffer is empty */
  assert(!queue_size(is->genome_node_buffer));

  for (;;) {
    /* open file if necessary */
    if (!is->file_is_open) {
      if (strarray_size(is->files) && is->next_file == strarray_size(is->files))
        break;
      if (strarray_size(is->files)) {
        if (strcmp(strarray_get(is->files, is->next_file), "-") == 0) {
          if (is->stdin_argument) {
            env_error_set(env,
                          "multiple specification of argument file \"-\"\n");
            had_err = -1;
            break;
          }
          is->fpin = NULL;
          is->file_is_open = true;
          is->stdin_argument = true;
        }
        else {
          is->fpin = genfile_xopen(strarray_get(is->files, is->next_file), "r",
                                   env);
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
        printf("processing file \"%s\"\n", strarray_size(is->files)
               ? strarray_get(is->files, is->next_file-1) : "stdin");
      }
      if (!had_err && is->fpin && is->be_verbose) {
        progressbar_start(&is->line_number,
                          file_number_of_lines(strarray_get(is->files,
                                                            is->next_file-1),
                                               env));
      }
    }

    assert(is->file_is_open);

    filenamestr = str_new_cstr(strarray_size(is->files)
                               ? strarray_get(is->files, is->next_file-1)
                               : "stdin", env);
    had_err = gff3parser_parse_genome_nodes(&status_code, is->gff3_parser,
                                            is->genome_node_buffer, filenamestr,
                                            &is->line_number, is->fpin, env);
    str_delete(filenamestr, env);
    if (had_err)
      break;

    if (status_code == EOF) {
      /* end of current file */
      if (is->be_verbose) progressbar_stop();
      genfile_xclose(is->fpin, env);
      is->fpin = NULL;
      is->file_is_open = false;
      gff3parser_reset(is->gff3_parser, env);
      if (!strarray_size(is->files))
        break;
      continue;
    }

    assert(queue_size(is->genome_node_buffer));

    /* make sure the parsed nodes are sorted */
    if (is->ensure_sorting) {
      for (i = 0; i < queue_size(is->genome_node_buffer); i++) {
        if (!genome_node_tree_is_sorted(&is->last_node, *(GenomeNode**)
                                        queue_get_elem(is->genome_node_buffer,
                                                       i), env)) {
          assert(is->last_node);
          /* a sorted stream can have at most one input file */
          assert(strarray_size(is->files) == 0 ||
                 strarray_size(is->files) == 1);
          env_error_set(env,
                    "the file %s is not sorted (example: line %lu and %lu)",
                    genome_node_get_filename(is->last_node),
                    genome_node_get_line_number(is->last_node),
                    genome_node_get_line_number(*(GenomeNode**)
                                    queue_get_elem(is->genome_node_buffer, i)));
          had_err = -1;
          break;
        }
      }
    }
    if (!had_err)
      *gn = *(GenomeNode**) queue_get(is->genome_node_buffer);
    return had_err;
  }
  *gn = NULL;
  return had_err;
}

static void gff3_in_stream_free(GenomeStream *gs, Env *env)
{
  GFF3InStream *gff3_in_stream = gff3_in_stream_cast(gs);
  strarray_delete(gff3_in_stream->files, env);
  while (queue_size(gff3_in_stream->genome_node_buffer)) {
    genome_node_rec_delete(*(GenomeNode**)
                           queue_get(gff3_in_stream->genome_node_buffer), env);
  }
  queue_delete(gff3_in_stream->genome_node_buffer, env);
  gff3parser_delete(gff3_in_stream->gff3_parser, env);
  genome_node_delete(gff3_in_stream->last_node, env);
  genfile_xclose(gff3_in_stream->fpin, env);
}

const GenomeStreamClass* gff3_in_stream_class(void)
{
  static const GenomeStreamClass gsc = { sizeof (GFF3InStream),
                                         gff3_in_stream_next_tree,
                                         gff3_in_stream_free };
  return &gsc;
}

static GenomeStream* gff3_in_stream_new(StrArray *files, /* takes ownership */
                                        bool ensure_sorting, bool be_verbose,
                                        Env *env)
{
  GenomeStream *gs = genome_stream_create(gff3_in_stream_class(),
                                          ensure_sorting, env);
  GFF3InStream *gff3_in_stream         = gff3_in_stream_cast(gs);
  gff3_in_stream->next_file              = 0;
  gff3_in_stream->files                  = files;
  gff3_in_stream->ensure_sorting         = ensure_sorting;
  gff3_in_stream->stdin_argument         = false;
  gff3_in_stream->file_is_open           = false;
  gff3_in_stream->fpin                   = NULL;
  gff3_in_stream->line_number            = 0;
  gff3_in_stream->genome_node_buffer     = queue_new(sizeof (GenomeNode*), env);
  gff3_in_stream->gff3_parser            = gff3parser_new(env);
  gff3_in_stream->last_node              = NULL;
  gff3_in_stream->be_verbose             = be_verbose;
  return gs;
}

void gff3_in_stream_set_offset(GenomeStream *gs, long offset)
{
  GFF3InStream *is = gff3_in_stream_cast(gs);
  gff3parser_set_offset(is->gff3_parser, offset);
}

int gff3_in_stream_set_offsetfile(GenomeStream *gs, Str *offsetfile, Env *env)
{
  GFF3InStream *is = gff3_in_stream_cast(gs);
  return gff3parser_set_offsetfile(is->gff3_parser, offsetfile, env);
}

int gff3_in_stream_set_chseqids(GenomeStream *gs, Str *chseqids, Env *env)
{
  GFF3InStream *is = gff3_in_stream_cast(gs);
  return gff3parser_set_chseqids(is->gff3_parser, chseqids, env);
}

GenomeStream* gff3_in_stream_new_unsorted(int num_of_files,
                                          const char **filenames,
                                          bool be_verbose, Env *env)
{
  int i;
  StrArray *files = strarray_new(env);
  for (i = 0; i < num_of_files; i++)
    strarray_add_cstr(files, filenames[i], env);
  return gff3_in_stream_new(files, false, be_verbose, env);
}

GenomeStream* gff3_in_stream_new_sorted(const char *filename, bool be_verbose,
                                        Env *env)
{
  StrArray *files = strarray_new(env);
  if (filename)
    strarray_add_cstr(files, filename, env);
  return gff3_in_stream_new(files, true, be_verbose, env);
}
