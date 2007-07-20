/*
  Copyright (c) 2006-2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2006-2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#include <assert.h>
#include <string.h>
#include "libgtcore/fileutils.h"
#include "libgtcore/progressbar.h"
#include "libgtext/genome_stream_rep.h"
#include "libgtext/gff3_in_stream.h"
#include "libgtext/gff3_parser.h"

struct GFF3InStream
{
  const GenomeStream parent_instance;
  int next_file;
  Array *files; /* contains char* to filenames */
  bool ensure_sorting,
       stdin_argument,
       non_stdin_file_is_open,
       be_verbose;
  FILE *fpin;
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
    if (is->fpin == NULL) {
      if (array_size(is->files) && is->next_file == array_size(is->files))
        break;
      assert(!is->non_stdin_file_is_open);
      if (array_size(is->files)) {
        if (strcmp(*(char**) array_get(is->files, is->next_file), "-") == 0) {
          if (is->stdin_argument) {
            env_error_set(env,
                          "multiple specification of argument file \"-\"\n");
            had_err = -1;
            break;
          }
          is->fpin = stdin;
          is->stdin_argument = true;
        }
        else {
          is->fpin = env_fa_xfopen(env, *(char**)
                                   array_get(is->files, is->next_file), "r");
          is->non_stdin_file_is_open = true;
        }
        is->next_file++;
      }
      else is->fpin = stdin;
      is->line_number = 0;

      if (is->be_verbose) {
        printf("processing file \"%s\"\n", array_size(is->files)
               ? *(char**) array_get(is->files, is->next_file-1) : "stdin");
      }
      if (is->non_stdin_file_is_open && is->be_verbose)
        progressbar_start(&is->line_number, file_number_of_lines(is->fpin));
    }

    assert(is->fpin); /* file is open */

    had_err = gff3parser_parse_genome_nodes(&status_code, is->gff3_parser,
                                            is->genome_node_buffer,
                                            array_size(is->files)
                                            ? *(char**)
                                              array_get(is->files,
                                                        is->next_file-1)
                                            : "stdin",
                                            &is->line_number, is->fpin, env);
    if (had_err)
      break;

    if (status_code == EOF) {
      /* end of current file */
      if (is->non_stdin_file_is_open) {
        assert(is->fpin);
        if (is->be_verbose) progressbar_stop();
        env_fa_xfclose(is->fpin, env);
        is->non_stdin_file_is_open = false;
      }
      is->fpin = NULL;
      gff3parser_reset(is->gff3_parser, env);
      if (!array_size(is->files))
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
          assert(array_size(is->files) == 0 || array_size(is->files) == 1);
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
  array_delete(gff3_in_stream->files, env);
  while (queue_size(gff3_in_stream->genome_node_buffer)) {
    genome_node_rec_delete(*(GenomeNode**)
                           queue_get(gff3_in_stream->genome_node_buffer), env);
  }
  queue_delete(gff3_in_stream->genome_node_buffer, env);
  gff3parser_delete(gff3_in_stream->gff3_parser, env);
  genome_node_delete(gff3_in_stream->last_node, env);
  if (gff3_in_stream->non_stdin_file_is_open)
    env_fa_xfclose(gff3_in_stream->fpin, env);
}

const GenomeStreamClass* gff3_in_stream_class(void)
{
  static const GenomeStreamClass gsc = { sizeof (GFF3InStream),
                                         gff3_in_stream_next_tree,
                                         gff3_in_stream_free };
  return &gsc;
}

static GenomeStream* gff3_in_stream_new(Array *files, /* takes ownership */
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
  gff3_in_stream->non_stdin_file_is_open = false;
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

GenomeStream* gff3_in_stream_new_unsorted(int num_of_files,
                                          const char **filenames,
                                          bool be_verbose, Env *env)
{
  int i;
  Array *files = array_new(sizeof (char*), env);
  for (i = 0; i < num_of_files; i++)
    array_add(files, filenames[i], env);
  return gff3_in_stream_new(files, false, be_verbose, env);
}

GenomeStream* gff3_in_stream_new_sorted(const char *filename, bool be_verbose,
                                        Env *env)
{
  Array *files = array_new(sizeof (char*), env);
  if (filename)
    array_add(files, filename, env);
  return gff3_in_stream_new(files, true, be_verbose, env);
}
