/*
  Copyright (c) 2006-2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2006-2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#include <assert.h>
#include <string.h>
#include "error.h"
#include "fileutils.h"
#include "genome_stream_rep.h"
#include "gff3_in_stream.h"
#include "gff3_parser.h"
#include "progressbar.h"
#include "queue.h"
#include "xansi.h"

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

static int gff3_in_stream_next_tree(GenomeStream *gs, GenomeNode **gn,
                                    /*@unused@*/ Log *l, Error *err)
{
  GFF3InStream *is = gff3_in_stream_cast(gs);
  unsigned long i;
  int has_err = 0, status_code;

  error_check(err);

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
            error_set(err, "multiple specification of argument file \"-\"\n");
            has_err = -1;
            break;
          }
          is->fpin = stdin;
          is->stdin_argument = true;
        }
        else {
          is->fpin = xfopen(*(char**) array_get(is->files, is->next_file), "r");
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

    has_err = gff3_parse_genome_nodes(&status_code, is->gff3_parser,
                                      is->genome_node_buffer,
                                      array_size(is->files)
                                      ? *(char**) array_get(is->files,
                                                            is->next_file-1)
                                      : "stdin",
                                      &is->line_number, is->fpin, err);
    if (has_err)
      break;

    if (status_code == EOF) {
      /* end of current file */
      if (is->non_stdin_file_is_open) {
        assert(is->fpin);
        if (is->be_verbose) progressbar_stop();
        xfclose(is->fpin);
        is->non_stdin_file_is_open = false;
      }
      is->fpin = NULL;
      gff3_reset(is->gff3_parser);
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
                                                       i))) {
          assert(is->last_node);
          /* a sorted stream can have at most one input file */
          assert(array_size(is->files) == 0 || array_size(is->files) == 1);
          error_set(err,
                    "the file %s is not sorted (example: line %lu and %lu)",
                    genome_node_get_filename(is->last_node),
                    genome_node_get_line_number(is->last_node),
                    genome_node_get_line_number(*(GenomeNode**)
                                    queue_get_elem(is->genome_node_buffer, i)));
          has_err = -1;
          break;
        }
      }
    }
    if (!has_err)
      *gn = *(GenomeNode**) queue_get(is->genome_node_buffer);
    return has_err;
  }
  *gn = NULL;
  return has_err;
}

static void gff3_in_stream_free(GenomeStream *gs)
{
  GFF3InStream *gff3_in_stream = gff3_in_stream_cast(gs);
  array_delete(gff3_in_stream->files);
  queue_delete(gff3_in_stream->genome_node_buffer);
  gff3_delete(gff3_in_stream->gff3_parser);
  genome_node_delete(gff3_in_stream->last_node);
}

const GenomeStreamClass* gff3_in_stream_class(void)
{
  static const GenomeStreamClass gsc = { sizeof (GFF3InStream),
                                         gff3_in_stream_next_tree,
                                         gff3_in_stream_free };
  return &gsc;
}

static GenomeStream* gff3_in_stream_new(Array *files, /* takes ownership */
                                        bool ensure_sorting, bool be_verbose)
{
  GenomeStream *gs = genome_stream_create(gff3_in_stream_class(),
                                          ensure_sorting);
  GFF3InStream *gff3_in_stream         = gff3_in_stream_cast(gs);
  gff3_in_stream->next_file              = 0;
  gff3_in_stream->files                  = files;
  gff3_in_stream->ensure_sorting         = ensure_sorting;
  gff3_in_stream->stdin_argument         = false;
  gff3_in_stream->non_stdin_file_is_open = false;
  gff3_in_stream->fpin                   = NULL;
  gff3_in_stream->line_number            = 0;
  gff3_in_stream->genome_node_buffer     = queue_new(sizeof (GenomeNode*));
  gff3_in_stream->gff3_parser            = gff3_new();
  gff3_in_stream->last_node              = NULL;
  gff3_in_stream->be_verbose             = be_verbose;
  return gs;
}

void gff3_in_stream_set_offset(GenomeStream *gs, long offset)
{
  GFF3InStream *is = gff3_in_stream_cast(gs);
  gff3_set_offset(is->gff3_parser, offset);
}

GenomeStream* gff3_in_stream_new_unsorted(int num_of_files,
                                           char **filenames,
                                           bool be_verbose)
{
  int i;
  Array *files = array_new(sizeof (char*));
  for (i = 0; i < num_of_files; i++)
    array_add(files, filenames[i]);
  return gff3_in_stream_new(files, false, be_verbose);
}

GenomeStream* gff3_in_stream_new_sorted(char *filename,
                                         bool be_verbose)
{
  Array *files = array_new(sizeof (char*));
  if (filename)
    array_add(files, filename);
  return gff3_in_stream_new(files, true, be_verbose);
}
