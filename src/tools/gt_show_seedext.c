/*
  Copyright (c) 2015 Joerg Winkler <joerg.winkler@studium.uni-hamburg.de>
  Copyright (c) 2015 Center for Bioinformatics, University of Hamburg

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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "core/ma_api.h"
#include "core/minmax.h"
#include "core/str_api.h"
#include "core/types_api.h"
#include "core/unused_api.h"
#include "core/encseq.h"
#include "match/revcompl.h"
#include "extended/linearalign.h"
#include "tools/gt_show_seedext.h"

typedef struct {
  bool show_alignment;
  bool seed_display;
  GtStr *filename;
} GtShowSeedextArguments;

typedef struct {
  GtUword alen, aseq, apos;
  char direction;
  GtUword blen, bseq, bpos, score, distance;
  double identity;
} GtShowSeedextCoords;

static void* gt_show_seedext_arguments_new(void)
{
  GtShowSeedextArguments *arguments = gt_calloc((size_t) 1, sizeof *arguments);
  arguments->filename = gt_str_new();
  return arguments;
}

static void gt_show_seedext_arguments_delete(void *tool_arguments)
{
  GtShowSeedextArguments *arguments = tool_arguments;
  if (arguments != NULL) {
    gt_str_delete(arguments->filename);
    gt_free(arguments);
  }
}

static GtOptionParser* gt_show_seedext_option_parser_new(void *tool_arguments)
{
  GtShowSeedextArguments *arguments = tool_arguments;
  GtOptionParser *op;
  GtOption *option;
  gt_assert(arguments);

  /* init */
  op = gt_option_parser_new("[-a] [-seed-display] -f <filename>",
                            "Parse output of a seed extension and show/verify "
                            "the alignment.");

  /* -a */
  option = gt_option_new_bool("a",
                              "show alignment",
                              &arguments->show_alignment,
                              false);
  gt_option_parser_add_option(op, option);

  /* -s */
  option = gt_option_new_bool("seed-display",
                              "display seeds if available",
                              &arguments->seed_display,
                              false);
  gt_option_parser_add_option(op, option);

  /* -f */
  option = gt_option_new_filename("f",
                                  "path to seed extension result file",
                                  arguments->filename);
  gt_option_is_mandatory(option);
  gt_option_parser_add_option(op, option);

  return op;
}

static int gt_show_seedext_arguments_check(GT_UNUSED int rest_argc,
                                       void *tool_arguments,
                                       GT_UNUSED GtError *err)
{
  GtShowSeedextArguments *arguments = tool_arguments;
  int had_err = 0;

  gt_error_check(err);
  gt_assert(arguments);

  if (arguments->filename == NULL || gt_str_length(arguments->filename) == 0) {
    gt_error_set(err,
                 "option -f requires a file name");
    had_err = -1;
  }
  return had_err;
}

/* Parse encseq input indices from -ii and -qii options in 1st line of file. */
static int gt_show_seedext_get_encseq_index(GtStr *ii,
                                            GtStr *qii,
                                            bool *mirror,
                                            const char *filename,
                                            GtError *err)
{
  FILE *file;
  int had_err = 0;

  gt_assert(ii != NULL && qii != NULL);
  file = fopen(filename, "r");
  if (file == NULL) {
    gt_error_set(err, "file %s does not exist", filename);
    had_err = -1;
  }
  if (!had_err) {
    GtStr *line_buffer = gt_str_new();
    /* read first line and evaluate tokens */
    if (gt_str_read_next_line(line_buffer,file) != EOF)
    {
      const char *tok = strtok(gt_str_get(line_buffer), " ");
      while (tok != NULL) {
        if (strcmp(tok, "-ii") == 0) {
          gt_str_set(ii, strtok(NULL, " ")); /* next token contains ii */
        } else if (strcmp(tok, "-qii") == 0) {
          gt_str_set(qii, strtok(NULL, " ")); /* next token contains qii */
        } else if (strcmp(tok, "-mirror") == 0) {
          *mirror = true; /* found -mirror option */
        }
        tok = strtok(NULL, " "); /* go to next token */
      }
      if (gt_str_length(ii) == 0UL) {
        gt_error_set(err, "need output of option string "
                     "(run gt seed_extend with -v or -verify)");
        had_err = -1;
      }
    } else
    {
      gt_error_set(err, "file %s is empty", filename);
      had_err = -1;
    }
    gt_str_delete(line_buffer);
  }
  if (file != NULL)
  {
    fclose(file);
  }
  return had_err;
}

static int gt_show_seedext_parse_extensions(const GtEncseq *aencseq,
                                            const GtEncseq *bencseq,
                                            const char *filename,
                                            bool show_alignment,
                                            bool seed_display,
                                            bool mirror,
                                            GtError *err)
{
  const GtUword width = 70;
  FILE *file;
  GtStr *line_buffer;
  GtUchar *asequence, *bsequence, *csequence;
  GtUword apos_ab = 0, bpos_ab = 0, maxseqlen = 0;
  int had_err = 0;
  GtShowSeedextCoords coords = {0, 0, 0, 'X', 0, 0, 0, 0, 0, 0.0};
  LinspaceManagement *linspace_spacemanager;
  GtScoreHandler *linspace_scorehandler;
  GtAlignment *alignment;
  GtUchar *alignment_show_buffer;
  const GtUchar *characters;
  GtUchar wildcardshow;

  gt_assert(aencseq && bencseq && filename);
  file = fopen(filename, "r");
  if (file == NULL) {
    gt_error_set(err, "file %s does not exist", filename);
    return -1;
  }
  linspace_spacemanager = gt_linspaceManagement_new();
  linspace_scorehandler = gt_scorehandler_new(0,1,0,1);;
  alignment = gt_alignment_new();
  alignment_show_buffer = gt_alignment_buffer_new(width);
  characters = gt_encseq_alphabetcharacters(aencseq);
  wildcardshow = gt_encseq_alphabetwildcardshow(aencseq);
  /* allocate buffers for alignment string and sequences */
  maxseqlen = MAX(gt_encseq_max_seq_length(aencseq),
                  gt_encseq_max_seq_length(bencseq)) + 1UL;
  line_buffer = gt_str_new();
  asequence = gt_malloc((mirror ? 3 : 2) * maxseqlen * sizeof *asequence);
  bsequence = asequence + maxseqlen;
  if (mirror) {
    csequence = bsequence + maxseqlen;
  }

  while (gt_str_read_next_line(line_buffer,file) != EOF)
  {
    const char *line_ptr = gt_str_get(line_buffer);
    /* ignore comment lines; but print seeds if -seed-display is set */
    if (line_ptr[0] != '\n')
    {
      if (line_ptr[0] == '#')
      {
        if (seed_display && strstr(line_ptr, "seed:") != NULL) {
          printf("%s", line_ptr);
        }
      } else
      {
        /* parse alignment string */
        int num = sscanf(line_ptr,
                         GT_WU " " GT_WU " " GT_WU " %c " GT_WU " " GT_WU " "
                         GT_WU " " GT_WU " " GT_WU " %lf",
                         &coords.alen, &coords.aseq, &coords.apos,
                         &coords.direction, &coords.blen, &coords.bseq,
                         &coords.bpos, &coords.score, &coords.distance,
                         &coords.identity);
        if (num == 10)
        {
          /* get sequences */
          apos_ab = gt_encseq_seqstartpos(aencseq, coords.aseq) + coords.apos;
          bpos_ab = gt_encseq_seqstartpos(bencseq, coords.bseq) + coords.bpos;
          gt_encseq_extract_encoded(aencseq,
                                    asequence,
                                    apos_ab,
                                    apos_ab + coords.alen - 1);
          gt_encseq_extract_encoded(bencseq,
                                    bsequence,
                                    bpos_ab,
                                    bpos_ab + coords.blen - 1);
          /* prepare reverse complement of 2nd sequence */
          if (mirror && coords.direction != 'F') {
            gt_copy_reverse_complement(csequence,bsequence,coords.blen);
            bsequence = csequence;
          }
          printf("%s\n",line_ptr);
          if (show_alignment) {
            GtUword edist = gt_computelinearspace_generic(linspace_spacemanager,
                                                          linspace_scorehandler,
                                                          alignment,
                                                          asequence, 0,
                                                          coords.alen,
                                                          bsequence, 0,
                                                          coords.blen);
            if (edist < coords.distance)
            {
              printf("# edist=" GT_WU "(smaller by " GT_WU ")\n",edist,
                                                   coords.distance - edist);
            }
            gt_assert(edist <= coords.distance);
            gt_alignment_show_generic(alignment_show_buffer,
                                      false,
                                      alignment,
                                      stdout,
                                      width,
                                      characters,
                                      wildcardshow);
            gt_alignment_reset(alignment);
          }
          if (mirror) {
            bsequence = asequence + maxseqlen;
          }
        }
      }
    }
    gt_str_reset(line_buffer);
  }
  fclose(file);
  gt_free(asequence);
  gt_str_delete(line_buffer);
  gt_alignment_delete(alignment);
  gt_free(alignment_show_buffer);
  gt_linspaceManagement_delete(linspace_spacemanager);
  gt_scorehandler_delete(linspace_scorehandler);
  return had_err;
}

static int gt_show_seedext_runner(GT_UNUSED int argc,
                                  GT_UNUSED const char **argv,
                                  GT_UNUSED int parsed_args,
                                  void *tool_arguments,
                                  GtError *err)
{
  GtShowSeedextArguments *arguments = tool_arguments;
  GtEncseq *aencseq = NULL, *bencseq = NULL;
  GtStr *ii = gt_str_new(),
        *qii = gt_str_new();
  bool mirror = false;
  int had_err = 0;

  gt_error_check(err);
  gt_assert(arguments);

  /* Parse option string in first line of file specified by filename. */
  had_err = gt_show_seedext_get_encseq_index(ii,
                                             qii,
                                             &mirror,
                                             gt_str_get(arguments->filename),
                                             err);
  printf("# file %s: mirror %sabled\n",
         gt_str_get(arguments->filename), mirror ? "en" : "dis");

  /* Load encseqs */
  if (!had_err) {
    GtEncseqLoader *encseq_loader = gt_encseq_loader_new();
    gt_encseq_loader_enable_autosupport(encseq_loader);
    aencseq = gt_encseq_loader_load(encseq_loader, gt_str_get(ii), err);
    if (aencseq == NULL) {
      had_err = -1;
    }

    if (!had_err) {
      if (gt_str_length(qii) != 0) {
        bencseq = gt_encseq_loader_load(encseq_loader, gt_str_get(qii), err);
      } else {
        bencseq = gt_encseq_ref(aencseq);
      }
      if (bencseq == NULL) {
        had_err = -1;
        gt_encseq_delete(aencseq);
      }
    }
    gt_encseq_loader_delete(encseq_loader);
  }
  gt_str_delete(ii);
  gt_str_delete(qii);

  /* Parse seed extensions. */
  if (!had_err) {
    had_err = gt_show_seedext_parse_extensions(aencseq,
                                               bencseq,
                                               gt_str_get(arguments->filename),
                                               arguments->show_alignment,
                                               arguments->seed_display,
                                               mirror,
                                               err);
    gt_encseq_delete(aencseq);
    gt_encseq_delete(bencseq);
  }
  return had_err;
}

GtTool* gt_show_seedext(void)
{
  return gt_tool_new(gt_show_seedext_arguments_new,
                     gt_show_seedext_arguments_delete,
                     gt_show_seedext_option_parser_new,
                     gt_show_seedext_arguments_check,
                     gt_show_seedext_runner);
}
