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
#include "core/complement.h"
#include "core/encseq_api.h"
#include "core/ma.h"
#include "core/ma_api.h"
#include "core/minmax.h"
#include "core/str_api.h"
#include "core/types_api.h"
#include "core/unused_api.h"
#include "extended/squarealign.h"
#include "match/test-pairwise.h"
#include "tools/gt_show_seedext.h"

typedef struct {
  bool show_alignment;
  bool seed_display;
  GtStr *filename;
} GtShowSeedextArguments;

typedef struct {
  GtUword alen, aseq, apos;
  char direction;
  GtUword blen, bseq, bpos;
  GtUword score, distance;
  double correlation;
} GtShowSeedextAlignment;

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
  const GtUword maxlinelength = 1000;
  FILE *file;
  int had_err = 0;

  gt_assert(ii && qii);

  file = fopen(filename, "r");
  if (file == NULL) {
    gt_error_set(err, "file %s does not exist", filename);
    had_err = -1;
  }
  if (!had_err) {
    char *buffer = gt_malloc(maxlinelength * sizeof *buffer);
    /* read first line and evaluate tokens */
    if (fgets(buffer, maxlinelength, file)) {
      char *tok = strtok(buffer, " ");
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
    } else {
      gt_error_set(err, "file %s is empty", filename);
      had_err = -1;
    }
    gt_free(buffer);
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
  const GtUword maxlinelength = 10000;
  FILE *file;
  char *buffer, *asequence, *bsequence, *csequence;
  GtUword apos_ab = 0, bpos_ab = 0, maxseqlen = 0, edist = 0;
  int num, had_err = 0;
  GtShowSeedextAlignment alignment = {0, 0, 0, 'X', 0, 0, 0, 0, 0, 0.0};

  gt_assert(aencseq && bencseq && filename);

  file = fopen(filename, "r");
  if (file == NULL) {
    gt_error_set(err, "file %s does not exist", filename);
    fclose(file);
    return -1;
  }

  /* allocate buffers for alignment string and sequences */
  maxseqlen = MAX(gt_encseq_max_seq_length(aencseq),
                  gt_encseq_max_seq_length(bencseq)) + 1UL;
  buffer = gt_malloc(maxlinelength * sizeof *buffer);
  asequence = gt_malloc((mirror ? 3 : 2) * maxseqlen * sizeof *asequence);
  bsequence = asequence + maxseqlen;
  if (mirror) {
    csequence = bsequence + maxseqlen;
  }

  while (fgets(buffer, maxlinelength, file)) {
    /* ignore comment lines; but print seeds if -seed-display is set */
    if (buffer[0] == '#') {
      if (seed_display && strstr(buffer, "seed:") != NULL) {
        printf("%s", buffer);
      }
      continue;
    }

    /* parse alignment string */
    num = sscanf(buffer,
                 GT_WU" "GT_WU" "GT_WU" %c "GT_WU" "GT_WU" "GT_WU" "GT_WU" "
                 GT_WU" %lf",
                 &alignment.alen, &alignment.aseq, &alignment.apos,
                 &alignment.direction, &alignment.blen, &alignment.bseq,
                 &alignment.bpos, &alignment.score, &alignment.distance,
                 &alignment.correlation);
    if (num != 10) {
      printf("alignment parse failed: %s", buffer);
      continue;
    }

    /* get sequences */
    apos_ab = alignment.apos + gt_encseq_seqstartpos(aencseq, alignment.aseq);
    bpos_ab = alignment.bpos + gt_encseq_seqstartpos(bencseq, alignment.bseq);
    gt_encseq_extract_decoded(aencseq,
                              asequence,
                              apos_ab,
                              apos_ab + alignment.alen - 1);
    gt_encseq_extract_decoded(bencseq,
                              bsequence,
                              bpos_ab,
                              bpos_ab + alignment.blen - 1);
    asequence[alignment.alen] = bsequence[alignment.blen] = '\0';
    /* prepare reverse complement of 2nd sequence */
    if (mirror && alignment.direction != 'F') {
      char *idx;
      for (idx = csequence; idx < csequence + alignment.blen; idx++) {
        gt_complement(idx, bsequence[csequence + alignment.blen -idx-1], NULL);
      }
      csequence[alignment.blen] = '\0';
      bsequence = csequence;
    }
    edist = gt_computegreedyunitedist((const GtUchar *)asequence,
                                      alignment.alen,
                                      (const GtUchar *)bsequence,
                                      alignment.blen);
    printf("# identity = %5.2lf\n",
           100.0 - 200.0 * edist / (alignment.alen + alignment.blen));

    if (show_alignment) {
      gt_print_edist_alignment((const GtUchar *)asequence, 0,
                               alignment.alen,
                               (const GtUchar *)bsequence, 0,
                               alignment.blen);
    }
    if (mirror) {
      bsequence = asequence + maxseqlen;
    }
  }
  fclose(file);
  gt_free(asequence);
  gt_free(buffer);
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
  printf("File %s: mirror %s\n",
         gt_str_get(arguments->filename), mirror ? "enabled" : "disabled");

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
