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
#include "core/bioseq.h"
#include "core/cstr.h"
#include "core/fasta.h"
#include "core/ma.h"
#include "core/option.h"
#include "core/splitter.h"
#include "core/unused_api.h"
#include "extended/genome_node.h"
#include "extended/genome_node_iterator.h"
#include "extended/gff3_escaping.h"
#include "extended/gff3_in_stream.h"
#include "extended/string_matching.h"
#include "tools/gt_extracttarget.h"

typedef struct {
  GT_StrArray *seqfiles;
} ExtractTargetArguments;

static void* gt_extracttarget_arguments_new(void)
{
  ExtractTargetArguments *arguments = gt_calloc(1, sizeof *arguments);
  arguments->seqfiles = gt_strarray_new();
  return arguments;
}

static void gt_extracttarget_arguments_delete(void *tool_arguments)
{
  ExtractTargetArguments *arguments = tool_arguments;
  if (!arguments) return;
  gt_strarray_delete(arguments->seqfiles);
  gt_free(arguments);
}

static OptionParser* gt_extracttarget_option_parser_new(void *tool_arguments)
{
  ExtractTargetArguments *arguments = tool_arguments;
  OptionParser *op;
  Option *o;
  assert(arguments);

  /* init */
  op = option_parser_new("[option ...] -seqfile sequence_file GFF3_file",
                         "Extract target sequences given in GFF3_file from "
                         "sequence_file.");

  /* -seqfile */
  o = option_new_filenamearray("seqfiles", "set the sequence file from which "
                               "to extract the features", arguments->seqfiles);
  option_is_mandatory(o);
  option_parser_add_option(op, o);

  option_parser_set_min_max_args(op, 1, 1);

  return op;
}

typedef struct {
  GT_Bioseq *bioseq;
  unsigned long seqnum;
} TargetInfo;

static bool show_target(GT_UNUSED unsigned long pos, void *data)
{
  TargetInfo *ti = data;
  assert(ti);
  fasta_show_entry(gt_bioseq_get_description(ti->bioseq, ti->seqnum),
                   gt_bioseq_get_sequence(ti->bioseq, ti->seqnum),
                   gt_bioseq_get_sequence_length(ti->bioseq, ti->seqnum), 0);
  return true;
}

static int extracttarget_from_seqfiles(const char *target,
                                       GT_StrArray *seqfiles,
                                       GT_Error *err)
{
  GT_Str *unescaped_target;
  char *escaped_target;
  Splitter *splitter;
  unsigned long i;
  int had_err = 0;
  gt_error_check(err);
  assert(target && seqfiles);
  splitter = splitter_new();
  unescaped_target = gt_str_new();
  escaped_target = gt_cstr_dup(target);
  splitter_split(splitter, escaped_target, strlen(escaped_target), ',');
  for (i = 0; !had_err && i < splitter_size(splitter); i++) {
    Splitter *blank_splitter;
    char *token = splitter_get_token(splitter, i);
    blank_splitter = splitter_new();
    splitter_split(blank_splitter, token, strlen(token), ' ');
    had_err = gff3_unescape(unescaped_target,
                            splitter_get_token(blank_splitter, 0),
                            strlen(splitter_get_token(blank_splitter, 0)), err);
    if (!had_err) {
      unsigned long j;
      for (j = 0; j < gt_strarray_size(seqfiles); j++) {
        unsigned long k;
        GT_Bioseq *bioseq;
        if (!(bioseq =  gt_bioseq_new(gt_strarray_get(seqfiles, j), err))) {
          had_err = -1;
          break;
        }
        for (k = 0; k < gt_bioseq_number_of_sequences(bioseq); k++) {
          TargetInfo target_info;
          const char *desc = gt_bioseq_get_description(bioseq, k);
          target_info.bioseq = bioseq;
          target_info.seqnum = k;
          string_matching_bmh(desc, strlen(desc), gt_str_get(unescaped_target),
                              gt_str_length(unescaped_target), show_target,
                              &target_info);
        }
        gt_bioseq_delete(bioseq);
      }
    }
    splitter_delete(blank_splitter);
  }
  gt_free(escaped_target);
  gt_str_delete(unescaped_target);
  splitter_delete(splitter);
  return had_err;
}

static int extracttarget_from_node(GT_GenomeNode *gn, GT_StrArray *seqfiles,
                                   GT_Error *err)
{
  GT_GenomeNodeIterator *gni;
  int had_err = 0;
  gt_error_check(err);
  assert(gn && seqfiles);
  if (gt_genome_node_cast(gt_genome_feature_class(), gn)) {
    const char *target;
    GT_GenomeNode *child;
    gni = gt_genome_node_iterator_new(gn);
    while (!had_err && (child = gt_genome_node_iterator_next(gni))) {
      if ((target = gt_genome_feature_get_attribute(child, "Target")))
        had_err = extracttarget_from_seqfiles(target, seqfiles, err);
    }
    gt_genome_node_iterator_delete(gni);
  }
  return had_err;
}

static int gt_extracttarget_runner(GT_UNUSED int argc, const char **argv,
                                   int parsed_args, void *tool_arguments,
                                   GT_Error *err)
{
  ExtractTargetArguments *arguments = tool_arguments;
  GenomeStream *gff3_in_stream;
  GT_GenomeNode *gn;
  int had_err;

  gt_error_check(err);
  assert(arguments);

  gff3_in_stream = gff3_in_stream_new_unsorted(1, argv + parsed_args, false,
                                               false);

  while (!(had_err = genome_stream_next_tree(gff3_in_stream, &gn, err)) && gn) {
    had_err = extracttarget_from_node(gn, arguments->seqfiles, err);
    gt_genome_node_rec_delete(gn);
  }

  genome_stream_delete(gff3_in_stream);

  return had_err;
}

Tool* gt_extracttarget(void)
{
  return tool_new(gt_extracttarget_arguments_new,
                  gt_extracttarget_arguments_delete,
                  gt_extracttarget_option_parser_new,
                  NULL,
                  gt_extracttarget_runner);
}
