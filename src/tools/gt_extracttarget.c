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
#include "core/cstr_api.h"
#include "core/fasta.h"
#include "core/ma.h"
#include "core/option_api.h"
#include "core/splitter.h"
#include "core/unused_api.h"
#include "extended/feature_node.h"
#include "extended/feature_node_iterator_api.h"
#include "extended/genome_node.h"
#include "extended/gff3_escaping.h"
#include "extended/gff3_in_stream.h"
#include "extended/string_matching.h"
#include "tools/gt_extracttarget.h"

typedef struct {
  GtStrArray *seqfiles;
} ExtractTargetArguments;

static void* gt_extracttarget_arguments_new(void)
{
  ExtractTargetArguments *arguments = gt_calloc(1, sizeof *arguments);
  arguments->seqfiles = gt_str_array_new();
  return arguments;
}

static void gt_extracttarget_arguments_delete(void *tool_arguments)
{
  ExtractTargetArguments *arguments = tool_arguments;
  if (!arguments) return;
  gt_str_array_delete(arguments->seqfiles);
  gt_free(arguments);
}

static GtOptionParser* gt_extracttarget_option_parser_new(void *tool_arguments)
{
  ExtractTargetArguments *arguments = tool_arguments;
  GtOptionParser *op;
  GtOption *o;
  gt_assert(arguments);

  /* init */
  op = gt_option_parser_new("[option ...] -seqfile sequence_file GFF3_file",
                         "Extract target sequences given in GFF3_file from "
                         "sequence_file.");

  /* -seqfile */
  o = gt_option_new_filename_array("seqfiles",
                                   "set the sequence file from which "
                                   "to extract the features",
                                   arguments->seqfiles);
  gt_option_is_mandatory(o);
  gt_option_parser_add_option(op, o);

  gt_option_parser_set_min_max_args(op, 1, 1);

  return op;
}

typedef struct {
  GtBioseq *bioseq;
  unsigned long seqnum;
} TargetInfo;

static bool show_target(GT_UNUSED unsigned long pos, void *data)
{
  TargetInfo *ti = data;
  char *seq;
  gt_assert(ti);
  seq = gt_bioseq_get_sequence(ti->bioseq, ti->seqnum);
  gt_fasta_show_entry(gt_bioseq_get_description(ti->bioseq, ti->seqnum),
                      seq,
                      gt_bioseq_get_sequence_length(ti->bioseq, ti->seqnum), 0,
                      NULL);
  gt_free(seq);
  return true;
}

static int extracttarget_from_seqfiles(const char *target,
                                       GtStrArray *seqfiles,
                                       GtError *err)
{
  GtStr *unescaped_target;
  char *escaped_target;
  GtSplitter *splitter;
  unsigned long i;
  int had_err = 0;
  gt_error_check(err);
  gt_assert(target && seqfiles);
  splitter = gt_splitter_new();
  unescaped_target = gt_str_new();
  escaped_target = gt_cstr_dup(target);
  gt_splitter_split(splitter, escaped_target, strlen(escaped_target), ',');
  for (i = 0; !had_err && i < gt_splitter_size(splitter); i++) {
    GtSplitter *blank_splitter;
    char *token = gt_splitter_get_token(splitter, i);
    blank_splitter = gt_splitter_new();
    gt_splitter_split(blank_splitter, token, strlen(token), ' ');
    had_err = gt_gff3_unescape(unescaped_target,
                               gt_splitter_get_token(blank_splitter, 0),
                               strlen(gt_splitter_get_token(blank_splitter, 0)),
                               err);
    if (!had_err) {
      unsigned long j;
      for (j = 0; j < gt_str_array_size(seqfiles); j++) {
        unsigned long k;
        GtBioseq *bioseq;
        if (!(bioseq =  gt_bioseq_new(gt_str_array_get(seqfiles, j), err))) {
          had_err = -1;
          break;
        }
        for (k = 0; k < gt_bioseq_number_of_sequences(bioseq); k++) {
          TargetInfo target_info;
          const char *desc = gt_bioseq_get_description(bioseq, k);
          target_info.bioseq = bioseq;
          target_info.seqnum = k;
          gt_string_matching_bmh(desc, strlen(desc),
                                 gt_str_get(unescaped_target),
                                 gt_str_length(unescaped_target), show_target,
                                 &target_info);
        }
        gt_bioseq_delete(bioseq);
      }
    }
    gt_splitter_delete(blank_splitter);
  }
  gt_free(escaped_target);
  gt_str_delete(unescaped_target);
  gt_splitter_delete(splitter);
  return had_err;
}

static int extracttarget_from_node(GtGenomeNode *gn, GtStrArray *seqfiles,
                                   GtError *err)
{
  GtFeatureNodeIterator *fni;
  int had_err = 0;
  gt_error_check(err);
  gt_assert(gn && seqfiles);
  /* XXX */
  if (gt_genome_node_cast(gt_feature_node_class(), gn)) {
    const char *target;
    GtFeatureNode *child;
    fni = gt_feature_node_iterator_new(gt_feature_node_cast(gn));
    while (!had_err && /* XXX remove cast */
           (child = (GtFeatureNode*) gt_feature_node_iterator_next(fni))) {
      if ((target = gt_feature_node_get_attribute(child, "Target")))
        had_err = extracttarget_from_seqfiles(target, seqfiles, err);
    }
    gt_feature_node_iterator_delete(fni);
  }
  return had_err;
}

static int gt_extracttarget_runner(GT_UNUSED int argc, const char **argv,
                                   int parsed_args, void *tool_arguments,
                                   GtError *err)
{
  ExtractTargetArguments *arguments = tool_arguments;
  GtNodeStream *gff3_in_stream;
  GtGenomeNode *gn;
  int had_err;

  gt_error_check(err);
  gt_assert(arguments);

  gff3_in_stream = gt_gff3_in_stream_new_unsorted(1, argv + parsed_args);

  while (!(had_err = gt_node_stream_next(gff3_in_stream, &gn, err)) && gn) {
    had_err = extracttarget_from_node(gn, arguments->seqfiles, err);
    gt_genome_node_delete(gn);
  }

  gt_node_stream_delete(gff3_in_stream);

  return had_err;
}

GtTool* gt_extracttarget(void)
{
  return gt_tool_new(gt_extracttarget_arguments_new,
                     gt_extracttarget_arguments_delete,
                     gt_extracttarget_option_parser_new,
                     NULL,
                     gt_extracttarget_runner);
}
