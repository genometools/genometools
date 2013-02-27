/*
  Copyright (c) 2011 Sascha Steinbiss <steinbiss@zbh.uni-hamburg.de>
  Copyright (c) 2011 Center for Bioinformatics, University of Hamburg

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

#include "core/codon_api.h"
#include "core/codon_iterator_simple_api.h"
#include "core/cstr_api.h"
#include "core/fasta.h"
#include "core/ma.h"
#include "core/output_file_api.h"
#include "core/seq_iterator_sequence_buffer.h"
#include "core/sequence_buffer.h"
#include "core/translator.h"
#include "core/unused_api.h"
#include "core/warning_api.h"
#include "extended/reverse_api.h"
#include "tools/gt_seqtranslate.h"

typedef struct {
  GtOutputFileInfo *ofi;
  GtFile *outfp;
  unsigned long fasta_width;
  bool reverse;
} GtTranslateArguments;

static void* gt_seqtranslate_arguments_new(void)
{
  GtTranslateArguments *arguments = gt_calloc(1, sizeof *arguments);
  arguments->ofi = gt_output_file_info_new();
  return arguments;
}

static void gt_seqtranslate_arguments_delete(void *tool_arguments)
{
  GtTranslateArguments *arguments = tool_arguments;
  if (!arguments) return;
  gt_output_file_info_delete(arguments->ofi);
  gt_file_delete(arguments->outfp);
  gt_free(arguments);
}

static GtOptionParser* gt_seqtranslate_option_parser_new(void *tool_arguments)
{
  GtTranslateArguments *arguments = tool_arguments;
  GtOptionParser *op;
  GtOption *option;
  gt_assert(arguments);

  /* init */
  op = gt_option_parser_new("[option ...] [files]",
                            "Translates a nucleotide sequence into a protein "
                            "sequence.");

  option = gt_option_new_bool("reverse", "also translate reverse complements",
                              &arguments->reverse, true);
  gt_option_parser_add_option(op, option);

  option = gt_option_new_ulong("fastawidth", "width of the FASTA output",
                              &arguments->fasta_width, 60);
  gt_option_parser_add_option(op, option);

  gt_output_file_register_options(op, &arguments->outfp, arguments->ofi);
  gt_option_parser_set_min_args(op, 1);

  return op;
}

static int gt_seqtranslate_arguments_check(GT_UNUSED int rest_argc,
                                        GT_UNUSED void *tool_arguments,
                                        GT_UNUSED GtError *err)
{
  GT_UNUSED GtTranslateArguments *arguments = tool_arguments;
  int had_err = 0;
  gt_error_check(err);
  gt_assert(arguments);

  return had_err;
}

static int gt_seqtranslate_do_translation(GtTranslateArguments *arguments,
                                       const char *sequence,
                                       unsigned long length,
                                       const char *desc,
                                       GtStr **translations,
                                       bool rev,
                                       GtError *err)
{
  GtTranslator *tr;
  GT_UNUSED GtTranslatorStatus trst;
  GtCodonIterator *ci;
  char translated;
  int had_err = 0;
  GtStr *str;
  unsigned int frame,
               i;

  ci = gt_codon_iterator_simple_new(sequence, length, err);
  tr = gt_translator_new(ci);
  trst = gt_translator_next(tr, &translated, &frame, err);
  while (trst == GT_TRANSLATOR_OK) {
    gt_str_append_char(translations[frame], translated);
    trst = gt_translator_next(tr, &translated, &frame, err);
  }
  gt_codon_iterator_delete(ci);
  gt_translator_delete(tr);
  if (trst == GT_TRANSLATOR_ERROR)
    return -1;
  str = gt_str_new();
  for (i = 0; i < 3; i++) {
    if (gt_str_length(translations[i]) > 0) {
      gt_str_append_cstr(str, desc);
      gt_str_append_cstr(str, " (");
      gt_str_append_ulong(str, i+1);
      gt_str_append_cstr(str, rev ? "-" : "+");
      gt_str_append_cstr(str, ")");
      gt_fasta_show_entry(gt_str_get(str), gt_str_get(translations[i]),
                          gt_str_length(translations[i]),
                          arguments->fasta_width, arguments->outfp);
      gt_str_reset(translations[i]);
      gt_str_reset(str);
    }
  }
  gt_str_delete(str);
  return had_err;
}

static int gt_seqtranslate_runner(int argc, const char **argv, int parsed_args,
                              void *tool_arguments, GT_UNUSED GtError *err)
{
  GtTranslateArguments *arguments = tool_arguments;
  GtSeqIterator *si = NULL;
  GtSequenceBuffer *sb = NULL;
  GtStrArray *infiles;
  int had_err = 0,
      rval,
      i;
  GtStr *translations[3];
  translations[0] = gt_str_new();
  translations[1] = gt_str_new();
  translations[2] = gt_str_new();

  gt_error_check(err);
  gt_assert(arguments);

  infiles = gt_str_array_new();
  for (i = parsed_args; i < argc; i++) {
    gt_str_array_add_cstr(infiles, argv[i]);
  }
  sb = gt_sequence_buffer_new_guess_type(infiles, err);
  if (!sb)
    had_err = -1;
  if (!had_err) {
    si = gt_seq_iterator_sequence_buffer_new_with_buffer(sb);
    if (!si)
      had_err = -1;
  }
  if (!had_err) {
    char *desc;
    const GtUchar *sequence;
    unsigned long len;
    while (!had_err && (rval = gt_seq_iterator_next(si,
                                                   &sequence,
                                                   &len, &desc, err))) {
      if (rval < 0) {
        had_err = -1;
        break;
      }
      if (len < GT_CODON_LENGTH) {
        gt_warning("sequence '%s' is shorter than codon length of %d, skipping",
                   desc, GT_CODON_LENGTH);
      } else {
        had_err = gt_seqtranslate_do_translation(arguments, (char*) sequence,
                                                 len, desc,
                                                 translations, false, err);
        if (!had_err && arguments->reverse) {
          char *revseq = gt_cstr_dup_nt((char*) sequence, len);
          had_err = gt_reverse_complement(revseq, len, err);
          if (!had_err) {
            had_err = gt_seqtranslate_do_translation(arguments, revseq, len,
                                                  desc, translations, true,
                                                  err);
          }
          gt_free(revseq);
        }
      }
    }
  }
  gt_str_delete(translations[0]);
  gt_str_delete(translations[1]);
  gt_str_delete(translations[2]);
  gt_str_array_delete(infiles);
  gt_seq_iterator_delete(si);
  gt_sequence_buffer_delete(sb);
  return had_err;
}

GtTool* gt_seqtranslate(void)
{
  return gt_tool_new(gt_seqtranslate_arguments_new,
                     gt_seqtranslate_arguments_delete,
                     gt_seqtranslate_option_parser_new,
                     gt_seqtranslate_arguments_check,
                     gt_seqtranslate_runner);
}
