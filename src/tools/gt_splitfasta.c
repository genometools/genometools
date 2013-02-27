/*
  Copyright (c) 2003-2010 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2003-2008 Center for Bioinformatics, University of Hamburg

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

#include "core/bioseq.h"
#include "core/fa.h"
#include "core/fasta.h"
#include "core/fileutils_api.h"
#include "core/ma.h"
#include "core/option_api.h"
#include "core/output_file.h"
#include "core/unused_api.h"
#include "core/xansi_api.h"
#include "tools/gt_splitfasta.h"

typedef struct {
  unsigned long max_filesize_in_MB, width;
  unsigned int num_files;
  GtStr *splitdesc;
  bool force;
} SplitfastaArguments;

static void* gt_splitfasta_arguments_new(void)
{
  SplitfastaArguments *arguments = gt_calloc(1, sizeof *arguments);
  arguments->splitdesc = gt_str_new();
  return arguments;
}

void gt_splitfasta_arguments_delete(void *tool_arguments)
{
  SplitfastaArguments *arguments = tool_arguments;
  if (!arguments) return;
  gt_str_delete(arguments->splitdesc);
  gt_free(arguments);
}

static GtOptionParser* gt_splitfasta_option_parser_new(void *tool_arguments)
{
  SplitfastaArguments *arguments = tool_arguments;
  GtOptionParser *op;
  GtOption *targetsize_option, *splitdesc_option, *o, *numfiles_option;
  gt_assert(arguments);
  op = gt_option_parser_new("[option ...] fastafile","Split the supplied fasta "
                         "file.");

  numfiles_option = gt_option_new_uint_min("numfiles",
                                           "set the number of target files "
                                           "",
                                           &arguments->num_files, 0,
                                           1);
  gt_option_parser_add_option(op, numfiles_option);

  targetsize_option = gt_option_new_ulong_min("targetsize",
                                              "set the target file "
                                              "size in MB",
                                              &arguments->max_filesize_in_MB,
                                              50, 1);
  gt_option_parser_add_option(op, targetsize_option);
  splitdesc_option = gt_option_new_string("splitdesc",
                                          "put every fasta entry in "
                                          "a separate file named by its "
                                          "description in the given directory",
                                          arguments->splitdesc, NULL);
  gt_option_parser_add_option(op, splitdesc_option);
  gt_option_exclude(targetsize_option, splitdesc_option);
  gt_option_exclude(numfiles_option, splitdesc_option);
  gt_option_exclude(numfiles_option, targetsize_option);
  o = gt_option_new_width(&arguments->width);
  gt_option_parser_add_option(op, o);
  o = gt_option_new_bool(GT_FORCE_OPT_CSTR, "force writing to output file",
                         &arguments->force, false);
  gt_option_parser_add_option(op, o);
  gt_option_parser_set_min_max_args(op, 1, 1);
  return op;
}

static unsigned long buf_contains_separator(char *buf, int offset,
                                            int read_bytes)
{
  char *cc;
  gt_assert(buf && offset < read_bytes);
  for (cc = buf + offset; cc < buf + read_bytes; cc++) {
    if (*cc == '>')
      return cc - buf + 1;
  }
  return 0;
}

static int split_description(const char *filename, GtStr *splitdesc,
                             unsigned long width, bool force, GtError *err)
{
  unsigned long i;
  GtBioseq *bioseq;
  GtStr *descname;
  int had_err = 0;
  gt_error_check(err);
  gt_assert(filename && splitdesc && gt_str_length(splitdesc));

  descname = gt_str_new();
  if (!(bioseq = gt_bioseq_new(filename, err)))
    had_err = -1;

  for (i = 0; !had_err && i < gt_bioseq_number_of_sequences(bioseq); i++) {
    GtFile *outfp;
    char *seq;
    gt_str_reset(descname);
    gt_str_append_str(descname, splitdesc);
    gt_str_append_char(descname, '/');
    gt_str_append_cstr(descname, gt_bioseq_get_description(bioseq, i));
    gt_str_append_cstr(descname, gt_file_suffix(filename));
    if (!(outfp = gt_output_file_xopen_forcecheck(gt_str_get(descname), "w",
                                                 force, err))) {
      had_err = -1;
      break;
    }
    seq = gt_bioseq_get_sequence(bioseq, i);
    gt_fasta_show_entry(gt_bioseq_get_description(bioseq, i), seq,
                        gt_bioseq_get_sequence_length(bioseq, i), width,
                        outfp);
    gt_free(seq);
    gt_file_delete(outfp);
  }

  gt_bioseq_delete(bioseq);
  gt_str_delete(descname);

  return had_err;
}

static int split_fasta_file(const char *filename, unsigned long max_filesize,
                            bool force, GtError *err)
{
  GtFile *srcfp = NULL, *destfp = NULL;
  GtStr *destfilename = NULL;
  unsigned long filenum = 0, bytecount = 0, separator_pos;
  int read_bytes, had_err = 0;
  char buf[BUFSIZ];

  gt_error_check(err);
  gt_assert(filename && max_filesize);

  /* open source file */
  srcfp = gt_file_xopen(filename, "r");
  gt_assert(srcfp);

  /* read start characters */
  if ((read_bytes = gt_file_xread(srcfp, buf, BUFSIZ)) == 0) {
    gt_error_set(err, "file \"%s\" is empty", filename);
    had_err = -1;
  }
  bytecount += read_bytes;

  /* make sure the file is in fasta format */
  if (!had_err && buf[0] != '>') {
    gt_error_set(err, "file is not in FASTA format");
    had_err = -1;
  }

  if (!had_err) {
    /* open destination file */
    destfilename = gt_str_new();
    gt_str_append_cstr_nt(destfilename, filename,
                          gt_file_basename_length(filename));
    gt_str_append_char(destfilename, '.');
    gt_str_append_ulong(destfilename, ++filenum);
    gt_str_append_cstr(destfilename,
                       gt_file_mode_suffix(gt_file_mode(srcfp)));
    if (!(destfp = gt_output_file_xopen_forcecheck(gt_str_get(destfilename),
                                                  "w",
                                                  force, err))) {
      had_err = -1;
    }
    if (!had_err)
      gt_file_xwrite(destfp, buf, read_bytes);

    while (!had_err &&
           (read_bytes = gt_file_xread(srcfp, buf, BUFSIZ)) != 0) {
      if (bytecount + read_bytes > max_filesize) {
        int offset = bytecount < max_filesize ? max_filesize - bytecount : 0;
        if ((separator_pos = buf_contains_separator(buf, offset, read_bytes))) {
          separator_pos--;
          gt_assert(separator_pos < read_bytes);
          if (separator_pos)
            gt_file_xwrite(destfp, buf, separator_pos);
          /* close current file */
          gt_file_delete(destfp);
          /* open new file */
          gt_str_reset(destfilename);
          gt_str_append_cstr_nt(destfilename, filename,
                                gt_file_basename_length(filename));
          gt_str_append_char(destfilename, '.');
          gt_str_append_ulong(destfilename, ++filenum);
          gt_str_append_cstr(destfilename,
                             gt_file_mode_suffix(gt_file_mode(srcfp)));
          if (!(destfp =
                  gt_output_file_xopen_forcecheck(gt_str_get(destfilename), "w",
                                                 force, err))) {
            had_err = -1;
            break;
          }
          bytecount = read_bytes - separator_pos; /* reset */
          gt_assert(buf[separator_pos] == '>');
          gt_file_xwrite(destfp, buf + separator_pos,
                         read_bytes - separator_pos);
          continue;
        }
      }
      bytecount += read_bytes;
      gt_file_xwrite(destfp, buf, read_bytes);
    }
  }

  /* free */
  gt_str_delete(destfilename);

  /* close current file */
  gt_file_delete(destfp);

  /* close source file */
  gt_file_delete(srcfp);

  return had_err;
}

static int gt_splitfasta_runner(GT_UNUSED int argc, const char **argv,
                                int parsed_args, void *tool_arguments,
                                GtError *err)
{
  SplitfastaArguments *arguments = tool_arguments;
  unsigned int num_files;
  int had_err;
  off_t file_size;
  const char* filename;
  gt_error_check(err);
  gt_assert(arguments);

  num_files = arguments->num_files;
  filename = argv[parsed_args];

  if (gt_str_length(arguments->splitdesc)) {
    had_err = split_description(filename, arguments->splitdesc,
                                arguments->width, arguments->force, err);
  }
  else {
    unsigned long max_filesize;
    if (num_files) {
      /* set the maxfile size based on requested number of files */
      file_size = gt_file_estimate_size(filename);
      max_filesize= file_size / num_files ;
    }
    else
      max_filesize= arguments->max_filesize_in_MB << 20;
    had_err = split_fasta_file(filename, max_filesize, arguments->force, err);
  }

  return had_err;
}

GtTool* gt_splitfasta(void)
{
  return gt_tool_new(gt_splitfasta_arguments_new,
                     gt_splitfasta_arguments_delete,
                     gt_splitfasta_option_parser_new,
                     NULL,
                     gt_splitfasta_runner);
}
