/*
  Copyright (c) 2003-2008 Gordon Gremme <gremme@zbh.uni-hamburg.de>
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

#include "libgtcore/bioseq.h"
#include "libgtcore/fa.h"
#include "libgtcore/fasta.h"
#include "libgtcore/fileutils.h"
#include "libgtcore/ma.h"
#include "libgtcore/option.h"
#include "libgtcore/outputfile.h"
#include "libgtcore/unused.h"
#include "libgtcore/xansi.h"
#include "tools/gt_splitfasta.h"

typedef struct {
  unsigned long max_filesize_in_MB;
  Str *splitdesc;
  bool force;
} SplitfastaArguments;

static void* gt_splitfasta_arguments_new(void)
{
  SplitfastaArguments *arguments = ma_calloc(1, sizeof *arguments);
  arguments->splitdesc = str_new();
  return arguments;
}

void gt_splitfasta_arguments_delete(void *tool_arguments)
{
  SplitfastaArguments *arguments = tool_arguments;
  if (!arguments) return;
  str_delete(arguments->splitdesc);
  ma_free(arguments);
}

static OptionParser* gt_splitfasta_option_parser_new(void *tool_arguments)
{
  SplitfastaArguments *arguments = tool_arguments;
  OptionParser *op;
  Option *targetsize_option, *splitdesc_option, *o;
  assert(arguments);
  op = option_parser_new("[option ...] fastafile","Split the supplied fasta "
                         "file.");
  targetsize_option = option_new_ulong_min("targetsize", "set the target file "
                                           "size in MB",
                                           &arguments->max_filesize_in_MB, 50,
                                           1);
  option_parser_add_option(op, targetsize_option);
  splitdesc_option = option_new_string("splitdesc", "put every fasta entry in "
                                       "a separate file named by its "
                                       "description in the given directory",
                                       arguments->splitdesc, NULL);
  option_parser_add_option(op, splitdesc_option);
  option_exclude(targetsize_option, splitdesc_option);
  o = option_new_bool(FORCE_OPT_CSTR, "force writing to output file",
                      &arguments->force, false);
  option_parser_add_option(op, o);
  option_parser_set_min_max_args(op, 1, 1);
  return op;
}

static unsigned long buf_contains_separator(char *buf)
{
  char *cc;
  assert(buf);
  for (cc = buf; cc < buf + BUFSIZ; cc++) {
    if (*cc == '>')
      return cc - buf + 1;
  }
  return 0;
}

static GenFile* genfile_xopen_forcecheck(const char *path, const char *mode,
                                         bool force, Error *err)
{
  if (!force && file_exists(path)) {
    error_set(err, "file \"%s\" exists already, use option -%s to overwrite",
              path, FORCE_OPT_CSTR);
    return NULL;
  }
  return genfile_xopen(path, mode);
}

static int split_description(const char *filename, Str *splitdesc, bool force,
                             Error *err)
{
  unsigned long i;
  Bioseq *bioseq;
  Str *descname;
  int had_err = 0;
  error_check(err);
  assert(filename && splitdesc && str_length(splitdesc));

  descname = str_new();
  if (!(bioseq = bioseq_new(filename, err)))
    had_err = -1;

  for (i = 0; !had_err && i < bioseq_number_of_sequences(bioseq); i++) {
    GenFile *outfp;
    str_reset(descname);
    str_append_str(descname, splitdesc);
    str_append_char(descname, '/');
    str_append_cstr(descname, bioseq_get_description(bioseq, i));
    str_append_cstr(descname, file_suffix(filename));
    if (!(outfp = genfile_xopen_forcecheck(str_get(descname), "w", force,
                                           err))) {
      had_err = -1;
      break;
    }
    fasta_show_entry_generic(bioseq_get_description(bioseq, i),
                             bioseq_get_sequence(bioseq, i),
                             bioseq_get_sequence_length(bioseq, i), 0, outfp);
    genfile_close(outfp);
  }

  bioseq_delete(bioseq);
  str_delete(descname);

  return had_err;
}

static int split_fasta_file(const char *filename,
                            unsigned long max_filesize_in_bytes, bool force,
                            Error *err)
{
  GenFile *srcfp = NULL, *destfp = NULL;
  Str *destfilename = NULL;
  unsigned long filenum = 0, bytecount = 0, separator_pos;
  int read_bytes, had_err = 0;
  char buf[BUFSIZ];
  error_check(err);
  assert(filename && max_filesize_in_bytes);

  /* open source file */
  srcfp = genfile_xopen(filename, "r");
  assert(srcfp);

  /* read start characters */
  if ((read_bytes = genfile_xread(srcfp, buf, BUFSIZ)) == 0) {
    error_set(err, "file \"%s\" is empty", filename);
    had_err = -1;
  }
  bytecount += read_bytes;

  /* make sure the file is in fasta format */
  if (!had_err && buf[0] != '>') {
    error_set(err, "file is not in FASTA format");
    had_err = -1;
  }

  if (!had_err) {
    /* open destination file */
    destfilename = str_new();
    str_append_cstr_nt(destfilename, filename,
                       genfile_basename_length(filename));
    str_append_char(destfilename, '.');
    str_append_ulong(destfilename, ++filenum);
    str_append_cstr(destfilename, genfilemode_suffix(genfile_mode(srcfp)));
    if (!(destfp = genfile_xopen_forcecheck(str_get(destfilename), "w", force,
                                            err))) {
      had_err = -1;
    }
    if (!had_err)
      genfile_xwrite(destfp, buf, read_bytes);

    while (!had_err && (read_bytes = genfile_xread(srcfp, buf, BUFSIZ)) != 0) {
      bytecount += read_bytes;
      if (bytecount > max_filesize_in_bytes &&
          (separator_pos = buf_contains_separator(buf))) {
        separator_pos--;
        assert(separator_pos < BUFSIZ);
        if (separator_pos)
          genfile_xwrite(destfp, buf, separator_pos);
        /* close current file */
        genfile_close(destfp);
        /* open new file */
        str_reset(destfilename);
        str_append_cstr_nt(destfilename, filename,
                           genfile_basename_length(filename));
        str_append_char(destfilename, '.');
        str_append_ulong(destfilename, ++filenum);
        str_append_cstr(destfilename, genfilemode_suffix(genfile_mode(srcfp)));
        if (!(destfp = genfile_xopen_forcecheck(str_get(destfilename), "w",
                                                force, err))) {
          had_err = -1;
          break;
        }
        bytecount = 0;
        assert(buf[separator_pos] == '>');
        genfile_xwrite(destfp, buf+separator_pos, read_bytes-separator_pos);
      }
      else
        genfile_xwrite(destfp, buf, read_bytes);
    }
  }

  /* free */
  str_delete(destfilename);

  /* close current file */
  genfile_close(destfp);

  /* close source file */
  genfile_close(srcfp);

  return had_err;
}

static int gt_splitfasta_runner(UNUSED int argc, const char **argv,
                                int parsed_args, void *tool_arguments,
                                Error *err)
{
  SplitfastaArguments *arguments = tool_arguments;
  unsigned long max_filesize_in_bytes;
  int had_err;
  error_check(err);
  assert(arguments);

  max_filesize_in_bytes = arguments->max_filesize_in_MB << 20;

  if (str_length(arguments->splitdesc)) {
    had_err = split_description(argv[parsed_args], arguments->splitdesc,
                                arguments->force, err);
  }
  else {
    had_err = split_fasta_file(argv[parsed_args], max_filesize_in_bytes,
                               arguments->force, err);
  }

  return had_err;
}

Tool* gt_splitfasta(void)
{
  return tool_new(gt_splitfasta_arguments_new,
                  gt_splitfasta_arguments_delete,
                  gt_splitfasta_option_parser_new,
                  NULL,
                  gt_splitfasta_runner);
}
