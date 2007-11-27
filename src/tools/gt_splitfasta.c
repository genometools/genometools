/*
  Copyright (c) 2003-2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2003-2007 Center for Bioinformatics, University of Hamburg

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

#include "libgtcore/fa.h"
#include "libgtcore/option.h"
#include "libgtcore/versionfunc.h"
#include "libgtcore/xansi.h"

static OPrval parse_options(int *parsed_args, unsigned long *max_filesize_in_MB,
                            int argc, const char **argv, Env *env)
{
  OptionParser *op;
  Option *o;
  OPrval oprval;
  env_error_check(env);
  op = option_parser_new("[option ...] fastafile","Split the supplied fasta "
                         "file.", env);
  o = option_new_ulong_min("targetsize", "set the target file size in MB",
                           max_filesize_in_MB, 50, 1, env);
  option_parser_add_option(op, o, env);
  oprval = option_parser_parse_min_max_args(op, parsed_args, argc, argv,
                                            versionfunc, 1, 1, env);
  option_parser_delete(op, env);
  return oprval;
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

int gt_splitfasta(int argc, const char **argv, Env *env)
{
  GenFile *srcfp = NULL;
  FILE *destfp = NULL;
  Str *destfilename = NULL;
  unsigned long filenum = 0, bytecount = 0, max_filesize_in_bytes,
                max_filesize_in_MB, separator_pos;
  int read_bytes, parsed_args, had_err = 0;
  char buf[BUFSIZ];

  /* option parsing */
  switch (parse_options(&parsed_args, &max_filesize_in_MB, argc, argv, env)) {
    case OPTIONPARSER_OK: break;
    case OPTIONPARSER_ERROR: return -1;
    case OPTIONPARSER_REQUESTS_EXIT: return 0;
  }
  assert(parsed_args + 1 == argc);
  max_filesize_in_bytes = max_filesize_in_MB << 20;

  /* open source file */
  srcfp = genfile_xopen(argv[parsed_args], "r", env);
  assert(srcfp);

  /* read start characters */
  if ((read_bytes = genfile_xread(srcfp, buf, BUFSIZ)) == 0) {
    env_error_set(env, "file \"%s\" is empty", argv[parsed_args]);
    had_err = -1;
  }
  bytecount += read_bytes;

  /* make sure the file is in fasta format */
  if (!had_err && buf[0] != '>') {
    env_error_set(env, "file is not in FASTA format");
    had_err = -1;
  }

  if (!had_err) {
    /* open destination file */
    destfilename = str_new();
    str_append_cstr_nt(destfilename, argv[parsed_args],
                       genfile_basename_length(argv[parsed_args]), env);
    str_append_char(destfilename, '.', env);
    str_append_ulong(destfilename, ++filenum, env);
    destfp = fa_xfopen(str_get(destfilename), "w");
    xfwrite(buf, 1, read_bytes, destfp);

    while ((read_bytes = genfile_xread(srcfp, buf, BUFSIZ)) != 0) {
      bytecount += read_bytes;
      if (bytecount > max_filesize_in_bytes &&
          (separator_pos = buf_contains_separator(buf))) {
        separator_pos--;
        assert(separator_pos < BUFSIZ);
        if (separator_pos)
          xfwrite(buf, 1, separator_pos, destfp);
        /* close current file */
        fa_xfclose(destfp);
        /* open new file */
        str_reset(destfilename);
        str_append_cstr_nt(destfilename, argv[parsed_args],
                           genfile_basename_length(argv[parsed_args]), env);
        str_append_char(destfilename, '.', env);
        str_append_ulong(destfilename, ++filenum, env);
        destfp = fa_xfopen(str_get(destfilename), "w");
        bytecount = 0;
        assert(buf[separator_pos] == '>');
        xfwrite(buf+separator_pos, 1, read_bytes-separator_pos, destfp);
      }
      else
        xfwrite(buf, 1, read_bytes, destfp);
    }
  }

  /* free */
  str_delete(destfilename);

  /* close current file */
  fa_xfclose(destfp);

  /* close source file */
  genfile_xclose(srcfp, env);

  return had_err;
}
