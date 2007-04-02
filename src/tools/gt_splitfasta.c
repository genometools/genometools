/*
  Copyright (c) 2003-2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2003-2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#include "gt.h"

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

int gt_splitfasta(int argc, const char **argv, Env *env)
{
  FILE *srcfp = NULL, *destfp = NULL;
  Str *destfilename = NULL;
  unsigned long filenum = 0, bytecount = 0, max_filesize_in_bytes,
                max_filesize_in_MB;
  int cc, parsed_args, has_err = 0;

  /* option parsing */
  switch (parse_options(&parsed_args, &max_filesize_in_MB, argc, argv, env)) {
    case OPTIONPARSER_OK: break;
    case OPTIONPARSER_ERROR: return -1;
    case OPTIONPARSER_REQUESTS_EXIT: return 0;
  }
  assert(parsed_args + 1 == argc);
  max_filesize_in_bytes = max_filesize_in_MB << 20;

  /* open source file */
  srcfp = env_fa_xfopen(env, argv[parsed_args], "r");
  assert(srcfp);

  /* read first character */
  if ((cc = getc(srcfp)) ==  EOF) {
    env_error_set(env, "file \"%s\" is empty", argv[parsed_args]);
    has_err = -1;
  }
  bytecount++;

  if (!has_err && cc != '>') {
    env_error_set(env, "file is not in FASTA format");
    has_err = -1;
  }

  if (!has_err) {
    /* open destination file */
    destfilename = str_new_cstr(argv[parsed_args], env);
    str_append_char(destfilename, '.', env);
    str_append_ulong(destfilename, ++filenum, env);
    destfp = env_fa_xfopen(env, str_get(destfilename), "w");
    xfputc(cc, destfp);

    /* XXX: this could me done more efficiently by reading and writing larger
       buffers instead of single characters */
    while ((cc = xfgetc(srcfp)) != EOF) {
      bytecount++;
      if (cc == '>' && bytecount > max_filesize_in_bytes) {
        /* close current file */
        env_fa_xfclose(destfp, env);
        /* open new file */
        str_reset(destfilename);
        str_append_cstr(destfilename, argv[parsed_args], env);
        str_append_char(destfilename, '.', env);
        str_append_ulong(destfilename, ++filenum, env);
        destfp = env_fa_xfopen(env, str_get(destfilename), "w");
        bytecount = 0;
      }
      xfputc(cc, destfp);
    }
  }

  /* free */
  str_delete(destfilename, env);

  /* close current file */
  env_fa_xfclose(destfp, env);

  /* close source file */
  env_fa_xfclose(srcfp, env);

  return has_err;
}
