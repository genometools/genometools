/*
  Copyright (c) 2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2007 Center for Bioinformatics, University of Hamburg

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

#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include "libgtcore/option.h"
#include "libgtcore/seqiterator.h"
#include "libgtcore/versionfunc.h"
#include "libgtcore/xposix.h"
#include "libgtcore/progressbar.h"

static OPrval parse_options(int *parsed_args, int argc, const char **argv,
                            Env *env)
{
  OptionParser *op;
  OPrval oprval;
  env_error_check(env);
  op = option_parser_new("file [...]", "Parse the supplied Fasta files.", env);
  oprval = option_parser_parse_min_args(op, parsed_args, argc, argv,
                                            versionfunc, 1, env);
  option_parser_delete(op, env);
  return oprval;
}

static off_t estimatetotalfilesizes(const StrArray *files)
{
  unsigned long filenum;
  off_t totalsize = 0;
  struct stat sb;
  GenFileMode gfm;
  int fd;

  for(filenum = 0; filenum < strarray_size(files); filenum++)
  {
    fd = xopen(strarray_get(files,filenum), O_RDONLY, 0);
    xfstat(fd, &sb);
    gfm = genfilemode_determine(strarray_get(files,filenum));
    if(gfm == GFM_UNCOMPRESSED)
    {
      totalsize += sb.st_size;
    } else
    {
      totalsize += (4*sb.st_size);
    }
    xclose(fd);
  }
  return totalsize;
}

int gt_seqiterator(int argc, const char **argv, Env *env)
{
  StrArray *files;
  SeqIterator *seqit;
  const Uchar *sequence;
  char *desc;
  unsigned long len;
  int i, parsed_args, had_err;
  off_t totalsize;

  env_error_check(env);

  /* option parsing */
  switch (parse_options(&parsed_args, argc, argv, env)) {
    case OPTIONPARSER_OK: break;
    case OPTIONPARSER_ERROR: return -1;
    case OPTIONPARSER_REQUESTS_EXIT: return 0;
  }

  files = strarray_new(env);
  for (i = parsed_args; i < argc; i++)
    strarray_add_cstr(files, argv[i], env);

  totalsize = estimatetotalfilesizes(files);
  printf("# estimated total size is %llu\n",(unsigned long long) totalsize);
  seqit = seqiterator_new(files, NULL, true, env);
  progressbar_start(getcurrentcounter(seqit), (unsigned long long) totalsize);
  while (true)
  {
    had_err = seqiterator_next(seqit, &sequence, &len, &desc, env);
    if(had_err != 1)
    {
      break;
    }
    env_ma_free(desc, env);
  }
  progressbar_stop();
  seqiterator_delete(seqit, env);
  strarray_delete(files, env);
  return had_err;
}
