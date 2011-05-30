/*
  Copyright (c) 2007-2008 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2007-2008 Center for Bioinformatics, University of Hamburg

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

#include <fcntl.h>
#include <sys/mman.h>
#include "core/error.h"
#include "core/option_api.h"
#include "core/progressbar.h"
#include "core/versionfunc.h"
#include "core/xposix.h"
#include "tools/gt_mmapandread.h"

static GtOPrval parse_options(int *parsed_args, int argc, const char **argv,
                              GtError *err)
{
  GtOptionParser *op;
  GtOPrval oprval;
  gt_error_check(err);
  op = gt_option_parser_new("file [...]",
                            "Map the supplied files into memory and "
                            "read them once.");
  gt_option_parser_set_min_args(op, 1);
  oprval = gt_option_parser_parse(op, parsed_args, argc, argv, gt_versionfunc,
                                  err);
  gt_option_parser_delete(op);
  return oprval;
}

int gt_mmapandread(int argc, const char **argv, GtError *err)
{
  int i, fd, parsed_args;
  void *map;
  struct stat sb;
  unsigned long long j;
  unsigned int byte = 0;
  gt_error_check(err);

  /* option parsing */
  switch (parse_options(&parsed_args, argc, argv, err)) {
    case GT_OPTION_PARSER_OK: break;
    case GT_OPTION_PARSER_ERROR: return -1;
    case GT_OPTION_PARSER_REQUESTS_EXIT: return 0;
  }

  /* iterate over all files */
  for (i = parsed_args; i < argc; i++) {
    /* open file */
    fd = gt_xopen(argv[i], O_RDONLY, 0);

    /* get file statistics */
    gt_xfstat(fd, &sb);

    if (sb.st_size == 0)
      printf("file \"%s\" is empty\n", argv[i]);
    else if (!(sb.st_mode & S_IFREG))
      printf("\"%s\" is not a regular file\n", argv[i]);
    else {
      /* map file */
      map = gt_xmmap(0, sb.st_size, PROT_READ, MAP_PRIVATE, fd, 0);

      /* read file */
      printf("reading file \"%s\"\n", argv[i]);
      j = 0;
      gt_progressbar_start(&j, (unsigned long long) sb.st_size);
      for (; j < (unsigned long long) sb.st_size; j++)
        byte |= (unsigned int) ((char*) map)[j];
      gt_progressbar_stop();

      /* unmap file */
      gt_xmunmap(map, sb.st_size);
    }

    /* close file */
    gt_xclose(fd);
  }

  if (!byte)
    printf("all read files contained only null characters\n");

  return 0;
}
