/*
  Copyright (c) 2007-2008, 2013 Gordon Gremme <gordon@gremme.org>
  Copyright (c) 2007-2008       Center for Bioinformatics, University of Hamburg

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
#ifndef _WIN32
#include <sys/mman.h>
#endif
#include "core/error.h"
#include "core/option_api.h"
#include "core/progressbar.h"
#include "core/unused_api.h"
#include "core/versionfunc.h"
#include "core/xposix.h"
#include "tools/gt_mmapandread.h"

static GtOptionParser* gt_mmapandread_option_parser_new(GT_UNUSED void
                                                        *tool_arguments)
{
  GtOptionParser *op;
  op = gt_option_parser_new("file [...]",
                            "Map the supplied files into memory and "
                            "read them once.");
  gt_option_parser_set_min_args(op, 1);
  return op;
}

static int gt_mmapandread_runner(int argc, const char **argv, int parsed_args,
                                 GT_UNUSED void *tool_arguments,
                                 GT_UNUSED GtError *err)
{
  int i, fd;
  void *map;
  struct stat sb;
  GtUint64 j;
  unsigned int byte = 0;
  gt_error_check(err);

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
#ifndef _WIN32
      map = gt_xmmap(0, sb.st_size, PROT_READ, MAP_PRIVATE, fd, 0);
#else
      /* XXX */
      fprintf(stderr, "mmapandread not implemented\n");
      exit(EXIT_FAILURE);
#endif

      /* read file */
      printf("reading file \"%s\"\n", argv[i]);
      j = 0;
      gt_progressbar_start(&j, (GtUint64) sb.st_size);
      for (; j < (GtUint64) sb.st_size; j++)
        byte |= (unsigned int) ((char*) map)[j];
      gt_progressbar_stop();

      /* unmap file */
#ifndef _WIN32
      gt_xmunmap(map, sb.st_size);
#else
      /* XXX */
#endif
    }

    /* close file */
    gt_xclose(fd);
  }

  if (!byte)
    printf("all read files contained only null characters\n");

  return 0;
}

GtTool* gt_mmapandread(void)
{
  return gt_tool_new(NULL,
                     NULL,
                     gt_mmapandread_option_parser_new,
                     NULL,
                     gt_mmapandread_runner);
}
