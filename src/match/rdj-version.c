/*
  Copyright (c) 2012 Giorgio Gonnella <gonnella@zbh.uni-hamburg.de>
  Copyright (c) 2012 Center for Bioinformatics, University of Hamburg

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
#include "core/ma.h"
#include "core/versionfunc.h"
#include "match/rdj-version.h"

void gt_readjoiner_show_version(const char *progname)
{
  char *clean_progname = NULL;
  size_t pnsize = strlen(progname), rdjsize = strlen(" "GT_READJOINER_CMD);
  printf("Readjoiner: a string graph-based sequence assembler\n\n");
  printf("version "GT_READJOINER_VERSION"\n\n");

  /* rm " readjoiner" from progname if possible, to avoid confusing the user */
  clean_progname = gt_malloc(pnsize);
  (void)strcpy(clean_progname, progname);
  if (pnsize > rdjsize && strcmp(clean_progname + (pnsize - rdjsize),
        " "GT_READJOINER_CMD) == 0)
  {
    clean_progname[pnsize - rdjsize] = '\0';
  }
  printf("GenomeTools version:\n");
  gt_versionfunc(clean_progname);
  gt_free(clean_progname);
}
