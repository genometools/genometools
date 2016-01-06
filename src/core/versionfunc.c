/*
  Copyright (c) 2006-2013 Gordon Gremme <gordon@gremme.org>
  Copyright (c) 2006-2008 Center for Bioinformatics, University of Hamburg

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

#include <stdio.h>
#include "gt_config.h"

void gt_versionfunc(const char *progname)
{
  printf("%s (GenomeTools) %s\n", progname, GT_VERSION);
  printf("Copyright (c) 2003-2016 G. Gremme, S. Steinbiss, S. Kurtz, and "
         "CONTRIBUTORS\n");
  printf("Copyright (c) 2003-2016 Center for Bioinformatics, University of "
         "Hamburg\n");
  printf("See LICENSE file or http://genometools.org/license.html for license "
         "details.\n\n");
  printf("Used compiler: %s\n", GT_CC);
  printf("Compile flags: %s\n", GT_CFLAGS);
}

void gt_showshortversion(const char *progname)
{
  printf("%s (GenomeTools) %s\n", progname, GT_VERSION);
}
