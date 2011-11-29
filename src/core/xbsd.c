/*
  Copyright (c) 2010 Gordon Gremme <gremme@zbh.uni-hamburg.de>

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
#include <string.h>
#include "core/xbsd.h"

static void gt_xflock_with_op(int fd, short l_type)
{
  struct flock f;
  memset(&f, 0, sizeof (f));
  f.l_type = l_type;
  if (fcntl(fd, F_SETLKW, &f)) {
    perror("cannot flock");
    exit(EXIT_FAILURE);
  }
}

void gt_xflock_shared(int fd)
{
  gt_xflock_with_op(fd, F_RDLCK);
}

void gt_xflock_exclusive(int fd)
{
  gt_xflock_with_op(fd, F_WRLCK);
}

void gt_xflock_unlock(int fd)
{
  gt_xflock_with_op(fd, F_UNLCK);
}
