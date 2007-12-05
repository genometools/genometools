/*
  Copyright (c) 2005-2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2005-2007 Center for Bioinformatics, University of Hamburg

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

#include <sys/stat.h>
#include <sys/mman.h>
#include <assert.h>
#include <fcntl.h>
#include "libgtcore/mmap.h"
#include "libgtcore/xposix.h"

static void* mmap_generic(const char *path, size_t *len, int prot)
{
  int fd;
  struct stat sb;
  void *map;
  assert(path && len);
  fd = open(path, O_RDONLY, 0);
  if (fd == -1)
    return NULL;
  xfstat(fd, &sb);
  if ((map = mmap(0, sb.st_size, prot, MAP_PRIVATE, fd, 0)) == MAP_FAILED)
    map = NULL;
  if (map)
    *len = sb.st_size;
  xclose(fd);
  return map;
}

void* mmap_read(const char *path, size_t *len)
{
  assert(path && len);
  return mmap_generic(path, len, PROT_READ);
}

void* mmap_write(const char *path, size_t *len)
{
  assert(path && len);
  return mmap_generic(path, len, PROT_WRITE);
}

static void* xmmap_generic(const char *path, size_t *len, int prot)
{
  int fd;
  struct stat sb;
  void *map;
  assert(path && len);
  fd = xopen(path, O_RDONLY, 0);
  xfstat(fd, &sb);
  map = xmmap(0, sb.st_size, prot, MAP_PRIVATE, fd, 0);
  *len = sb.st_size;
  xclose(fd);
  return map;
}

void* xmmap_read(const char *path, size_t *len)
{
  assert(path && len);
  return xmmap_generic(path, len, PROT_READ);
}

void* xmmap_write(const char *path, size_t *len)
{
  assert(path && len);
  return xmmap_generic(path, len, PROT_WRITE);
}
