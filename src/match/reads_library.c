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

#include "core/ma.h"
#include "core/file.h"
#include "match/reads_library.h"

int gt_reads_library_table_write(GtReadsLibrary *lib_table,
    unsigned long noflibs, const char *path, GtError *err)
{
  GtFile *file;
  file = gt_file_new(path, "w", err);
  if (file == NULL)
    return -1;
  else
  {
    gt_file_xwrite(file, &noflibs, sizeof (noflibs));
    gt_file_xwrite(file, lib_table, sizeof (*lib_table) * noflibs);
    gt_file_delete(file);
  }
  return 0;
}

GtReadsLibrary* gt_reads_library_table_read(const char *path, GtError *err,
    unsigned long *noflibs)
{
  GtFile *file;
  GtReadsLibrary *lib_table = NULL;
  file = gt_file_new(path, "r", err);
  if (file != NULL)
  {
    int freadretval;
    gt_assert(noflibs != NULL);
    freadretval = gt_file_xread(file, noflibs, sizeof (noflibs));
    if ((size_t)freadretval != sizeof (noflibs))
    {
      gt_error_set(err, "library file %s: "
          "error by reading number of libraries", path);
    }
    else if (*noflibs == 0)
    {
      gt_error_set(err, "library file %s: "
          "number of libraries is 0", path);
    }
    else
    {
      const size_t lib_table_size = sizeof (*lib_table) * (*noflibs);
      lib_table = gt_malloc(lib_table_size);
      freadretval = gt_file_xread(file, lib_table, lib_table_size);
      if ((size_t)freadretval != lib_table_size)
      {
        gt_error_set(err, "library file %s: "
            "error by reading libraries table "
            "(%lu bytes expected, %d bytes read)",
            path, (unsigned long)lib_table_size, freadretval);
        gt_free(lib_table);
        lib_table = NULL;
      }
    }
    gt_file_delete(file);
  }
  return lib_table;
}
