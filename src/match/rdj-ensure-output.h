/*
  Copyright (c) 2010-2011 Giorgio Gonnella <gonnella@zbh.uni-hamburg.de>
  Copyright (c) 2010-2011 Center for Bioinformatics, University of Hamburg

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

#ifndef RDJ_ENSURE_OUTPUT_H
#define RDJ_ENSURE_OUTPUT_H

#include <string.h>
#include "core/ensure.h"
#include "core/fa.h"
#include "core/xansi_api.h"

/* defines macros to make easier to test output to a temporary file */

#define GT_ENSURE_OUTPUT_DECLARE(MAXOUTSIZE)                               \
        GtStr *fnamestr;                                                   \
        char *fname;                                                       \
        FILE *handle;                                                      \
        GtFile *outfp;                                                     \
        size_t outsize;                                                    \
        char buffer[MAXOUTSIZE]

#define GT_ENSURE_OUTPUT_GETFILE(GTFILE)                                   \
        fnamestr = gt_str_new();                                           \
        handle = gt_xtmpfp(fnamestr);                                      \
        (GTFILE) = gt_file_new_from_fileptr(handle)

#define GT_ENSURE_OUTPUT_TESTFILE(GTFILE, EXPOUT)                          \
        (void)fflush(handle);                                              \
        gt_file_xrewind(GTFILE);                                           \
        outsize = (size_t)gt_file_xread((GTFILE), buffer, strlen(EXPOUT)); \
        gt_ensure(had_err, outsize == strlen(EXPOUT));                     \
        gt_ensure(had_err, memcmp((EXPOUT), buffer, outsize) == 0);        \
        if (had_err != 0)                                                  \
        {                                                                  \
          fprintf(stderr, "\nExpected file content: \n%s\n",               \
                  EXPOUT);                                                 \
          fprintf(stderr, "\nFile content actually read "                  \
                         "(up to %lu bytes):\n%s\n",                       \
               (unsigned long) strlen(EXPOUT), buffer);                    \
        }

#define GT_ENSURE_OUTPUT_RMFILE(GTFILE)                                    \
        gt_file_delete(GTFILE);                                            \
        fname = gt_str_get(fnamestr);                                      \
        gt_xremove(fname);                                                 \
        gt_str_delete(fnamestr)

#define GT_ENSURE_OUTPUT(FNCALL, EXPOUT)                                   \
        GT_ENSURE_OUTPUT_GETFILE(outfp);                                   \
        FNCALL;                                                            \
        GT_ENSURE_OUTPUT_TESTFILE(outfp, EXPOUT);                          \
        GT_ENSURE_OUTPUT_RMFILE(outfp)

/*
  Compares the output of the function call FNCALL
  to the char* EXPOUT (size and content).

  Usage: GT_ENSURE_OUTPUT_DECLARE in the declarations,
         then GT_ENSURE_OUTPUT to run the test (also multiple times);
         FNCALL should print something to the GtFile* outfp
         or to the FILE* handle

  Assumes:
  - int had_err is the return value variable of the unit test method
  - GtError *err (required by the gt_ensure macro)
  - reserved variables: fnamestr, fname, handle, outfp, outsize, buffer
*/

#endif
