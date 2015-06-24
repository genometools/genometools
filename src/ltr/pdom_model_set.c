/*
  Copyright (c) 2013-2014 Sascha Steinbiss <sascha@steinbiss.name>
  Copyright (c) 2013 Center for Bioinformatics, University of Hamburg

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

#include <stdlib.h>
#include <unistd.h>
#include <sys/types.h>
#ifndef _WIN32
#include <sys/wait.h>
#endif
#include "core/compat.h"
#include "core/error_api.h"
#include "core/fileutils_api.h"
#include "core/ma.h"
#include "core/md5_fingerprint_api.h"
#include "core/str_array_api.h"
#include "ltr/pdom_model_set.h"

#define GT_HMM_INDEX_SUFFIX "h3i"
#ifdef _WIN32
#ifndef WEXITSTATUS
#define WEXITSTATUS(x) (x)
#endif
#endif

struct GtPdomModelSet
{
  GtStr *filename;
};

#define PDOM_MODEL_SET_HMMER_NOT_FOUND "Please make sure that all HMMER " \
                                       "executables are in your PATH."

static void gt_pdom_model_set_delete_index(const char *idxname)
{
  char filename[BUFSIZ];
  gt_assert(idxname);
  (void) unlink(idxname);
  (void) snprintf(filename, BUFSIZ, "%s.h3f", idxname);
  (void) unlink(filename);
  (void) snprintf(filename, BUFSIZ, "%s.h3i", idxname);
  (void) unlink(filename);
  (void) snprintf(filename, BUFSIZ, "%s.h3m", idxname);
  (void) unlink(filename);
  (void) snprintf(filename, BUFSIZ, "%s.h3p", idxname);
  (void) unlink(filename);
}

GtPdomModelSet* gt_pdom_model_set_new(GtStrArray *hmmfiles, bool force,
                                      GtError *err)
{
  GtStr *concat_dbnames, *cmdline, *indexfilename = NULL;
  GtUword i;
  char *md5_hash, ch;
  const char *tmpdir;
  int had_err = 0, rval;
  FILE *dest;
  GtPdomModelSet *pdom_model_set;
  gt_assert(hmmfiles);
  gt_error_check(err);

  rval = system("hmmpress -h > /dev/null");
  if (WEXITSTATUS(rval) != 0) {
    gt_error_set(err, "Error running hmmpress. "
                      PDOM_MODEL_SET_HMMER_NOT_FOUND);
    return NULL;
  }
  rval = system("hmmconvert -h > /dev/null");
  if (WEXITSTATUS(rval) != 0) {
    gt_error_set(err, "Error running hmmconvert. "
                      PDOM_MODEL_SET_HMMER_NOT_FOUND);
    return NULL;
  }

  pdom_model_set = gt_calloc((size_t) 1, sizeof (GtPdomModelSet));
  concat_dbnames = gt_str_new();
  for (i = 0; !had_err && i < gt_str_array_size(hmmfiles); i++) {
    const char *filename = gt_str_array_get(hmmfiles, i);
    if (!gt_file_exists(filename)) {
      gt_error_set(err, "invalid HMM file: %s", filename);
      gt_str_delete(concat_dbnames);
      gt_free(pdom_model_set);
      return NULL;
    } else {
      gt_str_append_cstr(concat_dbnames, filename);
    }
  }
  if (!had_err) {
    pdom_model_set->filename = gt_str_new();
    if (!(tmpdir = getenv("TMPDIR")))
      tmpdir = "/tmp";
    gt_str_append_cstr(pdom_model_set->filename, tmpdir);
    gt_str_append_char(pdom_model_set->filename, GT_PATH_SEPARATOR);
    md5_hash = gt_md5_fingerprint(gt_str_get(concat_dbnames),
                                  gt_str_length(concat_dbnames));
    gt_str_append_cstr(pdom_model_set->filename, md5_hash);
    gt_free(md5_hash);
    indexfilename = gt_str_new_cstr(gt_str_get(pdom_model_set->filename));
    gt_str_append_cstr(indexfilename, GT_HMM_INDEX_SUFFIX);
  }

  if (!gt_file_exists(gt_str_get(indexfilename)) || force) {
    dest = fopen(gt_str_get(pdom_model_set->filename), "w+");
    if (!dest) {
      gt_error_set(err, "could not create file %s",
                 gt_str_get(pdom_model_set->filename));
      had_err = -1;
    }
    if (!had_err) {
      for (i = 0; !had_err && i < gt_str_array_size(hmmfiles); i++) {
        FILE *source;
        char cmd[BUFSIZ];
        (void) snprintf(cmd, BUFSIZ, "hmmconvert %s",
                        gt_str_array_get(hmmfiles, i));
        source = popen(cmd, "r");
        if (!source) {
          gt_error_set(err, "error opening/converting HMM file %s",
                            gt_str_array_get(hmmfiles, i));
          had_err = -1;
          gt_pdom_model_set_delete_index(gt_str_get(pdom_model_set->filename));
        }
        if (!had_err) {
          GtWord len = 0;
          while ((ch = fgetc(source)) != EOF) {
            (void) fputc(ch, dest);
            len++;
          }
          (void) pclose(source);
          if (len == 0) {
            gt_error_set(err, "invalid HMMER format encountered trying to "
                              "convert HMM file %s",
                              gt_str_array_get(hmmfiles, i));
            had_err = -1;
            gt_pdom_model_set_delete_index
                                         (gt_str_get(pdom_model_set->filename));
          }
        }
      }
      (void) fclose(dest);
    }
    if (!had_err) {
      cmdline = gt_str_new_cstr("hmmpress -f ");
      gt_str_append_str(cmdline, pdom_model_set->filename);
      gt_str_append_cstr(cmdline, "> /dev/null");   /* XXX: portability? */

      rval = system(gt_str_get(cmdline));
      gt_str_delete(cmdline);
      if (WEXITSTATUS(rval) != 0) {
        gt_error_set(err, "Error executing hmmpress. "
                          "Please make sure your pHMMs are valid and"
                          "share a common version (e.g. HMMER2 or HMMER3).");
        had_err = -1;
      }
    }
  }

  if (had_err) {
    gt_pdom_model_set_delete(pdom_model_set);
    pdom_model_set = NULL;
  }
  gt_str_delete(indexfilename);
  gt_str_delete(concat_dbnames);
  return pdom_model_set;
}

const char* gt_pdom_model_set_get_filename(GtPdomModelSet *set)
{
  gt_assert(set);
  return gt_str_get(set->filename);
}

void gt_pdom_model_set_delete(GtPdomModelSet *set)
{
  if (!set) return;
  gt_str_delete(set->filename);
  gt_free(set);
}
