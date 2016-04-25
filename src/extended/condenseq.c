/*
  Copyright (c) 2014 Dirk Willrodt <willrodt@zbh.uni-hamburg.de>
  Copyright (c) 2014 Center for Bioinformatics, University of Hamburg

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

#include <ctype.h>
#include <string.h>

#include "core/basename_api.h"
#include "core/chardef.h"
#include "core/cstr_api.h"
#include "core/divmodmul.h"
#include "core/fa.h"
#include "core/log_api.h"
#include "core/ma.h"
#include "core/safearith.h"
#include "core/undef_api.h"
#include "core/unused_api.h"
#include "core/warning_api.h"
#include "extended/condenseq.h"
#include "extended/condenseq_rep.h"
#include "extended/gff3_visitor_api.h"

#define GT_CONDENSEQ_BUFFSIZE = (size_t) 512

GtUword gt_condenseq_pos2seqnum(const GtCondenseq *condenseq, GtUword pos)
{
  GtUword ret = 0;
  if (condenseq->orig_num_seq > (GtUword) 1)
    ret = gt_intset_get_idx_smallest_geq(condenseq->ssptab, pos);
  return ret;
}

GtUword gt_condenseq_seqstartpos(const GtCondenseq *condenseq, GtUword seqnum)
{
  GtUword ret = 0;
  if (seqnum != 0)
    ret = gt_intset_get(condenseq->ssptab, seqnum - 1) + 1;
  return ret;
}

static inline GtUword condenseq_seqlength_help(const GtCondenseq *condenseq,
                                               GtUword seqnum,
                                               GtUword seqstart)
{
  GtUword end = condenseq->orig_len;
  if (seqnum < condenseq->orig_num_seq - 1)
    end = gt_intset_get(condenseq->ssptab, seqnum);
  return end - seqstart;
}

GtUword gt_condenseq_seqlength(const GtCondenseq *condenseq, GtUword seqnum)
{
  GtUword start = gt_condenseq_seqstartpos(condenseq, seqnum);
  return condenseq_seqlength_help(condenseq, seqnum, start);
}

static GtUword condenseq_next_sep(const GtCondenseq *ces, GtUword pos)
{
  GtUword seqidx = gt_condenseq_pos2seqnum(ces, pos);
  if (seqidx < ces->orig_num_seq - 1)
    return gt_intset_get(ces->ssptab, seqidx);
  return 0;
}

static GtIntset *condenseq_fill_tab(GtCondenseq *condenseq,
                                    const GtEncseq *orig_es)
{
  GtIntset *ssptab = NULL;
  GtUword max, idx;
  if (condenseq->orig_num_seq > (GtUword) 1) {
    max = gt_encseq_seqstartpos(orig_es, condenseq->orig_num_seq - 1);
    /* we store the internal separators, the end is explicit */
    ssptab = gt_intset_best_new(max - 1, condenseq->orig_num_seq - 1);
    for (idx = (GtUword) 1; idx < condenseq->orig_num_seq; ++idx) {
      GtUword pos = gt_encseq_seqstartpos(orig_es, idx) - 1;
      gt_assert(pos != 0);
      gt_intset_add(ssptab, pos);
    }
  }
  return ssptab;
}

static inline GtUword condenseq_idlen(const char *desc, GtUword desclen)
{
  GtUword idx;
  for (idx = 0; idx < desclen; ++idx) {
    if (isspace(desc[idx]) || desc[idx] == '\0')
      return idx;
  }
  return desclen;
}

static void condenseq_process_descriptions(GtCondenseq *condenseq,
                                           const GtEncseq *orig_es,
                                           GtLogger *logger)
{
  GtUword    *dist;
  const char *desc;
  char       *cur_id_startptr;
  GtUword     desclen,
              dist_idx,
              distsize = (GtUword) 128,
              idlen,
              idx,
              maxendidx = 0,
              maxlen = 0,
              minlen = GT_UWORD_MAX,
              wastedmem = 0,
              sdssize,
              cur_total_id_len = 0;
  bool        use_const_len;

  condenseq->ids_total_len = 0;
  dist = gt_calloc((size_t) distsize, sizeof (*dist));

  for (idx = 0; idx < condenseq->orig_num_seq; ++idx) {
    desc = gt_encseq_description(orig_es, &desclen, idx);
    idlen = condenseq_idlen(desc, desclen);
    if (distsize <= idlen) {
      dist = gt_realloc(dist, (size_t) (idlen + 1) * sizeof (*dist));
      for (dist_idx = distsize; dist_idx <= idlen; dist_idx++)
        dist[dist_idx] = 0;
      distsize = idlen + 1;
    }
    dist[idlen]++;
    if (idlen > maxlen)
      maxlen = idlen;
    if (idlen < minlen)
      minlen = idlen;
    maxendidx += idlen;
  }

  /* calculate memory we would waste if we assume equal length, and size if we
     store actual descriptions */
  for (dist_idx = minlen; dist_idx < maxlen; dist_idx++) {
    wastedmem += dist[dist_idx] * (maxlen - dist_idx);
    condenseq->ids_total_len += dist[dist_idx] * dist_idx;
  }
  condenseq->ids_total_len += dist_idx * dist[dist_idx];

  sdssize = (GtUword) gt_intset_best_memory_size(maxendidx,
                                                 condenseq->orig_num_seq);
  use_const_len = wastedmem < sdssize;

  if (use_const_len) {
    gt_logger_log(logger, "Condenseq descriptions will use const len, " GT_WU
                  ", \"wasting\" " GT_WU " bytes. SDS would use "
                  GT_WU " bytes",
                  maxlen, wastedmem, sdssize);
    condenseq->id_len = maxlen;
    condenseq->ids_total_len = maxlen * condenseq->orig_num_seq;
  }
  else {
    gt_logger_log(logger, "Condenseq descriptions will use sdstab with size "
                  GT_WU ". Const length would have wasted " GT_WU " bytes.",
                  sdssize, wastedmem);
    condenseq->sdstab = gt_intset_best_new(maxendidx, condenseq->orig_num_seq);
  }
  condenseq->orig_ids = gt_calloc((size_t) condenseq->ids_total_len,
                                  sizeof (*condenseq->orig_ids));

  cur_id_startptr = condenseq->orig_ids;
  for (idx = 0; idx < condenseq->orig_num_seq; ++idx) {
    desc = gt_encseq_description(orig_es, &desclen, idx);
    idlen = condenseq_idlen(desc, desclen);
    gt_assert(idlen <= maxlen);
    (void) memcpy(cur_id_startptr, desc, (size_t) idlen);
    if (use_const_len) {
      cur_id_startptr += maxlen;
      cur_total_id_len += maxlen;
    }
    else {
      cur_id_startptr += idlen;
      cur_total_id_len += idlen;
      gt_intset_add(condenseq->sdstab, cur_total_id_len);
    }
  }
  gt_assert(cur_total_id_len == condenseq->ids_total_len);
  gt_free(dist);
}

const char *gt_condenseq_description(const GtCondenseq *condenseq,
                                     GtUword *desclen,
                                     GtUword seqnum)
{
  gt_assert(condenseq != NULL);
  gt_assert(condenseq->orig_num_seq != 0);
  gt_assert(seqnum < condenseq->orig_num_seq);
  if (condenseq->id_len == GT_UNDEF_UWORD) {
    GtUword this = gt_intset_get(condenseq->sdstab, seqnum), previous;
    if (seqnum == 0) {
      *desclen = this;
      return condenseq->orig_ids;
    }
    previous = gt_intset_get(condenseq->sdstab, seqnum - 1);
    *desclen =  this - previous;
    return condenseq->orig_ids + previous;
  }
  else {
    const char *id = condenseq->orig_ids + seqnum * condenseq->id_len;
    *desclen = condenseq->id_len;
    while (id[*desclen - 1] == '\0')
      (*desclen)--;
    return condenseq->orig_ids + seqnum * condenseq->id_len;
  }
}

static GtCondenseq *condenseq_new_empty(const GtAlphabet *alph)
{
  GtCondenseq *condenseq = gt_malloc(sizeof (*condenseq));
  condenseq->alphabet = gt_alphabet_ref((GtAlphabet *) alph);

  condenseq->buffsize =
    condenseq->ldb_allocated =
    condenseq->ldb_nelems =
    condenseq->orig_len =
    condenseq->orig_num_seq =
    condenseq->ubuffsize =
    condenseq->udb_allocated =
    condenseq->udb_nelems = 0;

  condenseq->id_len = GT_UNDEF_UWORD;

  condenseq->buffer = NULL;
  condenseq->filename = NULL;
  condenseq->links = NULL;
  condenseq->orig_ids = NULL;
  condenseq->sdstab = NULL;
  condenseq->ssptab = NULL;
  condenseq->ubuffer = NULL;
  condenseq->unique_es = NULL;
  condenseq->uniques = NULL;

  return condenseq;
}

GtCondenseq *gt_condenseq_new(const GtEncseq *orig_es, GtLogger *logger)
{
  GtCondenseq *condenseq;
  condenseq = condenseq_new_empty(gt_encseq_alphabet(orig_es));

  condenseq->orig_num_seq = gt_encseq_num_of_sequences(orig_es);

  condenseq->ssptab = condenseq_fill_tab(condenseq, orig_es);
  condenseq->orig_len = gt_encseq_total_length(orig_es);

  condenseq_process_descriptions(condenseq, orig_es, logger);
  return condenseq;
}

GtUword gt_condenseq_num_links(const GtCondenseq *condenseq)
{
  return condenseq->ldb_nelems;
}

GtUword gt_condenseq_num_uniques(const GtCondenseq *condenseq)
{
  return condenseq->udb_nelems;
}

GtUword gt_condenseq_total_link_len(const GtCondenseq *condenseq)
{
  GtUword total = 0,
          i;

  for (i = 0; i < condenseq->ldb_nelems; i++) {
    total += condenseq->links[i].len;
  }

  return total;
}

GtUword gt_condenseq_total_unique_len(const GtCondenseq *condenseq)
{
  return gt_encseq_total_length(condenseq->unique_es);
}

GtUword gt_condenseq_total_length(const GtCondenseq *condenseq) {
  return condenseq->orig_len;
}

GtUword gt_condenseq_num_of_sequences(const GtCondenseq *condenseq)
{
  return condenseq->orig_num_seq;
}

#define GT_CONDENSEQ_ARRAY_RESIZE (GtUword) 10
#define GT_CONDENSEQ_ARRAY_INIT (GtUword) 100
GtUword gt_condenseq_array_size_increase(GtUword allocated)
{
  if (allocated != 0) {
    allocated *= 1.2;
  }
  else
    allocated = GT_CONDENSEQ_ARRAY_INIT;
  return allocated;
}

static inline void condenseq_udb_resize(GtCondenseq *condenseq)
{
  if (condenseq->udb_nelems == condenseq->udb_allocated) {
    condenseq->udb_allocated =
      gt_condenseq_array_size_increase(condenseq->udb_allocated);
    condenseq->uniques = gt_realloc(condenseq->uniques,
                              (size_t) condenseq->udb_allocated *
                              sizeof (*condenseq->uniques));
  }
}

static inline void condenseq_ldb_resize(GtCondenseq *condenseq)
{
  if (condenseq->ldb_nelems == condenseq->ldb_allocated) {
    condenseq->ldb_allocated =
      gt_condenseq_array_size_increase(condenseq->ldb_allocated);
    condenseq->links = gt_realloc(condenseq->links,
                            (size_t) condenseq->ldb_allocated *
                            sizeof (*condenseq->links));
  }
}

void gt_condenseq_add_unique_to_db(GtCondenseq *condenseq,
                                   GtUword orig_startpos,
                                   ces_unsigned len)
{
  gt_assert(len != 0);
  /* if previous unique and this one are not consecutive, add the new one */
  if (condenseq->udb_nelems == 0 ||
      condenseq->uniques[condenseq->udb_nelems - 1].orig_startpos +
      condenseq->uniques[condenseq->udb_nelems - 1].len != orig_startpos) {
    gt_assert(condenseq->udb_nelems == 0 ||
              condenseq->uniques[condenseq->udb_nelems - 1].orig_startpos +
              condenseq->uniques[condenseq->udb_nelems - 1].len <
              orig_startpos);
    gt_assert(condenseq->ldb_nelems == 0 ||
              condenseq->links[condenseq->ldb_nelems - 1].orig_startpos +
              condenseq->links[condenseq->ldb_nelems - 1].len <=
              orig_startpos);
    condenseq_udb_resize(condenseq);
    condenseq->uniques[condenseq->udb_nelems].orig_startpos = orig_startpos;
    condenseq->uniques[condenseq->udb_nelems].len = len;
    condenseq->uniques[condenseq->udb_nelems].links.spaceuint32_t = NULL;
    condenseq->udb_nelems++;
  }
  else {
    condenseq->uniques[condenseq->udb_nelems - 1].len += len;
  }
}

void gt_condenseq_add_link_to_db(GtCondenseq *condenseq, GtCondenseqLink link)
{
  condenseq_ldb_resize(condenseq);
  gt_assert(condenseq->links != NULL);
  gt_assert(condenseq->ldb_nelems == 0 ||
            condenseq->links[condenseq->ldb_nelems - 1].orig_startpos +
            condenseq->links[condenseq->ldb_nelems - 1].len <=
            link.orig_startpos);
  gt_assert(condenseq->udb_nelems == 0 ||
            condenseq->uniques[condenseq->udb_nelems - 1].orig_startpos +
            condenseq->uniques[condenseq->udb_nelems - 1].len <=
            link.orig_startpos);
  condenseq->links[condenseq->ldb_nelems] = link;
  condenseq->ldb_nelems++;
}

void gt_condenseq_delete(GtCondenseq *condenseq)
{
  if (condenseq != NULL) {
    GtUword i;
    for (i = 0; i < condenseq->ldb_nelems; i++) {
      gt_editscript_delete(condenseq->links[i].editscript);
    }
    for (i = 0; i < condenseq->udb_nelems; i++) {
      GT_FREEARRAY(&(condenseq->uniques[i].links), uint32_t);
    }
    gt_alphabet_delete(condenseq->alphabet);
    gt_encseq_delete(condenseq->unique_es);
    gt_free(condenseq->buffer);
    gt_free(condenseq->filename);
    gt_free(condenseq->links);
    gt_free(condenseq->orig_ids);
    gt_free(condenseq->ubuffer);
    gt_free(condenseq->uniques);
    gt_intset_delete(condenseq->sdstab);
    gt_intset_delete(condenseq->ssptab);

    gt_free(condenseq);
  }
}

#define gt_condenseq_io_one(elem) \
  io_func(&elem, sizeof (elem), (size_t) 1, fp, err)

/*generic IO function for one condenseq link entry*/
static int condenseq_linkentry_io(GtCondenseqLink* link,
                                  FILE* fp,
                                  GtIOFunc io_func,
                                  GtError *err)
{
  int had_err = 0;
  gt_assert (link != NULL);
  had_err = gt_condenseq_io_one(link->orig_startpos);
  if (!had_err)
    had_err = gt_condenseq_io_one(link->len);
  if (!had_err)
    had_err = gt_condenseq_io_one(link->unique_id);
  if (!had_err)
    had_err = gt_condenseq_io_one(link->unique_offset);
  if (!had_err) {
    link->editscript = gt_editscript_io(link->editscript, fp, err);
    if (link->editscript == NULL) {
      had_err = 1;
    }
  }
  return had_err;
}

/*generic IO function for one condenseq unique entry*/
static int condenseq_uniqueentry_io(GtCondenseqUnique* unique,
                                    FILE* fp,
                                    GtIOFunc io_func,
                                    GtError *err)
{
  int had_err = 0;
  gt_assert (unique != NULL);
  had_err = gt_condenseq_io_one(unique->orig_startpos);
  if (!had_err)
    had_err = gt_condenseq_io_one(unique->len);
  return had_err;
}

static int condenseq_io(GtCondenseq *condenseq,
                        FILE* fp,
                        GtIOFunc io_func,
                        GtError *err)
{
  int had_err = 0;
  int file_format = GT_CONDENSEQ_VERSION;
  GtUword idx;
  had_err = gt_condenseq_io_one(condenseq->orig_len);
  if (!had_err)
    had_err = gt_condenseq_io_one(file_format);
  if (!had_err && file_format != GT_CONDENSEQ_VERSION) {
    gt_error_set(err, "condenseq index is format version %d, current is "
                 "%d -- please re-encode",
                 file_format, GT_CONDENSEQ_VERSION);
    had_err = -1;
  }
  if (!had_err)
    had_err = gt_condenseq_io_one(condenseq->orig_num_seq);
  if (!had_err)
    had_err = gt_condenseq_io_one(condenseq->ldb_nelems);
  if (!had_err) {
    if (condenseq->ldb_nelems == 0) {
      gt_warning("compression of condenseq did not succeed in finding any "
                 "compressable similarities, maybe the input is to small or "
                 "the chosen parameters should be reconsidered.");
    }
    if (condenseq->links == NULL) {
      condenseq->links = gt_calloc((size_t) condenseq->ldb_nelems,
                                   sizeof (*condenseq->links));
      condenseq->ldb_allocated = condenseq->ldb_nelems;
    }

    had_err = gt_condenseq_io_one(condenseq->udb_nelems);
  }

  if (!had_err) {
    gt_assert(condenseq->udb_nelems > 0);

    if (condenseq->uniques == NULL) {
      condenseq->uniques = gt_malloc(sizeof (*condenseq->uniques) *
                                     condenseq->udb_nelems );
      condenseq->udb_allocated = condenseq->udb_nelems;
    }
  }

  for (idx = 0; !had_err && idx < condenseq->ldb_nelems; idx++) {
    had_err = condenseq_linkentry_io(&condenseq->links[idx], fp, io_func, err);
  }

  for (idx = 0; !had_err && idx < condenseq->udb_nelems; idx++) {
    had_err = condenseq_uniqueentry_io(&condenseq->uniques[idx], fp, io_func,
                                       err);
  }
  if (!had_err && condenseq->orig_num_seq > (GtUword) 1) {
    condenseq->ssptab = gt_intset_io(condenseq->ssptab, fp, err);
    if (condenseq->ssptab == NULL)
      had_err = 1;
  }
  if (!had_err)
    had_err = gt_condenseq_io_one(condenseq->id_len);
  if (!had_err) {
    if (condenseq->id_len == GT_UNDEF_UWORD) {
      condenseq->sdstab = gt_intset_io(condenseq->sdstab, fp, err);
      if (condenseq->sdstab == NULL)
        had_err = 1;
    }
  }
  if (!had_err)
    had_err = gt_condenseq_io_one(condenseq->ids_total_len);
  if (!had_err) {
    condenseq->orig_ids = gt_realloc(condenseq->orig_ids,
                                     (size_t) condenseq->ids_total_len);
    had_err = io_func(condenseq->orig_ids, sizeof (*condenseq->orig_ids),
                      (size_t) condenseq->ids_total_len, fp, err);
  }
  return had_err;
}

int gt_condenseq_write(GtCondenseq* condenseq, FILE* fp, GtError *err)
{
  gt_assert(condenseq != NULL && fp != NULL);
  return condenseq_io(condenseq, fp, gt_io_error_fwrite, err);
}

/*read condenseq data structure from file*/
GtCondenseq *gt_condenseq_new_from_file(const char *indexname,
                                        GtLogger *logger, GtError *err)
{
  int had_err = 0;
  FILE* fp;
  GtEncseqLoader *esl;
  GtEncseq *unique_es;
  GtCondenseq *condenseq = NULL;
  /*load unique_es*/
  esl = gt_encseq_loader_new();
  gt_encseq_loader_disable_autosupport(esl);
  gt_encseq_loader_drop_md5_support(esl);
  gt_encseq_loader_require_ssp_tab(esl);
  unique_es = gt_encseq_loader_load(esl, indexname, err);
  if (!unique_es)
    had_err = -1;
  if (!had_err) {
    gt_encseq_loader_delete(esl);
    condenseq = condenseq_new_empty(gt_encseq_alphabet(unique_es));
    condenseq->filename = gt_cstr_dup(indexname);
    condenseq->unique_es = unique_es;
    fp = gt_fa_fopen_with_suffix(indexname, GT_CONDENSEQ_FILE_SUFFIX,
                                 "rb", err);
    if (fp == NULL) {
      had_err = -1;
    }
    else {
      had_err = condenseq_io(condenseq, fp, gt_io_error_fread, err);
      if (!had_err) {
        GtUword i;
        gt_assert(condenseq->uniques);
        gt_assert(condenseq->links);
        gt_fa_fclose(fp);
        /*create link array for each unique entry*/
        for (i = 0; i < condenseq->udb_nelems; i++) {
          GT_INITARRAY(&(condenseq->uniques[i].links),uint32_t);
        }
        /* check for overflows */
        if (condenseq->ldb_nelems > (GtUword) ((uint32_t) 0 - (uint32_t) 1)) {
          gt_error_set(err, "Overflow, to many link-elements. Can't be stored");
          had_err = -1;
        }
        /* iterate through link entrys and store ids in corresponding unique
          entry array */
        for (i = 0; !had_err && (GtUword) i < condenseq->ldb_nelems; i++) {
          GtUword uid = condenseq->links[i].unique_id;
          gt_assert(uid < condenseq->udb_nelems);
          GT_STOREINARRAY(&(condenseq->uniques[uid].links),
                          uint32_t,
                          10,
                          (uint32_t) i);
        }
      }
    }
  }
  if (!had_err) {
    gt_assert(condenseq != NULL);
    if (condenseq->id_len != GT_UNDEF_UWORD)
      gt_logger_log(logger, "IDs const len: " GT_WU, condenseq->id_len);
    else
      gt_logger_log(logger, "using sdstab to access IDs");
  }
  if (had_err) {
    gt_condenseq_delete(condenseq);
    condenseq = NULL;
  }
  return (condenseq);
}

char *gt_condenseq_basefilename(const GtCondenseq *condenseq)
{
  char *basename = NULL,
       *suffix_ptr;
  if (condenseq->filename != NULL) {
    basename = gt_basename(condenseq->filename);
    if (strlen(basename) > (size_t) 1 &&
        (suffix_ptr = strrchr(basename + 1, '.')) != NULL) {
      /* remove suffix */
      *suffix_ptr = '\0';
    }
  }
  return basename;
}

GtStr *gt_condenseq_unique_fasta_file(const GtCondenseq *condenseq)
{
  GtStr *unique = NULL;
  gt_assert(condenseq != NULL);
  if (condenseq->filename != NULL) {
    unique = gt_str_new_cstr(condenseq->filename);
    gt_str_append_cstr(unique, ".fas");
  }
  return unique;
}

static GtUword condenseq_links_position_binsearch(const GtCondenseq *condenseq,
                                                  GtUword position)
{
  GtWord idx, low, high;
  gt_assert(condenseq && condenseq->ldb_nelems > 0);
  low = (GtWord) -1;
  gt_safe_assign(high, condenseq->ldb_nelems);
  idx = GT_DIV2(low + high);
  while (high - low > (GtWord) 1) {
    if (position < condenseq->links[idx].orig_startpos) {
      high = idx;
    }
    else {
      low = idx;
    }
    idx = GT_DIV2(low + high);
  }
  if (low > (GtWord) -1 && condenseq->links[idx].orig_startpos <= position)
    return (GtUword) idx;
  return condenseq->ldb_nelems;
}

GtUword gt_condenseq_uniques_position_binsearch(const GtCondenseq *condenseq,
                                                GtUword position)
{
  GtWord idx, low, high;
  gt_assert(condenseq && condenseq->udb_nelems > 0);
  low = (GtWord) -1;
  gt_safe_assign(high, condenseq->udb_nelems);
  idx = GT_DIV2(low + high);
  while (high - low > (GtWord) 1) {
    if (position < condenseq->uniques[idx].orig_startpos) {
      high = idx;
    }
    else {
      low = idx;
    }
    idx = GT_DIV2(low + high);
  }
  if (low > (GtWord) -1 && condenseq->uniques[idx].orig_startpos <= position)
    return (GtUword) idx;
  return condenseq->udb_nelems;
}

static GtUword condenseq_unique_extract_encoded(const GtCondenseq *cs,
                                                GtUword id,
                                                GtUchar *buffer,
                                                GtUword frompos,
                                                GtUword topos)
{
  GtCondenseqUnique unique = cs->uniques[id];
  GtUword startoffset,
          startpos,
          uniquelength,
          targetlength,
          endpos;
  gt_assert(unique.orig_startpos <= frompos);
  startoffset = frompos - unique.orig_startpos;
  gt_assert(startoffset < unique.len);
  startpos = gt_encseq_seqstartpos(cs->unique_es, id) + startoffset;
  uniquelength = unique.len - startoffset;
  targetlength = topos - frompos + 1;
  if (uniquelength < targetlength)
    endpos = startpos + uniquelength - 1;
  else
    endpos = startpos + targetlength - 1;

  gt_encseq_extract_encoded(cs->unique_es, buffer, startpos, endpos);
  return endpos - startpos + 1;
}

static GtUword condenseq_link_extract_encoded(const GtCondenseq *cs,
                                              GtUword id,
                                              GtUchar *buffer,
                                              GtUword frompos,
                                              GtUword topos)
{
  GtCondenseqLink link = cs->links[id];
  GtEditscript *editscript = link.editscript;
  GtUword unique_startpos,
          targetlength,
          startoffset,
          endpos,
          linklength,
          written;
  gt_assert(link.orig_startpos <= frompos);
  unique_startpos = gt_encseq_seqstartpos(cs->unique_es, link.unique_id);
  startoffset = frompos - link.orig_startpos;
  gt_assert(startoffset < link.len);
  linklength = link.len - startoffset;
  targetlength = topos - frompos + 1;
  if (linklength < targetlength)
    endpos = link.len - 1;
  else
    endpos = startoffset + targetlength - 1;
  written =
    gt_editscript_get_sub_sequence_v(editscript, cs->unique_es,
                                     unique_startpos + link.unique_offset,
                                     GT_READMODE_FORWARD, startoffset,
                                     endpos, buffer);
  gt_assert(written == endpos - startoffset + 1);
  return written;
}

const GtUchar *gt_condenseq_extract_encoded_range(GtCondenseq *condenseq,
                                                  GtRange range)
{
  GtUchar *buf;
  GtUword nextsep,
          linkid = 0,
          uniqueid,
          buffoffset = 0,
          length;
  GtCondenseqLink *link = NULL;
  GtCondenseqUnique *unique = NULL;

  gt_assert(condenseq && condenseq->udb_nelems != 0);
  gt_assert(condenseq->uniques[0].orig_startpos == 0);
  gt_assert(range.start <= range.end);
  gt_assert(range.end < condenseq->orig_len);

  nextsep = condenseq_next_sep(condenseq, range.start);
  uniqueid = gt_condenseq_uniques_position_binsearch(condenseq,
                                                     range.start);

  length = range.end - range.start + 1;

  /* TODO DW check if there is another way than using this buffer, so we could
     use const for condenseq here. */
  if (condenseq->ubuffer == NULL || condenseq->ubuffsize < length) {
    condenseq->ubuffer = gt_realloc(condenseq->ubuffer,
                                    sizeof (*condenseq->ubuffer) * length);
    condenseq->ubuffsize = length;
  }

  buf = condenseq->ubuffer;

  unique = &condenseq->uniques[uniqueid];

  if (unique->orig_startpos + unique->len <= range.start) {
    uniqueid++;
    if (uniqueid == condenseq->udb_nelems)
      unique = NULL;
    else {
      unique = &condenseq->uniques[uniqueid];
      gt_assert(unique->orig_startpos + unique->len > range.start);
    }
  }

  if (condenseq->ldb_nelems != 0) {
      linkid = condenseq_links_position_binsearch(condenseq, range.start);
    if (linkid == condenseq->ldb_nelems)
      linkid = 0;
    link = &condenseq->links[linkid];
  }
  if (link != NULL &&
      link->orig_startpos + link->len <= range.start) {
    linkid++;
    if (linkid == condenseq->ldb_nelems)
      link = NULL;
    else {
      link = &condenseq->links[linkid];
      gt_assert(link->orig_startpos + link->len > range.start);
    }
  }

  while (buffoffset < length) {
    gt_assert(nextsep == range.start + buffoffset ||
              unique != NULL ||
              link != NULL);
    if (nextsep != 0 && nextsep == range.start + buffoffset) {
      buf[buffoffset++] = SEPARATOR;
      nextsep = condenseq_next_sep(condenseq, range.start + buffoffset);
    }
    else if (unique != NULL &&
             (link == NULL || unique->orig_startpos < link->orig_startpos)) {
      buffoffset += condenseq_unique_extract_encoded(condenseq, uniqueid,
                                                     buf + buffoffset,
                                                     range.start + buffoffset,
                                                     range.end);
      if (++uniqueid == condenseq->udb_nelems)
        unique = NULL;
      else {
        unique = &condenseq->uniques[uniqueid];
        gt_assert(unique->orig_startpos + unique->len > range.start);
      }
    }
    else {
      gt_assert(link != NULL);
      buffoffset += condenseq_link_extract_encoded(condenseq, linkid,
                                                   buf + buffoffset,
                                                   range.start + buffoffset,
                                                   range.end);
      if (++linkid == condenseq->ldb_nelems)
        link = NULL;
      else {
        link = &condenseq->links[linkid];
        gt_assert(link->orig_startpos + link->len > range.start);
      }
    }
  }
  gt_assert(buffoffset == length);
  return buf;
}

const GtUchar *gt_condenseq_extract_encoded(GtCondenseq *condenseq,
                                            GtUword *length,
                                            GtUword id)
{
  GtRange range;
  range.start = gt_condenseq_seqstartpos(condenseq, id);
  if (id < condenseq->orig_num_seq - 1)
    /* -2 because of seperator */
    range.end = gt_condenseq_seqstartpos(condenseq, id + 1) - 2;
  else
    range.end = condenseq->orig_len - 1;
  *length = range.end - range.start + 1;
  return gt_condenseq_extract_encoded_range(condenseq, range);
}

const char *gt_condenseq_extract_decoded_range(GtCondenseq *condenseq,
                                               GtRange range,
                                               char separator)
{
  GtUword length = range.end - range.start + 1,
          idx;
  const GtUchar *ubuf;
  char *buf;
  gt_assert(range.start <= range.end);
  ubuf = gt_condenseq_extract_encoded_range(condenseq, range);
  if (condenseq->buffer == NULL || condenseq->buffsize < length) {
    condenseq->buffer = gt_realloc(condenseq->buffer,
                                   sizeof (*condenseq->buffer) * length);
    condenseq->buffsize = length;
  }
  buf = condenseq->buffer;
  for (idx = 0; idx < length; ++idx) {
    if (ubuf[idx] == SEPARATOR) {
      buf[idx] = separator;
    }
    else {
      buf[idx] = gt_alphabet_decode(condenseq->alphabet, ubuf[idx]);
    }
  }
  return buf;
}

const char *gt_condenseq_extract_decoded(GtCondenseq *condenseq,
                                         GtUword *length,
                                         GtUword id)
{
  GtRange range;
  range.start = gt_condenseq_seqstartpos(condenseq, id);
  if (id < condenseq->orig_num_seq - 1)
    /* -2 because of seperator */
    range.end = gt_condenseq_seqstartpos(condenseq, id + 1) - 2;
  else
    range.end = condenseq->orig_len - 1;
  *length = range.end - range.start + 1;
  return gt_condenseq_extract_decoded_range(condenseq, range, '\0');
}

/* static GtRange
gt_condenseq_convert_unique_range_to_global(const GtCondenseq *condenseq,
                                            GtUword unique_id,
                                            GtRange range)
{
  GtRange ret;
  GtCondenseqUnique *unique;
  gt_assert(condenseq != NULL);
  gt_assert(unique_id < condenseq->udb_nelems);
  unique = &condenseq->uniques[unique_id];
  ret.start = unique->orig_startpos + range.start;
  ret.end = unique->orig_startpos + range.end;
  return ret;
} */

GtUword gt_condenseq_each_redundant_seq(
                                       const GtCondenseq *condenseq,
                                       GtUword uid,
                                       GtCondenseqProcessExtractedSeqs callback,
                                       void *callback_data,
                                       GtError *err)
{
  int had_err = 0;
  GtUword num_seqs = (GtUword) 1, linkidx,
          orig_seqnum;
  const GtCondenseqUnique unique = condenseq->uniques[uid];

  orig_seqnum = gt_condenseq_pos2seqnum(condenseq, unique.orig_startpos);

  had_err = callback(callback_data, orig_seqnum, err);

  for (linkidx = 0;
       !had_err && linkidx < unique.links.nextfreeuint32_t;
       ++linkidx) {
    GtUword linkid = (GtUword) unique.links.spaceuint32_t[linkidx];
    const GtCondenseqLink link = condenseq->links[linkid];
    orig_seqnum = gt_condenseq_pos2seqnum(condenseq, link.orig_startpos);
    had_err = callback(callback_data, orig_seqnum, err);
    num_seqs++;
  }
  if (!had_err)
    return num_seqs;
  return (GtUword) had_err;
}

GtUword gt_condenseq_each_redundant_range(
                                      const GtCondenseq *condenseq,
                                      GtUword uid,
                                      GtRange urange,
                                      GtUword left_extend,
                                      GtUword right_extend,
                                      GtCondenseqProcessExtractedRange callback,
                                      void *callback_data,
                                      GtError *err)
{
  int had_err = 0;
  GtUword num_ranges = (GtUword) 1,
          linkidx,
          orig_seqnum,
          orig_seqstart,
          orig_seqend;
  const GtCondenseqUnique *unique;
  GtRange extract;

  gt_assert(condenseq != NULL);
  gt_assert(uid < condenseq->udb_nelems);
  gt_assert(urange.start <= urange.end);

  unique = &condenseq->uniques[uid];

  /* handle unique itself */
  orig_seqnum = gt_condenseq_pos2seqnum(condenseq, unique->orig_startpos);
  orig_seqstart = gt_condenseq_seqstartpos(condenseq, orig_seqnum);
  orig_seqend = orig_seqstart + condenseq_seqlength_help(condenseq, orig_seqnum,
                                                         orig_seqstart) - 1;
  extract.start = unique->orig_startpos + urange.start;
  extract.start = extract.start < left_extend ?
    0 : extract.start - left_extend;
  extract.start = extract.start < orig_seqstart ?
    orig_seqstart : extract.start;
  extract.end = unique->orig_startpos + urange.end + right_extend;
  extract.end = extract.end > orig_seqend ?
    orig_seqend : extract.end;

  gt_assert(extract.start <= extract.end);
  had_err = callback(callback_data, orig_seqnum, extract, err);

  for (linkidx = 0;
       !had_err && linkidx < unique->links.nextfreeuint32_t;
       ++linkidx) {
    const GtCondenseqLink *link =
      &condenseq->links[unique->links.spaceuint32_t[linkidx]];

    /* the second part is a little heuristic, the len of the link could be
       completely comprised of insertions, which would place it downstream of
       urange.start. But we assume ~ similar lengths for links and their unique
       counterpart */
    if (!(urange.end < link->unique_offset ||
        urange.start > link->unique_offset + link->len - 1)) {
      GtUword shift;
      orig_seqnum = gt_condenseq_pos2seqnum(condenseq, link->orig_startpos);
      orig_seqstart = gt_condenseq_seqstartpos(condenseq, orig_seqnum);
      orig_seqend = orig_seqstart + condenseq_seqlength_help(condenseq,
                                                             orig_seqnum,
                                                             orig_seqstart) - 1;
      extract.start = link->orig_startpos < left_extend ?
        0 : link->orig_startpos - left_extend;
      if (urange.start < link->unique_offset) {
        shift = link->unique_offset - urange.start;
        extract.start = extract.start < shift ?
          0 : extract.start - shift;
      }
      else {
        shift = urange.start - link->unique_offset;
        extract.start += shift;
      }
      extract.start = extract.start < orig_seqstart ?
        orig_seqstart : extract.start;
      /* see heuristic note above */
      extract.end = link->orig_startpos + right_extend + link->len;
      if (urange.end < link->unique_offset + link->len - 1) {
        shift = (link->unique_offset + link->len - 1) - urange.end;
        extract.end = extract.end < shift ?
          0 : extract.end - shift;
      }
      else {
        shift = urange.end - (link->unique_offset + link->len - 1);
        extract.end += shift;
      }
      extract.end = extract.end > orig_seqend ?
        orig_seqend : extract.end;
      gt_assert(extract.start <= extract.end);
      had_err = callback(callback_data, orig_seqnum, extract, err);
      num_ranges++;
    }
  }

  if (!had_err)
    return num_ranges;
  return 0;
}

int gt_condenseq_output_to_gff3(const GtCondenseq *condenseq,
                                GtError *err)
{
  int had_err = 0;
  GtUword idx,
          name_len,
          seqnum = 0, seqstart = 0, seqend = 0,
          desclen;
  GtStr *filename = NULL,
        *id = gt_str_new_cstr("U"),
        *name = gt_str_new_cstr("unique"),
        *parent_unique = gt_str_new_cstr("U"),
        *seqid = gt_str_new(),
        *source = gt_str_new_cstr("Condenseq");
  GtFile *outfile = NULL;
  GtGFF3Visitor *gffv = NULL;
  GtNodeVisitor *nodev = NULL;
  GtFeatureNode *fnode = NULL;
  GtGenomeNode *node = NULL;
  GtRange range;

  gt_assert(condenseq != NULL);

  filename = gt_str_new_cstr(gt_condenseq_basefilename(condenseq));

  name_len = gt_str_length(name);
  gt_str_append_cstr(filename, ".gff3");
  outfile = gt_file_new(gt_str_get(filename), "w", err);
  nodev = gt_gff3_visitor_new(outfile);
  gffv = (GtGFF3Visitor *) nodev;
  gt_gff3_visitor_retain_id_attributes(gffv);

  node = gt_feature_node_new(seqid, "experimental_feature", (GtUword) 1,
                             (GtUword) 1, GT_STRAND_BOTH);
  fnode = (GtFeatureNode*) node;
  gt_feature_node_set_source(fnode, source);
  for (idx = 0; !had_err && idx < condenseq->udb_nelems; ++idx) {
    GtCondenseqUnique uq = condenseq->uniques[idx];
    if (seqend <= uq.orig_startpos) {
      const char *desc;
      gt_genome_node_delete(node);
      seqnum = gt_condenseq_pos2seqnum(condenseq, uq.orig_startpos);
      seqstart = gt_condenseq_seqstartpos(condenseq, seqnum);
      seqend = seqstart + condenseq_seqlength_help(condenseq, seqnum, seqstart);
      desc = gt_condenseq_description(condenseq, &desclen, seqnum);
      gt_str_reset(seqid);
      gt_str_append_cstr_nt(seqid, desc, desclen);
      node = gt_feature_node_new(seqid, "experimental_feature", (GtUword) 1,
                                 (GtUword) 1, GT_STRAND_BOTH);
      fnode = (GtFeatureNode*) node;
      gt_feature_node_set_source(fnode, source);
    }
    gt_str_set_length(name, name_len);
    gt_str_append_uword(name, idx);
    gt_str_set_length(id, (GtUword) 1);
    gt_str_append_uword(id, idx);
    gt_feature_node_set_attribute(fnode, "Name", gt_str_get(name));
    gt_feature_node_set_attribute(fnode, "ID", gt_str_get(id));
    /* 1 Based coordinates! */
    range.start = uq.orig_startpos + 1 - seqstart;
    range.end = uq.orig_startpos + uq.len - seqstart;
    gt_genome_node_set_range(node, &range);
    had_err = gt_genome_node_accept(node, nodev, err);
  }
  gt_str_reset(name);
  gt_str_append_cstr(name, "link");
  gt_str_reset(id);
  gt_str_append_cstr(id, "L");
  name_len = gt_str_length(name);
  seqend = 0;
  for (idx = 0; !had_err && idx < condenseq->ldb_nelems; ++idx) {
    GtCondenseqLink link = condenseq->links[idx];
    if (seqend <= link.orig_startpos) {
      const char *desc;
      gt_genome_node_delete(node);
      seqnum = gt_condenseq_pos2seqnum(condenseq, link.orig_startpos);
      seqstart = gt_condenseq_seqstartpos(condenseq, seqnum);
      seqend = seqstart + condenseq_seqlength_help(condenseq, seqnum, seqstart);
      desc = gt_condenseq_description(condenseq, &desclen, seqnum);
      gt_str_reset(seqid);
      gt_str_append_cstr_nt(seqid, desc, desclen);
      node = gt_feature_node_new(seqid, "experimental_feature", (GtUword) 1,
                                 (GtUword) 1, GT_STRAND_BOTH);
      fnode = (GtFeatureNode*) node;
      gt_feature_node_set_source(fnode, source);
    }
    gt_str_set_length(name, name_len);
    gt_str_append_uword(name, idx);
    gt_str_set_length(id, (GtUword) 1);
    gt_str_append_uword(id, idx);
    gt_feature_node_set_attribute(fnode, "Name", gt_str_get(name));
    gt_feature_node_set_attribute(fnode, "ID", gt_str_get(id));
    gt_str_set_length(parent_unique, (GtUword) 1);
    gt_str_append_uword(parent_unique, link.unique_id);
    gt_feature_node_set_attribute(fnode, "Derives_from",
                                  gt_str_get(parent_unique));
    /* 1 Based coordinates! */
    range.start = link.orig_startpos + 1 - seqstart;
    range.end = link.orig_startpos + link.len - seqstart;
    gt_genome_node_set_range(node, &range);
    had_err = gt_genome_node_accept(node, nodev, err);
  }
  gt_file_delete(outfile);
  gt_genome_node_delete(node);
  gt_node_visitor_delete(nodev);
  gt_str_delete(filename);
  gt_str_delete(id);
  gt_str_delete(name);
  gt_str_delete(parent_unique);
  gt_str_delete(seqid);
  gt_str_delete(source);
  return had_err;
}

GtDiscDistri *gt_condenseq_unique_length_dist(const GtCondenseq *condenseq)
{
  GtUword idx;
  GtDiscDistri *res = gt_disc_distri_new();

  for (idx = 0; idx < condenseq->udb_nelems; idx++) {
    gt_disc_distri_add(res, condenseq->uniques[idx].len);
  }
  return res;
}

GtDiscDistri *gt_condenseq_link_length_dist(const GtCondenseq *condenseq)
{
  GtUword idx;
  GtDiscDistri *res = gt_disc_distri_new();

  for (idx = 0; idx < condenseq->ldb_nelems; idx++) {
    gt_disc_distri_add(res, condenseq->links[idx].len);
  }
  return res;
}

GtDiscDistri *gt_condenseq_link_comp_dist(const GtCondenseq *condenseq)
{
  GtUword idx;
  GtDiscDistri *res = gt_disc_distri_new();

  for (idx = 0; idx < condenseq->ldb_nelems; idx++) {
    GtEditscript *es = condenseq->links[idx].editscript;
    GtUword vlen;
    size_t size;
    vlen = gt_editscript_get_target_len(es);
    size = gt_editscript_size(es);
    gt_disc_distri_add(res, (GtUword) ((double) size/(double) vlen * 100));
  }
  return res;
}

GtUword gt_condenseq_unique_range_to_seqrange(const GtCondenseq *condenseq,
                                              GtUword uid,
                                              GtRange *urange)
{
  GtUword seqnum = 0, seqstart = 0;
  GtCondenseqUnique uq;

  gt_assert(condenseq != NULL);
  gt_assert(uid < condenseq->udb_nelems);

  uq = condenseq->uniques[uid];
  seqnum = gt_condenseq_pos2seqnum(condenseq, uq.orig_startpos);
  seqstart = gt_condenseq_seqstartpos(condenseq, seqnum);
  urange->start += uq.orig_startpos - seqstart;
  urange->end += uq.orig_startpos - seqstart;
  return seqnum;
}

/* static GtUword gt_condenseq_link_range_to_seqrange(GtCondenseq *condenseq,
                                                   GtUword lid,
                                                   GtRange *lrange)
{
  GtUword seqnum = 0, seqstart = 0;
  GtCondenseqLink lk;

  gt_assert(condenseq != NULL);
  gt_assert(lid < condenseq->ldb_nelems);

  lk = condenseq->links[lid];
  seqnum = gt_condenseq_pos2seqnum(condenseq, lk.orig_startpos);
  seqstart = gt_condenseq_seqstartpos(condenseq, seqnum);
  lrange->start = lk.orig_startpos - seqstart;
  lrange->end = lk.orig_startpos - seqstart;
  return seqnum;
} */

const GtEditscript *gt_condenseq_link_editscript(const GtCondenseq *condenseq,
                                                 GtUword lid)
{
  gt_assert(lid < condenseq->ldb_nelems);
  return condenseq->links[lid].editscript;
}

GtAlphabet *gt_condenseq_alphabet(const GtCondenseq *condenseq)
{
  return gt_alphabet_ref(condenseq->alphabet);
}

GtUword gt_condenseq_count_relevant_uniques(const GtCondenseq *condenseq,
                                            unsigned int min_align_len)
{
  GtUword idx, count = 0;
  for (idx = 0; idx < condenseq->udb_nelems; idx++) {
    if (condenseq->uniques[idx].len >= min_align_len)
      count++;
  }
  return count;
}

GtUword gt_condenseq_size(const GtCondenseq *condenseq,
                          GtUword *uniques,
                          GtUword *links,
                          GtUword *editscripts,
                          GtUword *descriptions,
                          GtUword *separators) {
  GtUword idx;
  *uniques = condenseq->udb_nelems * sizeof (*condenseq->uniques);
  for (idx = 0; idx < condenseq->udb_nelems; idx++) {
    *uniques += condenseq->uniques[idx].links.allocateduint32_t *
      sizeof (*condenseq->uniques[idx].links.spaceuint32_t);
  }
  *links = condenseq->ldb_nelems * sizeof (*condenseq->links);
  *editscripts = 0;
  for (idx = 0; idx < condenseq->ldb_nelems; idx++)
    *editscripts += gt_editscript_size(condenseq->links[idx].editscript);
  *descriptions = condenseq->ids_total_len;
  *descriptions += gt_intset_size_of_struct(condenseq->sdstab);
  *descriptions += gt_intset_size_of_rep(condenseq->sdstab);
  *separators = gt_intset_size_of_struct(condenseq->ssptab);
  *separators = gt_intset_size_of_rep(condenseq->ssptab);
  return *uniques + *links + *editscripts + *descriptions + *separators;
}
