/*
  Copyright (c) 2014 Dirk Willrodt <willrodt@zbh.uni-hamburg.de>
  Copyright (c) 2014 Florian Markowsky <1markows@informatik.uni-hamburg.de>
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

#ifndef S_SPLINT_S
#include <ctype.h>
#endif
#include <string.h>

#include "core/arraydef.h"
#include "core/divmodmul.h"
#include "core/encseq_api.h"
#include "core/ensure.h"
#include "core/error_api.h"
#include "core/fa.h"
#include "core/hashmap_api.h"
#include "core/log_api.h"
#include "core/logger_api.h"
#include "core/ma.h"
#include "core/qsort_r_api.h"
#include "core/range.h"
#include "core/readmode_api.h"
#include "core/safearith.h"
#include "core/str_array.h"
#include "core/str_array_api.h"
#include "core/types_api.h"
#include "core/undef_api.h"
#include "core/unused_api.h"
#include "core/xansi_api.h"
#include "extended/editscript.h"
#include "extended/intset.h"
#include "extended/n_r_encseq.h"
#include "extended/io_function_pointers.h"
#include "match/seqabstract.h"
#include "match/sfx-mappedstr.h"
#include "match/xdrop.h"

typedef struct GtNREncseqLink {
  GtEditscript *editscript;
  GtUword       len,
                orig_startpos,
                unique_id,
                unique_offset;
} GtNREncseqLink;

typedef struct GtNREncseqUnique {
  GtArrayuint32_t links;
  GtUword         len,
                  orig_startpos;
} GtNREncseqUnique;

struct GtNREncseq {
  GtAlphabet       *alphabet;
  GtEncseq         *unique_es;
  GtIntset         *sdstab;
  GtIntset         *ssptab;
  GtNREncseqLink   *links;
  GtNREncseqUnique *uniques;
  char             *orig_ids;
  GtUword           id_len, /* GT_UNDEF_UWORD if sdstab != NULL */
                    ids_total_len,
                    ldb_allocated,
                    ldb_nelems,
                    orig_length,
                    orig_num_seq,
                    udb_allocated,
                    udb_nelems;
};

static inline
GtUword gt_n_r_encseq_ssp_pos2seqnum(const GtNREncseq *nre, GtUword pos)
{
  GtUword ret = gt_intset_get_idx_smallest_geq(nre->ssptab, pos);
  return ret;
}

static inline GtUword gt_n_r_encseq_ssp_seqstartpos(const GtNREncseq *nre,
                                                    GtUword seqnum)
{
  GtUword ret;
  if (seqnum == 0)
    ret = 0;
  else
    ret = gt_intset_get(nre->ssptab, seqnum - 1) + 1;
  return ret;
}

static inline GtUword gt_n_r_encseq_ssp_seqlength(const GtNREncseq *nre,
                                                  GtUword seqnum)
{
  GtUword start = 0, end = nre->orig_length, ret;
  if (seqnum != 0)
    start = gt_intset_get(nre->ssptab, seqnum - 1) + 1;
  if (seqnum < nre->orig_num_seq - 1)
    end = gt_intset_get(nre->ssptab, seqnum);
  ret = end - start;
  return ret;
}

static GtIntset *gt_n_r_encseq_fill_ssp_tab(GtNREncseq *nre,
                                            const GtEncseq *orig_es)
{
  GtIntset *ssptab;
  GtUword max, idx;
  max = gt_encseq_seqstartpos(orig_es, nre->orig_num_seq - 1);
  /* we store the internal separators, the end is explicit */
  ssptab = gt_intset_best_new(max, nre->orig_num_seq - 1);
  for (idx = (GtUword) 1; idx < nre->orig_num_seq; ++idx) {
    GtUword pos = gt_encseq_seqstartpos(orig_es, idx) - 1;
    gt_assert(pos != 0);
    gt_intset_add(ssptab, pos);
  }
  return ssptab;
}

static inline GtUword gt_n_r_encseq_idlen(const char *desc, GtUword desclen)
{
  GtUword idx;
  for (idx = 0; idx < desclen; ++idx) {
    if (isspace(desc[idx]) || desc[idx] == '\0')
      return idx;
  }
  return desclen;
}

static void gt_n_r_encseq_process_ids(GtNREncseq *nre, const GtEncseq *orig_es,
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

  nre->ids_total_len = 0;
  dist = gt_calloc((size_t) distsize, sizeof (*dist));

  for (idx = 0; idx < nre->orig_num_seq; ++idx) {
    desc = gt_encseq_description(orig_es, &desclen, idx);
    idlen = gt_n_r_encseq_idlen(desc, desclen);
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
    nre->ids_total_len += dist[dist_idx] * dist_idx;
  }
  nre->ids_total_len += dist_idx * dist[dist_idx];

  sdssize = (GtUword) gt_intset_best_memory_size(maxendidx, nre->orig_num_seq);
  use_const_len = wastedmem < sdssize;

  if (use_const_len) {
    gt_logger_log(logger, "will use const len, " GT_WU ", \"wasting\" " GT_WU
                  " bytes. SDS would use " GT_WU " bytes",
                  maxlen, wastedmem, sdssize);
    nre->id_len = maxlen;
    nre->ids_total_len = maxlen * nre->orig_num_seq;
  }
  else {
    gt_logger_log(logger, "will use sdstab with size " GT_WU ". Would have "
                  "wasted " GT_WU " bytes.", sdssize, wastedmem);
    nre->sdstab = gt_intset_best_new(maxendidx, nre->orig_num_seq);
  }
  nre->orig_ids = gt_calloc((size_t) nre->ids_total_len,
                            sizeof (*nre->orig_ids));

  cur_id_startptr = nre->orig_ids;
  for (idx = 0; idx < nre->orig_num_seq; ++idx) {
    desc = gt_encseq_description(orig_es, &desclen, idx);
    idlen = gt_n_r_encseq_idlen(desc, desclen);
    gt_assert(idlen <= maxlen);
    (void) memcpy(cur_id_startptr, desc, (size_t) idlen);
    if (use_const_len) {
      cur_id_startptr += maxlen;
      cur_total_id_len += maxlen;
    }
    else {
      cur_id_startptr += idlen;
      cur_total_id_len += idlen;
      gt_intset_add(nre->sdstab, cur_total_id_len);
    }
  }
  gt_assert(cur_total_id_len == nre->ids_total_len);
  gt_free(dist);
}

static const char *gt_n_r_encseq_sdstab_get_id(const GtNREncseq *nre,
                                               GtUword *idlen,
                                               GtUword idx)
{
  gt_assert(idx < nre->orig_num_seq);
  if (nre->id_len == GT_UNDEF_UWORD) {
    GtUword this = gt_intset_get(nre->sdstab, idx), previous;
    if (idx == 0) {
      *idlen = this;
      return nre->orig_ids;
    }
    previous = gt_intset_get(nre->sdstab, idx - 1);
    *idlen =  this - previous;
    return nre->orig_ids + previous;
  }
  else {
    const char *id = nre->orig_ids + idx * nre->id_len;
    *idlen = nre->id_len;
    while (id[*idlen - 1] == '\0')
      (*idlen)--;
    return nre->orig_ids + idx * nre->id_len;
  }
}

static GtNREncseq *gt_n_r_encseq_new_empty(const GtAlphabet *alph)
{
  GtNREncseq *n_r_encseq = gt_malloc(sizeof (*n_r_encseq));
  n_r_encseq->alphabet = gt_alphabet_ref((GtAlphabet *) alph);

  n_r_encseq->ldb_allocated =
    n_r_encseq->ldb_nelems =
    n_r_encseq->orig_length =
    n_r_encseq->orig_num_seq =
    n_r_encseq->udb_allocated =
    n_r_encseq->udb_nelems = 0;
  n_r_encseq->id_len = GT_UNDEF_UWORD;

  n_r_encseq->links = NULL;
  n_r_encseq->orig_ids = NULL;
  n_r_encseq->sdstab = NULL;
  n_r_encseq->ssptab = NULL;
  n_r_encseq->unique_es = NULL;
  n_r_encseq->uniques = NULL;

  return n_r_encseq;
}

static GtNREncseq *gt_n_r_encseq_new(const GtEncseq *orig_es, GtLogger *logger)
{
  GtNREncseq *n_r_encseq;
  n_r_encseq = gt_n_r_encseq_new_empty(gt_encseq_alphabet(orig_es));

  n_r_encseq->orig_num_seq = gt_encseq_num_of_sequences(orig_es);

  n_r_encseq->ssptab = gt_n_r_encseq_fill_ssp_tab(n_r_encseq, orig_es);
  n_r_encseq->orig_length = gt_encseq_total_length(orig_es);

  gt_n_r_encseq_process_ids(n_r_encseq, orig_es, logger);
  return n_r_encseq;
}

void gt_n_r_encseq_print_info(const GtNREncseq *n_r_encseq)
{
  GtUword total_link_len = 0,
          av_uque_len = 0,
          i;
  for (i = 0; i < n_r_encseq->udb_nelems; i++) {
    av_uque_len += n_r_encseq->uniques[i].len;
    gt_log_log("unique - len: " GT_WU ", orig: " GT_WU,
               n_r_encseq->uniques[i].len,
               n_r_encseq->uniques[i].orig_startpos);
  }
  printf(GT_WU " unique entries\n", n_r_encseq->udb_nelems);
  if (n_r_encseq->udb_nelems > 0) {
    av_uque_len = av_uque_len / n_r_encseq->udb_nelems;
    printf("average length of unique fragments: " GT_WU "\n", av_uque_len);
  }
  for (i = 0; i < n_r_encseq->ldb_nelems; i++) {
    total_link_len += n_r_encseq->links[i].len;
    gt_log_log("link - len: " GT_WU ", orig: " GT_WU ", u id: " GT_WU,
               n_r_encseq->links[i].len, n_r_encseq->links[i].orig_startpos,
               n_r_encseq->links[i].unique_id);
  }
  printf(GT_WU " link entries\n",n_r_encseq->ldb_nelems);
  if (n_r_encseq->ldb_nelems > 0) {
    printf("total length of link fragments: " GT_WU "\n", total_link_len);
    printf("average length of link fragments: " GT_WU "\n",
           total_link_len / n_r_encseq->ldb_nelems);
  }
  if (n_r_encseq->unique_es != NULL)
    printf("total udb encseq length:" GT_WU "\n",
           gt_encseq_total_length(n_r_encseq->unique_es));
  printf("original db length:" GT_WU "\n",
         n_r_encseq->orig_length);
}

GtUword gt_n_r_encseq_get_orig_length(GtNREncseq *n_r_encseq) {
  return n_r_encseq->orig_length;
}

GtUword gt_n_r_encseq_get_unique_length(GtNREncseq *n_r_encseq) {
  return gt_encseq_total_length(n_r_encseq->unique_es);
}

GtUword gt_n_r_encseq_array_size_increase(GtUword allocated)
{
  if (allocated != 0) {
    allocated += (allocated) > (GtUword) 128 ?
      (GtUword) 128 :
      (allocated);
  }
  else
    allocated = (GtUword) 16;
  return allocated;
}

static inline void gt_n_r_encseq_udb_resize(GtNREncseq *nre)
{
  if (nre->udb_nelems == nre->udb_allocated) {
    nre->udb_allocated = gt_n_r_encseq_array_size_increase(nre->udb_allocated);
    nre->uniques = gt_realloc(nre->uniques,
                              (size_t) nre->udb_allocated *
                              sizeof (*nre->uniques));
  }
}

static inline void gt_n_r_encseq_ldb_resize(GtNREncseq *nre)
{
  if (nre->ldb_nelems == nre->ldb_allocated) {
    nre->ldb_allocated = gt_n_r_encseq_array_size_increase(nre->ldb_allocated);
    nre->links = gt_realloc(nre->links,
                            (size_t) nre->ldb_allocated *
                            sizeof (*nre->links));
  }
}

static void gt_n_r_encseq_add_unique_to_db(GtNREncseq *nre,
                                           GtUword orig_startpos,
                                           GtUword len)
{
  gt_assert(len != 0);
  /* if previous unique and this one are not consecutive, add the new one */
  if (nre->udb_nelems == 0 ||
      nre->uniques[nre->udb_nelems - 1].orig_startpos +
      nre->uniques[nre->udb_nelems - 1].len != orig_startpos) {
    gt_assert(nre->udb_nelems == 0 ||
              nre->uniques[nre->udb_nelems - 1].orig_startpos < orig_startpos);
    gt_n_r_encseq_udb_resize(nre);
    nre->uniques[nre->udb_nelems].orig_startpos = orig_startpos;
    nre->uniques[nre->udb_nelems].len = len;
    nre->uniques[nre->udb_nelems].links.spaceuint32_t = NULL;
    nre->udb_nelems++;
  }
  else {
    gt_log_log("found consecutive uniques");
    nre->uniques[nre->udb_nelems - 1].len += len;
  }
}

static void gt_n_r_encseq_add_link_to_db(GtNREncseq *nre,
                                         GtNREncseqLink link)
{
  gt_n_r_encseq_ldb_resize(nre);
  gt_assert(nre->links != NULL);
  gt_assert(nre->ldb_nelems == 0 ||
            nre->links[nre->ldb_nelems - 1].orig_startpos < link.orig_startpos);
  nre->links[nre->ldb_nelems].editscript    = link.editscript;
  nre->links[nre->ldb_nelems].len           = link.len;
  nre->links[nre->ldb_nelems].orig_startpos = link.orig_startpos;
  nre->links[nre->ldb_nelems].unique_id     = link.unique_id;
  nre->links[nre->ldb_nelems].unique_offset = link.unique_offset;
  nre->ldb_nelems++;
}

void gt_n_r_encseq_delete(GtNREncseq *n_r_encseq)
{
  if (n_r_encseq != NULL) {
    GtUword i;
    for (i = 0; i < n_r_encseq->ldb_nelems; i++) {
      gt_editscript_delete(n_r_encseq->links[i].editscript);
    }
    for (i = 0; i < n_r_encseq->udb_nelems; i++) {
      GT_FREEARRAY(&(n_r_encseq->uniques[i].links), uint32_t);
    }
    gt_alphabet_delete(n_r_encseq->alphabet);
    gt_encseq_delete(n_r_encseq->unique_es);
    gt_free(n_r_encseq->links);
    gt_free(n_r_encseq->orig_ids);
    gt_free(n_r_encseq->uniques);
    gt_intset_delete(n_r_encseq->sdstab);
    gt_intset_delete(n_r_encseq->ssptab);
    gt_free(n_r_encseq);
  }
}

#define gt_n_r_encseq_io_one(elem) \
  io_func(&elem, sizeof (elem), (size_t) 1, fp, err)

/*generic IO function for one n_r_encseq link entry*/
static int gt_n_r_encseq_linkentry_io(GtNREncseqLink* link,
                                      FILE* fp,
                                      GtIOFunc io_func,
                                      GtError *err)
{
  int had_err = 0;
  gt_assert (link != NULL);
  had_err = gt_n_r_encseq_io_one(link->orig_startpos);
  if (!had_err)
    had_err = gt_n_r_encseq_io_one(link->len);
  if (!had_err)
    had_err = gt_n_r_encseq_io_one(link->unique_id);
  if (!had_err)
    had_err = gt_n_r_encseq_io_one(link->unique_offset);
  if (!had_err) {
    link->editscript = gt_editscript_io(link->editscript, fp, err);
    if (link->editscript == NULL) {
      had_err = 1;
    }
  }
  return had_err;
}

/*generic IO function for one n_r_encseq unique entry*/
static int gt_n_r_encseq_uniqueentry_io(GtNREncseqUnique* unique,
                                        FILE* fp,
                                        GtIOFunc io_func,
                                        GtError *err)
{
  int had_err = 0;
  gt_assert (unique != NULL);
  had_err = gt_n_r_encseq_io_one(unique->orig_startpos);
  if (!had_err)
    had_err = gt_n_r_encseq_io_one(unique->len);
  return had_err;
}

/*generic IO function for n_r_encseq data structure*/
static int gt_n_r_encseq_io(GtNREncseq *nre,
                            FILE* fp,
                            GtIOFunc io_func,
                            GtError *err)
{
  int had_err = 0;
  GtUword idx;
  had_err = gt_n_r_encseq_io_one(nre->orig_length);
  if (!had_err)
    had_err = gt_n_r_encseq_io_one(nre->orig_num_seq);
  if (!had_err)
    had_err = gt_n_r_encseq_io_one(nre->ldb_nelems);
  if (!had_err) {
    gt_assert(nre->ldb_nelems > 0);
    if (nre->links == NULL) {
      nre->links = gt_calloc((size_t) nre->ldb_nelems, sizeof (*nre->links));
      nre->ldb_allocated = nre->ldb_nelems;
    }

    had_err = gt_n_r_encseq_io_one(nre->udb_nelems);
  }

  if (!had_err) {
    gt_assert(nre->udb_nelems > 0);

    if (nre->uniques == NULL) {
      nre->uniques = gt_malloc(sizeof (*nre->uniques) * nre->udb_nelems );
      nre->udb_allocated = nre->udb_nelems;
    }
  }

  for (idx = 0; !had_err && idx < nre->ldb_nelems; idx++) {
    had_err = gt_n_r_encseq_linkentry_io(&nre->links[idx], fp, io_func, err);
  }

  for (idx = 0; !had_err && idx < nre->udb_nelems; idx++) {
    had_err = gt_n_r_encseq_uniqueentry_io(&nre->uniques[idx], fp, io_func,
                                           err);
  }
  if (!had_err) {
    nre->ssptab = gt_intset_io(nre->ssptab, fp, err);
    if (nre->ssptab == NULL)
      had_err = 1;
  }
  if (!had_err)
    had_err = gt_n_r_encseq_io_one(nre->id_len);
  if (!had_err) {
    if (nre->id_len == GT_UNDEF_UWORD) {
      nre->sdstab = gt_intset_io(nre->sdstab, fp, err);
      if (nre->sdstab == NULL)
        had_err = 1;
    }
  }
  if (!had_err)
    had_err = gt_n_r_encseq_io_one(nre->ids_total_len);
  if (!had_err) {
    nre->orig_ids = gt_realloc(nre->orig_ids, (size_t) nre->ids_total_len);
    had_err = io_func(nre->orig_ids, sizeof (*nre->orig_ids),
                      (size_t) nre->ids_total_len, fp, err);
  }
  return had_err;
}

int gt_n_r_encseq_write(GtNREncseq* nre, FILE* fp, GtError *err)
{
  gt_assert(nre != NULL && fp != NULL);
  return gt_n_r_encseq_io(nre, fp, gt_io_error_fwrite, err);
}

/*read n_r_encseq data structure from file*/
GtNREncseq *gt_n_r_encseq_new_from_file(const char *basename_nre,
                                        GtLogger *logger, GtError *err)
{
  int had_err = 0;
  GtUword i;
  FILE* fp;
  GtEncseqLoader *esl;
  GtEncseq *unique_es;
  GtNREncseq *nre = NULL;
  /*load unique_es*/
  esl = gt_encseq_loader_new();
  unique_es = gt_encseq_loader_load(esl, basename_nre, err);
  if (!unique_es)
    had_err = -1;
  if (!had_err) {
    gt_encseq_loader_delete(esl);
    nre = gt_n_r_encseq_new_empty(gt_encseq_alphabet(unique_es));
    nre->unique_es = unique_es;
    fp = gt_fa_fopen_with_suffix(basename_nre, GT_NRENCSEQ_FILE_SUFFIX, "rb",
                                 err);
    if (fp == NULL) {
      had_err = -1;
    }
    else {
      had_err = gt_n_r_encseq_io(nre, fp, gt_io_error_fread, err);
      if (!had_err) {
        gt_assert(nre->uniques);
        gt_assert(nre->links);
        gt_fa_fclose(fp);
        /*create link array for eacht unique entry*/
        for (i = 0; i < nre->udb_nelems; i++) {
          GT_INITARRAY(&(nre->uniques[i].links),uint32_t);
        }
        /* iterate through link entrys and store ids in corresponding unique
          entry array */
        for (i = 0; i < nre->ldb_nelems; i++) {
          GtUword uid = nre->links[i].unique_id;
          gt_assert(uid < nre->udb_nelems);
          GT_STOREINARRAY(&(nre->uniques[uid].links),
                          uint32_t,
                          10,
                          (uint32_t) i);
        }
      }
    }
  }
  if (!had_err) {
    gt_assert(nre != NULL);
    if (nre->id_len != GT_UNDEF_UWORD)
      gt_logger_log(logger, "IDs const len: " GT_WU, nre->id_len);
    else
      gt_logger_log(logger, "using sdstab to access IDs");
  }
  if (had_err) {
    gt_n_r_encseq_delete(nre);
    nre = NULL;
  }
  return (nre);
}

/* TODO DW refactor and combine */
/* returns index of the link element with the biggest orig_startpos smaller than
   <position>
   if smallest is larger, return that. */
static GtUword gt_n_r_encseq_links_position_binsearch(GtNREncseq *nre,
                                                      GtUword position)
{
  GtWord idx, low, high;
  gt_assert(nre && nre->ldb_nelems > 0);
  low = (GtWord) -1;
  gt_safe_assign(high, nre->ldb_nelems);
  idx = GT_DIV2(low + high);
  while (high - low > (GtWord) 1) {
    if (position < nre->links[idx].orig_startpos) {
      high = idx;
    }
    else {
      low = idx;
    }
    idx = GT_DIV2(low + high);
  }
  if (low > (GtWord) -1 && nre->links[idx].orig_startpos <= position)
    return (GtUword) idx;
  return 0;
}

/* returns index of the unique element with the biggest orig_startpos smaller
   than <position>.
   if smallest is larger: return first. */
static GtUword gt_n_r_encseq_uniques_position_binsearch(GtNREncseq *nre,
                                                        GtUword position)
{
  GtWord idx, low, high;
  gt_assert(nre && nre->udb_nelems > 0);
  low = (GtWord) -1;
  gt_safe_assign(high, nre->udb_nelems);
  idx = GT_DIV2(low + high);

  while (high - low > (GtWord) 1) {
    if (position < nre->uniques[idx].orig_startpos) {
      high = idx;
    }
    else {
      low = idx;
    }
    idx = GT_DIV2(low + high);
  }
  if (low > (GtWord) -1 && nre->uniques[idx].orig_startpos <= position)
    return (GtUword) idx;
  return 0;
}

typedef struct GtNRECXdrop {
  GtXdropresources *left_xdrop_res,
                   *right_xdrop_res,
                   *best_left_res,
                   *best_right_res;
  GtSeqabstract    *current_seq_fwd,
                   *current_seq_bwd,
                   *unique_seq_fwd,
                   *unique_seq_bwd;
  GtWord            xdropscore;
} GtNRECXdrop;

/* circular storage for hits */
typedef struct GtNRECWindow {
  GtArrayGtUword **pos_arrs;
  GtUword *idxs;
  unsigned int next,
               count;
} GtNRECWindow;

typedef GtNREncseqLink
(*gt_n_r_e_compressor_extend_fkt)(GtNREncseqCompressor *nrec);

typedef struct GtNREncseqDiagonals {
  GtUword *diagonals;
} GtNREncseqDiagonals;

struct GtNREncseqCompressor {
  GtEncseq                       *input_es;
  GtHashmap                      *kmer_hash;
  GtKmercodeiterator             *adding_iter,
                                 *main_kmer_iter;
  GtLogger                       *logger;
  GtNREncseq                     *nre;
  GtNREncseqDiagonals            *diagonals;
  gt_n_r_e_compressor_extend_fkt  extend;
  GtNRECXdrop                     xdrop;
  GtNRECWindow                    window;
  GtUword                         current_orig_start,
                                  current_seq_len,
                                  current_seq_pos,
                                  current_seq_start,
                                  initsize,
                                  main_pos,
                                  main_seqnum,
                                  max_kmer_poss,
                                  minalignlen;
  unsigned int                    kmersize,
                                  windowsize;
  bool                            use_diagonals,
                                  extend_all_kmers;
};

static void gt_n_r_encseq_compressor_delete_hash_value(GtArrayGtUword *array)
{
  GT_FREEARRAY(array, GtUword);
  gt_free(array);
}

static void
gt_n_r_encseq_compressor_xdrop_init(GtXdropArbitraryscores *scores,
                                    GtWord xdropscore,
                                    GtNRECXdrop *xdrop)
{
  xdrop->left_xdrop_res = gt_xdrop_resources_new(scores);
  xdrop->right_xdrop_res = gt_xdrop_resources_new(scores);
  xdrop->best_left_res = gt_xdrop_resources_new(scores);
  xdrop->best_right_res = gt_xdrop_resources_new(scores);
  xdrop->current_seq_fwd = gt_seqabstract_new_empty();
  xdrop->current_seq_bwd = gt_seqabstract_new_empty();
  xdrop->unique_seq_fwd = gt_seqabstract_new_empty();
  xdrop->unique_seq_bwd = gt_seqabstract_new_empty();
  xdrop->xdropscore = xdropscore;
}

static void
gt_n_r_encseq_compressor_xdrop(GtNREncseqCompressor *nrec,
                               GtUword seed_pos,
                               GtUword match_pos,
                               GtRange current_bounds,
                               GtXdropbest *best_left_xdrop,
                               GtXdropbest *best_right_xdrop,
                               GtNREncseqLink *best_link,
                               GtUword *best_match)
{
  GtNREncseqUnique unique,
                   *unique_db = nrec->nre->uniques;
  GtXdropbest left_xdrop = {0,0,0,0,0}, right_xdrop = {0,0,0,0,0};
  GtRange unique_bounds;
  GtUword match_unique_id;
  GtNRECXdrop *xdrop = &nrec->xdrop;
  const bool forward = true;

  /* get bounds for this seeds uinque */
  match_unique_id = gt_n_r_encseq_uniques_position_binsearch(nrec->nre,
                                                             match_pos);
  unique = unique_db[match_unique_id];
  unique_bounds.start = unique.orig_startpos;
  unique_bounds.end = unique_bounds.start + unique.len;
  gt_log_log("unique: " GT_WU ":" GT_WU, unique_bounds.start,
             unique_bounds.end);
  gt_assert(unique_bounds.start <= match_pos);
  gt_assert(match_pos + nrec->kmersize <= unique_bounds.end);

  /* left xdrop */
  if (current_bounds.start < seed_pos &&
      unique_bounds.start < match_pos) {
    gt_seqabstract_reinit_encseq(xdrop->unique_seq_bwd,
                                 nrec->input_es,
                                 match_pos - unique_bounds.start,
                                 unique_bounds.start);
    gt_evalxdroparbitscoresextend(!forward,
                                  &left_xdrop,
                                  xdrop->left_xdrop_res,
                                  xdrop->unique_seq_bwd,
                                  xdrop->current_seq_bwd,
                                  xdrop->xdropscore);
  }
  /* right xdrop */
  if (seed_pos < current_bounds.end &&
      match_pos < unique_bounds.end) {
    gt_seqabstract_reinit_encseq(xdrop->unique_seq_fwd,
                                 nrec->input_es,
                                 unique_bounds.end - match_pos,
                                 match_pos);
    gt_evalxdroparbitscoresextend(forward,
                                  &right_xdrop,
                                  xdrop->right_xdrop_res,
                                  xdrop->unique_seq_fwd,
                                  xdrop->current_seq_fwd,
                                  xdrop->xdropscore);
  }

  /* ivalue corresponds to unique_seq and jvalue to current_seq */
  if (left_xdrop.jvalue + right_xdrop.jvalue >= nrec->minalignlen &&
      left_xdrop.score + right_xdrop.score >
      best_left_xdrop->score + best_right_xdrop->score) {
    GtXdropresources *swap = NULL;

    *best_left_xdrop = left_xdrop;
    *best_right_xdrop = right_xdrop;

    swap = xdrop->best_left_res;
    xdrop->best_left_res = xdrop->left_xdrop_res;
    xdrop->left_xdrop_res = swap;

    swap = xdrop->best_right_res;
    xdrop->best_right_res = xdrop->right_xdrop_res;
    xdrop->right_xdrop_res = swap;

    best_link->unique_id = match_unique_id;
    best_link->unique_offset = (match_pos - left_xdrop.ivalue) -
                               unique_bounds.start;
    best_link->len = best_left_xdrop->jvalue + best_right_xdrop->jvalue;
    best_link->orig_startpos = seed_pos;
    best_link->orig_startpos -= left_xdrop.jvalue;
    *best_match = match_pos;
  }
  gt_xdrop_resources_reset(xdrop->left_xdrop_res);
  gt_xdrop_resources_reset(xdrop->right_xdrop_res);
}

#define GT_NREC_WINDOWIDX(WIN,N) (WIN->next + N < WIN->count ? \
                                  WIN->next + N :              \
                                  WIN->next + N - WIN->count)
#define GT_NREC_LAST_WIN(WIN) (WIN->next == 0 ? \
                               WIN->count - 1 : \
                               WIN->next - 1)

static GtNREncseqLink
gt_n_r_e_compressor_extend_seeds(GtNREncseqCompressor *nrec)
{
  GtNREncseqLink best_link = {NULL, 0, 0, 0, 0};
  GtRange current_bounds;
  GtXdropbest best_left_xdrop = {0,0,0,0,0},
              best_right_xdrop = {0,0,0,0,0};
  GtArrayGtUword *positions;
  GtNRECWindow *win = &nrec->window;
  GtNRECXdrop *xdrop = &nrec->xdrop;
  GtUword best_match = GT_UNDEF_UWORD,
          idx_cur,
          l_match_pos,
          xdropstart = nrec->main_pos - nrec->windowsize + 1;
  const unsigned int max_win_idx = nrec->windowsize - 1;
  unsigned int idx_win;
  const bool forward = true;

  positions = win->pos_arrs[GT_NREC_WINDOWIDX(win, 0)];
  if (positions == NULL || nrec->window.count != nrec->windowsize)
    return best_link;

  /* get bounds for current */
  current_bounds.start = nrec->current_orig_start;
  current_bounds.end = gt_n_r_encseq_ssp_seqstartpos(nrec->nre,
                                                     nrec->main_seqnum) +
                       nrec->current_seq_len;
  gt_assert(nrec->current_seq_start + nrec->current_seq_len ==
            current_bounds.end);
  gt_assert(current_bounds.start <= xdropstart);
  gt_assert(nrec->main_pos + nrec->kmersize <= current_bounds.end);

  for (idx_win = 0; idx_win <= max_win_idx; idx_win++) {
    win->idxs[idx_win] = 0;
  }

  if (current_bounds.start < xdropstart) {
    gt_seqabstract_reinit_encseq(xdrop->current_seq_bwd,
                                 nrec->input_es,
                                 xdropstart -
                                 current_bounds.start,
                                 current_bounds.start);
  }
  if (xdropstart < current_bounds.end) {
    gt_seqabstract_reinit_encseq(xdrop->current_seq_fwd,
                                 nrec->input_es,
                                 current_bounds.end -
                                 xdropstart,
                                 xdropstart);
  }

  /* iterate over all known positions of left kmer */
  for (idx_cur = 0;
       idx_cur < positions->nextfreeGtUword;
       idx_cur++)
  {
    bool found = false;
    l_match_pos = positions->spaceGtUword[idx_cur];
    /* check if new position is already covered by current best alignment */
    if (best_match == GT_UNDEF_UWORD ||
        l_match_pos > best_match + best_right_xdrop.ivalue) {
      /* start with search for right hit at end of window */
      for (idx_win = nrec->windowsize - 1;
           !found && idx_win >= nrec->kmersize;
           idx_win--) {
        GtArrayGtUword *r_positions =
          win->pos_arrs[GT_NREC_WINDOWIDX(win, idx_win)];
        /* If NULL, there are no known positions for kmer at this window
           position */
        if (r_positions != NULL) {
          GtUword r_pos_idx;
          /* within each position array, remember last highest position, start
             there, because l_match_pos increases each iteration */
          for (r_pos_idx = win->idxs[idx_win];
               !found && r_pos_idx < r_positions->nextfreeGtUword;
               r_pos_idx++) {
            GtUword r_match_pos = r_positions->spaceGtUword[r_pos_idx];
            if (r_match_pos > l_match_pos + nrec->windowsize)
              break;
            if (r_match_pos > l_match_pos + nrec->kmersize) {
              found = true;
              gt_n_r_encseq_compressor_xdrop(nrec,
                                             xdropstart,
                                             l_match_pos,
                                             current_bounds,
                                             &best_left_xdrop,
                                             &best_right_xdrop,
                                             &best_link,
                                             &best_match);
            }
          }
          win->idxs[idx_win] = r_pos_idx;
        }
      }
    }
  }

  if (best_link.len > nrec->minalignlen) {
    GtMultieoplist *meops;
    gt_assert(best_link.orig_startpos >= current_bounds.start);
    gt_assert(best_link.orig_startpos + best_link.len <= current_bounds.end);
    if (best_right_xdrop.score > 0) {
      meops = gt_xdrop_backtrack(xdrop->best_right_res, &best_right_xdrop);
    }
    else
      meops = gt_multieoplist_new();
    if (best_left_xdrop.score > 0) {
      GtMultieoplist *meopsleft = gt_xdrop_backtrack(xdrop->best_left_res,
                                                     &best_left_xdrop);
      gt_multieoplist_combine(meops, meopsleft, !forward);
      gt_multieoplist_delete(meopsleft);
    }
    best_link.editscript =
      gt_editscript_new_with_sequences(nrec->input_es,
                                       meops,
                                       best_link.orig_startpos,
                                       GT_READMODE_FORWARD);
    gt_multieoplist_delete(meops);
  }
  else
    best_link.len = 0;
  return best_link;
}

static GtNREncseqLink
gt_n_r_e_compressor_extend_seeds_old(GtNREncseqCompressor *nrec)
{
  GtNREncseqLink best_link = {NULL, 0, 0, 0, 0};
  GtRange current_bounds;
  GtXdropbest best_left_xdrop = {0,0,0,0,0},
              best_right_xdrop = {0,0,0,0,0};
  GtArrayGtUword *positions;
  GtNRECWindow *win = &nrec->window;
  GtNRECXdrop *xdrop = &nrec->xdrop;
  GtUword best_match = GT_UNDEF_UWORD,
          idx_cur,
          match_pos,
          xdropstart = nrec->main_pos;
  const bool forward = true;

  positions =win->pos_arrs[GT_NREC_LAST_WIN(win)];
  if (positions == NULL)
    return best_link;

  /* get bounds for current */
  current_bounds.start = nrec->current_orig_start;
  current_bounds.end = gt_n_r_encseq_ssp_seqstartpos(nrec->nre,
                                                     nrec->main_seqnum) +
                       nrec->current_seq_len;
  gt_assert(nrec->current_seq_start + nrec->current_seq_len ==
            current_bounds.end);
  gt_assert(current_bounds.start <= xdropstart);
  gt_assert(nrec->main_pos + nrec->kmersize <= current_bounds.end);

  if (current_bounds.start < xdropstart) {
    gt_seqabstract_reinit_encseq(xdrop->current_seq_bwd,
                                 nrec->input_es,
                                 xdropstart -
                                 current_bounds.start,
                                 current_bounds.start);
  }
  if (xdropstart < current_bounds.end) {
    gt_seqabstract_reinit_encseq(xdrop->current_seq_fwd,
                                 nrec->input_es,
                                 current_bounds.end -
                                 xdropstart,
                                 xdropstart);
  }

  for (idx_cur = 0;
       idx_cur < positions->nextfreeGtUword;
       ++idx_cur) {
    match_pos = positions->spaceGtUword[idx_cur];
    gt_n_r_encseq_compressor_xdrop(nrec,
                                   xdropstart,
                                   match_pos,
                                   current_bounds,
                                   &best_left_xdrop,
                                   &best_right_xdrop,
                                   &best_link,
                                   &best_match);
  }

  if (best_link.len > nrec->minalignlen) {
    GtMultieoplist *meops;
    gt_assert(best_link.orig_startpos >= current_bounds.start);
    gt_assert(best_link.orig_startpos + best_link.len <= current_bounds.end);
    if (best_right_xdrop.score > 0) {
      meops = gt_xdrop_backtrack(xdrop->best_right_res, &best_right_xdrop);
    }
    else
      meops = gt_multieoplist_new();
    if (best_left_xdrop.score > 0) {
      GtMultieoplist *meopsleft = gt_xdrop_backtrack(xdrop->best_left_res,
                                                     &best_left_xdrop);
      gt_multieoplist_combine(meops, meopsleft, !forward);
      gt_multieoplist_delete(meopsleft);
    }
    best_link.editscript =
      gt_editscript_new_with_sequences(nrec->input_es,
                                       meops,
                                       best_link.orig_startpos,
                                       GT_READMODE_FORWARD);
    gt_multieoplist_delete(meops);
  }
  else
    best_link.len = 0;
  return best_link;
}

static GtNREncseqLink
gt_n_r_e_compressor_extend_diagonal_seeds(GtNREncseqCompressor *nrec)
{
  GtNREncseqLink best_link = {NULL, 0, 0, 0, 0};
  GtRange current_bounds;
  GtXdropbest best_left_xdrop = {0,0,0,0,0},
              best_right_xdrop = {0,0,0,0,0};
  GtArrayGtUword *match_positions;
  GtNRECWindow *win = &nrec->window;
  GtNRECXdrop *xdrop = &nrec->xdrop;
  GtNREncseqDiagonals *diags = nrec->diagonals;
  GtUword best_match = GT_UNDEF_UWORD,
          j_idx, i;
  const bool forward = true;

  match_positions = win->pos_arrs[GT_NREC_LAST_WIN(win)];
  if (match_positions == NULL)
    return best_link;

  i = nrec->main_pos;

  /* get bounds for current */
  current_bounds.start = nrec->current_orig_start;
  current_bounds.end = gt_n_r_encseq_ssp_seqstartpos(nrec->nre,
                                                     nrec->main_seqnum) +
                       nrec->current_seq_len;
  gt_log_log("current: " GT_WU ":" GT_WU, current_bounds.start,
             current_bounds.end);
  gt_assert(current_bounds.start <= nrec->main_pos);
  gt_assert(nrec->main_pos + nrec->kmersize <= current_bounds.end);

  gt_assert(match_positions != NULL);

  for (j_idx = 0;
       j_idx < match_positions->nextfreeGtUword;
       ++j_idx) {
    /* i = current position, j match positions, right part of seed pair */
    GtUword d,
            j = match_positions->spaceGtUword[j_idx];
    gt_assert(j <= i);
    d = i - j;

    /* check if we have already processed windowsize kmers of this sequence and
       for previous hit on diagonal. */
    if (nrec->main_pos - nrec->current_orig_start >=
          (GtUword) nrec->windowsize &&
        diags->diagonals[d] != GT_UNDEF_UWORD) {
      GtUword j_prime = diags->diagonals[d],
              distance;

      gt_assert(j_prime < j);
      distance = j - j_prime;

      if (distance > (GtUword) nrec->kmersize &&
          distance <= (GtUword) nrec->windowsize) {
        GtUword i_prime,
                midpoint_seed_i,
                midpoint_seed_j,
                midpoint_offset = GT_DIV2(distance + nrec->kmersize),
                j_unique = gt_n_r_encseq_uniques_position_binsearch(nrec->nre,
                                                                    j);

        /* as j >= j' and d = i - j = i' - j', i' = d + j' can not overflow */
        i_prime = d + j_prime;

        midpoint_seed_j = j_prime + midpoint_offset;
        midpoint_seed_i = i_prime + midpoint_offset;
        /* j and j_prime have to be from the same unique sequences AND
           midpoint_seed_j position has to be outside of the current
           best alignment. (only checks for '>' because the previous j was
           smaller, also note that i and j are reversed in xdrop) */
        if (j_prime >= nrec->nre->uniques[j_unique].orig_startpos &&
            (best_match == GT_UNDEF_UWORD ||
             midpoint_seed_j > best_match + best_right_xdrop.ivalue)) {
          gt_assert(midpoint_seed_i >= current_bounds.start);
          gt_assert(midpoint_seed_i <= current_bounds.end);
          gt_log_log("seed: " GT_WU " match: " GT_WU, midpoint_seed_i,
                     midpoint_seed_j);
          if (current_bounds.start < midpoint_seed_i) {
            gt_seqabstract_reinit_encseq(xdrop->current_seq_bwd,
                                         nrec->input_es,
                                         midpoint_seed_i -
                                         current_bounds.start,
                                         current_bounds.start);
          }
          if (midpoint_seed_i < current_bounds.end) {
            gt_seqabstract_reinit_encseq(xdrop->current_seq_fwd,
                                         nrec->input_es,
                                         current_bounds.end - midpoint_seed_i,
                                         midpoint_seed_i);
          }
          gt_n_r_encseq_compressor_xdrop(nrec,
                                         midpoint_seed_i,
                                         midpoint_seed_j,
                                         current_bounds,
                                         &best_left_xdrop,
                                         &best_right_xdrop,
                                         &best_link,
                                         &best_match);
        }
      }
      if (distance > (GtUword) nrec->kmersize)
        diags->diagonals[d] = j;
    }
    else
      diags->diagonals[d] = j;
  }

  if (best_link.len > nrec->minalignlen) {
    GtMultieoplist *meops;
    gt_log_log("orig_start: " GT_WU " len: " GT_WU, best_link.orig_startpos,
               best_link.len);
    gt_assert(best_link.orig_startpos >= current_bounds.start);
    gt_assert(best_link.orig_startpos + best_link.len <= current_bounds.end);
    if (best_right_xdrop.score > 0) {
      meops = gt_xdrop_backtrack(xdrop->best_right_res, &best_right_xdrop);
    }
    else
      meops = gt_multieoplist_new();
    if (best_left_xdrop.score > 0) {
      GtMultieoplist *meopsleft = gt_xdrop_backtrack(xdrop->best_left_res,
                                                     &best_left_xdrop);
      gt_multieoplist_combine(meops, meopsleft, !forward);
      gt_multieoplist_delete(meopsleft);
    }
    best_link.editscript =
      gt_editscript_new_with_sequences(nrec->input_es,
                                       meops,
                                       best_link.orig_startpos,
                                       GT_READMODE_FORWARD);
    gt_multieoplist_delete(meops);
  }
  else
    best_link.len = 0;
  return best_link;
}

GtNREncseqCompressor *gt_n_r_encseq_compressor_new(
                                                GtUword initsize,
                                                GtUword minalignlength,
                                                GtUword max_kmer_poss,
                                                GtWord xdropscore,
                                                GtXdropArbitraryscores *scores,
                                                unsigned int kmersize,
                                                unsigned int windowsize,
                                                GtLogger *logger)
{
  GtNREncseqCompressor *n_r_e_compressor;
  n_r_e_compressor = gt_malloc(sizeof (*n_r_e_compressor));
  n_r_e_compressor->kmer_hash =
    gt_hashmap_new(GT_HASH_DIRECT,
                   NULL,
                   (GtFree) gt_n_r_encseq_compressor_delete_hash_value);
  n_r_e_compressor->adding_iter = NULL;
  n_r_e_compressor->current_seq_pos = 0;
  n_r_e_compressor->current_orig_start = 0;
  n_r_e_compressor->initsize = initsize;
  n_r_e_compressor->kmersize = kmersize;
  n_r_e_compressor->logger = logger;
  n_r_e_compressor->main_kmer_iter = NULL;
  n_r_e_compressor->main_pos = 0;
  n_r_e_compressor->main_seqnum = 0;
  n_r_e_compressor->max_kmer_poss = max_kmer_poss;
  n_r_e_compressor->minalignlen = minalignlength;
  n_r_e_compressor->nre = NULL;
  n_r_e_compressor->windowsize = windowsize;
  n_r_e_compressor->diagonals = NULL;
  gt_n_r_encseq_compressor_xdrop_init(scores, xdropscore,
                                      &n_r_e_compressor->xdrop);
  n_r_e_compressor->window.next = 0;
  n_r_e_compressor->window.count = 0;
  n_r_e_compressor->window.pos_arrs =
    gt_calloc((size_t) windowsize,
              sizeof (*n_r_e_compressor->window.pos_arrs));
  n_r_e_compressor->window.idxs =
    gt_calloc((size_t) windowsize, sizeof (*n_r_e_compressor->window.idxs));
  gt_logger_log(logger, "Parameters: k: %u, win: %u, min algn: " GT_WU
                ", init: " GT_WU,
                kmersize, windowsize, minalignlength, initsize);
  n_r_e_compressor->extend = gt_n_r_e_compressor_extend_seeds;
  n_r_e_compressor->use_diagonals = false;
  n_r_e_compressor->extend_all_kmers = false;
  return n_r_e_compressor;
}

void
gt_n_r_encseq_compressor_disable_opt(GtNREncseqCompressor *n_r_e_compressor)
{
  gt_assert(n_r_e_compressor->use_diagonals == false);
  n_r_e_compressor->extend = gt_n_r_e_compressor_extend_seeds_old;
  n_r_e_compressor->extend_all_kmers = true;
}

void gt_n_r_encseq_compressor_enable_diagonal_filter(
                                        GtNREncseqCompressor *n_r_e_compressor)
{
  gt_assert(n_r_e_compressor->extend_all_kmers == false);
  n_r_e_compressor->extend = gt_n_r_e_compressor_extend_diagonal_seeds;
  n_r_e_compressor->use_diagonals = true;
}

void gt_n_r_encseq_compressor_delete(GtNREncseqCompressor *n_r_e_compressor)
{
  if (n_r_e_compressor != NULL) {
    gt_hashmap_delete(n_r_e_compressor->kmer_hash);
    gt_xdrop_resources_delete(n_r_e_compressor->xdrop.left_xdrop_res);
    gt_xdrop_resources_delete(n_r_e_compressor->xdrop.right_xdrop_res);
    gt_xdrop_resources_delete(n_r_e_compressor->xdrop.best_left_res);
    gt_xdrop_resources_delete(n_r_e_compressor->xdrop.best_right_res);
    gt_seqabstract_delete(n_r_e_compressor->xdrop.unique_seq_fwd);
    gt_seqabstract_delete(n_r_e_compressor->xdrop.unique_seq_bwd);
    gt_seqabstract_delete(n_r_e_compressor->xdrop.current_seq_fwd);
    gt_seqabstract_delete(n_r_e_compressor->xdrop.current_seq_bwd);
    gt_free(n_r_e_compressor->window.pos_arrs);
    gt_free(n_r_e_compressor->window.idxs);
    gt_free(n_r_e_compressor);
  }
}

static void gt_n_r_e_compressor_add_kmer(GtNREncseqCompressor *n_r_e_compressor,
                                         GtCodetype kmercode,
                                         GtUword position)
{
  GtArrayGtUword *arr =
    (GtArrayGtUword *) gt_hashmap_get(n_r_e_compressor->kmer_hash,
                                      (void *) kmercode);
  if (arr == NULL) {
    arr = gt_malloc(sizeof (*arr));
    GT_INITARRAY(arr, GtUword);
    gt_hashmap_add(n_r_e_compressor->kmer_hash,
                   (void *) kmercode,
                   (void *) arr);
  }
  /* we only store the first max_kmer_poss in the array */
  /* TODO DW check other ways, ringbuffer, resetting the array, maximum total of
     positions? */
  if (arr->nextfreeGtUword < n_r_e_compressor->max_kmer_poss) {
    GtUword extend = GT_DIV8(n_r_e_compressor->max_kmer_poss) <
      (GtUword) GT_NRENCSEQ_MIN_KMER_POS ?
      (GtUword) GT_NRENCSEQ_MIN_KMER_POS :
      GT_DIV8(n_r_e_compressor->max_kmer_poss);
    GT_STOREINARRAY(arr,
                    GtUword,
                    extend,
                    position);
  }
}

typedef enum {
  GT_NREC_CONT,
  GT_NREC_EOD,
  GT_NREC_RESET,
} GtNRECState ;

static GtNRECState gt_n_r_e_compressor_reset_pos_and_main_iter_to_pos(
                                                     GtNREncseqCompressor *nrec,
                                                     GtUword pos)
{
  unsigned int idx;
  if (pos >= nrec->nre->orig_length) {
    return GT_NREC_EOD;
  }
  nrec->current_orig_start =
    nrec->main_pos = pos;
  nrec->current_seq_pos = nrec->main_pos - nrec->current_seq_start;
  nrec->window.count = 0;
  nrec->window.next = 0;
  for (idx = 0; idx < nrec->windowsize; idx++)
    nrec->window.pos_arrs[idx] = NULL;
  gt_kmercodeiterator_reset(nrec->main_kmer_iter,
                            GT_READMODE_FORWARD,
                            nrec->main_pos);
  return GT_NREC_RESET;
}

static GtNRECState gt_n_r_e_compressor_reset_pos_and_main_iter_to_current_seq(
                                                     GtNREncseqCompressor *nrec)
{
  if (nrec->main_seqnum >= nrec->nre->orig_num_seq) {
    return GT_NREC_EOD;
  }
  nrec->current_seq_start = gt_n_r_encseq_ssp_seqstartpos(nrec->nre,
                                                          nrec->main_seqnum);
  return gt_n_r_e_compressor_reset_pos_and_main_iter_to_pos(
                                                       nrec,
                                                       nrec->current_seq_start);
}

  static GtNRECState
gt_n_r_e_compressor_skip_short_seqs(GtNREncseqCompressor *nrec)
{
  GtUword start;

  while (nrec->main_seqnum < nrec->nre->orig_num_seq &&
         (nrec->current_seq_len =
          gt_n_r_encseq_ssp_seqlength(nrec->nre, nrec->main_seqnum)) <
         nrec->minalignlen) {
    start = gt_n_r_encseq_ssp_seqstartpos(nrec->nre, nrec->main_seqnum);
    gt_n_r_encseq_add_unique_to_db(nrec->nre, start, nrec->current_seq_len);
    nrec->main_seqnum++;
  }
  return nrec->main_seqnum >= nrec->nre->orig_num_seq ?
    GT_NREC_EOD : GT_NREC_CONT;
}

static GtNRECState gt_n_r_e_compressor_handle_seqend(GtNREncseqCompressor *nrec)
{
  GtNRECState state;
  /* rest of sequence length */
  GtUword length =
    nrec->current_seq_len - nrec->current_seq_pos;
  /* add length of unique befor this pos */
  length += nrec->main_pos - nrec->current_orig_start;
  if (length != 0) {
    gt_n_r_encseq_add_unique_to_db(nrec->nre,
                                   nrec->current_orig_start,
                                   length);
  }
  nrec->main_seqnum++;
  state = gt_n_r_e_compressor_skip_short_seqs(nrec);
  if (state == GT_NREC_CONT) {
    state = gt_n_r_e_compressor_reset_pos_and_main_iter_to_current_seq(nrec);
  }
  return state;
}

/* <end> beeing the first position NOT to include in the stretch of kmers */
static void gt_n_r_e_compressor_add_kmers(GtNREncseqCompressor *nrec,
                                          GtUword start,
                                          GtUword end)
{
  const GtKmercode *addcode;
  gt_assert(start < end);
  /* not long enough for an alignment, no need to store kmers */
  if (start > end - nrec->minalignlen) {
    return;
  }
  if (nrec->adding_iter == NULL) {
    nrec->adding_iter = gt_kmercodeiterator_encseq_new(nrec->input_es,
                                                       GT_READMODE_FORWARD,
                                                       nrec->kmersize,
                                                       start);
  }
  else
    gt_kmercodeiterator_reset(nrec->adding_iter, GT_READMODE_FORWARD, start);

  if (!gt_kmercodeiterator_inputexhausted(nrec->adding_iter)) {
    /* get first one */
    addcode = gt_kmercodeiterator_encseq_next(nrec->adding_iter);
    if (!addcode->definedspecialposition) {
      gt_n_r_e_compressor_add_kmer(nrec, addcode->code, start);
    }
    while (start < end - nrec->kmersize  &&
           (addcode =
            gt_kmercodeiterator_encseq_next(nrec->adding_iter)) != NULL) {
      start++;
      if (!addcode->definedspecialposition) {
        gt_n_r_e_compressor_add_kmer(nrec, addcode->code, start);
      }
    }
  }
}

static void gt_n_r_e_compressor_add_current_unique_kmers(
                                           GT_UNUSED GtNREncseqCompressor *nrec)
{
  GtUword addpos = nrec->current_orig_start;
  gt_n_r_e_compressor_add_kmers(nrec, addpos, nrec->main_pos);
}

  static GtNRECState
gt_n_r_e_compressor_extend_seed_kmer(GtNREncseqCompressor *nrec)
{
  GtNRECState state = GT_NREC_CONT;
  GtNREncseqLink link = nrec->extend(nrec);

  if (link.len >= nrec->minalignlen) {
    GtUword remaining,
            unique_len = link.orig_startpos - nrec->current_orig_start;

    gt_log_log("found link at " GT_WU ", of length " GT_WU,
               link.orig_startpos, link.len);
    /* gt_editscript_show(link.editscript, gt_encseq_alphabet(nrec->input_es));
    */
    gt_n_r_encseq_add_link_to_db(nrec->nre, link);

    if (nrec->current_orig_start < link.orig_startpos) {
      /* TODO DW check if I add unnecessary kmers to the DB */
      gt_n_r_e_compressor_add_kmers(nrec, nrec->current_orig_start,
                                    link.orig_startpos);
      gt_n_r_encseq_add_unique_to_db(nrec->nre,
                                     nrec->current_orig_start,
                                     unique_len);
    }

    if (nrec->main_pos <= link.orig_startpos + link.len)
      state = gt_n_r_e_compressor_reset_pos_and_main_iter_to_pos(
                                                          nrec,
                                                          link.orig_startpos +
                                                          link.len);
    else {
      nrec->current_orig_start = link.orig_startpos + link.len;
      nrec->window.count = 0;
      nrec->window.next = 0;
    }

    remaining = nrec->current_seq_len - nrec->current_seq_pos;
    if (state != GT_NREC_EOD &&
        remaining < nrec->minalignlen) {
      state = gt_n_r_e_compressor_handle_seqend(nrec);
    }
  }
  return state;
}

static void  gt_n_r_encseq_compressor_advance_win(GtNREncseqCompressor *nrec,
                                                  GtArrayGtUword *positions)
{
  GtNRECWindow *win = &nrec->window;
  gt_assert(win->next != 0 ||
            (win->count == 0 || win->count == nrec->windowsize));
  win->pos_arrs[win->next] = positions;
  win->next++;
  if (win->next == nrec->windowsize)
    win->next = 0;
  if (win->count < nrec->windowsize)
    win->count++;
}

  static GtNRECState
gt_n_r_e_compressor_process_kmer(GtNREncseqCompressor *nrec,
                                 const GtKmercode *main_kmercode)
{
  GtNRECState state = GT_NREC_CONT;
  if (!main_kmercode->definedspecialposition) {
    GtArrayGtUword *positions = gt_hashmap_get(nrec->kmer_hash,
                                               (void *) main_kmercode->code);
    gt_n_r_encseq_compressor_advance_win(nrec, positions);
    state = gt_n_r_e_compressor_extend_seed_kmer(nrec);
  }
  /* check if special in kmer is end of sequence, we can add previous kmers and
     the current unique to the databese */
  else if (nrec->current_seq_pos + nrec->kmersize >= nrec->current_seq_len) {
    gt_n_r_e_compressor_add_current_unique_kmers(nrec);
    state = gt_n_r_e_compressor_handle_seqend(nrec);
    nrec->window.count = 0;
    gt_assert(state == GT_NREC_RESET || state == GT_NREC_EOD);
  }
  return state;
}

  static GtNRECState
gt_n_r_e_compressor_process_init_kmer_position(GtNREncseqCompressor *nrec,
                                               const GtKmercode *kcode)
{
  GtNRECState state = GT_NREC_CONT;
  if (!kcode->definedspecialposition) {
    gt_n_r_e_compressor_add_kmer(nrec, kcode->code, nrec->main_pos);
    if (nrec->initsize != 0)
      nrec->initsize--;
  }
  /* check if special in kmer is end of sequence, we can stop iterating on the
     current piece of sequence and add it completely to the database*/
  else if (nrec->current_seq_pos + nrec->kmersize >= nrec->current_seq_len) {
    state = gt_n_r_e_compressor_handle_seqend(nrec);
    if (state == GT_NREC_RESET) {
      gt_assert(nrec->main_pos == nrec->current_orig_start);
      gt_assert(nrec->main_pos ==
                gt_n_r_encseq_ssp_seqstartpos(nrec->nre,
                                              nrec->main_seqnum));
      gt_assert(nrec->current_seq_pos == 0);
    }
  }
  return state;
}

  static GtNRECState
gt_n_r_e_compressor_init_add_seqend_to_unique(GtNREncseqCompressor *nrec,
                                              const GtKmercode *main_kmercode)
{
  GtUword seqnum = nrec->main_seqnum;
  GtNRECState state = GT_NREC_CONT;
  while (state == GT_NREC_CONT &&
         seqnum == nrec->main_seqnum &&
         (main_kmercode =
          gt_kmercodeiterator_encseq_next(nrec->main_kmer_iter)) != NULL) {
    nrec->main_pos++;
    nrec->current_seq_pos++;
    state = gt_n_r_e_compressor_process_init_kmer_position(nrec, main_kmercode);
  }
  return state;
}

static int gt_n_r_e_compressor_init_kmerhash(GtNREncseqCompressor *nrec,
                                             GtError *err)
{
  int had_err = 0;
  const GtKmercode *main_kmercode = NULL;
  GtNRECState state;

  state = gt_n_r_e_compressor_skip_short_seqs(nrec);
  if (state == GT_NREC_CONT)
    state = gt_n_r_e_compressor_reset_pos_and_main_iter_to_current_seq(nrec);
  if (state != GT_NREC_EOD &&
      !gt_kmercodeiterator_inputexhausted(nrec->main_kmer_iter)) {
    /* handle first kmer */
    if ((main_kmercode =
         gt_kmercodeiterator_encseq_next(nrec->main_kmer_iter)) != NULL) {
      state = gt_n_r_e_compressor_process_init_kmer_position(nrec,
                                                             main_kmercode);
    }
    while (state == GT_NREC_CONT &&
           nrec->initsize != 0 &&
           (main_kmercode =
            gt_kmercodeiterator_encseq_next(nrec->main_kmer_iter)) != NULL) {
      nrec->main_pos++;
      nrec->current_seq_pos++;
      state = gt_n_r_e_compressor_process_init_kmer_position(nrec,
                                                             main_kmercode);
      /* handle first kmer after reset of position, state will either be CONT or
         EOD afterwards. */
      if (state == GT_NREC_RESET &&
          (main_kmercode =
           gt_kmercodeiterator_encseq_next(nrec->main_kmer_iter)) != NULL) {
        state = gt_n_r_e_compressor_process_init_kmer_position(nrec,
                                                               main_kmercode);
        gt_assert(state == GT_NREC_CONT || state == GT_NREC_EOD);
      }
    }
  }
  else {
    gt_error_set(err, "Sequence seems to contain no kmers or only too short "
                 "sequences, check input data.");
    had_err = -1;
  }
  if (state == GT_NREC_CONT && !had_err) {
    GtUword remaining = nrec->current_seq_len -
      (nrec->current_seq_pos + nrec->kmersize);
    if (remaining < nrec->minalignlen) {
      state = gt_n_r_e_compressor_init_add_seqend_to_unique(nrec,
                                                            main_kmercode);
    }
    else {
      gt_n_r_encseq_add_unique_to_db(nrec->nre, nrec->current_seq_start,
                                     nrec->current_seq_pos + nrec->kmersize);
      state =
        gt_n_r_e_compressor_reset_pos_and_main_iter_to_pos(nrec,
                                                           nrec->main_pos +
                                                           nrec->kmersize);
    }
  }
  if (!had_err && (state == GT_NREC_EOD || nrec->initsize != 0)) {
    gt_error_set(err, "Reached end of input befor initial k-mer hash was "
                 "filed, suggest smaller initial");
    had_err = -1;
  }
  return had_err;
}

/* scan the seq and fill tables */
static int gt_n_r_e_compressor_analyse(GtNREncseqCompressor *n_r_e_compressor,
                                       GtError *err)
{
  const GtKmercode *main_kmercode = NULL;
  GtNRECState state = GT_NREC_RESET;
  int had_err = 0;

  n_r_e_compressor->main_kmer_iter =
    gt_kmercodeiterator_encseq_new(n_r_e_compressor->input_es,
                                   GT_READMODE_FORWARD,
                                   n_r_e_compressor->kmersize,
                                   n_r_e_compressor->main_pos);
  had_err = gt_n_r_e_compressor_init_kmerhash(n_r_e_compressor, err);
  /* we are now within one sequence, and the rest of it is long enough, or we
     are at the beginning of a sequence that is long enough */
  if (!had_err &&
      !gt_kmercodeiterator_inputexhausted(n_r_e_compressor->main_kmer_iter)) {
    /* handle first kmer (because the iter will always be reset after init) */
    while ((main_kmercode =
            gt_kmercodeiterator_encseq_next(n_r_e_compressor->main_kmer_iter)
           ) != NULL &&
           state == GT_NREC_RESET) {
      state = gt_n_r_e_compressor_process_kmer(n_r_e_compressor, main_kmercode);
    }
    while (state == GT_NREC_CONT &&
           (main_kmercode =
            gt_kmercodeiterator_encseq_next(n_r_e_compressor->main_kmer_iter)
           ) != NULL) {
      n_r_e_compressor->main_pos++;
      n_r_e_compressor->current_seq_pos++;
      state = gt_n_r_e_compressor_process_kmer(n_r_e_compressor, main_kmercode);
      /* handle first kmer after reset of position, state will either be CONT or
         EOD afterwards. */
      while (state == GT_NREC_RESET &&
             (main_kmercode =
              gt_kmercodeiterator_encseq_next(n_r_e_compressor->main_kmer_iter)
             ) != NULL) {
        state = gt_n_r_e_compressor_process_kmer(n_r_e_compressor,
                                                 main_kmercode);
      }
    }
    if (state != GT_NREC_EOD) {
      had_err = -1;
      gt_error_set(err, "Processing of kmers stopped, "
                   "but end of data not reached");
    }
  }
  gt_kmercodeiterator_delete(n_r_e_compressor->main_kmer_iter);
  gt_kmercodeiterator_delete(n_r_e_compressor->adding_iter);
  n_r_e_compressor->main_kmer_iter = NULL;
  n_r_e_compressor->adding_iter = NULL;
  return had_err;
}

  static void
gt_n_r_encseq_compressor_write_unique_fasta(GtNREncseqCompressor *nrec,
                                            FILE *fp)
{
  GtUword idx;
  char *buffer = NULL;
  unsigned int buffsize = 0;
  GtNREncseqUnique current;
  for (idx = 0; idx < nrec->nre->udb_nelems; ++idx) {
    current = nrec->nre->uniques[idx];
    gt_assert(current.len != 0);
    fprintf(fp, ">unique" GT_WU ", start: " GT_WU ", len: " GT_WU "\n",
            idx, current.orig_startpos, current.len);
    if (buffsize < (unsigned int) current.len) {
      gt_safe_assign(buffsize, current.len + 1);
      buffer = gt_realloc(buffer, buffsize * sizeof (*buffer));
    }
    gt_encseq_extract_decoded(nrec->input_es, buffer, current.orig_startpos,
                              current.orig_startpos + current.len - 1);
    fprintf(fp, "%.*s\n", (int) current.len, buffer);
  }
  gt_free(buffer);
}

  static GtNREncseqDiagonals *
gt_n_r_encseq_compressor_diagonals_new(GtUword length)
{
  GtUword i;
  GtNREncseqDiagonals *diagonals = gt_malloc(sizeof (*diagonals));
  /* we can't store more than GT_UWORD_MAX diagonals */
  diagonals->diagonals = gt_malloc((size_t) length *
                                   sizeof (*diagonals->diagonals));
  for (i = 0; i < length; ++i) {
    diagonals->diagonals[i] = GT_UWORD_MAX;
  }
  return diagonals;
}

  static void
gt_n_r_encseq_compressor_diagonals_delete(GtNREncseqDiagonals *diagonals)
{
  if (diagonals != NULL) {
    gt_free(diagonals->diagonals);
    gt_free(diagonals);
  }
}

int gt_n_r_encseq_compressor_compress(GtNREncseqCompressor *n_r_e_compressor,
                                      GtStr *basename,
                                      GtEncseq *encseq,
                                      GtLogger *logger,
                                      GtError *err)
{
  int had_err = 0;
  GtNREncseq *nre;
  FILE *fp = NULL;
  gt_assert(n_r_e_compressor != NULL);
  gt_assert(encseq != NULL);
  n_r_e_compressor->input_es = encseq;
  nre = gt_n_r_encseq_new(encseq, logger);
  n_r_e_compressor->nre = nre;
  if (n_r_e_compressor->use_diagonals)
    n_r_e_compressor->diagonals =
      gt_n_r_encseq_compressor_diagonals_new(nre->orig_length);

  had_err = gt_n_r_e_compressor_analyse(n_r_e_compressor, err);

  if (!had_err) {
    had_err = gt_alphabet_to_file(nre->alphabet, gt_str_get(basename), err);
  }
  if (!had_err) {
    fp = gt_fa_fopen_with_suffix(gt_str_get(basename),
                                 GT_NRENCSEQ_FILE_SUFFIX,
                                 "w", err);
    if (fp == NULL)
      had_err = -1;
  }
  if (!had_err) {
    had_err = gt_n_r_encseq_write(n_r_e_compressor->nre, fp, err);
    gt_fa_xfclose(fp);
    fp = NULL;
  }
  if (!had_err) {
    fp = gt_fa_fopen_with_suffix(gt_str_get(basename), ".fas", "w", err);
    if (fp == NULL)
      had_err = -1;
  }
  if (!had_err) {
    GtEncseqEncoder *esenc;
    GtStrArray *toencode;
    gt_n_r_encseq_compressor_write_unique_fasta(n_r_e_compressor, fp);
    gt_fa_xfclose(fp);
    /*encoding of unique_es*/
    esenc = gt_encseq_encoder_new();
    gt_encseq_encoder_disable_description_support(esenc);
    gt_encseq_encoder_do_not_create_des_tab(esenc);
    gt_encseq_encoder_do_not_create_sds_tab(esenc);
    gt_encseq_encoder_do_not_create_md5_tab(esenc);
    toencode = gt_str_array_new();
    gt_str_array_add(toencode, basename);
    gt_str_append_cstr(gt_str_array_get_str(toencode, 0), ".fas");
    had_err = gt_encseq_encoder_encode(esenc,
                                       toencode,
                                       gt_str_get(basename),
                                       err);
    gt_encseq_encoder_delete(esenc);
    gt_str_array_delete(toencode);
  }
  n_r_e_compressor->input_es = NULL;
  gt_n_r_encseq_delete(nre);
  n_r_e_compressor->nre = NULL;
  gt_n_r_encseq_compressor_diagonals_delete(n_r_e_compressor->diagonals);
  n_r_e_compressor->diagonals = NULL;
  return had_err;
}

/*BEGIN simple GtUword binary tree implementation*/
typedef struct simpleBinaryTreeNode {
  GtUword value;
  struct simpleBinaryTreeNode *left,
                              *right;
} simpleBinaryTreeNode;

static simpleBinaryTreeNode* binary_tree_new_node(GtUword elem)
{
  simpleBinaryTreeNode *bt = gt_malloc(sizeof(*bt));
  bt->value = elem;
  bt->left = NULL;
  bt->right = NULL;
  return bt;
}

static void binary_tree_insert(simpleBinaryTreeNode *bt, GtUword elem)
{
  simpleBinaryTreeNode *current = bt;
  while (true) {
    if (elem < current->value) {
      if (current->left == NULL) {
        current->left = binary_tree_new_node(elem);
        break;
      } else {
        current = current->left;
      }
    } else {
      if (current->right == NULL) {
        current->right = binary_tree_new_node(elem);
        break;
      } else {
        current = current->right;
      }
    }
  }
}

static int binary_tree_search(simpleBinaryTreeNode *bt, GtUword elem)
{
  simpleBinaryTreeNode *current = bt;
  while (current != NULL) {
    if (elem < current->value) {
      current = current->left;
    } else if (elem > current->value) {
      current = current->right;
    } else {
      return 1;
    }
  }
  return 0;
}

static void binary_tree_delete(simpleBinaryTreeNode *bt)
{
  if (bt != NULL) {
    binary_tree_delete(bt->left);
    binary_tree_delete(bt->right);
  }
  gt_free(bt);
}
/*END simple GtUword binary tree implementation*/

struct GtNREncseqDecompressor
{
  GtNREncseq *nre;
  GtUword    *extraction_uids,
             exuid_space_mallocd,
             num_of_seqs_to_extr;
};

GtNREncseqDecompressor *gt_n_r_encseq_decompressor_new(GtNREncseq *nre)
{
  GtNREncseqDecompressor * n_r_e_decompressor;
  n_r_e_decompressor = gt_malloc(sizeof (*n_r_e_decompressor));
  n_r_e_decompressor->nre = nre;
  n_r_e_decompressor->num_of_seqs_to_extr =
    n_r_e_decompressor->exuid_space_mallocd = 0;
  n_r_e_decompressor->extraction_uids = NULL;
  return n_r_e_decompressor;
}

  void
gt_n_r_encseq_decompressor_delete(GtNREncseqDecompressor *n_r_e_decompressor)
{
  if (n_r_e_decompressor != NULL) {
    gt_free(n_r_e_decompressor->extraction_uids);
    gt_free(n_r_e_decompressor);
  }
}

/*adds uentry id of unique to decompressor->extraction_uids*/
void gt_n_r_encseq_decompressor_add_unique_idx_to_extract(
                                                  GtNREncseqDecompressor *nred,
                                                  GtUword uentry_id)
{
  if (nred->extraction_uids == NULL) {
    nred->exuid_space_mallocd = (GtUword) 10;
    nred->extraction_uids = gt_malloc(sizeof(*nred->extraction_uids) *
                                      nred->exuid_space_mallocd);
  } else if (nred->num_of_seqs_to_extr == nred->exuid_space_mallocd) {
    nred->exuid_space_mallocd += 10;
    nred->extraction_uids = gt_realloc(nred->extraction_uids,
                                       sizeof (*nred->extraction_uids) *
                                       nred->exuid_space_mallocd);
  }
  nred->extraction_uids[nred->num_of_seqs_to_extr] = uentry_id;
  nred->num_of_seqs_to_extr++;
}

typedef struct{
  GtFile             *fp;
  GtNREncseqUnique unique;
  GtNREncseqLink   link;
  GtRange          *extraction_range;
  char             *buffer;
  GtUword          unique_id,
                   num_printed_chars,
                   buffsize;
  bool             fas_header;
} GtNREncseqPrintState;

#define gt_nre_write_one(gt_file,byte) gt_file_xwrite(gt_file, byte,(size_t) 1)

/*checks if, according to current printing position derived from <printstate>, a
  separator has to be printed. Prints either a fasta header or '|' symbol,
  depending on <printstate>->fas_header.*/
static void gt_n_r_encseq_check_separator(const GtNREncseq *nre,
                                          GtNREncseqPrintState *printstate)
{
  GtUword position = printstate->extraction_range->start +
    printstate->num_printed_chars;
  if (printstate->num_printed_chars == 0 && printstate->fas_header) {
    const char* desc;
    GtUword desclen;
    GtUword seqnum = gt_n_r_encseq_ssp_pos2seqnum(nre, position);
    gt_nre_write_one(printstate->fp,">");
    desc = gt_n_r_encseq_sdstab_get_id(nre,
                                       &desclen,
                                       seqnum);
    gt_file_xwrite(printstate->fp, (void*) desc, (size_t) desclen);
    gt_nre_write_one(printstate->fp, "\n");
  } else if (gt_intset_is_member(nre->ssptab, position)) {
    if (printstate->fas_header) {
      const char* desc;
      GtUword desclen;
      GtUword seqnum = gt_n_r_encseq_ssp_pos2seqnum(nre, position + 1);
      gt_file_xwrite(printstate->fp, "\n>", (size_t) 2);
      printstate->num_printed_chars++;
      desc = gt_n_r_encseq_sdstab_get_id(nre,
                                         &desclen,
                                         seqnum);
      gt_file_xwrite(printstate->fp, (void*) desc, (size_t) desclen);
      gt_nre_write_one(printstate->fp,"\n");
    } else {
      gt_nre_write_one(printstate->fp,"|");
      printstate->num_printed_chars++;
    }
  }
}

/*sets offset and length according to <printstate>->extraction_range and
  current unique, loads sequence in printstate->buffer and prints it*/
static void gt_n_r_encseq_process_unique_print(const GtNREncseq *nre,
                                               GtNREncseqPrintState *printstate)
{
  GtUword offset = 0,
          len = printstate->unique.len;
  /*set offset if first element*/
  if (printstate->num_printed_chars == (GtUword) 0) {
    offset = printstate->extraction_range->start -
      printstate->unique.orig_startpos;
    gt_assert(offset < len);
    len -= offset;
  }
  /*set real length if last element*/
  if (printstate->extraction_range->end < (printstate->unique.orig_startpos +
                                           printstate->unique.len)) {
    gt_assert(printstate->extraction_range->end >
              printstate->unique.orig_startpos);
    len = printstate->extraction_range->end - offset -
      printstate->unique.orig_startpos + 1;
  }
  if (len > printstate->buffsize) {
    gt_safe_assign(printstate->buffsize, len);
    printstate->buffer = gt_realloc(printstate->buffer,
                          sizeof (*printstate->buffer) * printstate->buffsize);
  }
  gt_encseq_extract_decoded(nre->unique_es,
         printstate->buffer,
         gt_encseq_seqstartpos(nre->unique_es,
                               printstate->unique_id) + offset,
         gt_encseq_seqstartpos(nre->unique_es, printstate->unique_id) +
                                                              offset + len - 1);
  gt_file_xwrite(printstate->fp, printstate->buffer, (size_t)(len));
  printstate->num_printed_chars += len;
}

static void gt_n_r_encseq_process_link_print(const GtNREncseq *nre,
                                             GtNREncseqPrintState *printstate)
{
  GtUword offset = 0,
          len = printstate->link.len,
          i;
  GtAlphabet *alph = nre->alphabet;
  char *offset_buff;
  /*set offset if first element*/
  if (printstate->num_printed_chars == (GtUword) 0) {
    offset = printstate->extraction_range->start -
      printstate->link.orig_startpos;
    gt_assert(offset < len);
    len -= offset;
  }
  /*set real length if last element*/
  if (printstate->extraction_range->end < (printstate->link.orig_startpos +
                                           printstate->link.len)) {
    len = printstate->extraction_range->end - offset -
      printstate->link.orig_startpos + 1;
  }
  if (len + offset >= printstate->buffsize) {
    gt_safe_assign(printstate->buffsize, len + offset + 1);
    printstate->buffer = gt_realloc(printstate->buffer,
                          sizeof (*printstate->buffer) * printstate->buffsize);
  }
  (void) gt_editscript_get_sequence(
                           printstate->link.editscript,
                           nre->unique_es,
                           gt_encseq_seqstartpos(nre->unique_es,
                                                 printstate->link.unique_id) +
                             printstate->link.unique_offset,
                           GT_READMODE_FORWARD,
                           (GtUchar *)printstate->buffer);
  gt_alphabet_decode_seq_to_cstr(alph,
                                 printstate->buffer,
                                 (GtUchar *) printstate->buffer,
                                 len + offset);
  offset_buff = gt_malloc(sizeof(*offset_buff) * len);
  for (i=0; i < len; i++) {
    offset_buff[i] = printstate->buffer[offset + i];
  }
  gt_file_xwrite(printstate->fp, offset_buff, (size_t)(len));
  printstate->num_printed_chars += len;
  gt_free(offset_buff);
}

/*decompresses segment of original sequence given by <range>*/
int gt_n_r_encseq_decompressor_extract_originrange(GtFile* fp,
                                                   GtNREncseqDecompressor *nred,
                                                   GtRange *range,
                                                   bool fasta,
                                                   GtError *err)
{
  int had_err = 0;
  GtNREncseqPrintState *printstate;
  GtUword uidx, lidx;
  GtNREncseqUnique *uniques = nred->nre->uniques;
  GtNREncseqLink *links = nred->nre->links;
  gt_assert(nred->nre->ldb_nelems != 0);
  /*verify range*/
  if (range->start < uniques[0].orig_startpos ||
      range->end >= nred->nre->orig_length) {
    gt_error_set(err,"Range out of bounds.");
    had_err = -1;
  }
  if (!had_err) {
    gt_error_check(err);
    /*search index of unique fragment containing or preceding range->start*/
    uidx =  gt_n_r_encseq_uniques_position_binsearch(nred->nre,
                                                     range->start);
    gt_assert(uidx < nred->nre->udb_nelems);
    /*search index of link fragment containing or preceding range->start*/
    lidx =  gt_n_r_encseq_links_position_binsearch(nred->nre,
                                                   range->start);
    gt_assert(lidx < nred->nre->ldb_nelems);

    printstate = gt_malloc(sizeof(*printstate));
    printstate->fp = fp;
    printstate->fas_header = fasta;
    printstate->unique = uniques[uidx];
    printstate->unique_id = uidx;
    printstate->link = links[lidx];
    printstate->extraction_range = range;
    printstate->buffsize = printstate->unique.len;
    printstate->buffer = gt_malloc(sizeof (*printstate->buffer) *
                                   printstate->buffsize);
    printstate->num_printed_chars = 0;

    if ((uniques[uidx].orig_startpos + uniques[uidx].len) <= range->start) {
      uidx++;
      if (uidx == nred->nre->udb_nelems) {
        printstate->unique.orig_startpos = GT_UWORD_MAX;
      }
      else {
        printstate->unique_id = uidx;
        printstate->unique = uniques[uidx];
      }
    }
    if ((links[lidx].orig_startpos + links[lidx].len) <= range->start) {
      lidx++;
      if (lidx == nred->nre->ldb_nelems) {
        printstate->link.orig_startpos = GT_UWORD_MAX;
      }
      else
        printstate->link = links[lidx];
    }

    while ((printstate->num_printed_chars) < (range->end - range->start)) {
      gt_n_r_encseq_check_separator(nred->nre, printstate);
      if (printstate->unique.orig_startpos < printstate->link.orig_startpos) {
        gt_n_r_encseq_process_unique_print(nred->nre, printstate);
        uidx++;
        if (uidx == nred->nre->udb_nelems) {
          printstate->unique.orig_startpos = GT_UWORD_MAX;
        } else {
          printstate->unique_id = uidx;
          printstate->unique = uniques[printstate->unique_id];
        }
      } else {
        gt_n_r_encseq_process_link_print(nred->nre, printstate);
        lidx++;
        if (lidx == nred->nre->ldb_nelems) {
          printstate->link.orig_startpos = GT_UWORD_MAX;
        } else {
          printstate->link = links[lidx];
        }
      }
    }
    /*
    gt_assert(printstate->num_printed_chars == range->end - range->start + 1);
    */
    gt_free(printstate->buffer);
    gt_free(printstate);
  }
  return (had_err);
}

/*calls gt_n_r_encseq_decompressor_extract_originrange with range of whole
  original sequence*/
int gt_n_r_encseq_decompressor_extract_origin_complete(
                                                   GtFile* fp,
                                                   GtNREncseqDecompressor *nred,
                                                   bool fasta,
                                                   GtError *err)
{
  int had_err = 0;
  GtRange *range = gt_malloc(sizeof(*range));
  range->start = 0;
  range->end = nred->nre->orig_length - 1;
  had_err = gt_n_r_encseq_decompressor_extract_originrange(fp, nred, range,
                                                           fasta, err);
  gt_free(range);
  return had_err;
}

/*prints complete protein printstate->unique originates from*/
static int gt_n_r_encseq_print_complete_unique_origin_prot(
                                               GtNREncseqDecompressor *nred,
                                               GtNREncseqPrintState *printstate,
                                               GtError *err)
{
  int had_err = 0;
  GtUword seqnum =
    gt_n_r_encseq_ssp_pos2seqnum(nred->nre, printstate->unique.orig_startpos);
  GtUword seqlength = gt_n_r_encseq_ssp_seqlength(nred->nre, seqnum);
  GtRange *range = gt_malloc( sizeof (*range));
  range->start = gt_n_r_encseq_ssp_seqstartpos(nred->nre, seqnum);
  range->end = range->start + seqlength - 1;
  had_err = gt_n_r_encseq_decompressor_extract_originrange(printstate->fp,
                                                           nred,
                                                           range,
                                                           true,
                                                           err);
  gt_free(range);
  printstate->num_printed_chars += seqlength;
  gt_assert(printstate->num_printed_chars != 0);
  return had_err;
}

/*prints complete protein printstate->link originates from*/
static int gt_n_r_encseq_print_complete_link_origin_prot(
                                            GtNREncseqDecompressor *nred,
                                            GtNREncseqPrintState *printstate,
                                            GtError *err)
{
  int had_err = 0;
  GtUword seqnum = gt_n_r_encseq_ssp_pos2seqnum(nred->nre,
                                                printstate->link.orig_startpos);
  GtUword seqlength = gt_n_r_encseq_ssp_seqlength(nred->nre, seqnum);
  GtRange *range = gt_malloc( sizeof (*range));
  range->start = gt_n_r_encseq_ssp_seqstartpos(nred->nre, seqnum);
  range->end = range->start + seqlength - 1;
  had_err = gt_n_r_encseq_decompressor_extract_originrange(printstate->fp,
                                                           nred,
                                                           range,
                                                           true,
                                                           err);
  gt_free(range);
  printstate->num_printed_chars += seqlength;
  gt_assert(printstate->num_printed_chars != 0);
  return had_err;
}

/* decompresses all sequences pointing to given relative unique range in unique
   entry <uentry_id> of the non redundant database.
   returns !0 on error.
   */
static int gt_n_r_encseq_decompressor_extract_from_uniqueid(
                                                 GtFile *fp,
                                                 simpleBinaryTreeNode **visited,
                                                 GtNREncseqDecompressor *nred,
                                                 GtUword uentry_id,
                                                 GtUword *chars_printed,
                                                 GtError *err)
{
  int had_err = 0;
  int i;
  GtUword cur_link_id, seqnum;
  GtNREncseqUnique *uniques = nred->nre->uniques;
  GtNREncseqLink *links = nred->nre->links;
  GtNREncseqPrintState *printstate;
  printstate = gt_malloc(sizeof(*printstate));
  printstate->fp = fp;
  printstate->fas_header = true; /*TODO DW should this be hard coded?*/
  printstate->unique = uniques[uentry_id];
  printstate->unique_id = uentry_id;
  printstate->buffsize = printstate->unique.len;
  printstate->buffer = gt_malloc(sizeof (*printstate->buffer) *
                                 printstate->buffsize);
  printstate->num_printed_chars = 0;
  seqnum = gt_n_r_encseq_ssp_pos2seqnum(
                                   nred->nre,
                                   nred->nre->uniques[uentry_id].orig_startpos);
  if (!*visited) {
    *visited = binary_tree_new_node(seqnum);
  } else if (!binary_tree_search(*visited, seqnum)) {
    binary_tree_insert(*visited, seqnum);
    /*write unique entry*/
    had_err = gt_n_r_encseq_print_complete_unique_origin_prot(nred, printstate,
                                                              err);
    if (!had_err)
      gt_nre_write_one(fp,"\n");
  }
  /*write each link entry*/
  for (i = 0;
       !had_err && i < (int) uniques[uentry_id].links.nextfreeuint32_t;
       i++)
  {
    cur_link_id = (GtUword) uniques[uentry_id].links.spaceuint32_t[i];
    printstate->link = links[cur_link_id];
    seqnum = gt_n_r_encseq_ssp_pos2seqnum(nred->nre,
                                          printstate->link.orig_startpos);
    if (!binary_tree_search(*visited, seqnum)) {
      binary_tree_insert(*visited, seqnum);
      had_err =  gt_n_r_encseq_print_complete_link_origin_prot(nred,
                                                               printstate,
                                                               err);
      if (!had_err)
        gt_nre_write_one(fp,"\n");
    }
  }
  if (!had_err) {
    *chars_printed += printstate->num_printed_chars;
  }
  gt_free(printstate->buffer);
  gt_free(printstate);

  return had_err;
}

int gt_n_r_encseq_decompressor_start_unique_extraction(
                                                   GtFile *fp,
                                                   GtNREncseqDecompressor *nred,
                                                   GtUword *chars_printed,
                                                   GtError *err)
{
  int had_err = 0;
  GtUword uentry_id, i;
  simpleBinaryTreeNode *visited = NULL;
  gt_assert(*chars_printed == 0);
  for (i=0; !had_err && i < nred->num_of_seqs_to_extr; i++) {
    uentry_id = nred->extraction_uids[i];
    had_err = gt_n_r_encseq_decompressor_extract_from_uniqueid(fp,
                                                               &visited,
                                                               nred,
                                                               uentry_id,
                                                               chars_printed,
                                                               err);
    gt_error_check(err);
  }
  binary_tree_delete(visited);
  return had_err;
}

int gt_n_r_encseq_unit_test(GT_UNUSED GtError *err)
{
  int had_err = 0;
  gt_error_check(err);
  return had_err;
}
