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
#include "extended/n_r_encseq.h"
#include "extended/io_function_pointers.h"
#include "match/seqabstract.h"
#include "match/sfx-mappedstr.h"
#include "match/xdrop.h"

typedef struct {
  GtEditscript *editscript;
  GtUword       orig_startpos,
                unique_id,
                unique_offset,
                len;
} GtNREncseqLink;

typedef struct {
  GtUword orig_startpos,
          len;
  GtArrayuint32_t links;
} GtNREncseqUnique;

struct GtNREncseq
{
  GtEncseq         *unique_es;
  const GtEncseq   /* TODO: remove orig_es when the rest works.
                      we need:
                      desctab,
                      sdstab,
                      ssptab
                    */
                   *orig_es;
  GtNREncseqLink   *links;
  GtNREncseqUnique *uniques;
  GtAlphabet       *alphabet;
  GtUword           ldb_allocated,
                    ldb_nelems,
                    udb_allocated,
                    udb_nelems,
                    orig_length;
};

static GtNREncseq *gt_n_r_encseq_new_empty(const GtEncseq *orig_es)
{
  GtNREncseq *n_r_encseq;
  n_r_encseq = gt_malloc(sizeof (*n_r_encseq));
  n_r_encseq->alphabet = gt_alphabet_ref(gt_encseq_alphabet(orig_es));
  n_r_encseq->ldb_allocated = n_r_encseq->udb_allocated = 0;
  n_r_encseq->ldb_nelems = n_r_encseq->udb_nelems = 0;
  n_r_encseq->links = NULL;
  n_r_encseq->orig_es = orig_es;
  n_r_encseq->unique_es = NULL;
  n_r_encseq->uniques = NULL;
  return n_r_encseq;
}

static GtNREncseq *gt_n_r_encseq_new(const GtEncseq *orig_es)
{
  GtNREncseq *n_r_encseq;
  n_r_encseq = gt_n_r_encseq_new_empty(orig_es);
  n_r_encseq->uniques = NULL;
  n_r_encseq->links = NULL;
  n_r_encseq->ldb_allocated = n_r_encseq->udb_allocated = 0;
  return n_r_encseq;
}

void gt_n_r_encseq_print_info(const GtNREncseq *n_r_encseq)
{
  GtUword av_link_len = 0,
          av_uque_len = 0,
          i;
  gt_log_log("length,uniqueid");
  for (i = 0; gt_log_enabled() && i < n_r_encseq->ldb_nelems; i++) {
    av_link_len += n_r_encseq->links[i].len;
    gt_log_log("link: " GT_WU "," GT_WU, n_r_encseq->links[i].len,
        n_r_encseq->links[i].unique_id);
  }
  for (i = 0; i < n_r_encseq->udb_nelems; i++) {
    av_uque_len += n_r_encseq->uniques[i].len;
    gt_log_log("unique: " GT_WU, n_r_encseq->uniques[i].len);
  }
  printf(GT_WU " link entries\n",n_r_encseq->ldb_nelems);
  if (n_r_encseq->ldb_nelems > 0) {
    av_link_len = av_link_len / n_r_encseq->ldb_nelems;
    printf("average length of link fragments: " GT_WU "\n", av_link_len);
  }
  printf(GT_WU " unique entries\n", n_r_encseq->udb_nelems);
  if (n_r_encseq->udb_nelems > 0) {
    av_uque_len = av_uque_len / n_r_encseq->udb_nelems;
    printf("average length of unique fragments: " GT_WU "\n", av_uque_len);
  }
  if (n_r_encseq->unique_es != NULL)
    printf("total udb length:" GT_WU "\n",
           gt_encseq_total_length(n_r_encseq->unique_es));
  if (n_r_encseq->orig_es != NULL)
    printf("original db length:" GT_WU "\n",
           gt_encseq_total_length(n_r_encseq->orig_es));
}

GtUword gt_n_r_encseq_get_orig_length(GtNREncseq *n_r_encseq) {
  return n_r_encseq->orig_length;
}

GtUword gt_n_r_encseq_get_unique_length(GtNREncseq *n_r_encseq) {
  return gt_encseq_total_length(n_r_encseq->unique_es);
}

GtUword gt_n_r_encseq_resize(GtUword allocated) {
  if (allocated != 0) {
    allocated += (allocated) > (GtUword) 128 ?
      (GtUword) 128 :
      (allocated);
  }
  else
    allocated = (GtUword) 16;
  return allocated;
}

static void gt_n_r_encseq_check_udb_size(GtNREncseq *nre)
{
  if (nre->udb_nelems == nre->udb_allocated) {
    nre->udb_allocated = gt_n_r_encseq_resize(nre->udb_allocated);
    nre->uniques = gt_realloc(nre->uniques,
                              (size_t) nre->udb_allocated *
                              sizeof (*nre->uniques));
  }
}

static void gt_n_r_encseq_check_ldb_size(GtNREncseq *nre)
{
  if (nre->ldb_nelems == nre->ldb_allocated) {
    nre->ldb_allocated = gt_n_r_encseq_resize(nre->ldb_allocated);
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
    gt_n_r_encseq_check_udb_size(nre);
    nre->uniques[nre->udb_nelems].orig_startpos = orig_startpos;
    nre->uniques[nre->udb_nelems].len = len;
    nre->uniques[nre->udb_nelems].links.spaceuint32_t = NULL;
    nre->udb_nelems++;
  }
  else {
    nre->uniques[nre->udb_nelems - 1].len += len;
  }
}

static void gt_n_r_encseq_add_link_to_db(GtNREncseq *nre,
                                         GtNREncseqLink link)
{
  gt_n_r_encseq_check_ldb_size(nre);
  gt_assert(nre->links != NULL); /* why does scan build need this here, but not
                                    in udb? */
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
    gt_free(n_r_encseq->uniques);
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
                            GT_UNUSED GtError *err)
{
  int had_err = 0;
  GtUword idx;
  had_err = gt_n_r_encseq_io_one(nre->orig_length);
  if (!had_err)
    had_err = gt_n_r_encseq_io_one(nre->ldb_allocated);
  if (!had_err)
    had_err = gt_n_r_encseq_io_one(nre->ldb_nelems);
  if (!had_err)
    had_err = gt_n_r_encseq_io_one(nre->udb_allocated);
  if (!had_err) {
    gt_assert(nre->udb_allocated > 0);
    had_err = gt_n_r_encseq_io_one(nre->udb_nelems);
  }

  if (!had_err) {
    gt_assert(nre->udb_nelems > 0);

    if (nre->links == NULL) {
      nre->links = gt_calloc((size_t) nre->ldb_nelems, sizeof (*nre->links));
    }
  }
  for (idx = 0; !had_err && idx < nre->ldb_nelems; idx++) {
    had_err = gt_n_r_encseq_linkentry_io(&nre->links[idx], fp, io_func, err);
  }
  if (!had_err && nre->uniques == NULL) {
    nre->uniques = gt_malloc(sizeof (*nre->uniques) * nre->udb_nelems );
  }
  for (idx = 0; !had_err && idx < nre->udb_nelems; idx++) {
    had_err = gt_n_r_encseq_uniqueentry_io(&nre->uniques[idx], fp, io_func,
                                           err);
  }
  return had_err;
}

int gt_n_r_encseq_write_nre(GtNREncseq* nre, FILE* fp, GtError *err)
{
  gt_assert(nre != NULL && fp != NULL);
  return gt_n_r_encseq_io(nre, fp, gt_io_error_fwrite, err);
}

void gt_n_r_encseq_write_unique_fasta(GtNREncseq *nre,
                                      FILE *fp)
{
  GtUword idx;
  char *buffer = NULL;
  unsigned int buffsize = 0;
  GtNREncseqUnique current;
  for (idx = 0; idx < nre->udb_nelems; ++idx) {
    current = nre->uniques[idx];
    gt_assert(current.len != 0);
    fprintf(fp, ">unique" GT_WU ", start: " GT_WU ", len: " GT_WU "\n",
            idx, current.orig_startpos, current.len);
    if (buffsize < (unsigned int) current.len) {
      gt_safe_assign(buffsize, current.len + 1);
      buffer = gt_realloc(buffer, buffsize * sizeof (*buffer));
    }
    gt_encseq_extract_decoded(nre->orig_es, buffer, current.orig_startpos,
                              current.orig_startpos + current.len - 1);
    fprintf(fp, "%.*s\n", (int) current.len, buffer);
  }
  gt_free(buffer);
}

/*read n_r_encseq data structure from file*/
GtNREncseq *gt_n_r_encseq_new_from_file(const char *basename_nre,
                                        const GtEncseq *orig_es,
                                        GtError *err)
{
  int had_err = 0;
  GtUword i;
  FILE* fp;
  GtEncseqLoader *esl;
  GtEncseq *unique_es;
  GtNREncseq *nre = gt_n_r_encseq_new_empty(orig_es);
  /*load unique_es*/
  esl = gt_encseq_loader_new();
  unique_es = gt_encseq_loader_load(esl, basename_nre, err);
  if (!unique_es)
    had_err = -1;
  nre->unique_es = unique_es;
  gt_encseq_loader_delete(esl);
  if (!had_err) {
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
  if (had_err) {
    gt_n_r_encseq_delete(nre);
    return NULL;
  }
  return (nre);
}

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

typedef struct GtNRECXdrop
{
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
typedef struct GtNRECWindow
{
  GtArrayGtUword **pos_arrs;
  GtUword *idxs;
  unsigned int next,
               count;
} GtNRECWindow;

typedef GtNREncseqLink
(*gt_n_r_e_compressor_extend_fkt)(GtNREncseqCompressor *nrec);

struct GtNREncseqCompressor
{
  GtEncseq                      *input_es;
  GtHashmap                     *kmer_hash;
  GtKmercodeiterator            *main_kmer_iter,
                                *adding_iter;
  GtLogger                      *logger;
  GtNREncseq                    *nre;
  gt_n_r_e_compressor_extend_fkt extend;
  GtNRECXdrop                    xdrop;
  GtNRECWindow                   window;
  GtUword                        current_seq_len,
                                 current_seq_pos,
                                 current_seq_start,
                                 current_orig_start,
                                 initsize,
                                 input_len,
                                 main_pos,
                                 main_seqnum,
                                 max_kmer_poss,
                                 minalignlen,
                                 input_num_seq;
  unsigned int                   kmersize,
                                 windowsize;
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

static GtNREncseqLink
gt_n_r_e_compressor_extend_seeds(GtNREncseqCompressor *nrec);

static GtNREncseqLink
gt_n_r_e_compressor_extend_seeds_old(GtNREncseqCompressor *nrec);

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
  n_r_e_compressor->input_num_seq = 0;
  n_r_e_compressor->windowsize = windowsize;
  gt_n_r_encseq_compressor_xdrop_init(scores, xdropscore,
                                      &n_r_e_compressor->xdrop);
  n_r_e_compressor->window.next = 0;
  n_r_e_compressor->window.count = 0;
  n_r_e_compressor->window.pos_arrs =
    gt_calloc((size_t) windowsize, sizeof (*n_r_e_compressor->window.pos_arrs));
  n_r_e_compressor->window.idxs =
    gt_calloc((size_t) windowsize, sizeof (*n_r_e_compressor->window.idxs));
  gt_logger_log(logger, "Parameters: k: %u, win: %u, min algn: " GT_WU
                ", init: " GT_WU,
                kmersize, windowsize, minalignlength, initsize);
  n_r_e_compressor->extend = gt_n_r_e_compressor_extend_seeds;
  return n_r_e_compressor;
}

void
gt_n_r_encseq_compressor_disable_opt(GtNREncseqCompressor *n_r_e_compressor)
{
  n_r_e_compressor->extend = gt_n_r_e_compressor_extend_seeds_old;
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
  /* TODO: check other ways, ringbuffer, resetting the array, maximum total of
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
  if (pos >= nrec->input_len) {
    return GT_NREC_EOD;
  }
  nrec->current_orig_start =
    nrec->main_pos = pos;
  nrec->current_seq_pos = nrec->main_pos - nrec->current_seq_start;
  gt_kmercodeiterator_reset(nrec->main_kmer_iter,
                            GT_READMODE_FORWARD,
                            nrec->main_pos);
  return GT_NREC_RESET;
}

static GtNRECState gt_n_r_e_compressor_reset_pos_and_main_iter_to_current_seq(
                                                     GtNREncseqCompressor *nrec)
{
  if (nrec->main_seqnum >= nrec->input_num_seq) {
    return GT_NREC_EOD;
  }
  nrec->current_seq_start = gt_encseq_seqstartpos(nrec->input_es,
                                                  nrec->main_seqnum);
  return gt_n_r_e_compressor_reset_pos_and_main_iter_to_pos(
                                                       nrec,
                                                       nrec->current_seq_start);
}

static GtNRECState
gt_n_r_e_compressor_skip_short_seqs(GtNREncseqCompressor *nrec)
{
  GtUword start;

  while (nrec->main_seqnum < nrec->input_num_seq &&
         (nrec->current_seq_len = gt_encseq_seqlength(nrec->input_es,
                                                      nrec->main_seqnum)) <
         nrec->minalignlen) {
    start = gt_encseq_seqstartpos(nrec->input_es, nrec->main_seqnum);
    gt_n_r_encseq_add_unique_to_db(nrec->nre, start, nrec->current_seq_len);
    nrec->main_seqnum++;
  }
  return nrec->main_seqnum >= nrec->input_num_seq ? GT_NREC_EOD : GT_NREC_CONT;
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
  if (start > end - nrec->kmersize) {
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

static void
gt_n_r_encseq_compressor_xdrop(GtNREncseqCompressor *nrec,
                               GtUword seed_pos,
                               GtRange current_bounds,
                               GtXdropbest *best_left_xdrop,
                               GtXdropbest *best_right_xdrop,
                               GtNREncseqLink *best_link,
                               GtUword *best_seed)
{
  GtNREncseqUnique unique,
                   *unique_db = nrec->nre->uniques;
  GtXdropbest left_xdrop = {0,0,0,0,0}, right_xdrop = {0,0,0,0,0};
  GtRange unique_bounds;
  GtUword seed_unique_id;
  GtNRECXdrop *xdrop = &nrec->xdrop;
  const bool forward = true;

  /* get bounds for this seeds uinque */
  seed_unique_id = gt_n_r_encseq_uniques_position_binsearch(nrec->nre,
                                                            seed_pos);
  unique = unique_db[seed_unique_id];
  unique_bounds.start = unique.orig_startpos;
  unique_bounds.end = unique_bounds.start + unique.len;
  gt_assert(unique_bounds.start <= seed_pos);
  gt_assert(seed_pos + nrec->kmersize <= unique_bounds.end);

  /* left xdrop */
  if (current_bounds.start < nrec->main_pos - nrec->windowsize &&
      unique_bounds.start < seed_pos) {
    gt_seqabstract_reinit_encseq(xdrop->unique_seq_bwd,
                                 nrec->input_es,
                                 seed_pos - unique_bounds.start,
                                 unique_bounds.start);
    gt_evalxdroparbitscoresextend(!forward,
                                  &left_xdrop,
                                  xdrop->left_xdrop_res,
                                  xdrop->unique_seq_bwd,
                                  xdrop->current_seq_bwd,
                                  xdrop->xdropscore);
  }
  /* right xdrop */
  if (nrec->main_pos - nrec->windowsize < current_bounds.end &&
      seed_pos < unique_bounds.end) {
    gt_seqabstract_reinit_encseq(xdrop->unique_seq_fwd,
                                 nrec->input_es,
                                 unique_bounds.end -
                                 seed_pos,
                                 seed_pos);
    gt_evalxdroparbitscoresextend(forward,
                                  &right_xdrop,
                                  xdrop->right_xdrop_res,
                                  xdrop->unique_seq_fwd,
                                  xdrop->current_seq_fwd,
                                  xdrop->xdropscore);
  }

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
    best_link->unique_id = seed_unique_id;
    best_link->unique_offset = (seed_pos - left_xdrop.ivalue) -
      unique_bounds.start;
    best_link->len = best_left_xdrop->jvalue + best_right_xdrop->jvalue;
    *best_seed = seed_pos;
  }
  gt_xdrop_resources_reset(xdrop->left_xdrop_res);
  gt_xdrop_resources_reset(xdrop->right_xdrop_res);
}

#define GT_NREC_WINDOWIDX(WIN,N) (WIN->next + N < WIN->count ? \
                                  WIN->next + N :              \
                                  WIN->next + N - WIN->count)

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
  GtUword best_seed = GT_UNDEF_UWORD,
          idx_cur,
          l_hit_pos,
          xdropstart = nrec->main_pos - nrec->windowsize + 1;
  const unsigned int max_win_idx = nrec->windowsize - 1;
  unsigned int idx_win;
  const bool forward = true;

  /* get bounds for current */
  current_bounds.start = nrec->current_orig_start;
  current_bounds.end = gt_encseq_seqstartpos(nrec->input_es,
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

  positions =win->pos_arrs[GT_NREC_WINDOWIDX(win, 0)];

  for (idx_cur = 0;
       idx_cur < positions->nextfreeGtUword;
       idx_cur++)
  {
    bool found = false;
    l_hit_pos = positions->spaceGtUword[idx_cur];
    if (best_seed == GT_UWORD_MAX ||
        l_hit_pos > best_seed + best_right_xdrop.ivalue) {
      /* start with search for right hit at end of window */
      for (idx_win = nrec->windowsize - 1;
           !found && idx_win >= nrec->kmersize;
           idx_win--) {
        GtArrayGtUword *r_positions =
          win->pos_arrs[GT_NREC_WINDOWIDX(win, idx_win)];
        if (r_positions != NULL) {
          GtUword r_pos_idx;
          /* within each position array, remember last highest position, start
             there, because l_hit_pos increases each iteration */
          for (r_pos_idx = win->idxs[idx_win];
               !found && r_pos_idx < r_positions->nextfreeGtUword;
               r_pos_idx++) {
            GtUword r_hit_pos = r_positions->spaceGtUword[r_pos_idx];
            if (r_hit_pos > l_hit_pos + nrec->windowsize)
              break;
            if (r_hit_pos > l_hit_pos + nrec->kmersize) {
              found = true;
              gt_n_r_encseq_compressor_xdrop(nrec,
                                             l_hit_pos,
                                             current_bounds,
                                             &best_left_xdrop,
                                             &best_right_xdrop,
                                             &best_link,
                                             &best_seed);
            }
          }
          win->idxs[idx_win] = r_pos_idx;
        }
      }
    }
  }

  if (best_link.len > nrec->minalignlen) {
    GtMultieoplist *meops;
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
    best_link.orig_startpos = xdropstart;
    best_link.orig_startpos -= best_left_xdrop.jvalue;
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
  GtUword best_seed = GT_UNDEF_UWORD,
          idx_cur,
          seed_pos,
          xdropstart = nrec->main_pos - nrec->windowsize + 1;
  const bool forward = true;

  /* get bounds for current */
  current_bounds.start = nrec->current_orig_start;
  current_bounds.end = gt_encseq_seqstartpos(nrec->input_es,
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

  positions =win->pos_arrs[GT_NREC_WINDOWIDX(win, 0)];

  for (idx_cur = 0;
       idx_cur < positions->nextfreeGtUword;
       ++idx_cur) {
    seed_pos = positions->spaceGtUword[idx_cur];
    gt_n_r_encseq_compressor_xdrop(nrec,
                                   seed_pos,
                                   current_bounds,
                                   &best_left_xdrop,
                                   &best_right_xdrop,
                                   &best_link,
                                   &best_seed);
  }

  if (best_link.len > nrec->minalignlen) {
    GtMultieoplist *meops;
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
    best_link.orig_startpos = xdropstart;
    best_link.orig_startpos -= best_left_xdrop.jvalue;
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

static GtNRECState
gt_n_r_e_compressor_extend_seed_kmer(GtNREncseqCompressor *nrec)
{
  GtNRECState state = GT_NREC_CONT;
  GtNREncseqLink link = nrec->extend(nrec);

  if (link.len >= nrec->minalignlen) {
    GtUword remaining,
            unique_len = link.orig_startpos - nrec->current_orig_start;

    gt_n_r_encseq_add_link_to_db(nrec->nre, link);

    if (nrec->current_orig_start < link.orig_startpos) {
      gt_n_r_e_compressor_add_kmers(nrec, nrec->current_orig_start,
                                    link.orig_startpos);
      gt_n_r_encseq_add_unique_to_db(nrec->nre,
                                     nrec->current_orig_start,
                                     unique_len);
    }

    state = gt_n_r_e_compressor_reset_pos_and_main_iter_to_pos(
                                                        nrec,
                                                        link.orig_startpos +
                                                        link.len);

    nrec->window.count = 0;

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
  win->pos_arrs[win->next++] = positions;
  if (win->next == nrec->windowsize)
    win->next = 0;
  if (win->count < nrec->windowsize)
    win->count++;
}

static GtNRECState gt_n_r_e_compressor_process_kmer(GtNREncseqCompressor *nrec,
                                                const GtKmercode *main_kmercode)
{
  GtNRECState state = GT_NREC_CONT;
  if (!main_kmercode->definedspecialposition) {
    GtArrayGtUword *positions = gt_hashmap_get(nrec->kmer_hash,
                                               (void *) main_kmercode->code);
    gt_n_r_encseq_compressor_advance_win(nrec, positions);
    if (nrec->window.count == nrec->windowsize &&
        nrec->window.pos_arrs[GT_NREC_WINDOWIDX((&nrec->window), 0)] != NULL) {
      state = gt_n_r_e_compressor_extend_seed_kmer(nrec);
    }
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
      gt_assert(nrec->main_pos == gt_encseq_seqstartpos(nrec->input_es,
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
           gt_kmercodeiterator_encseq_next(n_r_e_compressor->main_kmer_iter))
        != NULL &&
        state == GT_NREC_RESET) {
      state = gt_n_r_e_compressor_process_kmer(n_r_e_compressor, main_kmercode);
    }
    while (state == GT_NREC_CONT &&
           (main_kmercode =
              gt_kmercodeiterator_encseq_next(n_r_e_compressor->main_kmer_iter))
           != NULL) {
      n_r_e_compressor->main_pos++;
      n_r_e_compressor->current_seq_pos++;
      state = gt_n_r_e_compressor_process_kmer(n_r_e_compressor, main_kmercode);
      /* handle first kmer after reset of position, state will either be CONT or
         EOD afterwards. */
      while (state == GT_NREC_RESET &&
          (main_kmercode =
             gt_kmercodeiterator_encseq_next(n_r_e_compressor->main_kmer_iter))
          != NULL) {
        state = gt_n_r_e_compressor_process_kmer(n_r_e_compressor,
                                                 main_kmercode);
      }
    }
    if (state != GT_NREC_EOD) {
      had_err = -1;
      gt_error_set(err, "end of data not reached");
    }
  }
  gt_kmercodeiterator_delete(n_r_e_compressor->main_kmer_iter);
  gt_kmercodeiterator_delete(n_r_e_compressor->adding_iter);
  n_r_e_compressor->main_kmer_iter = NULL;
  n_r_e_compressor->adding_iter = NULL;
  return had_err;
}

int gt_n_r_encseq_compressor_compress(GtNREncseqCompressor *n_r_e_compressor,
                                      GtStr *basename,
                                      GtEncseq *encseq,
                                      GtError *err)
{
  int had_err = 0;
  GtNREncseq *nre;
  FILE *fp = NULL;
  gt_assert(n_r_e_compressor != NULL);
  gt_assert(encseq != NULL);
  n_r_e_compressor->input_es = encseq;
  n_r_e_compressor->input_len = gt_encseq_total_length(encseq);
  n_r_e_compressor->input_num_seq = gt_encseq_num_of_sequences(encseq);
  nre = gt_n_r_encseq_new(encseq);
  nre->orig_length = gt_encseq_total_length(encseq);
  n_r_e_compressor->nre = nre;

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
    had_err = gt_n_r_encseq_write_nre(n_r_e_compressor->nre, fp, err);
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
    gt_n_r_encseq_write_unique_fasta(n_r_e_compressor->nre, fp);
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
  return had_err;
}

/*BEGIN simple GtUword binary tree implementation*/
typedef struct simpleBinaryTreeNode {
  GtUword value;
  struct simpleBinaryTreeNode *left,
                              *right;
} simpleBinaryTreeNode;

static simpleBinaryTreeNode* binary_tree_new_node(GtUword elem) {
  simpleBinaryTreeNode *bt = gt_malloc(sizeof(*bt));
  bt->value = elem;
  bt->left = NULL;
  bt->right = NULL;
  return bt;
}

static void binary_tree_insert(simpleBinaryTreeNode *bt, GtUword elem) {
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

static int binary_tree_search(simpleBinaryTreeNode *bt, GtUword elem) {
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

static void binary_tree_delete(simpleBinaryTreeNode *bt) {
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
  FILE             *fp;
  GtNREncseqUnique unique;
  GtNREncseqLink   link;
  GtRange          *extraction_range;
  char             *buffer;
  GtUword          unique_id,
                   num_printed_chars,
                   buffsize;
  bool             fas_header;
} GtNREncseqPrintState;

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
      GtUword seqnum = gt_encseq_seqnum(nre->orig_es, position);
      gt_xfwrite_one(">",printstate->fp);
      desc = gt_encseq_description(nre->orig_es,
                                   &desclen,
                                   seqnum);
      gt_xfwrite(desc, sizeof (char), (size_t) desclen, printstate->fp);
      gt_xfwrite_one("\n",printstate->fp);
   } else if (gt_encseq_position_is_separator(nre->orig_es, position,
                                      GT_READMODE_FORWARD)) {
    if (printstate->fas_header) {
      const char* desc;
      GtUword desclen;
      GtUword seqnum = gt_encseq_seqnum(nre->orig_es, position + 1);
      gt_xfwrite("\n>", sizeof (char), (size_t) 2, printstate->fp);
      printstate->num_printed_chars++;
      desc = gt_encseq_description(nre->orig_es,
                                 &desclen,
                                 seqnum);
      gt_xfwrite(desc, sizeof (char), (size_t) desclen, printstate->fp);
      gt_xfwrite_one("\n",printstate->fp);
    } else {
      gt_xfwrite_one("|",printstate->fp);
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
  gt_xfwrite(printstate->buffer, sizeof (*printstate->buffer), (size_t)(len),
             printstate->fp);
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
  (void) gt_editscript_get_sequence(printstate->link.editscript,
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
  gt_xfwrite(offset_buff, sizeof (*offset_buff), (size_t)(len),
             printstate->fp);
  printstate->num_printed_chars += len;
  gt_free(offset_buff);
}

/*decompresses segment of original sequence given by <range>*/
int gt_n_r_encseq_decompressor_extract_originrange(FILE* fp,
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
      printstate->unique_id = uidx;
      printstate->unique = uniques[uidx];
      if (uidx == nred->nre->udb_nelems) {
          printstate->unique.orig_startpos = GT_UWORD_MAX;
      }
    }
    if ((links[lidx].orig_startpos + links[lidx].len) <= range->start) {
      lidx++;
      printstate->link = links[lidx];
      if (lidx == nred->nre->ldb_nelems) {
        printstate->link.orig_startpos = GT_UWORD_MAX;
      }
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
    gt_free(printstate->buffer);
    gt_free(printstate);
  }
  return (had_err);
}

/*calls gt_n_r_encseq_decompressor_extract_originrange with range of whole
  original sequence*/
int gt_n_r_encseq_decompressor_extract_origin_complete(FILE* fp,
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
  GtUword seqnum = gt_encseq_seqnum(nred->nre->orig_es,
                            printstate->unique.orig_startpos);
  GtRange *range = gt_malloc( sizeof (*range));
  range->start = gt_encseq_seqstartpos(nred->nre->orig_es, seqnum);
  range->end = range->start + gt_encseq_seqlength(nred->nre->orig_es, seqnum)
                            - 1;
  had_err = gt_n_r_encseq_decompressor_extract_originrange(printstate->fp,
                                                    nred,
                                                    range,
                                                    true,
                                                    err);
  gt_free(range);
  printstate->num_printed_chars += gt_encseq_seqlength(nred->nre->orig_es,
                                                       seqnum);
  return had_err;
}

/*prints complete protein printstate->link originates from*/
static int gt_n_r_encseq_print_complete_link_origin_prot(
                                            GtNREncseqDecompressor *nred,
                                            GtNREncseqPrintState *printstate,
                                            GtError *err)
{
  int had_err = 0;
  GtUword seqnum = gt_encseq_seqnum(nred->nre->orig_es,
                            printstate->link.orig_startpos);
  GtRange *range = gt_malloc( sizeof (*range));
  range->start = gt_encseq_seqstartpos(nred->nre->orig_es, seqnum);
  range->end = range->start + gt_encseq_seqlength(nred->nre->orig_es, seqnum)
                            - 1;
  had_err = gt_n_r_encseq_decompressor_extract_originrange(printstate->fp,
                                                    nred,
                                                    range,
                                                    true,
                                                    err);
  gt_free(range);
  printstate->num_printed_chars += gt_encseq_seqlength(nred->nre->orig_es,
                                                       seqnum);
  return had_err;
}

/*decompresses all sequences pointing to given relative unique range in unique
  entry <uentry_id> of the non redundant database.*/
static GtUword gt_n_r_encseq_decompressor_extract_from_uniqueid(
                                    FILE *fp,
                                    simpleBinaryTreeNode **visited,
                                    GtNREncseqDecompressor *nred,
                                    GtUword uentry_id,
                                    GtError *err)
{
  int i;
  GtUword chars_printed;
  GtUword cur_link_id, seqnum;
  GtNREncseqUnique *uniques = nred->nre->uniques;
  GtNREncseqLink *links = nred->nre->links;
  GtNREncseqPrintState *printstate;
  printstate = gt_malloc(sizeof(*printstate));
  printstate->fp = fp;
  printstate->fas_header = true; /*TODO should this be hard coded?*/
  printstate->unique = uniques[uentry_id];
  printstate->unique_id = uentry_id;
  printstate->buffsize = printstate->unique.len;
  printstate->buffer = gt_malloc(sizeof (*printstate->buffer) *
                                          printstate->buffsize);
  printstate->num_printed_chars = 0;
  seqnum = gt_encseq_seqnum(nred->nre->orig_es,
                            nred->nre->uniques[uentry_id].orig_startpos);
  if (!*visited) {
    *visited = binary_tree_new_node(seqnum);
  } else if (!binary_tree_search(*visited, seqnum)) {
    binary_tree_insert(*visited, seqnum);
    /*write unique entry*/
    (void) gt_n_r_encseq_print_complete_unique_origin_prot(nred, printstate,
                                                           err);
    gt_xfwrite_one("\n",fp);
    gt_error_check(err);
  }
  /*write each link entry*/
  for (i = 0; i < (int) uniques[uentry_id].links.nextfreeuint32_t; i++)
  {
    cur_link_id = (GtUword) uniques[uentry_id].links.spaceuint32_t[i];
    printstate->link = links[cur_link_id];
    seqnum = gt_encseq_seqnum(nred->nre->orig_es,
                                     printstate->link.orig_startpos);
     if (!binary_tree_search(*visited, seqnum)) {
      binary_tree_insert(*visited, seqnum);
        (void) gt_n_r_encseq_print_complete_link_origin_prot(nred,
                                                             printstate,
                                                             err);
        gt_xfwrite_one("\n",fp);
        gt_error_check(err);
    }
  }
  chars_printed = printstate->num_printed_chars;
  printstate->num_printed_chars = 0;
  gt_free(printstate->buffer);
  gt_free(printstate);

  return chars_printed;
}

/*decompresses all sequences pointing to given absolute unique ranges of the non
  redundant database.*/
GtUword gt_n_r_encseq_decompressor_start_unique_extraction(FILE *fp,
                                                  GtNREncseqDecompressor *nred,
                                                  GtError *err)
{
  GtUword uentry_id, i, chars_printed = 0;
  simpleBinaryTreeNode *visited = NULL;
  for (i=0; i < nred->num_of_seqs_to_extr; i++) {
    uentry_id = nred->extraction_uids[i];
    chars_printed += gt_n_r_encseq_decompressor_extract_from_uniqueid(fp,
                                                             &visited,
                                                             nred,
                                                             uentry_id,
                                                             err);
    gt_error_check(err);
  }
  binary_tree_delete(visited);
  return chars_printed;
}

int gt_n_r_encseq_unit_test(GT_UNUSED GtError *err)
{
  int had_err = 0;
  gt_error_check(err);
  return had_err;
}
