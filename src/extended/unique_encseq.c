/*
 Copyright (c) 2013 Ole Eigenbrod <ole.eigenbrod@gmx.de>
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

#include <string.h>

#include "core/array_api.h"
#include "core/ensure.h"
#include "core/error_api.h"
#include "core/fa.h"
#include "core/fileutils_api.h"
#include "core/log_api.h"
#include "core/logger_api.h"
#include "core/ma_api.h"
#include "core/minmax.h"
#include "core/types_api.h"
#include "core/undef_api.h"
#include "core/unused_api.h"
#include "core/xansi_api.h"
#include "extended/editscript.h"
#include "extended/io_function_pointers.h"
#include "extended/unique_encseq.h"
#include "extended/unique_encseq_rep.h"
#include "match/echoseq.h"
#include "match/kmer2string.h"

struct GtUniqueEncseq {
};

typedef struct {
  GtXdropbest      *xdropbest;
  GtXdropresources *res;
  GtUword           position1,
                    position2,
                    seed_seqidx;
} GtBestAlignment;

typedef struct {
  GtBestAlignment    *bestAlignment;
  GtKmercodeiterator *kmercodeit1;
  GtXdropbest        *xdropbest;
  GtUword             maxalignlength,
                      position1,
                      position2;
} GtKmerhitInfo;

GtUniqueEncseq *gt_unique_encseq_new(void)
{
  GtUniqueEncseq *unique_encseq;
  unique_encseq = gt_malloc(sizeof (GtUniqueEncseq));
  return unique_encseq;
}

void gt_unique_encseq_delete(GtUniqueEncseq *unique_encseq)
{
  gt_free(unique_encseq);
}

GtUniqueEncseqDB *gt_unique_encseq_new_db(GtEncseq *encseq)
{
  GtUword idx;
  GtUniqueEncseqDB *uedb = gt_malloc(sizeof (GtUniqueEncseqDB));
  gt_assert(encseq);
  uedb->fragmentdb = gt_malloc(1000 *
      sizeof (GtUniqueEncseqDBentry));
  uedb->nelems = 0;
  uedb->maxelems = 1000UL;
  uedb->cumulength = 0;
  uedb->nseq = gt_encseq_num_of_sequences(encseq);
  uedb->ssp = gt_calloc((size_t) uedb->nseq, sizeof (*uedb->ssp));
  uedb->sde = gt_calloc((size_t) uedb->nseq, sizeof (*uedb->sde));
  uedb->desc = NULL;
  for (idx = 0; idx < uedb->nseq; idx++) {
    uedb->ssp[idx] = gt_encseq_seqstartpos(encseq, idx);
  }
  return (uedb);
}

void gt_unique_encseq_delete_db(GtUniqueEncseqDB *uedb)
{
  if (uedb != NULL) {
    GtUword i;
    for (i = 0; i < uedb->nelems; i++) {
      if (uedb->fragmentdb[i].type == Unique) {
        gt_free(uedb->fragmentdb[i].entry.unique);
      }
      else if (uedb->fragmentdb[i].type == Link) {
        /* gt_assert(uedb->fragmentdb[i].entry.link != NULL); */
        /* gt_assert(uedb->fragmentdb[i].entry.link->editscript != NULL); */
        gt_editscript_delete(uedb->fragmentdb[i].entry.link->editscript);
        gt_free(uedb->fragmentdb[i].entry.link);
      }
    }
    gt_free(uedb->sde);
    gt_free(uedb->desc);
    gt_free(uedb->fragmentdb);
    gt_free(uedb->ssp);
    gt_free(uedb);
  }
}

/* finds the unique encseq database entry index given a query position */
static GtUword
gt_unique_encseq_uniquedb_binsearch(GtUniqueEncseqDB *uedb,
                                    GtUword querypos,
                                    GtUword idxmin,
                                    GtUword idxmax)
{
  GtUword idx;
  gt_assert(uedb && uedb->nelems > 0 && querypos <=
            uedb->fragmentdb[uedb->nelems - 1].orig_endpos);
  idx = idxmin + (idxmax - idxmin + 1) / 2;
  if (querypos < uedb->fragmentdb[idx].orig_startpos) {
    return (gt_unique_encseq_uniquedb_binsearch(uedb,
                                                querypos,
                                                idxmin,
                                                idx - 1));
  }
  else if (querypos > uedb->fragmentdb[idx].orig_endpos) {
    return (gt_unique_encseq_uniquedb_binsearch(uedb,
                                                querypos,
                                                idx + 1,
                                                idxmax));
  }
  else {
    return (idx);
  }
}

/* finds the sequence start position (original encseq) given a query position */
static bool gt_unique_encseq_ssp_binsearch(GtUword *ssp,
                                           GtUword querypos,
                                           GtUword idxmin,
                                           GtUword idxmax)
{
  GtUword idx;
  if (idxmin > idxmax) {
    return (false);
  }
  idx = idxmin + (idxmax - idxmin + 1) / 2;
  if (querypos < ssp[idx]) {
    return (gt_unique_encseq_ssp_binsearch(ssp, querypos, idxmin, idx - 1));
  }
  else if (querypos > ssp[idx]) {
    return (gt_unique_encseq_ssp_binsearch(ssp, querypos, idx + 1, idxmax));
  }
  else if (querypos == ssp[idx]) {
    return (true);
  }
  else {
    return (false);
  }
}

/* prints coarse stats about the unique encseq database*/
void gt_unique_encseq_database_stats_coarse(GtUniqueEncseqDB *uedb,
                                            GtFile *outfp)
{
  GtUword idx, nuniques = 0, nlinks = 0, cumulen_uniq = 0,
      cumulen_link = 0, len;
  GtUniqueEncseqDBentry entry;
  for (idx = 0; idx < uedb->nelems; idx++) {
    entry = uedb->fragmentdb[idx];
    len = entry.orig_endpos - entry.orig_startpos + 1;
    if (entry.type == Link) {
      cumulen_uniq += len;
      nlinks++;
    }
    else if (entry.type == Unique) {
      cumulen_link += len;
      nuniques++;
    }
  }
  gt_file_xprintf(outfp, "ORIGINAL DATABASE STATS\n");
  gt_file_xprintf(outfp, "Number of sequences in original DB: " GT_WU "\n",
                  uedb->nseq);
  gt_file_xprintf(outfp,
                  "Total length of database sequences: " GT_WU "\n",
                  uedb->fragmentdb[uedb->nelems - 1].orig_endpos);
  gt_file_xprintf(outfp, "Size of stored descriptions: " GT_WU " bytes\n",
                  uedb->sde[uedb->nseq - 1]);

  gt_file_xprintf(outfp, "\nCOMPRESSED DATABASE STATS\n");
  gt_file_xprintf(outfp, "Number of entries in compressed db: " GT_WU "\n",
                  uedb->nelems);
  gt_file_xprintf(outfp, "Number of unique entries in compressed db: " GT_WU
                  "\n", nuniques);
  gt_file_xprintf(outfp, "Number of linked entries in compressed db: " GT_WU
                  "\n", nlinks);
}

void gt_unique_encseq_database_stats_fine(GtUniqueEncseqDB *uedb,
                                          GtFile *outfp)
{
  GtUword *len,
          *linkedidx,
          *linkoffset,
          *ndelperlink,
          *ninsperlink,
          *nlinksperunique,
          *nmatchperlink,
          *nmmperlink,
          *scriptsize,
          idx,
          seqnum = 0,
          tmplen;
  GtUniqueEncseqDBentry entry;

  len = gt_calloc((size_t) uedb->nelems * 9, sizeof (GtUword));
  linkedidx = len + uedb->nelems;
  linkoffset = linkedidx + uedb->nelems;
  ndelperlink = linkoffset + uedb->nelems;
  ninsperlink = ndelperlink + uedb->nelems;
  nlinksperunique = ninsperlink + uedb->nelems;
  nmatchperlink = nlinksperunique + uedb->nelems;
  nmmperlink = nmatchperlink + uedb->nelems;
  scriptsize = nmmperlink + uedb->nelems;

  for (idx = 0; idx < uedb->nelems; idx++) {
    entry = uedb->fragmentdb[idx];
    tmplen = entry.orig_endpos - entry.orig_startpos + 1;
    if (entry.type == Link) {
      GtUword match, mismatch, insertion, deletion;
      gt_editscript_get_stats(entry.entry.link->editscript,
                              &match, &mismatch, &insertion, &deletion);
      nmatchperlink[idx] += match;
      nmmperlink[idx] += mismatch;
      ninsperlink[idx] += insertion;
      ndelperlink[idx] += deletion;

      nlinksperunique[entry.entry.link->aligned_seq_idx]++;
      linkedidx[idx] = entry.entry.link->aligned_seq_idx;
      linkoffset[idx] = entry.entry.link->offset;
      scriptsize[idx] =
        (GtUword) gt_editscript_size(entry.entry.link->editscript);
    }
    len[idx] = tmplen;
  }

  gt_file_xprintf(outfp, "idx\ttype\tstart\tend\tlength\tnlinks\t"
                  "nmatch\tnmismatch\tnins\tndel\tseqidx\tlinkedidx\t"
                  "linkoffset\tscriptsize\n");
  for (idx = 0; idx < uedb->nelems; idx++) {
    entry = uedb->fragmentdb[idx];
    if (entry.type == Link) {
      gt_file_xprintf(outfp, GT_WU "\tlink\t", idx);
    }
    else if (entry.type == Unique) {
      gt_file_xprintf(outfp, GT_WU "\tunique\t", idx);
    }
    if ((seqnum < (uedb->nseq - 1))
        && entry.orig_startpos >= uedb->ssp[seqnum + 1]) {
      seqnum++;
    }
    gt_file_xprintf(outfp,
            GT_WU "\t" GT_WU "\t" GT_WU "\t" GT_WU "\t" GT_WU "\t" GT_WU "\t"
            GT_WU "\t" GT_WU "\t" GT_WU "\t" GT_WU "\t" GT_WU "\t" GT_WU "\n",
            entry.orig_startpos,
            entry.orig_endpos,
            len[idx],
            nlinksperunique[idx],
            nmatchperlink[idx],
            nmmperlink[idx],
            ninsperlink[idx],
            ndelperlink[idx],
            seqnum,
            linkedidx[idx],
            linkoffset[idx],
            scriptsize[idx]);
  }

  gt_free(len);
}

/* given a genomic range, the sequence is printed to a file pointer.
 unique sequence parts are read from the encseq, while linked parts
 are reconstructed with the help of edit scripts. */
int gt_unique_encseq_get_sequence_from_range(GtRange *range,
                                             GtEncseq *unique_encseq,
                                             GtUniqueEncseqDB *uedb,
                                             FILE *fp,
                                             GtError *err)
{
  char c;
  int had_err = 0;
  GtUword aligned_seq_idx,
          buffer_len = 0,
          comp_start,
          decode_end,
          decode_end2,
          decode_start,
          end_offset,
          endidx = 0,
          fragment_len,
          i, j,
          offset,
          start_offset,
          startidx = 0,
          unique_len;
  GtUword editscript_seqlen;
  GtAlphabet *alphabet = gt_encseq_alphabet(unique_encseq);
  GtEncseqReader *esr =
      gt_encseq_create_reader_with_readmode(unique_encseq,
                                            GT_READMODE_FORWARD,
                                            0);
  GtUchar *buffer = NULL;
  GtUniqueEncseqDBentry entry;
  GtUniqueEncseqUniqueEntry *unique_entry;
  GtUniqueEncseqLinkEntry *link_entry;

  if (range->start < uedb->fragmentdb[0].orig_startpos ||
      range->end > uedb->fragmentdb[uedb->nelems - 1].orig_endpos) {
    had_err = -1;
    gt_error_set(err, "range query is out of bounds");
  }

  if (!had_err) {
    startidx = gt_unique_encseq_uniquedb_binsearch(uedb,
                                                   range->start,
                                                   0,
                                                   uedb->nelems - 1);
    endidx = gt_unique_encseq_uniquedb_binsearch(uedb,
                                                 range->end,
                                                 0,
                                                 uedb->nelems - 1);
    start_offset = range->start - uedb->fragmentdb[startidx].orig_startpos;
    end_offset = range->end - uedb->fragmentdb[endidx].orig_startpos;
  }
  for (i = startidx; i <= endidx && !had_err; i++) {
    entry = uedb->fragmentdb[i];
    fragment_len = entry.orig_endpos - entry.orig_startpos + 1;
    if (startidx == endidx) {
      decode_start = start_offset;
      decode_end = end_offset;
    }
    else if (i == startidx) {
      decode_start = start_offset;
      decode_end = fragment_len - 1;
    }
    else if (i == endidx) {
      decode_start = 0;
      decode_end = end_offset;
    }
    else {
      decode_start = 0;
      decode_end = fragment_len - 1;
    }
    if (entry.type == Link) {
      link_entry = entry.entry.link;
      if (buffer_len < fragment_len) {
        buffer = gt_realloc(buffer, sizeof (*buffer) * (fragment_len));
        buffer_len = fragment_len;
      }
      gt_assert(buffer != NULL);
      aligned_seq_idx = link_entry->aligned_seq_idx;
      offset = link_entry->offset;
      comp_start =
          uedb->fragmentdb[aligned_seq_idx].entry.unique->compressed_start;
      editscript_seqlen = gt_editscript_get_sequence(link_entry->editscript,
                                                     unique_encseq,
                                                     comp_start + offset,
                                                     GT_READMODE_FORWARD,
                                                     buffer);
      decode_end2 = MIN(decode_end, editscript_seqlen - 1);
      for (j = decode_start; j <= decode_end2; j++) {
        gt_xfputc(gt_alphabet_decode(alphabet, buffer[j]), fp);
      }
    }
    else if (uedb->fragmentdb[i].type == Unique) {
      unique_entry = uedb->fragmentdb[i].entry.unique;
      comp_start = unique_entry->compressed_start;
      unique_len = MIN(decode_end - decode_start + 1,
          gt_encseq_total_length(unique_encseq) - comp_start);
      gt_encseq_reader_reinit_with_readmode(esr,
                                            unique_encseq,
                                            GT_READMODE_FORWARD,
                                            comp_start + decode_start);
      for (j = 0; j < unique_len; j++) {
        c = gt_encseq_reader_next_decoded_char(esr);
        if (c == (char) SEPARATOR) {
          continue;
        }
        else {
          gt_xfputc(c, fp);
        }
      }
    }
    else {
      gt_error_set(err, "entry type in unique encseq database is invalid");
      had_err = -1;
    }
    if (!had_err
        && gt_unique_encseq_ssp_binsearch(uedb->ssp,
                                          entry.orig_startpos + decode_end + 1,
                                          0,
                                          uedb->nseq - 1)) {
      gt_xfputc('|', fp);
    }
  }
  gt_xfputc('\n', fp);
  gt_encseq_reader_delete(esr);
  gt_free(buffer);
  return (had_err);
}

/* given an index (original encseq), the sequence description is printed
 to a file pointer */
static void gt_unique_encseq_print_seq_desc(GtUniqueEncseqDB *uedb,
                                            GtUword seqnum,
                                            FILE *fp)
{
  GtUword sds, idx;
  sds = seqnum == 0 ? 0 : uedb->sde[seqnum - 1] + 1;
  gt_xfputc('>', fp);
  for (idx = sds; idx <= uedb->sde[seqnum]; idx++) {
    if (uedb->desc[idx] == '\0') {
      continue;
    }
    gt_xfputc(uedb->desc[idx], fp);
  }
  if (seqnum == uedb->nseq - 1) {
    gt_xfputc('\n', fp);
  }
}

/* prints the sequence of a given index (of the original encseq)
 to a file pointer. just determines the range and calls the the
 get_sequence_from_range function*/
int gt_unique_encseq_get_sequence_from_idx(GtUword idx,
                                           GtEncseq *unique_encseq,
                                           GtUniqueEncseqDB *uedb,
                                           FILE *fp,
                                           GtError *err)
{
  GtRange range;
  int had_err = 0;
  if (idx >= uedb->nseq) {
    gt_error_set(err,"range end is out of bounds");
    had_err = -1;
  }
  if (!had_err) {
    range.start = uedb->ssp[idx];
    if (idx == uedb->nseq - 1) {
      range.end = uedb->cumulength - 1;
    }
    else {
      range.end = uedb->ssp[idx + 1] - 2;
    }
    gt_unique_encseq_print_seq_desc(uedb, idx, fp);
    had_err = gt_unique_encseq_get_sequence_from_range(&range,
                                                       unique_encseq,
                                                       uedb,
                                                       fp,
                                                       err);
  }
  return (had_err);
}

/* returns the entry index of the unique encseq for a given query position*/
static GtUword gt_unique_encseq_get_seqnum_binsearch(
    GtUniqueEncseqDB *uedb,
    GtUword querypos,
    GtUword idxmin,
    GtUword idxmax)
{
  GtUword idx;
  gt_assert(uedb && uedb->nseq > 0 && querypos <=
            uedb->fragmentdb[uedb->nelems - 1].orig_endpos);
  idx = idxmin + (idxmax - idxmin + 1) / 2;
  if (querypos < uedb->ssp[idx]) {
    return (gt_unique_encseq_get_seqnum_binsearch(uedb,
                                                  querypos,
                                                  idxmin,
                                                  idx - 1));
  }
  else if ((idx < (uedb->nseq - 1)) && (querypos >= uedb->ssp[idx + 1])) {
    return (gt_unique_encseq_get_seqnum_binsearch(uedb,
                                                  querypos,
                                                  idx + 1,
                                                  idxmax));
  }
  else {
    return (idx);
  }
}

/* stores all unique sequence fragments from the original encseq in a
 FASTA file (that should later be converted to an encseq again) */
static void gt_unique_encseq_udb2fasta(GtUniqueEncseqDB *uedb,
                                       const GtEncseq *encseq,
                                       FILE *fp)
{
  char *desc = NULL;
  GtUword i,
          len,
          unique_count = 0,
          endpoint,
          seqnum,
          max = GT_UWORD_MAX;
  short max_desc_len = 0;
  while (max > (GtUword) 9) {
    max_desc_len++;
    max /= 10;
  }
  desc = gt_calloc((size_t)(max_desc_len + 1), sizeof (*desc));
  for (i = 0; i < uedb->nelems; i++) {
    if (uedb->fragmentdb[i].type == Unique) {
      sprintf(desc, GT_WU, unique_count);
      len = uedb->fragmentdb[i].orig_endpos - uedb->fragmentdb[i].orig_startpos
            + 1;
      seqnum = gt_unique_encseq_get_seqnum_binsearch(
                                             uedb,
                                             uedb->fragmentdb[i].orig_startpos,
                                             0,
                                             uedb->nseq - 1);
      endpoint = gt_encseq_seqstartpos(encseq, seqnum)
                 + gt_encseq_seqlength(encseq, seqnum);
      if (uedb->fragmentdb[i].orig_endpos == endpoint) {
        gt_assert(len > 0);
        len--;
      }
      gt_assert(len > 0);
      gt_encseq2fastaoutput(fp,
                            desc,
                            encseq,
                            GT_READMODE_FORWARD,
                            uedb->fragmentdb[i].orig_startpos,
                            len,
                            60UL);
      unique_count++;
    }
  }
  gt_free(desc);
}

typedef void (*UniqueEncseqIOFunc)(void *ptr,
                                   size_t size,
                                   size_t nmemb,
                                   FILE *stream);

/* generic IO function for an unique encseq unique entry*/
static GtUniqueEncseqUniqueEntry *gt_unique_encseq_unique_io(
                                              GtUniqueEncseqUniqueEntry *unique,
                                              FILE *fp,
                                              GtIOFunc io_func,
                                              GtError *err)
{
  if (unique == NULL ) {
    unique = gt_calloc((size_t) 1, sizeof (*unique));
  }
  if (io_func(&unique->compressed_start,
              sizeof (unique->compressed_start),
              (size_t) 1,
              fp, err)) {
    gt_free(unique);
    unique = NULL;
  }
  return (unique);
}

/* generic IO function for an unique encseq link entry */
static GtUniqueEncseqLinkEntry *gt_unique_encseq_link_io(
                                                  GtUniqueEncseqLinkEntry *link,
                                                  FILE *fp,
                                                  GtIOFunc io_func,
                                                  GtError *err)
{
  int had_err = 0;
  if (link == NULL ) {
    link = gt_calloc((size_t) 1, sizeof (*link));
  }
  had_err = io_func(&link->aligned_seq_idx, sizeof (link->aligned_seq_idx),
                    (size_t) 1, fp, err);
  if (!had_err)
    had_err = io_func(&link->offset, sizeof (link->offset), (size_t) 1,
                      fp, err);
  if (!had_err) {
    link->editscript = gt_editscript_io(link->editscript, fp, err);
    if (link->editscript == NULL) {
      had_err = 1;
      gt_assert(gt_error_is_set(err));
    }
  }
  if (had_err) {
    gt_free(link);
    link = NULL;
  }
  return (link);
}

/* generic IO function for an unique encseq database entry */
static int gt_unique_encseq_uedbentry_io(GtUniqueEncseqDBentry *entry,
                                         FILE *fp,
                                         GtIOFunc io_func,
                                         GtError *err)
{
  int had_err = 0;
  had_err = io_func(&entry->orig_startpos, sizeof (entry->orig_startpos),
                    (size_t) 1, fp, err);
  if (!had_err)
    had_err = io_func(&entry->orig_endpos, sizeof (entry->orig_endpos),
                      (size_t) 1, fp, err);
  if (!had_err) {
    gt_assert(entry->orig_endpos >= entry->orig_startpos);
    had_err = io_func(&entry->type, sizeof (entry->type), (size_t) 1, fp, err);
  }
  if (!had_err) {
    gt_assert(entry->type == Unique || entry->type == Link);
    if (entry->type == Link) {
      entry->entry.link = gt_unique_encseq_link_io(entry->entry.link,
                                                   fp,
                                                   io_func, err);
      if (entry->entry.link == NULL) {
        had_err = 1;
        gt_assert(gt_error_is_set(err));
      }

    }
    else if (entry->type == Unique) {
      entry->entry.unique = gt_unique_encseq_unique_io(entry->entry.unique,
                                                       fp,
                                                       io_func, err);
      if (entry->entry.unique == NULL) {
        had_err = 1;
        gt_assert(gt_error_is_set(err));
      }
    }
  }
  return had_err;
}

/* generic IO function for the unique encseq data structure */
static int gt_unique_encseq_uedb_io(GtUniqueEncseqDB *uedb,
                                    FILE *fp,
                                    GtIOFunc io_func,
                                    GtError *err)
{
  int had_err = 0;
  GtUword idx;
  had_err = io_func(&uedb->nelems, sizeof (uedb->nelems), (size_t) 1, fp, err);
  if (!had_err) {
    gt_assert(uedb->nelems > 0);
    had_err = io_func(&uedb->maxelems, sizeof (uedb->maxelems), (size_t) 1, fp,
                      err);
  }
  if (!had_err) {
    gt_assert(uedb->maxelems >= uedb->nelems);
    had_err = io_func(&uedb->cumulength, sizeof (uedb->cumulength), (size_t) 1,
                      fp, err);
  }
  if (!had_err) {
    gt_assert(uedb->cumulength > 0);
    had_err = io_func(&uedb->nseq, sizeof (uedb->nseq), (size_t) 1, fp, err);
  }
  if (!had_err) {
    gt_assert(uedb->nseq > 0);
    if (uedb->ssp == NULL ) {
      uedb->ssp = gt_calloc((size_t) uedb->nseq ,sizeof (*uedb->ssp));
    }
    had_err = io_func(uedb->ssp, sizeof (*uedb->ssp), (size_t) uedb->nseq, fp,
                      err);
  }
  if (!had_err) {
    if (uedb->sde == NULL ) {
      uedb->sde = gt_calloc((size_t) uedb->nseq ,sizeof (*uedb->sde));
    }
    had_err = io_func(uedb->sde, sizeof (*uedb->sde), (size_t) uedb->nseq, fp,
                      err);
  }
  if (!had_err) {
    if (uedb->desc == NULL ) {
      uedb->desc = gt_calloc((size_t) uedb->sde[uedb->nseq - 1] + 1,
                             sizeof (*uedb->desc));
    }
    had_err = io_func(uedb->desc,
                      sizeof (*uedb->desc),
                      (size_t) uedb->sde[uedb->nseq - 1],
                      fp, err);
  }
  if (!had_err) {
    if (uedb->fragmentdb == NULL ) {
      uedb->fragmentdb = gt_calloc((size_t) uedb->nelems ,
                                   sizeof (*uedb->fragmentdb));
    }
    for (idx = 0; !had_err && idx < uedb->nelems; idx++) {
      had_err = gt_unique_encseq_uedbentry_io(&uedb->fragmentdb[idx], fp,
                                              io_func, err);
      if (had_err)
        gt_assert(gt_error_is_set(err));
    }
  }
  return had_err;
}

/* function to write a unique encseq data structure to a file*/
static int gt_unique_encseq_uedb_write(GtUniqueEncseqDB *uedb,
                                       FILE *fp, GtError *err)
{
  gt_assert(uedb != NULL && fp != NULL);
  return gt_unique_encseq_uedb_io(uedb, fp, gt_io_error_fwrite, err);
}

GtUniqueEncseqDB *gt_unique_encseq_uedb_read(const char *basename, GtError *err)
{
  int had_err = 0;
  FILE *fp;
  GtUniqueEncseqDB *uedb = NULL;

  uedb = gt_calloc((size_t) 1, sizeof (GtUniqueEncseqDB));
  fp = gt_fa_fopen_with_suffix(basename, GT_UNIQUE_ENCSEQ_FILE_SUFFIX, "rb",
                               err);
  if (fp == NULL) {
    return uedb;
  }
  had_err = gt_unique_encseq_uedb_io(uedb, fp, gt_io_error_fread, err);
  if (had_err) {
    gt_unique_encseq_delete_db(uedb);
    uedb = NULL;
  }
  gt_fa_fclose(fp);
  return (uedb);
}

/* function to copy the original sequence descriptions into the unique encseq*/
static void gt_unique_encseq_copy_desc(GtUniqueEncseqDB *uedb,
                                       const GtEncseq *encseq)
{
  GtUword idx, desclen, cumulen = 0;
  char *descstart;
  for (idx = 0; idx < uedb->nseq; idx++) {
    (void) gt_encseq_description(encseq, &desclen, idx);
    cumulen += desclen + 1;
    uedb->sde[idx] = cumulen - 1;
  }
  uedb->desc = gt_calloc((size_t) cumulen, sizeof (*uedb->desc));
  descstart = (char *) gt_encseq_description(encseq, &desclen, 0);
  memcpy(uedb->desc, descstart, (size_t) cumulen);
}

/* function to generate an encseq containing all unique sequence fragments
 of the original encseq. therefore, a fasta file is generated which is then
 encoded as an encseq. unfortunately, it is not possible to make an
 encseq-to-encseq-transformation.*/
int gt_unique_encseq_encseq2uniqueencseq(GtUniqueEncseqDB *uedb,
                                         const GtEncseq *encseq,
                                         const char *indexname,
                                         GtError *err)
{
  GtEncseqEncoder *encoder;
  int had_err = 0;
  GtStr *tmp;
  GtStrArray *seqfiles;
  FILE *fp, *fp2;
  char *filename = NULL;
  GT_UNUSED int nchars;

  fp = fopen("TMP.FASTA", "w");
  gt_unique_encseq_udb2fasta(uedb, encseq, fp);
  (void) fclose(fp);
  tmp = gt_str_new_cstr("TMP.FASTA");
  seqfiles = gt_str_array_new();
  gt_str_array_add(seqfiles, tmp);
  encoder = gt_encseq_encoder_new();
  gt_encseq_encoder_disable_description_support(encoder);
  had_err = gt_encseq_encoder_encode(encoder, seqfiles, indexname, err);
  if (had_err)
    gt_assert(gt_error_is_set(err));
  if (!had_err) {
    gt_unique_encseq_copy_desc(uedb, encseq);
    filename = gt_calloc(strlen(indexname) + 6, sizeof (char));
    nchars = snprintf(filename, strlen(indexname) + 6, "%s.uedb", indexname);
    gt_assert(nchars == (int) (strlen(indexname) + 5));
    fp2 = fopen(filename, "wb");
    had_err = gt_unique_encseq_uedb_write(uedb, fp2, err);
    if (had_err)
      gt_assert(gt_error_is_set(err));
    (void) fclose(fp2);

    if (!had_err && !gt_log_enabled()) {
      had_err = remove("TMP.FASTA");
      if (had_err != 0) {
        gt_error_set(err, "TMP file could not be deleted successfully");
      }
    }
  }
  gt_str_delete(tmp);
  gt_str_array_delete(seqfiles);
  gt_encseq_encoder_delete(encoder);
  gt_free(filename);
  return (had_err);
}

/* prints the reconstructed sequences. just for debugging */
void gt_unique_encseq_print_editscripts(GtUniqueEncseqDB *uedb,
                                        GtEncseq *unique_encseq,
                                        FILE *fp)
{
  GtUword i, align_seq_idx, fragment_len, buffer_len = 0, unique_start;
  GT_UNUSED GtUword testlen;
  GtUchar *buffer = NULL;
  GtUniqueEncseqDBentry entry;

  for (i = 0; i < uedb->nelems; i++) {
    entry = uedb->fragmentdb[i];
    if (entry.type == Link) {
      fragment_len = entry.orig_endpos - entry.orig_startpos + 1;
      if (buffer_len < fragment_len) {
        buffer = gt_realloc(buffer, sizeof (*buffer) * fragment_len);
        buffer_len = fragment_len;
      }
      align_seq_idx = entry.entry.link->aligned_seq_idx;
      unique_start =
          uedb->fragmentdb[align_seq_idx].entry.unique->compressed_start;
      testlen = gt_editscript_get_sequence(entry.entry.link->editscript,
                                           unique_encseq,
                                           unique_start +
                                             entry.entry.link->offset,
                                           GT_READMODE_FORWARD,
                                           buffer);
      gt_assert(testlen == fragment_len);
      gt_alphabet_decode_seq_to_fp(gt_encseq_alphabet(unique_encseq),
                                   fp,
                                   buffer,
                                   fragment_len);
      gt_xfputc('\n', fp);
    }
  }
  gt_free(buffer);
}

/* function to check whether all sequence positions are covered by entries
 of the unique encseq */
int gt_unique_encseq_check_db(GtUniqueEncseqDB *uedb,
                              GtLogger *logger,
                              GtError *err)
{
  GtUword i, cumulen = 0;
  int had_err = 0;
  GtUniqueEncseqDBentry entry;
  for (i = 0; i < uedb->nelems && !had_err; i++) {
    entry = uedb->fragmentdb[i];
    if (entry.orig_startpos == cumulen) {
      cumulen += entry.orig_endpos - entry.orig_startpos + 1;
    }
    else {
      had_err = -1;
      gt_error_set(err, "cumulative database length is wrong!");
    }
  }
  if (!had_err && cumulen == uedb->cumulength) {
    gt_logger_log(logger,
                  "database has covered all sequence positions!\n");
  }
  else {
    gt_error_set(err, "database does not cover all sequence positions!\n"
                 "cumulen " GT_WU " total " GT_WU,
                 uedb->cumulength, cumulen);
    had_err = -1;
  }
  return (had_err);
}

/* advances the unique encseq and reallocs space for 1000 more entries*/
static void gt_unique_encseq_extend_db(GtUniqueEncseqDB *uedb)
{
  GtUniqueEncseqDBentry *fragmentdb;
  gt_assert(uedb->nelems == uedb->maxelems);
  fragmentdb = gt_realloc(uedb->fragmentdb,
      sizeof (GtUniqueEncseqDBentry) * (uedb->maxelems + 1000));
  uedb->fragmentdb = fragmentdb;
  uedb->maxelems += 1000;
}

/* resets the kmeriterator within the encseq. also resets seqnum, seqlen and
 seqstartpos, if the last inserted segment ends at a sequence end*/
static void gt_unique_encseq_reset_kmeriterator(GtUniqueEncseqInfo *ueinfo,
                                                GtUword endpos)
{
  gt_kmercodeiterator_reset(ueinfo->kmercodeitMain,
                            GT_READMODE_FORWARD,
                            MIN(ueinfo->totallength - 1,
                                ueinfo->nextUniqueStartpos));
  if (endpos == ueinfo->seqstartpos + ueinfo->seqlen) {
    if (ueinfo->seqnum == ueinfo->nSequences - 1)
      return;
    else {
      ueinfo->seqnum++;
      ueinfo->seqlen = gt_encseq_seqlength(ueinfo->encseq, ueinfo->seqnum);
      ueinfo->seqstartpos = gt_encseq_seqstartpos(ueinfo->encseq,
                                                  ueinfo->seqnum);
    }
  }
}

/* inserts the boundaries of a unique sequence into the unique encseq */
static void gt_unique_encseq_uniquedb_insert(GtUniqueEncseqInfo *ueinfo,
                                             GtUword startpos,
                                             GtUword endpos,
                                             GtUword seqnum)
{
  GtUniqueEncseqDB *uedb = ueinfo->uniqueencseqdb;
  if (startpos > endpos) {
    return;
  }
  if (endpos == ueinfo->totallength) {
    endpos--;
  }
  if (uedb->nelems == uedb->maxelems) {
    gt_unique_encseq_extend_db(uedb);
  }

  uedb->fragmentdb[uedb->nelems].orig_startpos = startpos;
  uedb->fragmentdb[uedb->nelems].orig_endpos = endpos;
  uedb->fragmentdb[uedb->nelems].type = Unique;
  uedb->fragmentdb[uedb->nelems].entry.unique =
      gt_malloc(sizeof (GtUniqueEncseqUniqueEntry));
  uedb->fragmentdb[uedb->nelems].entry.unique->compressed_start =
      ueinfo->unique_cumulen;
  uedb->nelems++;
  uedb->cumulength += (endpos - startpos + 1);
  if (endpos + ueinfo->kmersize
      >= gt_encseq_seqstartpos(ueinfo->encseq, seqnum)
         + gt_encseq_seqlength(ueinfo->encseq, seqnum)) {
    ueinfo->unique_cumulen += (endpos - startpos + 1);
  }
  else {
    ueinfo->unique_cumulen += (endpos - startpos + 2);
  }

  ueinfo->nextUniqueStartpos = endpos + 1;
  ueinfo->nPosSinceInsert = 0;
  if (ueinfo->nextUniqueStartpos < ueinfo->totallength - ueinfo->kmersize + 1) {
    gt_unique_encseq_reset_kmeriterator(ueinfo, endpos);
  }
  else {
    gt_kmercodeiterator_encseq_setexhausted(ueinfo->kmercodeitMain, true);
  }
}

/* inserts an editscript into the unique encseq*/
static void gt_unique_encseq_alignment_insert(GtUniqueEncseqInfo *ueinfo,
                                              GtKmerhitInfo *kmerhitinfo,
                                              GtEditscript *editscript)
{
  GtUword aligned_seqidx, startpos, endpos;

  GtUniqueEncseqDB *uedb = ueinfo->uniqueencseqdb;

  startpos = kmerhitinfo->position1;
  endpos = kmerhitinfo->position1
           + kmerhitinfo->bestAlignment->xdropbest->jvalue;
  aligned_seqidx = kmerhitinfo->bestAlignment->seed_seqidx;
  if ((endpos == ueinfo->totallength)
      || !(endpos == ueinfo->seqstartpos + ueinfo->seqlen)) {
    endpos--;
  }
  gt_assert(startpos <= endpos);

  if (uedb->nelems == uedb->maxelems) {
    gt_unique_encseq_extend_db(uedb);
  }

  uedb->fragmentdb[uedb->nelems].orig_startpos = startpos;
  uedb->fragmentdb[uedb->nelems].orig_endpos = endpos;
  uedb->fragmentdb[uedb->nelems].type = Link;
  uedb->fragmentdb[uedb->nelems].entry.link =
      gt_calloc((size_t) 1, sizeof (GtUniqueEncseqLinkEntry));

  uedb->fragmentdb[uedb->nelems].entry.link->aligned_seq_idx = aligned_seqidx;
  uedb->fragmentdb[uedb->nelems].entry.link->offset =
      kmerhitinfo->bestAlignment->position2 -
        ueinfo->uniqueencseqdb->fragmentdb[aligned_seqidx].orig_startpos;
  uedb->fragmentdb[uedb->nelems].entry.link->editscript = editscript;
  uedb->cumulength += (endpos - startpos + 1);
  uedb->nelems++;

  ueinfo->nextUniqueStartpos = endpos + 1;
  ueinfo->nPosSinceInsert = 0;

  if (ueinfo->nextUniqueStartpos < ueinfo->totallength - ueinfo->kmersize + 1) {
    gt_unique_encseq_reset_kmeriterator(ueinfo, endpos);
  }
  else if (ueinfo->nextUniqueStartpos == ueinfo->totallength) {
    gt_kmercodeiterator_encseq_setexhausted(ueinfo->kmercodeitMain, true);
  }
  else {
    gt_unique_encseq_uniquedb_insert(ueinfo,
                                     ueinfo->nextUniqueStartpos,
                                     ueinfo->totallength - 1,
                                     ueinfo->seqnum);
  }
}

/* searches the kmer-hash for the position of a given kmer. */
static bool gt_unique_encseq_array_binsearch(GtArray *arr,
                                             GtUword querypos,
                                             GtUword idxmin,
                                             GtUword idxmax)
{
  GtUword idx;
  gt_assert(gt_array_size(arr) > 0);
  if (idxmax < idxmin)
    return false;
  idx = idxmin + (idxmax - idxmin + 1) / 2;
  gt_assert(gt_array_size(arr) > idx);
  if (*(GtUword *) gt_array_get(arr, idx) > querypos)
    if (idx == 0)
      return false;
    else
      return gt_unique_encseq_array_binsearch(arr, querypos, idxmin, idx - 1);
  else if (*(GtUword *) gt_array_get(arr, idx) < querypos)
    return gt_unique_encseq_array_binsearch(arr, querypos, idx + 1, idxmax);
  else
    return true;
}

/*performs the xdrop-extension using the given scoring scheme.*/
static void gt_unique_encseq_xdrop_extension(GtUniqueEncseqInfo *ueinfo,
                                             GtKmerhitInfo *kmerhitinfo)
{
  const GtEncseq *encseq = ueinfo->encseq;
  GtUword len, pos1, pos2;
  pos1 = kmerhitinfo->position1;
  pos2 = kmerhitinfo->position2;
  len = kmerhitinfo->maxalignlength;
  gt_assert(pos1 > pos2);

  gt_seqabstract_reinit_encseq(ueinfo->useq, encseq, len, pos2);
  gt_seqabstract_reinit_encseq(ueinfo->vseq, encseq, len, pos1);

  gt_evalxdroparbitscoresextend(true,
                                kmerhitinfo->xdropbest,
                                ueinfo->res,
                                ueinfo->useq,
                                ueinfo->vseq,
                                ueinfo->xdropbelowscore);

  gt_logger_log(ueinfo->logger,
                "xdrop alignment score: " GT_WD ", ival: " GT_WU ", jval "
                GT_WU " at: " GT_WU " and " GT_WU,
                (GtWord) kmerhitinfo->xdropbest->score,
                kmerhitinfo->xdropbest->ivalue,
                kmerhitinfo->xdropbest->jvalue,
                pos2, pos1);
}

/* performs the non-overlapping kmer extension and calls the xdrop extension
 if the kmer extension is successful.*/
static void gt_unique_encseq_process_kmer_hit(GtUniqueEncseqInfo *ueinfo,
                                              GtKmerhitInfo *kmerhitinfo)
{
  unsigned int kmersize = ueinfo->kmersize,
      windowsize = ueinfo->windowsize, i,
      hitcount = 1U, nohitcount = 0;
  bool hitfound;
  GtArray *arr;
  kmerhitinfo->xdropbest->score = 0;
  kmerhitinfo->xdropbest->ivalue = 0;
  kmerhitinfo->xdropbest->jvalue = 0;

  gt_kmercodeiterator_reset(kmerhitinfo->kmercodeit1,
                            GT_READMODE_FORWARD,
                            kmerhitinfo->position1 + kmersize);

  for (i = 0; i < windowsize - kmersize; i++) {
    ueinfo->kmercodeExt =
      gt_kmercodeiterator_encseq_next(kmerhitinfo->kmercodeit1);
    if (i % kmersize == 0) {
      if (ueinfo->kmercodeExt->definedspecialposition) {
        nohitcount++;
      }
      else {
        arr = gt_hashmap_get(ueinfo->hashmap,
                             (const void *) ueinfo->kmercodeExt->code);
        if (arr == NULL || gt_array_size(arr) == 0) {
          nohitcount++;
        }
        else {
          gt_assert(gt_array_size(arr) > 0);
          hitfound = gt_unique_encseq_array_binsearch(arr,
                                                      kmerhitinfo->position2
                                                      + kmersize
                                                      + i,
                                                      0,
                                                      gt_array_size(arr) - 1);
          if (!hitfound) {
            nohitcount++;
          }
          else {
            hitcount++;
            if (hitcount == ueinfo->minkmerhits) {
              gt_logger_log(ueinfo->logger,
                            "kmer extension qualifies for xdrop extension");
              ueinfo->kmerhitcount++;
              break;
            }
          }
        }
      }
      if (nohitcount > (ueinfo->maxkmerhits - ueinfo->minkmerhits)) {
        return;
      }
    }
  }
  if (hitcount >= ueinfo->minkmerhits) {
    gt_unique_encseq_xdrop_extension(ueinfo, kmerhitinfo);
  }
  return;
}

/*simply adds the kmer to the kmer-hash. the kmercode is the key, the position
 is the value. */
static void gt_unique_encseq_add_kmer_to_kmerhash(GtUniqueEncseqInfo *ueinfo,
                                                  GtCodetype kmercode,
                                                  GtUword position)
{
  GtArray *arr =
      (GtArray *) gt_hashmap_get(ueinfo->hashmap,
                                 (void *) kmercode);
  if (arr != NULL ) {
    if (gt_array_size(arr) < ueinfo->nkmers) {
      gt_array_add(arr, position);
    }
    ueinfo->kmercount++;
  }
  else {
    arr = gt_array_new(sizeof (GtUword));
    gt_array_add(arr, position);
    gt_hashmap_add(ueinfo->hashmap,
                   (void *) kmercode,
                   (void *) arr);
    ueinfo->kmercount++;
  }
}

/* this function adds the last #(minalignlength) kmers into the kmer-hash.
  this is neccessary to only store kmers in the hash which can be used for a
  valid alignment. */
static void gt_unique_encseq_insert_kmer_block(GtUniqueEncseqInfo *ueinfo)
{
  GtUword idx;
  gt_kmercodeiterator_reset(ueinfo->kmercodeitAdd,
                            GT_READMODE_FORWARD,
                            ueinfo->nextUniqueStartpos);
 for (idx = 0; idx < ueinfo->minalignlength;
      idx++) {
   ueinfo->kmercodeAdd = gt_kmercodeiterator_encseq_next(ueinfo->kmercodeitAdd);
   gt_unique_encseq_add_kmer_to_kmerhash(ueinfo,
                                         ueinfo->kmercodeAdd->code,
                                         ueinfo->nextUniqueStartpos + idx);
 }
}

/* if the processed kmer has not yet been seen, it is inserted into the
 kmer-hash (if more than #(minalignlength) kmers have been read).
 otherwise, it is checked for each position, where the kmer occurs,
 whether a successful alignment can be computed.
 if the remaining sequence length of position2 is too short for an alignment,
 the kmer is removed from the hash, as it can never be the start of a valid
 alignment. the process_kmer_hit function is called and all possible
 alignments are computed.
 the editscript from the best alignment is inserted into the links table.
 the sequence part before the alignment is inserted into the uniqueDB. */
static void gt_unique_encseq_process_kmer_seed(GtUniqueEncseqInfo *ueinfo,
                                               GtUword position,
                                               const GtKmercode *kmercode,
                                               GtKmerhitInfo *kmerhitinfo)
{
  GtUword i,
          matchlen,
          maxlen1,
          maxlen2,
          maxlen,
          minalignlen,
          offset,
          position2,
          seqend2,
          seed_seqidx;
  bool foundalignment = false;
  GtLogger *logger = ueinfo->logger;
  GtArray * arr = (GtArray *) gt_hashmap_get(ueinfo->hashmap,
                                             (void *) kmercode->code);
  GtXdropresources *res_tmp = NULL;
  GtMultieoplist *multieops = NULL;
  GtEditscript *editscript = NULL;
  GtUniqueEncseqDB *uedb = ueinfo->uniqueencseqdb;

  char *buffer = NULL;
  GtUword buflen = 0;

  minalignlen = ueinfo->minalignlength;

  gt_assert(ueinfo->currentposition == position);
  if (arr != NULL ) {

    kmerhitinfo->bestAlignment->xdropbest->score = 0;
    kmerhitinfo->bestAlignment->xdropbest->best_d = GT_UNDEF_LONG;
    kmerhitinfo->bestAlignment->xdropbest->best_k = GT_UNDEF_LONG;
    kmerhitinfo->bestAlignment->xdropbest->ivalue = GT_UNDEF_ULONG;
    kmerhitinfo->bestAlignment->xdropbest->jvalue = GT_UNDEF_ULONG;
    /* TODO compute possible seed alignments in parallel*/
    for (i = 0; i < gt_array_size(arr); i++) {
      position2 = *(GtUword *) gt_array_get(arr, i);
      if ((uedb->nelems == 0) ||
          position2 > uedb->fragmentdb[uedb->nelems - 1].orig_endpos) {
        seqend2 = position - 1;
        seed_seqidx = uedb->nelems;
      }
      else {
        seed_seqidx = gt_unique_encseq_uniquedb_binsearch(
                                            ueinfo->uniqueencseqdb,
                                            position2,
                                            0,
                                            ueinfo->uniqueencseqdb->nelems - 1);
        seqend2 = uedb->fragmentdb[seed_seqidx].orig_endpos;
      }
      maxlen2 = seqend2 - position2;
      maxlen1 = MIN(ueinfo->seqstartpos + ueinfo->seqlen - position,
          ueinfo->totallength - position);
      offset = position - position2;
      kmerhitinfo->maxalignlength = MIN(maxlen1, maxlen2);
      if (offset < minalignlen || kmerhitinfo->maxalignlength < minalignlen) {
        if (maxlen2 <
            (minalignlen - ueinfo->windowsize + ueinfo->kmersize - 1)) {
          gt_array_rem(arr, i--);
        }
        continue;
      }
       /* XXX to much redundancies why do we have all of these variables in 3
          places? */
      kmerhitinfo->position1 = position;
      kmerhitinfo->position2 = position2;
      gt_unique_encseq_process_kmer_hit(ueinfo, kmerhitinfo);
      matchlen = kmerhitinfo->xdropbest->jvalue;
      if (matchlen >= ueinfo->minalignlength
          && kmerhitinfo->xdropbest->score
             > kmerhitinfo->bestAlignment->xdropbest->score) {
        /*store the best alignment*/
        kmerhitinfo->bestAlignment->position1 = position;
        kmerhitinfo->bestAlignment->position2 = kmerhitinfo->position2;
        kmerhitinfo->bestAlignment->seed_seqidx = seed_seqidx;
        kmerhitinfo->bestAlignment->xdropbest->best_d =
          kmerhitinfo->xdropbest->best_d;
        kmerhitinfo->bestAlignment->xdropbest->best_k =
          kmerhitinfo->xdropbest->best_k;
        kmerhitinfo->bestAlignment->xdropbest->ivalue =
          kmerhitinfo->xdropbest->ivalue;
        kmerhitinfo->bestAlignment->xdropbest->jvalue =
          kmerhitinfo->xdropbest->jvalue;
        kmerhitinfo->bestAlignment->xdropbest->score =
          kmerhitinfo->xdropbest->score;
        res_tmp = kmerhitinfo->bestAlignment->res;
        kmerhitinfo->bestAlignment->res = ueinfo->res;
        ueinfo->res = res_tmp;
        foundalignment = true;
      }
    }
    if (!foundalignment) {
      if (ueinfo->nPosSinceInsert ==
          ueinfo->minalignlength) {
        gt_unique_encseq_insert_kmer_block(ueinfo);
      }
      else if (ueinfo->nPosSinceInsert >
               ueinfo->minalignlength) {
        gt_unique_encseq_add_kmer_to_kmerhash(ueinfo,
                                              ueinfo->kmercodeMain->code,
                                              ueinfo->currentposition);
      }
    }
    else {
      GtBestAlignment *ba = kmerhitinfo->bestAlignment;
      /*store best alignment in links table and maybe update the unique db*/
      gt_logger_log(logger, "match against unique db found!");

      gt_unique_encseq_uniquedb_insert(
                                  ueinfo,
                                  ueinfo->nextUniqueStartpos,
                                  position - 1,
                                  gt_encseq_seqnum(ueinfo->encseq,
                                                   ueinfo->nextUniqueStartpos));

      maxlen = MAX(ba->xdropbest->ivalue,
          ba->xdropbest->jvalue);
      if (buflen < maxlen) {
        buffer = gt_realloc(buffer, (size_t) (maxlen + 1));
      }
      gt_log_log("from: " GT_WU " to: " GT_WU,
                 ba->position2,
                 ba->position2 +
                   ba->xdropbest->ivalue - 1);
      gt_encseq_extract_decoded(ueinfo->encseq, buffer,
                                ba->position2,
                                ba->position2 +
                                  ba->xdropbest->ivalue - 1);
      gt_log_log("U: %.*s", (int) ba->xdropbest->ivalue, buffer);
      gt_log_log("from: " GT_WU " to: " GT_WU,
                 ba->position1,
                 ba->position1 +
                   ba->xdropbest->jvalue - 1);
      gt_encseq_extract_decoded(ueinfo->encseq, buffer,
                                ba->position1,
                                ba->position1 +
                                  ba->xdropbest->jvalue - 1);
      gt_log_log("V: %.*s", (int) ba->xdropbest->jvalue, buffer);
      multieops = gt_xdrop_backtrack(ba->res,
                                     ba->xdropbest);
      gt_multieoplist_show(multieops, stdout);
      editscript = gt_editscript_new_with_sequences(
                                          ueinfo->encseq,
                                          multieops,
                                          ba->position1,
                                          GT_READMODE_FORWARD);
      gt_editscript_show(editscript, ueinfo->alphabet);
      gt_multieoplist_delete(multieops);
      multieops = NULL;
      gt_unique_encseq_alignment_insert(ueinfo, kmerhitinfo, editscript);
      ueinfo->alignmentcount++;
    }
  }
  else {
    gt_unique_encseq_add_kmer_to_kmerhash(ueinfo,
                                          kmercode->code,
                                          ueinfo->currentposition);
  }
  gt_free(buffer);
}

/*adds the kmer to the kmer-hash of the uniqueDB. inserts an entry into the
 uniqueDB, if the end of a sequence or the end of the init-phase is reached.*/
static void gt_unique_encseq_init_uniquedb(GtUniqueEncseqInfo *ueinfo)
{
  gt_logger_log(ueinfo->logger, "initializing unique DB");
  gt_unique_encseq_add_kmer_to_kmerhash(ueinfo,
                                        ueinfo->kmercodeMain->code,
                                        ueinfo->currentposition);
  if (ueinfo->maxlen == (GtUword) ueinfo->kmersize) {
    gt_unique_encseq_uniquedb_insert(ueinfo,
                                     ueinfo->nextUniqueStartpos,
                                     ueinfo->seqlen + ueinfo->seqstartpos,
                                     ueinfo->seqnum);
  }
  else if (ueinfo->initKmerCounter + ueinfo->kmersize
           == ueinfo->initUniqueDBsize) {
    gt_unique_encseq_uniquedb_insert(ueinfo,
                                     ueinfo->nextUniqueStartpos,
                                     ueinfo->currentposition + ueinfo->kmersize
                                     - 1,
                                     ueinfo->seqnum);
  }
  ueinfo->initKmerCounter++;
}

/* processes kmers with special characters. the function simply inserts
 the remaining sequence to the uniqueDB, if the remaining sequence length
 is smaller than the given minimal alignment size.*/
static void gt_unique_encseq_process_special_kmer(GtUniqueEncseqInfo *ueinfo)
{
  if (ueinfo->currentposition
      >= ueinfo->seqstartpos + ueinfo->seqlen - ueinfo->kmersize) {
    gt_unique_encseq_uniquedb_insert(ueinfo,
                                     ueinfo->nextUniqueStartpos,
                                     ueinfo->seqstartpos + ueinfo->seqlen,
                                     ueinfo->seqnum);
  }
  return;
}

/*function to process a kmer depending its sequence position and
 whether it contains a special character. special kmers are processed
 in a seperate function. if the initialization of the uniqueDB is not
 finished, the kmer is added to the unique db. if the remaining sequence
 length is smaller than the given minimal alignment size, the remaining
 sequence is added to the uniqueDB without processing the kmers.
 */
static void gt_unique_encseq_processkmercode(GtUniqueEncseqInfo *ueinfo,
                                             GtKmerhitInfo *kmerhitinfo)
{
  GtUword position;
  const GtKmercode *kmercode = ueinfo->kmercodeMain;
  GtLogger * logger = ueinfo->logger;
  unsigned int kmersize = ueinfo->kmersize;
  ueinfo->currentposition =
    position = gt_kmercodeiterator_encseq_get_currentpos(ueinfo->kmercodeitMain)
      - ueinfo->kmersize;
  ueinfo->nPosSinceInsert++;

  if (kmercode->definedspecialposition) {
    gt_logger_log(ueinfo->logger,
                  "special found at position " GT_WU,
                  position);
    gt_unique_encseq_process_special_kmer(ueinfo);
  }
  else {
    ueinfo->maxlen = ueinfo->seqlen + ueinfo->seqstartpos - position;
    gt_assert(ueinfo->maxlen >= (GtUword) kmersize);

    if (ueinfo->initKmerCounter < ueinfo->initUniqueDBsize - kmersize + 1) {
      gt_unique_encseq_init_uniquedb(ueinfo);
    }
    else if ((ueinfo->maxlen <
        ueinfo->minalignlength)) {
      gt_logger_log(logger,
                    "wont proceed matching, remaining sequence is "
                    "shorter than windowsize (or min alignment length)");
      if (position == ueinfo->seqstartpos + ueinfo->seqlen - kmersize) {
        gt_unique_encseq_uniquedb_insert(ueinfo,
                                         ueinfo->nextUniqueStartpos,
                                         ueinfo->seqstartpos + ueinfo->seqlen,
                                         ueinfo->seqnum);
      }
    }
    else {
      gt_unique_encseq_process_kmer_seed(ueinfo,
                                         position,
                                         kmercode,
                                         kmerhitinfo);
    }
  }
}

static GtKmerhitInfo *gt_unique_encseq_kmerhitinfo_new(
    GtUniqueEncseqInfo *ueinfo)
{
  GtKmerhitInfo *kmerhitinfo = gt_malloc(sizeof (GtKmerhitInfo));
  GtKmercodeiterator *kmercodeit1;
  GtXdropbest *xdropbest1 = gt_calloc((size_t) 1, sizeof (GtXdropbest));
  GtXdropbest *xdropbest2 = gt_calloc((size_t) 1, sizeof (GtXdropbest));
  GtBestAlignment * bestAlignment = gt_calloc((size_t) 1,
      sizeof (GtBestAlignment));

  kmercodeit1 = gt_kmercodeiterator_encseq_new(ueinfo->encseq,
                                               GT_READMODE_FORWARD,
                                               ueinfo->kmersize,
                                               0);
  kmerhitinfo->kmercodeit1 = kmercodeit1;
  kmerhitinfo->xdropbest = xdropbest1;
  kmerhitinfo->bestAlignment = bestAlignment;
  kmerhitinfo->bestAlignment->xdropbest = xdropbest2;
  kmerhitinfo->bestAlignment->res = gt_xdrop_resources_new(ueinfo->arbitscores);
  return (kmerhitinfo);
}

static void gt_unique_encseq_kmerhitinfo_delete(GtKmerhitInfo *kmerhitinfo)
{
  gt_kmercodeiterator_delete(kmerhitinfo->kmercodeit1);
  gt_free(kmerhitinfo->xdropbest);
  gt_free(kmerhitinfo->bestAlignment->xdropbest);
  gt_xdrop_resources_delete(kmerhitinfo->bestAlignment->res);
  gt_free(kmerhitinfo->bestAlignment);
  gt_free(kmerhitinfo);
}

/*function to iterate over all kmers of an encseq*/
void gt_unique_encseq_process_all_kmers(GtUniqueEncseqInfo *ueinfo)
{
  GtKmerhitInfo *kmerhitinfo = gt_unique_encseq_kmerhitinfo_new(ueinfo);

  if (ueinfo->totallength < (GtUword) ueinfo->kmersize) {
    return;
  }
  while (!gt_kmercodeiterator_inputexhausted(ueinfo->kmercodeitMain)) {
    ueinfo->kmercodeMain =
        gt_kmercodeiterator_encseq_next(ueinfo->kmercodeitMain);
    if (ueinfo->kmercodeMain != NULL )
      gt_unique_encseq_processkmercode(ueinfo, kmerhitinfo);
    else
      break;
  }
  gt_unique_encseq_kmerhitinfo_delete(kmerhitinfo);
}

int gt_unique_encseq_unit_test(GtError *err)
{
  int had_err = 0;

  gt_ensure(0 == 0);
  return (had_err);
}
