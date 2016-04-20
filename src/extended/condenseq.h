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

#ifndef CONDENSEQ_H
#define CONDENSEQ_H

#include "core/file_api.h"
#include "core/logger_api.h"
#include "core/range_api.h"
#include "core/str_api.h"
#include "core/types_api.h"
#include "core/disc_distri_api.h"
#include "extended/editscript.h"

#define GT_CONDENSEQ_FILE_SUFFIX ".cse"
#define GT_CONDENSEQ_VERSION 1

/* The <GtCondenseq> class efficiently stores Sequences, either DNA or protein
   by finding redundancies and storing only references with editscripts of
   those.
   */
typedef struct GtCondenseq GtCondenseq;

/* Return new <GtCondenseq> object filled with data read from file with basename
   <indexname>.
   Returns <NULL> on error, or fails if read/file access failes. */
GtCondenseq*       gt_condenseq_new_from_file(const char *indexname,
                                              GtLogger *logger, GtError *err);
/* Write <condenseq> to File <fp>, fails hard on io-error, returns error value
   and sets <err> accordingly on data errors. */
int                gt_condenseq_write(GtCondenseq *condenseq, FILE* fp,
                                      GtError *err);
/* Return the number of sequences contained in the original sequence collection.
 */
GtUword            gt_condenseq_num_of_sequences(const GtCondenseq *condenseq);

/* Return the total length of the original sequence collection, including
   separators (like <GtEncseq>). */
GtUword            gt_condenseq_total_length(const GtCondenseq *condenseq);

/* Return the number of links of <condenseq> */
GtUword            gt_condenseq_num_links(const GtCondenseq *condenseq);

/* Return the number of unique elements of <condenseq> */
GtUword            gt_condenseq_num_uniques(const GtCondenseq *condenseq);

/* Return the total length of all links in <condenseq> */
GtUword            gt_condenseq_total_link_len(const GtCondenseq *condenseq);

/* Return the total length of all uniques in <condenseq> */
GtUword            gt_condenseq_total_unique_len(const GtCondenseq *condenseq);

/* Returns the number of the sequence corresponding to <pos> */
GtUword            gt_condenseq_pos2seqnum(const GtCondenseq *condenseq,
                                           GtUword pos);
/* Returns the starting position of the <seqnum>th sequence in <condenseq> */
GtUword            gt_condenseq_seqstartpos(const GtCondenseq *condenseq,
                                            GtUword seqnum);
/* Returns the id of the <seqnum>th sequence in <condenseq>, sets <*desclen> to
   the length of that string. */
const char*        gt_condenseq_description(const GtCondenseq *condenseq,
                                            GtUword *desclen,
                                            GtUword seqnum);
/* Returns the lenght of sequence <seqnum> from <condenseq> */
GtUword            gt_condenseq_seqlength(const GtCondenseq *condenseq,
                                          GtUword seqnum);
/* Returns the encoded representation of the <id>s sequence of
   <condenseq>. <length> will be set to the length of that sequence. Fails for
   <id>s out of range. */
const GtUchar*     gt_condenseq_extract_encoded(GtCondenseq *condenseq,
                                                GtUword *length,
                                                GtUword id);
/* Returns the encoded representation of the substring defined by (inclusive)
   range <range> of <condenseq>.
   If the positions include a sequence separator this will be set to the non
   printable SEPARATOR (core/chardef.h). */
const GtUchar*     gt_condenseq_extract_encoded_range(GtCondenseq *condenseq,
                                                      GtRange range);
/* Returns the decoded representation of the <id>s sequence of
   <condenseq>. <length> will be set to the length of that sequence. Fails for
   <id>s out of range. */
const char*        gt_condenseq_extract_decoded(GtCondenseq *condenseq,
                                                GtUword *length,
                                                GtUword id);
/* Return the decoded version of the substring defined by (inclusive)
   range <range> of <condenseq>.
   If the positions include a sequence separator <separator> will be used in
   <*buffer>. */
const char*        gt_condenseq_extract_decoded_range(GtCondenseq *condenseq,
                                                      GtRange range,
                                                      char separator);
/* Function type used to process redundant seqs, should return != 0 on error
   and set <err> accordingly. */
typedef int (GtCondenseqProcessExtractedSeqs)(void *data,
                                              GtUword seqid,
                                              GtError *err);
/* Returns the number of sequences within <condenseq> that are similar to unique
   given by <uid>.
   Each such sequence, including the sequence of <uid>, will be passed to
   <callback> along with <callback_data>, its number and <err>.
   Returns 0 on error as there is at least one similar range, the input range
   itself. */
GtUword            gt_condenseq_each_redundant_seq(
                                       const GtCondenseq *condenseq,
                                       GtUword uid,
                                       GtCondenseqProcessExtractedSeqs callback,
                                       void *callback_data,
                                       GtError *err);
/* Function type used to process redundant ranges, should return != 0 on error
   and set <err> accordingly. */
typedef int        (GtCondenseqProcessExtractedRange)(void *data,
                                                      GtUword seqid,
                                                      GtRange seqrange,
                                                      GtError *err);
/* Returns the number of ranges within <condenseq> that are similar to range
   given by <uid> and <urange> which should be relative to the start of the
   unique range.
   Each such range will be passed to <callback> along with <callback_data>, its
   number and <err>. If extend is > 0 the resulting ranges will be extended in
   the obvious direction by that amount without crossing sequence boundaries.
   Returns 0 on error as there is at least one similar range, the input range
   itself. */
GtUword             gt_condenseq_each_redundant_range(
                                      const GtCondenseq *condenseq,
                                      GtUword uid,
                                      GtRange urange,
                                      GtUword left_extend,
                                      GtUword right_extend,
                                      GtCondenseqProcessExtractedRange callback,
                                      void *callback_data,
                                      GtError *err);
/* Returns a pointer to a string comprised of the basefilename of the archive,
   excluding any suffix or path, caller is responsible for freeing the received
   pointer.
   Returns NULL if no such file exists. */
char*               gt_condenseq_basefilename(const GtCondenseq *condenseq);

/* Returns a new <GtStr> object containing the path to the fasta file containing
   the unique part of the archive in <condenseq>.
   Returns NULL if no such file exists. */
GtStr*              gt_condenseq_unique_fasta_file(
                                                  const GtCondenseq *condenseq);
/* Returns the number of uniques with length at least <min_align_len> */
GtUword             gt_condenseq_count_relevant_uniques(
                                                   const GtCondenseq *condenseq,
                                                   unsigned int min_align_len);
/* return size in bytes of the <condenseq>, other parameters will be set to the
   sise in bytes of the various elements. where <links> does not contain the
   size of the included editscript, but only the pointer to it. */
GtUword             gt_condenseq_size(const GtCondenseq *condenseq,
                                      GtUword *uniques,
                                      GtUword *links,
                                      GtUword *editscripts,
                                      GtUword *descriptions,
                                      GtUword *separators);
/* Creates an gff3-File with the basename of the index containing the unique and
   link ranges as experimental_features */
int                 gt_condenseq_output_to_gff3(const GtCondenseq *condenseq,
                                                GtError *err);
/* Returns a distripution of the lengths of the unique elements. */
GtDiscDistri*       gt_condenseq_unique_length_dist(
                                                  const GtCondenseq *condenseq);
/* Returns a distripution of the lengths of the link elements. */
GtDiscDistri*       gt_condenseq_link_length_dist(const GtCondenseq *condenseq);

/* Returns a distripution of the compression ratios of all editscripts. */
GtDiscDistri*       gt_condenseq_link_comp_dist(const GtCondenseq *condenseq);

/* Returns the original seqnum from which the unique with <uid> derives, changes
   urange (which has coordinate relative to that unique) so it represents the
   same range but relative to the sequence collection. */
GtUword             gt_condenseq_unique_range_to_seqrange(
                                                   const GtCondenseq *condenseq,
                                                   GtUword uid,
                                                   GtRange *urange);
/* Returns a pointer to the editscript corresponding to <lid> in <condenseq>,
   <condenseq> retains ownership of that editscript. */
const GtEditscript* gt_condenseq_link_editscript(const GtCondenseq *condenseq,
                                                 GtUword lid);
/* Returns a reference to the <GtAlphabet> on which the sequences within
   <condenseq> are based. */
GtAlphabet*         gt_condenseq_alphabet(const GtCondenseq *condenseq);

/* Free space for <condenseq> */
void                gt_condenseq_delete(GtCondenseq *condenseq);
#endif
