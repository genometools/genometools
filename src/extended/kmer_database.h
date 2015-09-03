/*
  Copyright (c) 2014 Andreas Blaufelder <9blaufel@informatik.uni-hamburg.de>

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

#ifndef KMER_DATABASE_H
#define KMER_DATABASE_H

#include "core/arraydef.h"
#include "core/codetype.h"
#include "core/encseq_api.h"
#include "core/error_api.h"
#include "core/logger_api.h"

/* TODO DW:
   update this, comment and actual functionality are not the same */
/* The <GtKmerDatabase> class stores all kmers occuring in a set of given
   files. Size of k and of the alphabet need to be known beforehand */
typedef struct GtKmerDatabase GtKmerDatabase;

/* An instance of the <GtKmerStartpos> class stores every occuring
   startposition for a given kmer */
typedef struct {
  GtUword *startpos,
          *unique_ids,
          no_positions;
} GtKmerStartpos;

/* Returns new <GtKmerDatabase> object for a given <encseq> file.
   <alphabet_size> specifies the size of the alphabet used in the <encseq>
   files to be read, <kmer_size> specifies k and <sb_man_nu_kmers> gives the
   number of kmers to be buffered before they are inserted into the Database. */
GtKmerDatabase* gt_kmer_database_new(unsigned int alpabet_size,
                                     unsigned int kmer_size,
                                     GtUword sb_max_nu_kmers,
                                     GtEncseq *encseq);

/* Frees space for <GtKmerDatabase>. */
void            gt_kmer_database_delete(GtKmerDatabase *kdb);

/* Adds a single kmer to the Database, used for testing purposes. */
void            gt_kmer_database_add_kmer(GtKmerDatabase *kdb,
                                          GtCodetype kmercode,
                                          GtUword startpos,
                                          GtUword id);

/* Writes the current content of the internal buffer into the Database. */
void            gt_kmer_database_flush(GtKmerDatabase *kdb);

/* Fills the internal buffer with kmers from an interval from <encseq>.
   The interval is specified by <start> and <end> (including <end>). When
   the buffer is full all of it's content will be written to the Database.
   <id> is expected to be unique for each interval and increasing in value. */
void            gt_kmer_database_add_interval(GtKmerDatabase *kdb,
                                              GtUword start,
                                              GtUword end,
                                              GtUword id);

/* Returns an <GtKmerStartpos> object, which returns all startpositions
   of kmer specified through <kmercode>. */
GtKmerStartpos  gt_kmer_database_get_startpos(GtKmerDatabase *kdb,
                                              GtCodetype kmercode);

/* If a kmer occurs more than <cutoff> times it won't be included in the
   <GtKmerDatabase>. */
void            gt_kmer_database_set_cutoff(GtKmerDatabase *kdb,
                                            GtUword cutoff);

/* Disables the use of a cutoff. */
void            gt_kmer_database_disable_cutoff(GtKmerDatabase *kdb);

/* If a kmer occurs more often than the mean/mean_fraction times it won't be
   included in the <GtKmerDatabase>, unless it only occurs min_cutoff times. */
void            gt_kmer_database_use_mean_cutoff(GtKmerDatabase *kdb,
                                                 GtUword mean_fraction,
                                                 GtUword min_cutoff);

/* If <cutoff> and this option are set, kmers occuring more than <cutoff> times
   will be removed during the process. Function fails if no cutoff is set. */
void            gt_kmer_database_set_prune(GtKmerDatabase *kdb);

/* Disables pruning within kmer-database */
void            gt_kmer_database_disable_prune(GtKmerDatabase *kdb);

/* Returns the number of kmers inserted in the <GtKmerDatabase>. */
GtUword         gt_kmer_database_get_kmer_count(GtKmerDatabase *kdb);

/* Returns the arithmetic mean for all kmers inserted in the <GtKmerDatabase>,
   based on the number of different kmers already inserted. If a cutoff is set
   this returns the mean for all kmers which would have been inserted without
   a cutoff. */
GtUword         gt_kmer_database_get_mean_nu_of_occ(GtKmerDatabase *kdb);

/* Returns the minimum occurrence of all occuring kmers in the
   <GtKmerDatabase>. */
GtUword         gt_kmer_database_get_min_nu_of_occ(GtKmerDatabase *kdb);

/* Compares two Databases <a> and <b>. If they differ it returns -1 and an
   appropriate error message is set. */
int             gt_kmer_database_compare(GtKmerDatabase *a,
                                         GtKmerDatabase *b,
                                         GtError *err);

/* Checks if a given Database satisfies it's constraints. If not it returns -1
   and an appropriate error message is set. */
int             gt_kmer_database_check_consistency(GtKmerDatabase *kdb,
                                                   GtError *err);

/* Returns the number of bytes allocated for the <GtKmerDatabase> structure
   (excluding the internal buffer and constants). */
GtUword         gt_kmer_database_get_byte_size(GtKmerDatabase *kdb);

/* Returns the number of bytes the <GtKmerDatabase> structure currently uses
   (excluding the internal buffer and constants). */
GtUword         gt_kmer_database_get_used_size(GtKmerDatabase *kdb);
/* Prints out the content of a Database to the given <GtLogger> object. For
   every kmer in the Database it's code and it's startpositions are printed. */
void            gt_kmer_database_print(GtKmerDatabase *kdb,
                                       GtLogger *logger,
                                       bool verbose);

/* Prints of the content of the internal buffer to the given <GtLogger> object.
   For every kmer in the buffer it's code ans startpositions are printed. */
void            gt_kmer_database_print_buffer(GtKmerDatabase *kdb,
                                              GtLogger *logger);

/* A function to test basic funtionality of the <GtKmerDatabase> class. */
int             gt_kmer_database_unit_test(GtError *err);

#endif
