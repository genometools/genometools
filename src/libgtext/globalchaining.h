/*
  Copyright (c)       2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2004       Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
  Copyright (c) 2004, 2007 Center for Bioinformatics, University of Hamburg

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

#ifndef GLOBALCHAINING_H
#define GLOBALCHAINING_H

#include "libgtcore/error.h"
#include "libgtext/chain.h"

typedef struct {
  unsigned long startpos1, /* start of fragment in first sequence */
                endpos1,   /* end of fragment in first sequence */
                startpos2, /* start of fragment in second sequence */
                endpos2;   /* end of fragment in second sequence */
  long weight;             /* weight of fragment */
} Fragment;

typedef void (*ChainProc)(Chain*, Fragment*, void*);

/* Perform global chaining with overlaps of <num_of_fragments> many <fragments>
   in O(<num_of_fragments>^2) time.
   Two fragments can maximally be <max_gap_width> many bases away.
   For all global chains of maximal score, the ChainProc function is called.
   Thereby, ChainProc does not get the ownership of the Chain. */
void globalchaining_max(Fragment *fragments, unsigned long num_of_fragments,
                        unsigned long max_gap_width, ChainProc, void *cpinfo);

/* Perform global chaining with overlaps of <num_of_fragments> many <fragments>
   in O(<num_of_fragments>^2) time.
   Two fragments can maximally be <max_gap_width> many bases away.
   For all non-overlapping global chains with a coverage of more then
   <mincoverage> of the sequence in dimension 1 (with length <seqlen1>), the
   ChainProc function is called. Thereby, ChainProc does not get the ownership
   of the Chain. */
void globalchaining_coverage(Fragment *fragments,
                             unsigned long num_of_fragments,
                             unsigned long max_gap_width, unsigned long seqlen1,
                             double mincoverage, ChainProc, void *cpinfo);

#endif
