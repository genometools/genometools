/*
  Copyright (c) 2011 Giorgio Gonnella <gonnella@zbh.uni-hamburg.de>
  Copyright (c) 2011 Center for Bioinformatics, University of Hamburg

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

#ifndef MD5SET_H
#define MD5SET_H

/* A set which represents a (possibly large) set of sequences, allowing to
   check whether a sequence or its reverse complement belongs to the set;
   only 128-bit MD5 hashes are stored, not the sequences themselves
   (a low probability of failure due to collisions exists). */
typedef struct GtMD5Set GtMD5Set;

typedef enum {
  GT_MD5SET_ERROR = -1,
  GT_MD5SET_NOT_FOUND,
  GT_MD5SET_FOUND,
  GT_MD5SET_RC_FOUND
} GtMD5SetStatus;

/* Create a new <GtMD5Set> with <nof_elements> sequences. If the number of
   sequences is not known, set <nof_elements> to 0. */
GtMD5Set* gt_md5set_new(unsigned long nof_elements);

/* Deletes a <GtMD5Set> and frees all associated memory. */
void      gt_md5set_delete(GtMD5Set *set);

/* Calculates the MD5 hash of an upper case copy of <seq> and adds it to
   <set> if not present. If <both_strands> is true, the MD5 hash of the reverse
   complement is calculated as well and the forward direction MD5 of <seq> is
   added only if both forward and reverse complement MD5s are not present in
   <set>.
   Returns a negative value on error, <err> is set accordingly. Otherwise,
   returns <GT_MD5SET_NOT_FOUND> if <seq> could be added to the set,
   <GT_MD5SET_FOUND> if <seq> was already present in the forward direction, and
   <GT_MD5SET_RC_FOUND> if <seq> was already present in the reverse complement
   direction. */
GtMD5SetStatus gt_md5set_add_sequence(GtMD5Set *set, const char* seq,
                                 unsigned long seqlen, bool both_strands,
                                 GtError *err);

#endif
