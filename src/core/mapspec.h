/*
  Copyright (c) 2007      Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
  Copyright (c)      2011 Sascha Steinbiss <steinbiss@zbh.uni-hamburg.de>
  Copyright (c) 2007-2011 Center for Bioinformatics, University of Hamburg

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

#ifndef MAPSPEC_H
#define MAPSPEC_H

#include "core/arraydef.h"
#include "core/bitpackarray.h"
#include "core/chardef.h"
#include "core/error.h"
#include "core/intbits.h"
#include "core/filelengthvalues.h"
#include "core/pairbwtidx.h"
#include "core/str_api.h"
#include "core/types_api.h"
#include "core/ulongbound.h"

/* The <GtMapspec> class contains the structure of a memory region mapped to
   and from a file in terms of the types and lengths of its serialized
   components. */
typedef struct GtMapspec GtMapspec;

/* This function is used to describe the structure of the serialized
   data by iteratively adding pointers to the sources (for writing) or
   destinations (for reading) of the data of a given type and length to the
   <mapspec>. The <data> pointer can be used to pass arbitrary data into this
   function. If this function is called from a writing context, <readwrite> is
   <TRUE>, otherwise it is <FALSE>. */
typedef void (*GtMapspecSetupFunc)(GtMapspec *mapspec, void *data,
                                   bool readwrite);

/* Adds a C string of length <n>, to be read from or written to <ptr> to the
   <mapspec>. */
#define gt_mapspec_add_char(MAPSPEC, PTR, N)\
        gt_mapspec_add_char_ptr(MAPSPEC, &(PTR), N);
/* Adds a C string of length <n>, to be read from or written to the address to
   which <ptr> points to the <mapspec>. Usually, this method is not used
   directly, but <gt_mapspec_add_char()> is used instead. */
void gt_mapspec_add_char_ptr(GtMapspec *mapspec, char **ptr, unsigned long n);
/* Adds a unsigned character array of length <n>, to be read from or written
   to <ptr> to the <mapspec>. */
#define gt_mapspec_add_uchar(MAPSPEC, PTR, N)\
        gt_mapspec_add_uchar_ptr(MAPSPEC, &(PTR), N);
/* Adds a unsigned character array of length <n>, to be read from or written
   to <ptr> to the <mapspec>. */
void gt_mapspec_add_uchar_ptr(GtMapspec *mapspec, GtUchar **ptr,
                              unsigned long n);
/* Adds a uint16_t array of length <n>, to be read from or written
   to <ptr> to the <mapspec>. */
#define gt_mapspec_add_uint16(MAPSPEC, PTR, N)\
        gt_mapspec_add_uint16_ptr(MAPSPEC, &(PTR), N);
/* Adds a uint16_t array of length <n>, to be read from or written
   to <ptr> to the <mapspec>. */
void gt_mapspec_add_uint16_ptr(GtMapspec *mapspec, uint16_t **ptr,
                               unsigned long n);
/* Adds a unsigned long array of length <n>, to be read from or written
   to <ptr> to the <mapspec>. */
#define gt_mapspec_add_ulong(MAPSPEC, PTR, N)\
        gt_mapspec_add_ulong_ptr(MAPSPEC, &(PTR), N);
/* Adds a unsigned long array of length <n>, to be read from or written
   to <ptr> to the <mapspec>. */
void gt_mapspec_add_ulong_ptr(GtMapspec *mapspec, unsigned long **ptr,
                              unsigned long n);
/* Adds a <GtUlongBound> array of length <n>, to be read from or written
   to <ptr> to the <mapspec>. */
#define gt_mapspec_add_ulongbound(MAPSPEC, PTR, N)\
        gt_mapspec_add_ulongbound_ptr(MAPSPEC, &(PTR), N);
/* Adds a <GtUlongBound> array of length <n>, to be read from or written
   to <ptr> to the <mapspec>. */
void gt_mapspec_add_ulongbound_ptr(GtMapspec *mapspec, GtUlongBound **ptr,
                                   unsigned long n);
/* Adds a uint32_t array of length <n>, to be read from or written
   to <ptr> to the <mapspec>. */
#define gt_mapspec_add_uint32(MAPSPEC, PTR, N)\
        gt_mapspec_add_uint32_ptr(MAPSPEC, &(PTR), N);
/* Adds a uint32_t array of length <n>, to be read from or written
   to <ptr> to the <mapspec>. */
void gt_mapspec_add_uint32_ptr(GtMapspec *mapspec, uint32_t **ptr,
                               unsigned long n);
/* Adds a uint64_t array of length <n>, to be read from or written
   to <ptr> to the <mapspec>. */
#define gt_mapspec_add_uint64(MAPSPEC, PTR, N)\
        gt_mapspec_add_uint64_ptr(MAPSPEC, &(PTR), N);
/* Adds a uint64_t array of length <n>, to be read from or written
   to <ptr> to the <mapspec>. */
void gt_mapspec_add_uint64_ptr(GtMapspec *mapspec, uint64_t **ptr,
                               unsigned long n);
/* Adds a <GtBitsequence> array of length <n>, to be read from or written
   to <ptr> to the <mapspec>. */
#define gt_mapspec_add_bitsequence(MAPSPEC, PTR, N)\
        gt_mapspec_add_bitsequence_ptr(MAPSPEC, &(PTR), N);
/* Adds a <GtBitsequence> array of length <n>, to be read from or written
   to <ptr> to the <mapspec>. */
void gt_mapspec_add_bitsequence_ptr(GtMapspec *mapspec, GtBitsequence **ptr,
                                    unsigned long n);
/* Adds a <GtTwobitencoding> array of length <n>, to be read from or written
   to <ptr> to the <mapspec>. */
#define gt_mapspec_add_twobitencoding(MAPSPEC, PTR, N)\
        gt_mapspec_add_twobitencoding_ptr(MAPSPEC, &(PTR), N);
/* Adds a <GtTwobitencoding> array of length <n>, to be read from or written
   to <ptr> to the <mapspec>. */
void gt_mapspec_add_twobitencoding_ptr(GtMapspec *mapspec,
                                       GtTwobitencoding **ptr, unsigned long n);
/* Adds a <GtSpecialcharinfo> array of length <n>, to be read from or written
   to <ptr> to the <mapspec>. */
#define gt_mapspec_add_specialcharinfo(MAPSPEC, PTR, N)\
        gt_mapspec_add_specialcharinfo_ptr(MAPSPEC, &(PTR), N);
/* Adds a <GtSpecialcharinfo> array of length <n>, to be read from or written
   to <ptr> to the <mapspec>. */
void gt_mapspec_add_specialcharinfo_ptr(GtMapspec *mapspec,
                                        GtSpecialcharinfo **ptr,
                                        unsigned long n);
/* Adds a <BitElem> array of length <n>, to be read from or written
   to <ptr> to the <mapspec>. */
#define gt_mapspec_add_bitelem(MAPSPEC, PTR, N)\
        gt_mapspec_add_bitelem_ptr(MAPSPEC, &(PTR), N);
/* Adds a <BitElem> array of length <n>, to be read from or written
   to <ptr> to the <mapspec>. */
void gt_mapspec_add_bitelem_ptr(GtMapspec *mapspec, BitElem **ptr,
                                unsigned long n);
/* Adds a <GtFilelengthvalues> array of length <n>, to be read from or written
   to <ptr> to the <mapspec>. */
#define gt_mapspec_add_filelengthvalues(MAPSPEC, PTR, N)\
        gt_mapspec_add_filelengthvalues_ptr(MAPSPEC, &(PTR), N);
/* Adds a <GtFilelengthvalues> array of length <n>, to be read from or written
   to <ptr> to the <mapspec>. */
void gt_mapspec_add_filelengthvalues_ptr(GtMapspec *mapspec,
                                         GtFilelengthvalues **ptr,
                                         unsigned long n);
/* Adds a <GtPairBwtidx> array of length <n>, to be read from or written
   to <ptr> to the <mapspec>. */
#define gt_mapspec_add_pairbwtindex(MAPSPEC, PTR, N)\
        gt_mapspec_add_pairbwtindex_ptr(MAPSPEC, &(PTR), N);
/* Adds a <GtPairBwtidx> array of length <n>, to be read from or written
   to <ptr> to the <mapspec>. */
void gt_mapspec_add_pairbwtindex_ptr(GtMapspec *mapspec, GtPairBwtidx **ptr,
                                     unsigned long n);

/* Runs <setup> to build the map specification using <data> if given,
   then maps the file specified by <filename> with the expected size
   <expectedsize> and fills the pointers in the map specification with their
   corresponding values. The beginning of the mapped area is written to
   <mapped>. Returns 0 on success, -1 otherwise. <err> is set accordingly. */
int  gt_mapspec_read(GtMapspecSetupFunc setup, void *data,
                     const GtStr *filename, unsigned long expectedsize,
                     void **mapped, GtError *err);
/* Runs <setup> to build the map specification using <data> if given,
   then writes the data at the pointers given in the map specification to the
   file specified by <fp> with the expected file size <expectedsize>.
   Returns 0 on success, -1 otherwise. <err> is set accordingly. */
int  gt_mapspec_write(GtMapspecSetupFunc setup, FILE *fp, void *data,
                      unsigned long expectedsize, GtError *err);
/* Pads file specified by <fp> at position <byteoffset> with zero bytes up to
   the next word boundary. The amount of padding in bytes is written to
   <bytes_written>. Returns 0 on success, -1 otherwise. <err> is set
   accordingly.*/
int  gt_mapspec_pad(FILE *fp, unsigned long *bytes_written,
                    unsigned long byteoffset, GtError *err);

#endif
