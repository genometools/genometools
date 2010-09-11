/*
  Copyright (c) 2009      Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
  Copyright (c)      2010 Sascha Steinbiss <steinbiss@zbh.uni-hamburg.de>
  Copyright (c)      2010 Dirk Willrodt <dwillrodt@zbh.uni-hamburg.de>
  Copyright (c) 2009-2010 Center for Bioinformatics, University of Hamburg

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

#ifndef ENCSEQ_REP_H
#define ENCSEQ_REP_H

/*
  The contents of this file is to be considered private
  implementation detail but, whenever the code is compiled with option
  GT_INLINEDENCSEQ, is exposed to the compiler solely for performance
  optimization. So we can compare the time overhead of a bytearray
  implementation of strings to all other representations implemented in
  encseq.c.
*/

#include "core/alphabet.h"
#include "core/bitpackarray.h"
#include "core/chardef.h"
#include "core/encseq_access_type.h"
#include "core/filelengthvalues.h"
#include "core/intbits.h"
#include "core/types_api.h"
#include "core/str_array_api.h"
#include "core/defined-types.h"
#include "core/types_api.h"
#include "core/thread.h"

typedef uint32_t Uint32;

struct GtEncseq
{
  /* Common part */
  unsigned long *satcharptr; /* need for writing char */
  GtEncseqAccessType sat;
  void *mappedptr; /* NULL or pointer to the mapped space block */
  unsigned long numofspecialstostore;
  unsigned long *totallengthptr,
                totallength;
  unsigned long numofdbsequences,
                *numofdbsequencesptr; /* need for writing numofdbsequences */
  unsigned long numofdbfiles, *numofdbfilesptr;
  unsigned long lengthofdbfilenames, *lengthofdbfilenamesptr;
  unsigned long sizeofrep;
  const char *name;
  GtUchar(*deliverchar)(const GtEncseq *,unsigned long);
  const char *delivercharname;
  GtUchar(*delivercharnospecial)(const GtEncseq *,unsigned long);
  const char *delivercharnospecialname;
  GtUchar(*seqdeliverchar)(GtEncseqReader *);
  const char *seqdelivercharname;
  bool(*delivercontainsspecial)(const GtEncseq *,
                                GtReadmode,
                                GtEncseqReader *,
                                unsigned long,
                                unsigned long);
  const char *delivercontainsspecialname;
  unsigned long numofspecialcells; 
  /* encseq->totallength/encseq->maxspecialtype + 1;*/
  unsigned int maxspecialtype;  /* maximal value of special type */
  unsigned long *characterdistribution;
  GtSpecialcharinfo *specialcharinfoptr, /* need for writing specialcharinfo */
                  specialcharinfo; /* information about specialcharacters */
  Definedunsignedlong equallength;
  GtStrArray *filenametab;    /* table of filenames */
  char *firstfilename;
  GtFilelengthvalues *filelengthtab;  /* table of length of files */

  char *destab;
  bool hasallocateddestab,
       hasallocatedssptab,
       hasallocatedsdstab;
  unsigned long destablength, *sdstab;

  GtAlphabet *alpha;   /* alphabet representation */

  unsigned long *ssptab; /* (if numofdbsequences = 1 then NULL  else
                                                           numofdbsequences  -1)
                                                           entries */
  unsigned long *fsptab; /* (if numofdbfiles = 1 then NULL  else
                                                      numofdbfiles  -1
                                                      entries */

  /* only for GT_ACCESS_TYPE_BITACCESS,
              GT_ACCESS_TYPE_UCHARTABLES,
              GT_ACCESS_TYPE_USHORTTABLES,
              GT_ACCESS_TYPE_UINT32TABLES */

  GtTwobitencoding *twobitencoding;
  unsigned long unitsoftwobitencoding;

  /* only for GT_ACCESS_TYPE_UCHARTABLES,
              GT_ACCESS_TYPE_USHORTTABLES,
              GT_ACCESS_TYPE_UINT32TABLES */

  /* only for  GT_ACCESS_TYPE_DIRECTACCESS */
  GtUchar *plainseq;
  bool hasplainseqptr;

  /* only for GT_ACCESS_TYPE_BYTECOMPRESS */
  BitPackArray *bitpackarray;
  unsigned int numofchars; /* used to have faster access in getencodedchar */

  /* only for GT_ACCESS_TYPE_BITACCESS */
  GtBitsequence *specialbits;

  /* only for GT_ACCESS_TYPE_UCHARTABLES */
  GtUchar *ucharspecialpositions,
          *ucharspecialrangelength;
  unsigned long *ucharendspecialsubsUint;

  /* only for GT_ACCESS_TYPE_USHORTTABLES */
  GtUshort *ushortspecialpositions,
           *ushortspecialrangelength;
  unsigned long *ushortendspecialsubsUint;

  /* only for GT_ACCESS_TYPE_UINT32TABLES */
  Uint32 *uint32specialpositions,
         *uint32specialrangelength;
  unsigned long *uint32endspecialsubsUint;

  unsigned long reference_count;
  GtMutex *refcount_lock;
};
#endif
