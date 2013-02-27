/*
  Copyright (c) 2009-2011 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
  Copyright (c)      2011 Sascha Steinbiss <steinbiss@zbh.uni-hamburg.de>
  Copyright (c)      2010 Dirk Willrodt <willrodt@zbh.uni-hamburg.de>
  Copyright (c) 2009-2011 Center for Bioinformatics, University of Hamburg

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
  implementation detail.
*/

#include "core/alphabet.h"
#include "core/bitpackarray.h"
#include "core/chardef.h"
#include "core/encseq_access_type.h"
#include "core/encseq_api.h"
#include "core/filelengthvalues.h"
#include "core/intbits.h"
#include "core/md5_tab.h"
#include "core/types_api.h"
#include "core/str_array_api.h"
#include "core/defined-types.h"
#include "core/types_api.h"
#include "core/thread_api.h"

typedef struct
{
  GtUchar *positions,
          *rangelengths;
  unsigned long *endidxinpage;
  unsigned long *mappositions;
  unsigned long numofpages;
  unsigned long numofpositionstostore;
  unsigned long maxrangevalue;  /* maximal value of special type */
} GtSWtable_uchar;

typedef struct
{
  uint16_t *positions,
           *rangelengths;
  unsigned long *endidxinpage;
  unsigned long *mappositions;
  unsigned long numofpages;
  unsigned long numofpositionstostore;
  unsigned long maxrangevalue;  /* maximal value of special type */
} GtSWtable_uint16;

typedef struct
{
  uint32_t *positions,
           *rangelengths;
  unsigned long *endidxinpage;
  unsigned long *mappositions;
  unsigned long numofpages;
  unsigned long numofpositionstostore;
  unsigned long maxrangevalue;  /* maximal value of special type */
} GtSWtable_uint32;

typedef union
{
  GtSWtable_uchar st_uchar;
  GtSWtable_uint16 st_uint16;
  GtSWtable_uint32 st_uint32;
} GtSWtable;

typedef struct
{
  GtUchar *is64bitptr;
  unsigned long *versionptr,
                *satcharptr,
                *totallengthptr,
                *numofdbsequencesptr,
                *numofdbfilesptr,
                *lengthofdbfilenamesptr,
                *minseqlenptr,
                *maxseqlenptr,
                *alphatypeptr,
                *lengthofalphadefptr,
                *numofallcharsptr;
  GtUchar *maxsubalphasizeptr;
  GtSpecialcharinfo *specialcharinfoptr;
  char *firstfilename,
       *alphadef;
  GtFilelengthvalues *filelengthtab;
  unsigned long *characterdistribution;
} GtEncseqHeaderPtr;

typedef struct
{
  unsigned long *classstartpositionsptr;
  char *maxcharsptr,
       *allcharsptr;
  GtUchar *subsymbolmapptr;
} GtExceptionTablePtr;

struct GtEncseq
{
  /* Common part */
  GtEncseqAccessType sat, satsep;
  const char *satname;
  void *mappedptr, /* NULL or pointer to the mapped space block */
       *ssptabmappedptr, /* NULL or pointer to the mapped space block */
       *oistabmappedptr;
  bool has_specialranges,
       has_wildcardranges,
       has_ssptab;
  unsigned long totallength,
                logicaltotallength,
                numofdbsequences,
                logicalnumofdbsequences,
                numofdbfiles,
                lengthofdbfilenames,
                sizeofrep;
  char *indexname;

  GtEncseqHeaderPtr headerptr;

  unsigned long version;
  GtUchar is64bit;

  GtUchar(*seqdeliverchar)(GtEncseqReader *);
  const char *seqdelivercharname;
  bool(*delivercontainsspecial)(const GtEncseq *,
                                GtReadmode,
                                GtEncseqReader *,
                                unsigned long,
                                unsigned long);
  const char *delivercontainsspecialname;
  bool(*issinglepositioninspecialrange)(const GtEncseq *,unsigned long);
  const char *issinglepositioninspecialrangename;
  bool(*issinglepositioninwildcardrange)(const GtEncseq *,unsigned long);
  const char *issinglepositioninwildcardrangename;
  bool(*getexceptionmapping)(const GtEncseq *,unsigned long*, unsigned long);
  const char *getexceptionmappingname;
  bool(*issinglepositionseparator)(const GtEncseq *,unsigned long);
  const char *issinglepositionseparatorname;

  unsigned int leastprobablecharacter;
  GtSpecialcharinfo specialcharinfo; /* information about specialcharacters */
  Definedunsignedlong equallength;
  GtStrArray *filenametab;    /* table of filenames */

  char *destab;
  bool hasallocateddestab,
       hasallocatedssptab,
       hasallocatedsdstab;
  unsigned long destablength, *sdstab;

  /* alphabet representation */
  GtAlphabet *alpha;
  char *alphadef;
  unsigned long lengthofalphadef,
                alphatype;

  /* separator index structure */
  GtSWtable ssptab;

  /* file start position table */
  unsigned long *fsptab; /* is NULL when numofdbfiles is 1
                            otherwise has numofdbfiles - 1 entries */

  /* only for GT_ACCESS_TYPE_EQUALLENGTH,
              GT_ACCESS_TYPE_BITACCESS,
              GT_ACCESS_TYPE_UCHARTABLES,
              GT_ACCESS_TYPE_USHORTTABLES,
              GT_ACCESS_TYPE_UINT32TABLES */
  GtTwobitencoding *twobitencoding;
  unsigned long unitsoftwobitencoding;

  /* only for  GT_ACCESS_TYPE_DIRECTACCESS */
  GtUchar *plainseq;
  bool hasplainseqptr;

  /* only for GT_ACCESS_TYPE_BYTECOMPRESS */
  BitPackArray *bitpackarray;
  unsigned int numofchars;

  /* only for GT_ACCESS_TYPE_BITACCESS */
  GtBitsequence *specialbits;

  /* only for GT_ACCESS_TYPE_UCHARTABLES,
              GT_ACCESS_TYPE_USHORTTABLES,
              GT_ACCESS_TYPE_UINT32TABLES */
  GtSWtable wildcardrangetable;

  /* for lossless reproduction of original sequences */
  bool has_exceptiontable;
  GtSWtable exceptiontable;
  BitPackArray *exceptions;
  unsigned long *classstartpositions,
                numofallchars;
  char *maxchars,
       *allchars;
  unsigned char *subsymbolmap,
                maxsubalphasize;
  GtExceptionTablePtr exceptionheaderptr;
  GtEncseqAccessType oissat;

  /* reference counting */
  unsigned long reference_count;
  GtMutex *refcount_lock;

  /* MD5 sums */
  GtMD5Tab *md5_tab;

  bool hasmirror,
       accesstype_via_utables;

  unsigned long minseqlen,
                maxseqlen;
};
#endif
