/*
  Copyright (c) 2009      Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
  Copyright (c)      2010 Sascha Steinbiss <steinbiss@zbh.uni-hamburg.de>
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

#ifndef ENCODEDSEQUENCE_REP_H
#define ENCODEDSEQUENCE_REP_H

/*
  The contents of this file is to be considered private
  implementation detail but, whenever the code is compiled with option
  INLINEDENCSEQ, is exposed to the compiler solely for performance
  optimization. So we can compare the time overhead of a bytearray
  implementation of strings to all other representations implemented in
  encodedseq.c.
*/

#include "core/alphabet.h"
#include "core/bitpackarray.h"
#include "core/filelengthvalues.h"
#include "core/intbits.h"
#include "core/intdef.h"

#include "core/str_array_api.h"
#include "core/symboldef.h"

typedef enum
{
  Viadirectaccess,
  Viabytecompress,
  Viabitaccess,
  Viauchartables,
  Viaushorttables,
  Viauint32tables,
  Undefpositionaccesstype
} GtPositionaccesstype;

typedef uint32_t Uint32;

struct GtEncodedsequence
{
  /* Common part */
  unsigned long *satcharptr; /* need for writing char */
  GtPositionaccesstype sat;
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
  GtUchar(*deliverchar)(const GtEncodedsequence *,unsigned long);
  const char *delivercharname;
  GtUchar(*delivercharnospecial)(const GtEncodedsequence *,unsigned long);
  const char *delivercharnospecialname;
  GtUchar(*seqdeliverchar)(const GtEncodedsequence *,
                         GtEncodedsequenceScanstate *,unsigned long);
  const char *seqdelivercharname;
  bool(*delivercontainsspecial)(const GtEncodedsequence *,
                                bool,
                                GtEncodedsequenceScanstate *,
                                unsigned long,
                                unsigned long);
  const char *delivercontainsspecialname;
  unsigned int maxspecialtype;  /* maximal value of special type */
  unsigned long *characterdistribution;
  GtSpecialcharinfo *specialcharinfoptr, /* need for writing specialcharinfo */
                  specialcharinfo; /* information about specialcharacters */
  GtStrArray *filenametab;    /* table of filenames */
  char *firstfilename;
  GtFilelengthvalues *filelengthtab;  /* table of length of files */

  const char *destab;
  unsigned long destablength, *sdstab;

  const GtAlphabet *alpha;   /* alphabet representation */

  const unsigned long *ssptab; /* (if numofdbsequences = 1 then NULL  else
                                                         numofdbsequences  -1)
                           entries */

  /* only for Viabitaccess,
              Viauchartables,
              Viaushorttables,
              Viauint32tables */

  GtTwobitencoding *twobitencoding;
  unsigned long unitsoftwobitencoding;

  /* only for Viauchartables,
              Viaushorttables,
              Viauint32tables */

  /* only for Viadirectaccess */
  GtUchar *plainseq;
  bool hasplainseqptr;

  /* only for Viabytecompress */
  BitPackArray *bitpackarray;
  unsigned int numofchars; /* used to have faster access in getencodedchar */

  /* only for Viabitaccess */
  GtBitsequence *specialbits;

  /* only for Viauchartables */
  GtUchar *ucharspecialpositions,
        *ucharspecialrangelength;
  unsigned long *ucharendspecialsubsUint;

  /* only for Viaushorttables */
  GtUshort *ushortspecialpositions,
         *ushortspecialrangelength;
  unsigned long *ushortendspecialsubsUint;

  /* only for Viauint32tables */
  Uint32 *uint32specialpositions,
         *uint32specialrangelength;
  unsigned long *uint32endspecialsubsUint;
};
#endif
