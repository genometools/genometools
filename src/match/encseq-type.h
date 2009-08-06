/*
  Copyright (c) 2009 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
  Copyright (c) 2009 Center for Bioinformatics, University of Hamburg

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
#ifndef ENCSEQ_TYPE_H
#define ENCSEQ_TYPE_H

/*
  The contents of this file is to be considered private
  implementation detail but, whenever the code is compiled with option
  INLINEDENCSEQ, is exposed to the compiler solely for performance
  optimization. So we can compare the time overhead of a bytearray
  implementation of strings to all other represenations implemented in
  encodedseq.c.
*/

#include "core/alphabet.h"
#include "core/symboldef.h"
#include "core/str_array_api.h"
#include "core/filelengthvalues.h"
#include "bitpack-itf.h"
#include "seqpos-def.h"
#include "intbits.h"
#include "ushort-def.h"

typedef enum
{
  Viadirectaccess,
  Viabytecompress,
  Viabitaccess,
  Viauchartables,
  Viaushorttables,
  Viauint32tables,
  Undefpositionaccesstype
} Positionaccesstype;

typedef uint32_t Uint32;

struct Encodedsequence
{
  /* Common part */
  unsigned long *satcharptr; /* need for writing char */
  Positionaccesstype sat;
  void *mappedptr; /* NULL or pointer to the mapped space block */
  unsigned long numofspecialstostore;
  Seqpos *totallengthptr,
         totallength;
  unsigned long numofdbsequences,
                *numofdbsequencesptr; /* need for writing numofdbsequences */
  unsigned long sizeofrep;
  const char *name;
  GtUchar(*deliverchar)(const Encodedsequence *,Seqpos);
  const char *delivercharname;
  GtUchar(*delivercharnospecial)(const Encodedsequence *,Seqpos);
  const char *delivercharnospecialname;
  GtUchar(*seqdeliverchar)(const Encodedsequence *,
                         Encodedsequencescanstate *,Seqpos);
  const char *seqdelivercharname;
  bool(*delivercontainsspecial)(const Encodedsequence *,
                                bool,Encodedsequencescanstate *,Seqpos,Seqpos);
  const char *delivercontainsspecialname;
  unsigned int maxspecialtype;  /* maximal value of special type */
  unsigned long *characterdistribution;
  Specialcharinfo *specialcharinfoptr, /* need for writing specialcharinfo */
                  specialcharinfo; /* information about specialcharacters */

  const GtStrArray *filenametab;    /* table of filenames */
  const Filelengthvalues *filelengthtab;  /* table of length of files */

  const char *destab;
  unsigned long destablength, *sdstab;

  const GtAlphabet *alpha;   /* alphabet representation */

  const Seqpos *ssptab; /* (if numofdbsequences = 1 then NULL  else
                                                         numofdbsequences  -1)
                           entries */

  /* only for Viabitaccess,
              Viauchartables,
              Viaushorttables,
              Viauint32tables */

  Twobitencoding *twobitencoding;
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
  Bitsequence *specialbits;

  /* only for Viauchartables */
  GtUchar *ucharspecialpositions,
        *ucharspecialrangelength;
  unsigned long *ucharendspecialsubsUint;

  /* only for Viaushorttables */
  Ushort *ushortspecialpositions,
         *ushortspecialrangelength;
  unsigned long *ushortendspecialsubsUint;

  /* only for Viauint32tables */
  Uint32 *uint32specialpositions,
         *uint32specialrangelength;
  unsigned long *uint32endspecialsubsUint;
};
#endif
