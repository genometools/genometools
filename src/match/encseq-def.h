/*
  Copyright (c) 2007 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
  Copyright (c) 2007 Center for Bioinformatics, University of Hamburg

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

#ifndef ENCSEQ_DEF_H
#define ENCSEQ_DEF_H

#include "core/alphabet.h"
#include "core/chardef.h"
#include "core/str.h"
#include "core/str_array.h"
#include "core/symboldef.h"
#include "core/filelengthvalues.h"
#include "core/disc_distri.h"
#include "seqpos-def.h"
#include "intcode-def.h"
#include "intbits.h"
#include "readmode-def.h"
#include "verbose-def.h"

#define DBFILEKEY "dbfile="

#ifdef SKDEBUG
#define CHECKENCCHAR(CC,ENCSEQ,POS,READMODE)\
        {\
          GtUchar cctmp = getencodedchar(ENCSEQ,POS,READMODE);\
          if ((CC) != cctmp)\
          {\
            printf("file %s, line %d: pos = %lu:cc = %u != %u = ccreal\n",\
                   __FILE__,__LINE__,\
                   (unsigned long) (POS),\
                   (unsigned int) (CC),\
                   (unsigned int) cctmp);\
            exit(GT_EXIT_PROGRAMMING_ERROR);\
          }\
        }
#else
#define CHECKENCCHAR(CC,ENCSEQ,POS,READMODE)
#endif

typedef struct
{
  Seqpos leftpos,
         rightpos;
} Sequencerange;          /* \Typedef{Sequencerange} */

typedef struct
{
  Seqpos seqstartpos,  /* the position of the first character in the encseq */
         seqlength;    /* the length of the sequence */
} Seqinfo;             /* \Typedef{Seqinfo} */

typedef struct
{
  Twobitencoding tbe;           /* two bit encoding */
  unsigned int unitsnotspecial; /* units which are not special */
  Seqpos position;
} EndofTwobitencoding;

typedef struct Encodedsequence Encodedsequence;
typedef struct Encodedsequencescanstate Encodedsequencescanstate;
typedef struct Specialrangeiterator Specialrangeiterator;

#undef INLINEDENCSEQ
#ifdef INLINEDENCSEQ
#include "encseq-type.h"
#endif

#ifdef INLINEDENCSEQ
#define getencseqtotallength(ENCSEQ) ((ENCSEQ)->totallength)
#else
Seqpos getencseqtotallength(const Encodedsequence *encseq);
#endif

#ifdef INLINEDENCSEQ
#define getencseqnumofdbsequences(ENCSEQ) ((ENCSEQ)->numofdbsequences)
#else
unsigned long getencseqnumofdbsequences(const Encodedsequence *encseq);
#endif

#define REVERSEPOS(TOTALLENGTH,POS) ((TOTALLENGTH) - 1 - (POS))

#ifdef INLINEDENCSEQ
#define MAKECOMPL(CC)\
        (ISSPECIAL(CC) ? (CC) : (GtUchar) 3 - (CC))
/*@unused@*/ static inline GtUchar getencodedchar(const Encodedsequence *encseq,
                                                Seqpos pos,
                                                Readmode readmode)
{
  return (readmode == Forwardmode)
          ? encseq->plainseq[pos]
          : ((readmode == Reversemode)
            ? encseq->plainseq[REVERSEPOS(encseq->totallength,pos)]
            : ((readmode == Complementmode)
              ? MAKECOMPL(encseq->plainseq[pos])
              : MAKECOMPL(encseq->plainseq[
                           REVERSEPOS(encseq->totallength,pos)])
              )
            )
         ;
}

#define extractencodedchar(ENCSEQ,POS,RM)\
        getencodedchar(ENCSEQ,POS,RM)
#else
GtUchar getencodedchar(const Encodedsequence *encseq,Seqpos pos,
                     Readmode readmode);
GtUchar extractencodedchar(const Encodedsequence *encseq,
                           Seqpos pos,
                           Readmode readmode);
#endif

#ifdef INLINEDENCSEQ
#define getencodedcharnospecial(ENCSEQ,POS,RM)\
        getencodedchar(ENCSEQ,POS,RM)
#else
GtUchar getencodedcharnospecial(const Encodedsequence *encseq,
                              Seqpos pos,
                              Readmode readmode);
#endif

#ifdef INLINEDENCSEQ
#define sequentialgetencodedchar(ENCSEQ,ENCSEQSTATE,POS,READMODE)\
        getencodedchar(ENCSEQ,POS,READMODE)
#else
GtUchar sequentialgetencodedchar(const Encodedsequence *encseq,
                               Encodedsequencescanstate *esr,
                               Seqpos pos,
                               Readmode readmode);
#endif

void extract2bitenc(bool fwd,
                    EndofTwobitencoding *ptbe,
                    const Encodedsequence *encseq,
                    Encodedsequencescanstate *esr,
                    Seqpos startpos);

typedef struct
{
  unsigned int common;
  bool leftspecial, rightspecial;
  Seqpos finaldepth;
} GtCommonunits;

int compareTwobitencodings(bool fwd,
                           bool complement,
                           GtCommonunits *commonunits,
                           const EndofTwobitencoding *ptbe1,
                           const EndofTwobitencoding *ptbe2);

uint64_t detencseqofsatviatables(int kind,
                                 Seqpos totallength,
                                 unsigned long numofdbfiles,
                                 unsigned long lengthofdbfilenames,
                                 Seqpos specialranges,
                                 unsigned int numofchars);

void plainseq2bytecode(GtUchar *bytecode,const GtUchar *seq,unsigned long len);

void sequence2bytecode(GtUchar *dest,const Encodedsequence *encseq,
                       Seqpos startindex,Seqpos len);

int flushencseqfile(const GtStr *indexname,Encodedsequence *encseq,GtError*);

Encodedsequencescanstate *newEncodedsequencescanstate(void);

void encodedsequence_free(Encodedsequence **encseqptr);

void initEncodedsequencescanstate(Encodedsequencescanstate *esr,
                                  const Encodedsequence *encseq,
                                  Readmode readmode,
                                  Seqpos startpos);

void initEncodedsequencescanstategeneric(Encodedsequencescanstate *esr,
                                         const Encodedsequence *encseq,
                                         bool moveforward,
                                         Seqpos startpos);

void freeEncodedsequencescanstate(Encodedsequencescanstate **esr);

/*@null@*/ Encodedsequence *files2encodedsequence(
                                    bool withrange,
                                    const GtStrArray *filenametab,
                                    const Filelengthvalues *filelengthtab,
                                    bool plainformat,
                                    Seqpos totallength,
                                    unsigned long numofsequences,
                                    const Seqpos *specialrangestab,
                                    const GtAlphabet *alphabet,
                                    const char *str_sat,
                                    unsigned long *characterdistribution,
                                    const Specialcharinfo *specialcharinfo,
                                    Verboseinfo *verboseinfo,
                                    GtError *err);

/*@null@*/ Encodedsequence *mapencodedsequence(bool withrange,
                                               const GtStr *indexname,
                                               bool withesqtab,
                                               bool withdestab,
                                               bool withsdstab,
                                               bool withssptab,
                                               Verboseinfo *verboseinfo,
                                               GtError *err);

void checkallsequencedescriptions(const Encodedsequence *encseq);

Encodedsequence *plain2encodedsequence(bool withrange,
                                       const GtUchar *seq1,
                                       Seqpos len1,
                                       const GtUchar *seq2,
                                       unsigned long len2,
                                       const GtAlphabet *alpha,
                                       Verboseinfo *verboseinfo);

Specialrangeiterator *newspecialrangeiterator(const Encodedsequence *encseq,
                                              bool moveforward);

bool hasspecialranges(const Encodedsequence *encseq);

bool hasfastspecialrangeenumerator(const Encodedsequence *encseq);

bool possibletocmpbitwise(const Encodedsequence *encseq);

bool nextspecialrangeiterator(Sequencerange *range,Specialrangeiterator *sri);

void freespecialrangeiterator(Specialrangeiterator **sri);

/*@null@*/ const char *encseqaccessname(const Encodedsequence *encseq);

void encseqextract(GtUchar *buffer,
                   const Encodedsequence *encseq,
                   Seqpos frompos,
                   Seqpos topos);

Codetype extractprefixcode(unsigned int *unitsnotspecial,
                           const Encodedsequence *encseq,
                           const Codetype *filltable,
                           Readmode readmode,
                           Encodedsequencescanstate *esr,
                           const Codetype **multimappower,
                           Seqpos frompos,
                           unsigned int prefixlength);

int comparewithonespecial(bool *leftspecial,
                          bool *rightspecial,
                          const Encodedsequence *encseq,
                          bool fwd,
                          bool complement,
                          Seqpos pos1,
                          Seqpos pos2,
                          Seqpos depth,
                          Seqpos maxdepth);

int compareEncseqsequences(GtCommonunits *commonunits,
                           const Encodedsequence *encseq,
                           bool fwd,
                           bool complement,
                           Encodedsequencescanstate *esr1,
                           Encodedsequencescanstate *esr2,
                           Seqpos pos1,
                           Seqpos pos2,
                           Seqpos depth);

int compareEncseqsequencesmaxdepth(GtCommonunits *commonunits,
                                   const Encodedsequence *encseq,
                                   bool fwd,
                                   bool complement,
                                   Encodedsequencescanstate *esr1,
                                   Encodedsequencescanstate *esr2,
                                   Seqpos pos1,
                                   Seqpos pos2,
                                   Seqpos depth,
                                   Seqpos maxdepth);

/* some check functions called in test-encseq.c */

int multicharactercompare(const Encodedsequence *encseq,
                          bool fwd,
                          bool complement,
                          Encodedsequencescanstate *esr1,
                          Seqpos pos1,
                          Encodedsequencescanstate *esr2,
                          Seqpos pos2);

void checkextractunitatpos(const Encodedsequence *encseq,
                           bool fwd,bool complement);

void checkextractspecialbits(const Encodedsequence *encseq,bool fwd);

void multicharactercompare_withtest(const Encodedsequence *encseq,
                                    bool fwd,
                                    bool complement,
                                    Encodedsequencescanstate *esr1,
                                    Seqpos pos1,
                                    Encodedsequencescanstate *esr2,
                                    Seqpos pos2);

void showsequenceatstartpos(FILE *fp,
                            bool fwd,
                            bool complement,
                            const Encodedsequence *encseq,
                            Seqpos startpos);

bool containsspecial(const Encodedsequence *encseq,
                     bool moveforward,
                     Encodedsequencescanstate *esrspace,
                     Seqpos startpos,
                     Seqpos len);

int getsatforcevalue(const char *str,GtError *err);

/* check if the marked positions are correct */

void checkmarkpos(const Encodedsequence *encseq);

/* for a array of recordseparator, obtain the sequence
 * number from the given position */

unsigned long getrecordnumSeqpos(const Seqpos *recordseps,
                                 unsigned long numofrecords,
                                 Seqpos totalwidth,
                                 Seqpos position);

/* for a given Encodedsequence mapped with withssptab=true, obtain the sequence
 * number from the given position */

unsigned long getencseqfrompos2seqnum(const Encodedsequence *encseq,
                                      Seqpos position);

/* for a given Encodedsequence and a sequencenumber, fill the Seqinfo
 * structure */

void getencseqSeqinfo(Seqinfo *seqinfo,
                      const Encodedsequence *encseq,
                      unsigned long seqnum);

/* for a give  Encodedsequence and a sequencenumber return a pointer to
   the description of the sequence and store the length of the description
   in desclen */

const char *retrievesequencedescription(unsigned long *desclen,
                                        const Encodedsequence *encseq,
                                        unsigned long seqnum);

/* here are some functions to extract the different components of the
 * specialcharinfo included in encseq */

Seqpos getencseqspecialcharacters(const Encodedsequence *encseq);

Seqpos getencseqspecialranges(const Encodedsequence *encseq);

Seqpos getencseqrealspecialranges(const Encodedsequence *encseq);

Seqpos getencseqlengthofspecialprefix(const Encodedsequence *encseq);

Seqpos getencseqlengthofspecialsuffix(const Encodedsequence *encseq);

/* In case an Encodedsequence is not mapped, we still need to obtain the
   Specialcharainfo. This is done by the following function */

int readSpecialcharinfo(Specialcharinfo *specialcharinfo,
                        const GtStr *indexname,GtError *err);

/* some functions to obtain some components from the Alphabet pointed to
   by encseq->alpha */

unsigned int getencseqAlphabetnumofchars(const Encodedsequence *encseq);

const GtUchar *getencseqAlphabetsymbolmap(const Encodedsequence *encseq);

const GtAlphabet *getencseqAlphabet(const Encodedsequence *encseq);

const GtUchar *getencseqAlphabetcharacters(const Encodedsequence *encseq);

GtUchar getencseqAlphabetwildcardshow(const Encodedsequence *encseq);

/* Obtain the filenametable and the filelengthtable from the
   Encodedsequence */

const GtStrArray *getencseqfilenametab(const Encodedsequence *encseq);

unsigned long getencseqcharactercount(const Encodedsequence *encseq,GtUchar cc);

/* some function to remove reference from an Encodedsequence to prevent that
   the referenced alphabet or filenametab are freed */

void removealpharef(Encodedsequence *encseq);

void removefilenametabref(Encodedsequence *encseq);

void showgetencodedcharcounters(void);

void gt_showsequencefeatures(Verboseinfo *verboseinfo,
                             const Encodedsequence *encseq,bool withfilenames);

unsigned long determinelengthofdbfilenames(const GtStrArray *filenametab);

int gt_inputfiles2sequencekeyvalues(
        const GtStr *indexname,
        Seqpos *totallength,
        Specialcharinfo *specialcharinfo,
        unsigned int forcetable,
        Seqpos *specialrangestab,
        const GtStrArray *filenametab,
        Filelengthvalues **filelengthtab,
        const GtAlphabet *alpha,
        bool plainformat,
        bool outdestab,
        bool outsdstab,
        bool outkystab,
        bool outkyssort,
        unsigned long *characterdistribution,
        bool outssptab,
        ArraySeqpos *sequenceseppos,
        Verboseinfo *verboseinfo,
        GtError *err);

FILE *opendestabfile(const GtStr *indexname,const char *mode,GtError *err);

FILE *openssptabfile(const GtStr *indexname,const char *mode,GtError *err);

int comparetwosuffixes(const Encodedsequence *encseq,
                       Readmode readmode,
                       Seqpos *maxlcp,
                       bool specialsareequal,
                       bool specialsareequalatdepth0,
                       Seqpos maxdepth,
                       Seqpos start1,
                       Seqpos start2,
                       Encodedsequencescanstate *esr1,
                       Encodedsequencescanstate *esr2);

int comparetwostrings(const Encodedsequence *encseq,
                      bool fwd,
                      bool complement,
                      Seqpos *maxcommon,
                      Seqpos pos1,
                      Seqpos pos2,
                      Seqpos maxdepth);

int comparetwostringsgeneric(const Encodedsequence *encseq,
                             bool fwd,
                             bool complement,
                             Seqpos *maxcommon,
                             Seqpos pos1,
                             Seqpos pos2,
                             Seqpos depth,
                             Seqpos maxdepth);

#endif
