/*
  Copyright (c) 2003-2009 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2000-2004 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
  Copyright (c) 2003-2008 Center for Bioinformatics, University of Hamburg

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

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <math.h>
#include "core/chardef.h"
#include "core/codon_api.h"
#include "core/minmax.h"
#include "core/translator_api.h"
#include "core/unused_api.h"
#include "gth/gthalignment.h"

#define OUTCHAR(C)         gt_file_xfputc((C),outfp)

#define OUTSTRINGNUM(N,S)  gt_file_xprintf(outfp,"%*lu",(int) (N), (S))

#define CHECK_FOR_IDENTICAL_SUBSTRING_OF_LENGTH_0(EOPPTR)\
        if ((EOPPTR) == 0)\
        {\
          fprintf(stderr,"identical substring of length 0 not allowed\n");\
          exit(EXIT_FAILURE);\
        }

#define GENOMICDNAPTR_STILL_VALID\
        gt_assert(genomicdnaptr < aftergenomicdnaline)

#define GENOMICPROTEINPTR_STILL_VALID\
        gt_assert(genomicproteinptr < aftergenomicproteinline)

#define REFERENCEPROTEINPTR_STILL_VALID\
        gt_assert(referenceproteinptr < afterreferenceproteinline)

#define EQUAL_AMINO_ACID_CHAR  '|'
#define POSITIVE_SCORE_CHAR    '+'
#define ZERO_SCORE_CHAR        '.'
#define NEGATIVE_SCORE_CHAR    ' '

/*EE
  This file implements functions to formatting and showing an alignment.
  The alignment is formatted in two steps. In the first step we
  take the list of edit operations and compute the two lines of the alignment
  with inserted abstract gap symbols for insertions and deletions.
  The characters are however stored as integers in the range [0..alphasize-1].
  In the second phase the two lines are formatted such that they fit
  on the given linewidth. Additionally, the position offset for the
  line are shown on the right. Moreover, if the given linewidth contains
  at least one mismatch or indel, this is shown as an extra line
  marking the corresponding columns by the symbol
  \texttt{\symbol{34}\symbol{33}\symbol{34}}.
*/

#define NUMWIDTH  12  /* width of position right of alignment */

static unsigned long calctotallengthofintron(unsigned long *numofeopsinintron,
                                    Editoperation *intronstart,
                                    Editoperation *alignmentstart)
{
  Editoperation *eopptr;
  unsigned long totalintronlength = 0;
  bool breakforloop = false;

  *numofeopsinintron = 0;

  for (eopptr = intronstart; eopptr >= alignmentstart; eopptr--) {
    switch (*eopptr) {
      case MISMATCHEOP:
      case DELETIONEOP:
      case INSERTIONEOP:
        breakforloop = true;
        break;
      default:
        switch (*eopptr & ~MAXIDENTICALLENGTH) {
          case 0:                               /* match */
            breakforloop = true;
            break;
          case DELETIONEOP:                     /* intron */
            (*numofeopsinintron)++;
            totalintronlength += *eopptr & MAXIDENTICALLENGTH;
            break;
        }
    }
    if (breakforloop)
      break;
  }

  return totalintronlength;
}

/*
  The following function implements the first step mentioned above, i.e.\
  it fills the \texttt{firstlinecompare} and the \texttt{secondlinecompare}
  of the alignment.
*/

static unsigned long fillthelines(GtUchar *firstline,GtUchar *secondline,
                                  const GtUchar *useq, const GtUchar *vseq,
                                  Editoperation *alignment,
                                  unsigned long lenalg, unsigned long linewidth,
                                  unsigned long showintronmaxlen,
                                  GtArrayShortIntronInfo *shortintroninfo)
{
  Editoperation *eopptr;
  unsigned long l, intronlength, totalintronlength, numofeopsinintron,
                restlength, end, completeshortintronlen = 0, i = 0, j = 0;
  GtUchar *fptr = firstline, *sptr = secondline;
  ShortIntronInfo oneshortintron;

  /* first phase: construct first line of alignment */
  for (eopptr = alignment+lenalg-1; eopptr >= alignment; eopptr--) {
    switch (*eopptr)
    {
      case MISMATCHEOP:
      case DELETIONEOP:
        *fptr++ = useq[i];
        i++;
        break;
      case INSERTIONEOP:
        *fptr++ = (GtUchar) ABSTRACTGAPSYMBOL;
        break;
      default:
        switch (*eopptr & ~MAXIDENTICALLENGTH) {
        case 0:                               /* match */
          CHECK_FOR_IDENTICAL_SUBSTRING_OF_LENGTH_0(*eopptr);
          for (l = 0; l < (unsigned long) *eopptr; l++) {
            *fptr++ = useq[i];
            i++;
          }
          break;
        case DELETIONEOP:                     /* intron */
          if (showintronmaxlen > 0) {
            /* compute total intron length */
            totalintronlength = calctotallengthofintron(&numofeopsinintron,
                                                        eopptr, alignment);

            /* compute the rest length which is necessary to fill the current
               line */
            restlength = (unsigned long) (fptr - firstline) % linewidth;
            if (restlength != 0) {
              /* restlength is strictly smaller than linewidth" */
              gt_assert(restlength < linewidth);
              restlength = linewidth - restlength;
            }

            /* check if introducing a ``short'' intron here is necessary */
            if (totalintronlength >=
                showintronmaxlen + restlength + 2 * linewidth ) {
               /* introduce ``short'' intron */

               /* write the beginning of the intron in the conventional way */
               for (l = 0; l < restlength + linewidth; l++) {
                 *fptr++ = useq[i];
                 i++;
               }

               /* store short intron information if necessary */
               if (shortintroninfo) {
                 oneshortintron.start  = (unsigned long) (fptr - firstline +
                                                 completeshortintronlen);
                 end = oneshortintron.start +
                       ((totalintronlength - restlength - 2 * linewidth)
                        / linewidth) * linewidth - 1;
                 oneshortintron.length  = end - oneshortintron.start + 1;
                 completeshortintronlen += oneshortintron.length;

                 GT_STOREINARRAY(shortintroninfo, ShortIntronInfo,
                                 8, oneshortintron);
               }

               /* add length of short intron to i */
               i += ((totalintronlength - restlength - 2 * linewidth)
                     / linewidth) * linewidth;

               /* write the end  of the intron in the conventional way */
               for (l = 0;
                    l < linewidth +
                        ((totalintronlength - restlength - 2 * linewidth) %
                         linewidth);
                    l++) {
                 *fptr++ = useq[i];
                 i++;
               }

               /* place eopptr to the last multi editoperation of this intron
                  because the intron has already been completely processed */
               eopptr -= (numofeopsinintron - 1);

               /* intron is already processed, therefore continue */
               continue;
            }
          }

          intronlength = *eopptr & MAXIDENTICALLENGTH;
          for (l = 0; l < intronlength; l++) {
            *fptr++ = useq[i];
            i++;
          }
        break;
        default: gt_assert(0);
        }
    }
  }
  /* first phase: construct second line of alignment */
  for (eopptr = alignment+lenalg-1; eopptr >= alignment; eopptr--) {
    switch (*eopptr) {
      case MISMATCHEOP:
      case INSERTIONEOP:
        *sptr++ = vseq[j];
        j++;
        break;
      case DELETIONEOP:
        *sptr++ = (GtUchar) ABSTRACTGAPSYMBOL;
        break;
      default:
        switch (*eopptr & ~MAXIDENTICALLENGTH) {
        case 0:                               /* match */
          CHECK_FOR_IDENTICAL_SUBSTRING_OF_LENGTH_0(*eopptr);
          for (l = 0; l < (unsigned long) *eopptr; l++) {
            *sptr++ = vseq[j];
            j++;
          }
          break;
        case DELETIONEOP:                     /* intron */
          if (showintronmaxlen > 0) {
            /* compute total intron length */
            totalintronlength = calctotallengthofintron(&numofeopsinintron,
                                                        eopptr, alignment);

            /* compute the rest length which is necessary to fill the current
               line */
            restlength = (unsigned long) (sptr - secondline) % linewidth;
            if (restlength != 0) {
              /* restlength is strictly smaller than linewidth */
              gt_assert(restlength < linewidth);
              restlength = linewidth - restlength;
            }

            /* check if introducing a ``short'' intron here is necessary */
            if (totalintronlength >=
                showintronmaxlen + restlength + 2 * linewidth ) {
               /* write the rest of the intron in the conventional way */
               for (l = 0;
                    l < restlength + 2 * linewidth +
                        ((totalintronlength - restlength - 2 * linewidth) %
                         linewidth);
                    l++) {
                 *sptr++ = (GtUchar) ABSTRACTINTRONSYMBOL;
               }

               /* place eopptr to the last multi editoperation of this intron
                  because the intron has already been completely processed */
               eopptr -= (numofeopsinintron - 1);

               /* intron is already processed, therefore continue */
               continue;
            }
          }

          intronlength = *eopptr & MAXIDENTICALLENGTH;
          for (l = 0; l < intronlength; l++) {
            *sptr++ = (GtUchar) ABSTRACTINTRONSYMBOL;
          }
          break;
        default: gt_assert(0);
        }
    }
  }

  return (unsigned long) (sptr-secondline);
}

static void match_mismatch_genomicdnaline(GtUchar **genomicdnaptr,
                                          bool
                                          *processing_intron_with_1_base_left,
                                          bool
                                          *processing_intron_with_2_bases_left,
                                          const GtUchar *genseqorig,
                                          unsigned long *genseqindex)
{
  if (*processing_intron_with_2_bases_left) {
    /* this means we are after an intron with 2 bases left
       this bases have already been shown, therefore we only show 1 more */
    *processing_intron_with_2_bases_left = false;
    **genomicdnaptr = genseqorig[(*genseqindex)++];
    (*genomicdnaptr)++;
  }
  else {
    if (*processing_intron_with_1_base_left) {
      /* this means we are after an intron with 1 base left
         this base has already been shown, therefore we only show 2 more */
      *processing_intron_with_1_base_left = false;
    }
    else {
      **genomicdnaptr = genseqorig[(*genseqindex)++];
      (*genomicdnaptr)++;
    }
    **genomicdnaptr = genseqorig[(*genseqindex)++];
    (*genomicdnaptr)++;
    **genomicdnaptr = genseqorig[(*genseqindex)++];
    (*genomicdnaptr)++;
  }
}

/*
  The following function is used to construct the genomic DNA line of
  protein alignments.
*/

static unsigned long construct_genomic_dna_line(GtUchar *genomicdnaline,
                                       GT_UNUSED unsigned long
                                       lengthofgenomicdnaline,
                                       const GtUchar *genseqorig,
                                       Editoperation *alignment,
                                       unsigned long lenalg)
{
  Editoperation *eopptr;
  Eoptype eoptype;
  unsigned long eoplength, l,
       i = 0;
  GtUchar *genomicdnaptr = genomicdnaline;
#ifndef NDEBUG
  GtUchar *aftergenomicdnaline = genomicdnaline + lengthofgenomicdnaline;
#endif
  bool processing_intron_with_1_base_left  = false,
       processing_intron_with_2_bases_left = false;

  for (eopptr = alignment + lenalg - 1; eopptr >= alignment; eopptr--) {
    eoptype   = gt_editoperation_type(*eopptr, true);
    eoplength = gt_editoperation_length(*eopptr, true);

    /* we are not processing two intron types at the same time */
    gt_assert(!(processing_intron_with_1_base_left &&
             processing_intron_with_2_bases_left));

    switch (eoptype) {
      case EOP_TYPE_DELETION:
        gt_assert(!processing_intron_with_1_base_left);
        gt_assert(!processing_intron_with_2_bases_left);
        /*@fallthrough@*/
      case EOP_TYPE_MISMATCH:
        gt_assert(eoplength == 1);
        match_mismatch_genomicdnaline(&genomicdnaptr,
                                      &processing_intron_with_1_base_left,
                                      &processing_intron_with_2_bases_left,
                                      genseqorig,
                                      &i);
        break;
      case EOP_TYPE_INSERTION:
        gt_assert(eoplength == 1);
        gt_assert(!processing_intron_with_1_base_left);
        gt_assert(!processing_intron_with_2_bases_left);
        GENOMICDNAPTR_STILL_VALID;
        *genomicdnaptr++ = (GtUchar) ABSTRACTGAPSYMBOL;
        GENOMICDNAPTR_STILL_VALID;
        *genomicdnaptr++ = (GtUchar) ABSTRACTGAPSYMBOL;
        GENOMICDNAPTR_STILL_VALID;
        *genomicdnaptr++ = (GtUchar) ABSTRACTGAPSYMBOL;
        break;
      case EOP_TYPE_MISMATCH_WITH_1_GAP:
      case EOP_TYPE_DELETION_WITH_1_GAP:
        gt_assert(eoplength == 1);
        gt_assert(!processing_intron_with_1_base_left);
        gt_assert(!processing_intron_with_2_bases_left);
        GENOMICDNAPTR_STILL_VALID;
        *genomicdnaptr++ = genseqorig[i++];
        GENOMICDNAPTR_STILL_VALID;
        *genomicdnaptr++ = (GtUchar) ABSTRACTGAPSYMBOL;
        GENOMICDNAPTR_STILL_VALID;
        *genomicdnaptr++ = genseqorig[i++];
        break;
      case EOP_TYPE_MISMATCH_WITH_2_GAPS:
      case EOP_TYPE_DELETION_WITH_2_GAPS:
        gt_assert(eoplength == 1);
        gt_assert(!processing_intron_with_1_base_left);
        gt_assert(!processing_intron_with_2_bases_left);
        GENOMICDNAPTR_STILL_VALID;
        *genomicdnaptr++ = (GtUchar) ABSTRACTGAPSYMBOL;
        GENOMICDNAPTR_STILL_VALID;
        *genomicdnaptr++ = genseqorig[i++];
        GENOMICDNAPTR_STILL_VALID;
        *genomicdnaptr++ = (GtUchar) ABSTRACTGAPSYMBOL;
        break;
      case EOP_TYPE_MATCH:
        for (l = 0; l < eoplength; l++) {
          match_mismatch_genomicdnaline(&genomicdnaptr,
                                        &processing_intron_with_1_base_left,
                                        &processing_intron_with_2_bases_left,
                                        genseqorig, &i);
        }
        break;
      case EOP_TYPE_INTRON_WITH_2_BASES_LEFT:
        if (!processing_intron_with_2_bases_left) {
          processing_intron_with_2_bases_left = true;
          GENOMICDNAPTR_STILL_VALID;
          *genomicdnaptr++ = genseqorig[i++];
          GENOMICDNAPTR_STILL_VALID;
          *genomicdnaptr++ = genseqorig[i++];
        }
        /* skip the next case statement and process this intron */
        goto process_intron;
      case EOP_TYPE_INTRON_WITH_1_BASE_LEFT:
        if (!processing_intron_with_1_base_left) {
          processing_intron_with_1_base_left = true;
          GENOMICDNAPTR_STILL_VALID;
          *genomicdnaptr++ = genseqorig[i++];
        }
        /*@fallthrough@*/
      case EOP_TYPE_INTRON:
      process_intron:
        for (l = 0; l < eoplength; l++) {
          GENOMICDNAPTR_STILL_VALID;
          *genomicdnaptr++ = genseqorig[i++];
        }
        break;
      default: gt_assert(0);
    }
  }

  gt_assert(!processing_intron_with_1_base_left);
  gt_assert(!processing_intron_with_2_bases_left);

  return (unsigned long) (genomicdnaptr - genomicdnaline);
}

static void match_mismatch_genomicproteinline(GtUchar **genomicproteinptr,
                                          bool
                                          *processing_intron_with_1_base_left,
                                          bool
                                          *processing_intron_with_2_bases_left,
                                          const GtUchar *genseqorig,
                                          GtUchar *first_base_left,
                                          GtUchar *second_base_left,
                                          GtUchar *dummyptr,
                                          unsigned long *genseqindex,
                                          unsigned long translationschemenumber)
{
  GtUchar dna[GT_CODON_LENGTH];
  GtTransTable *transtable;
  char codon;
  GT_UNUSED int rval;

  transtable = gt_trans_table_new(translationschemenumber, NULL);
  /* XXX: the validity of the translation table has to be checked before */
  gt_assert(transtable);

  if (*processing_intron_with_2_bases_left) {
    *processing_intron_with_2_bases_left = false;
    gt_assert(*first_base_left != UNDEFCHAR);
    gt_assert(*second_base_left != UNDEFCHAR);
    gt_assert(dummyptr != NULL);
    dna[0] = *first_base_left;
    dna[1] = *second_base_left;
    dna[2] = genseqorig[(*genseqindex)++];
    *first_base_left  = (GtUchar) UNDEFCHAR;
    *second_base_left = (GtUchar) UNDEFCHAR;

    rval = gt_trans_table_translate_codon(transtable, dna[0], dna[1], dna[2],
                                          &codon, NULL);
    /* since the sequence has been preprocessed before, the codon translation
       should not fail */
    gt_assert(!rval);

    /* set dummy pointer */
    *dummyptr = codon;
    /* show second blank here */
    **genomicproteinptr = (GtUchar) ' ';
    (*genomicproteinptr)++;
  }
  else {
    if (*processing_intron_with_1_base_left) {
      *processing_intron_with_1_base_left = false;
      gt_assert(*first_base_left != UNDEFCHAR);
      dna[0] = *first_base_left;
      *first_base_left = (GtUchar) UNDEFCHAR;
    }
    else {
      dna[0] = genseqorig[(*genseqindex)++];
    }
    dna[1] = genseqorig[(*genseqindex)++];
    dna[2] = genseqorig[(*genseqindex)++];

    rval = gt_trans_table_translate_codon(transtable, dna[0], dna[1], dna[2],
                                          &codon, NULL);
    /* since the sequence has been preprocessed before, the codon translation
       should not fail */
    gt_assert(!rval);

    (**genomicproteinptr) = (GtUchar) ' ';
    (*genomicproteinptr)++;
    (**genomicproteinptr) = codon;
    (*genomicproteinptr)++;
    (**genomicproteinptr) = (GtUchar) ' ';
    (*genomicproteinptr)++;
  }

  gt_trans_table_delete(transtable);
}

/*
  The following function is used to construct the genomic protein line of
  protein alignments.
*/

static unsigned long construct_genomic_protein_line(GtUchar *genomicproteinline,
                                                    GT_UNUSED unsigned long
                                                    lengthofgenomicproteinline,
                                                    const GtUchar *genseqorig,
                                                    Editoperation *alignment,
                                                    unsigned long lenalg,
                                                    unsigned long
                                                    translationschemenumber)
{
  Editoperation *eopptr;
  Eoptype eoptype;
  unsigned long eoplength, l,
       i = 0;
  GtUchar first_base_left    = (GtUchar) UNDEFCHAR,
        second_base_left   = (GtUchar) UNDEFCHAR,
        *dummyptr          = NULL,
        *genomicproteinptr = genomicproteinline;
#ifndef NDEBUG
  GtUchar *aftergenomicproteinline   = genomicproteinline +
                                     lengthofgenomicproteinline;
#endif
  bool processing_intron_with_1_base_left  = false,
       processing_intron_with_2_bases_left = false;

  for (eopptr = alignment + lenalg - 1; eopptr >= alignment; eopptr--) {
    eoptype   = gt_editoperation_type(*eopptr, true);
    eoplength = gt_editoperation_length(*eopptr, true);

    gt_assert(!(processing_intron_with_1_base_left &&
             processing_intron_with_2_bases_left));

    switch (eoptype) {
      case EOP_TYPE_MISMATCH:
        gt_assert(eoplength == 1);
        match_mismatch_genomicproteinline(&genomicproteinptr,
                                          &processing_intron_with_1_base_left,
                                          &processing_intron_with_2_bases_left,
                                          genseqorig,
                                          &first_base_left,
                                          &second_base_left,
                                          dummyptr,
                                          &i,
                                          translationschemenumber);
        break;
      case EOP_TYPE_DELETION:
        /* skip three characters in genomic dna (fall through case!) */
        i++;
        /*@fallthrough@*/
      case EOP_TYPE_MISMATCH_WITH_1_GAP:
      case EOP_TYPE_DELETION_WITH_1_GAP:
        /* skip two characters in genomic dna (fall through case!) */
        i++;
        /*@fallthrough@*/
      case EOP_TYPE_MISMATCH_WITH_2_GAPS:
      case EOP_TYPE_DELETION_WITH_2_GAPS:
        /* skip one characters in genomic dna (fall through case!) */
        i++;
        /*@fallthrough@*/
      case EOP_TYPE_INSERTION:
        gt_assert(eoplength == 1);
        gt_assert(!processing_intron_with_1_base_left);
        gt_assert(!processing_intron_with_2_bases_left);
        GENOMICPROTEINPTR_STILL_VALID;
        *genomicproteinptr++ = (GtUchar) ' ';
        GENOMICPROTEINPTR_STILL_VALID;
        *genomicproteinptr++ = (GtUchar) ' ';
        GENOMICPROTEINPTR_STILL_VALID;
        *genomicproteinptr++ = (GtUchar) ' ';
        break;
      case EOP_TYPE_MATCH:
        for (l = 0; l < eoplength; l++) {
          match_mismatch_genomicproteinline(&genomicproteinptr,
                                            &processing_intron_with_1_base_left,
                                           &processing_intron_with_2_bases_left,
                                            genseqorig,
                                            &first_base_left,
                                            &second_base_left,
                                            dummyptr,
                                            &i,
                                            translationschemenumber);
        }
        break;
      case EOP_TYPE_INTRON_WITH_2_BASES_LEFT:
        if (!processing_intron_with_2_bases_left) {
          processing_intron_with_2_bases_left = true;
          /* save the first two bases */
          first_base_left  = genseqorig[i++];
          second_base_left = genseqorig[i++];
          /* print first blank here already */
          GENOMICPROTEINPTR_STILL_VALID;
          *genomicproteinptr++ = (GtUchar) ' ';
          /* save dummy pointer */
          GENOMICPROTEINPTR_STILL_VALID;
          dummyptr = genomicproteinptr++;
        }
        /* skip the next case statement and process this intron */
        goto process_intron;
      case EOP_TYPE_INTRON_WITH_1_BASE_LEFT:
        if (!processing_intron_with_1_base_left)
        {
          processing_intron_with_1_base_left = true;
          first_base_left = genseqorig[i++];
        }
        /*@fallthrough@*/
      case EOP_TYPE_INTRON:
      process_intron:
        for (l = 0; l < eoplength; l++) {
          GENOMICPROTEINPTR_STILL_VALID;
          *genomicproteinptr++ = (GtUchar) ' ';
          i++;
        }
        break;
      default: gt_assert(0);
    }
  }

  gt_assert(!processing_intron_with_1_base_left);
  gt_assert(!processing_intron_with_2_bases_left);

  return (unsigned long) (genomicproteinptr - genomicproteinline);
}

static void match_mismatch_referenceproteinline(GtUchar **referenceproteinptr,
                                          bool
                                          *processing_intron_with_1_base_left,
                                          bool
                                          *processing_intron_with_2_bases_left,
                                          const GtUchar *refseqorig,
                                          unsigned long *refseqindex)
{
  if (*processing_intron_with_2_bases_left) {
    *processing_intron_with_2_bases_left = false;
    **referenceproteinptr = (GtUchar) ' ';
    (*referenceproteinptr)++;
  }
  else {
    if (*processing_intron_with_1_base_left) {
      *processing_intron_with_1_base_left = false;
    }
    else {
      **referenceproteinptr = (GtUchar) ' ';
      (*referenceproteinptr)++;
    }
    **referenceproteinptr = refseqorig[(*refseqindex)++];
    (*referenceproteinptr)++;
    **referenceproteinptr = (GtUchar) ' ';
    (*referenceproteinptr)++;
  }
}

/*
  The following function is used to construct the genomic protein line of
  protein alignments.
*/

static unsigned long construct_reference_protein_line(GtUchar
                                                      *referenceproteinline,
                                                      GT_UNUSED unsigned long
                                                   lengthofreferenceproteinline,
                                                      const GtUchar *refseqorig,
                                                      Editoperation *alignment,
                                                      unsigned long lenalg)
{
  Editoperation *eopptr;
  Eoptype eoptype;
  unsigned long eoplength, l,
       i = 0;
  GtUchar *referenceproteinptr = referenceproteinline;
#ifndef NDEBUG
  GtUchar *afterreferenceproteinline = referenceproteinline +
                                     lengthofreferenceproteinline;
#endif
  bool processing_intron_with_1_base_left  = false,
       processing_intron_with_2_bases_left = false;

  for (eopptr = alignment + lenalg - 1; eopptr >= alignment; eopptr--) {
    eoptype   = gt_editoperation_type(*eopptr, true);
    eoplength = gt_editoperation_length(*eopptr, true);

    gt_assert(!(processing_intron_with_1_base_left &&
             processing_intron_with_2_bases_left));

    switch (eoptype) {
      case EOP_TYPE_MISMATCH_WITH_1_GAP:
      case EOP_TYPE_MISMATCH_WITH_2_GAPS:
      case EOP_TYPE_INSERTION:
        gt_assert(!processing_intron_with_1_base_left);
        gt_assert(!processing_intron_with_2_bases_left);
        /*@fallthrough@*/
      case EOP_TYPE_MISMATCH:
        gt_assert(eoplength == 1);
        match_mismatch_referenceproteinline(&referenceproteinptr,
                                           &processing_intron_with_1_base_left,
                                           &processing_intron_with_2_bases_left,
                                           refseqorig,
                                           &i);
        break;
      case EOP_TYPE_DELETION:
      case EOP_TYPE_DELETION_WITH_1_GAP:
      case EOP_TYPE_DELETION_WITH_2_GAPS:
        gt_assert(eoplength == 1);
        gt_assert(!processing_intron_with_1_base_left);
        gt_assert(!processing_intron_with_2_bases_left);
        REFERENCEPROTEINPTR_STILL_VALID;
        *referenceproteinptr++ = (GtUchar) ' ';
        REFERENCEPROTEINPTR_STILL_VALID;
        *referenceproteinptr++ = (GtUchar) ABSTRACTGAPSYMBOL;
        REFERENCEPROTEINPTR_STILL_VALID;
        *referenceproteinptr++ = (GtUchar) ' ';
        break;
      case EOP_TYPE_MATCH:
        for (l = 0; l < eoplength; l++)
        {
          match_mismatch_referenceproteinline(&referenceproteinptr,
                                           &processing_intron_with_1_base_left,
                                           &processing_intron_with_2_bases_left,
                                           refseqorig,
                                           &i);
        }
        break;
      case EOP_TYPE_INTRON_WITH_2_BASES_LEFT:
        if (!processing_intron_with_2_bases_left)
        {
          processing_intron_with_2_bases_left = true;
          REFERENCEPROTEINPTR_STILL_VALID;
          *referenceproteinptr++ = (GtUchar) ' ';
          REFERENCEPROTEINPTR_STILL_VALID;
          *referenceproteinptr++ = refseqorig[i++];
        }
        /* skip the next case statement and process this intron */
        goto process_intron;
      case EOP_TYPE_INTRON_WITH_1_BASE_LEFT:
        if (!processing_intron_with_1_base_left)
        {
          processing_intron_with_1_base_left = true;
          REFERENCEPROTEINPTR_STILL_VALID;
          *referenceproteinptr++ = (GtUchar) ' ';
        }
        /*@fallthrough@*/
      case EOP_TYPE_INTRON:
      process_intron:
        for (l = 0; l < eoplength; l++)
        {
          REFERENCEPROTEINPTR_STILL_VALID;
          *referenceproteinptr++ = (GtUchar) ABSTRACTINTRONSYMBOL;
        }
        break;
      default: gt_assert(0);
    }
  }

  gt_assert(!processing_intron_with_1_base_left);
  gt_assert(!processing_intron_with_2_bases_left);

  return (unsigned long) (referenceproteinptr - referenceproteinline);
}

static unsigned long filltheproteinlines(GtUchar *genomicdnaline,
                                GtUchar *genomicproteinline,
                                GtUchar *referenceproteinline,
                                unsigned long lengthofgenomicdnaline,
                                unsigned long lengthofgenomicproteinline,
                                unsigned long lengthofreferenceproteinline,
                                const GtUchar *genseqorig,
                                const GtUchar *refseqorig,
                                Editoperation *alignment,
                                unsigned long lenalg,
                                unsigned long translationschemenumber)
{
  unsigned long GT_UNUSED genomicdnalinelen,
       GT_UNUSED genomicproteinlinelen,
       referenceproteinlinelen;

  /* construct genomic DNA line of alignment */
  genomicdnalinelen = construct_genomic_dna_line(genomicdnaline,
                                                 lengthofgenomicdnaline,
                                                 genseqorig,
                                                 alignment,
                                                 lenalg);

  /* construct genomic protein line of alignment */
  genomicproteinlinelen =
    construct_genomic_protein_line(genomicproteinline,
                                   lengthofgenomicproteinline,
                                   genseqorig,
                                   alignment,
                                   lenalg,
                                   translationschemenumber);

  /* construct reference protein line of alignment */
  referenceproteinlinelen =
    construct_reference_protein_line(referenceproteinline,
                                     lengthofreferenceproteinline,
                                     refseqorig,
                                     alignment,
                                     lenalg);

  /* all three lines have the same length */
  gt_assert(genomicdnalinelen == genomicproteinlinelen);
  gt_assert(genomicproteinlinelen == referenceproteinlinelen);

  return referenceproteinlinelen;
}

static GtUchar implodewildcard(GtUchar c, bool wildcardimplosion,
                             GtAlphabet *alphabet)
{
  if (wildcardimplosion && (gt_alphabet_encode(alphabet, c) == WILDCARD)) {
    if (islower(c))
      return (GtUchar) tolower(gt_alphabet_wildcard_show(alphabet));
    else
      return (GtUchar) toupper(gt_alphabet_wildcard_show(alphabet));
  }
  else
    return c;
}

/*
  The following function formats a sequence \texttt{s} of length \texttt{len}.
  Each character code is shown as the corresponding character w.r.t.
  the given alphabet. Each abstract gap symbol is shown as a
  \texttt{\symbol{34}-\symbol{34}}. The output goes to the file pointer
  \texttt{outfp}.
*/

static void formatseqwithgaps(GtFile *outfp, GtUchar *sorig, unsigned long len,
                              unsigned long *insertioncount,
                              bool countproteininsertions, GtAlphabet *alphabet,
                              bool wildcardimplosion)
{
  unsigned long i;

  for (i=0; i<len; i++)
  {
    if ((i != 0) && (i%10 == 0))
    {
      OUTCHAR(' ');
    }
    if (sorig[i] == (GtUchar) ABSTRACTGAPSYMBOL)
    {
      OUTCHAR(CONCRETEGAPSYMBOL);
      if (insertioncount)
      {
        if (countproteininsertions)
        {
          *insertioncount += GT_CODON_LENGTH;
        }
        else
        {
          (*insertioncount)++;
        }
      }
    }
    else if (sorig[i] == (GtUchar) ABSTRACTINTRONSYMBOL)
    {
      OUTCHAR(CONCRETEINTRONSYMBOL);
      if (insertioncount)
      {
        (*insertioncount)++;
      }
    }
    else
    {
       gt_file_xfputc(implodewildcard(toupper(sorig[i]), wildcardimplosion,
                                      alphabet), outfp);
    }
  }
}

/*
  For the alignment lines \texttt{firstlinecompare} and
  \texttt{secondlinecompare}, both of length \texttt{len}, the following
  function shows a line
  with the symbol \texttt{\symbol{34}\symbol{33}\symbol{34}} in each column
  corresponding to a mismatch or an indel. The output goes to the file
  pointer \texttt{outfp}.
*/

#define CHECKORIGEQUAL(OKAY,A,B)\
        if ((A) == (B))\
        {\
          OKAY = true;\
        } else\
        {\
          if (islower((A)))\
          { \
            OKAY = ((A) == (GtUchar) tolower(B)) ? true : false;\
          } else\
          {\
            OKAY = ((A) == (GtUchar) toupper((B))) ? true : false;\
          }\
        }

static void showeditopline(GtFile *outfp,
                           GtUchar *firstlineorig,
                           GtUchar *secondlineorig,
                           unsigned long len,
                           GtAlphabet *alphabet)
{
  unsigned long i;
  GtUchar acompare, bcompare, aorig, borig;
  bool charequal, charinline = false;

  for (i=0; i<len; i++)
  {
    acompare = tolower(firstlineorig[i]);
    bcompare = tolower(secondlineorig[i]);
    if (acompare == bcompare && acompare != (GtUchar) CONCRETEGAPSYMBOL &&
        gt_alphabet_encode(alphabet, acompare) != (GtUchar) WILDCARD)
    {
      charinline = true;
      break;
    }
    aorig = firstlineorig[i];
    borig = secondlineorig[i];
    CHECKORIGEQUAL(charequal,aorig,borig);
    if (!charequal)
    {
      charinline = true;
      break;
    }
  }
  if (charinline)
  {
    for (i=0; i<len; i++)
    {
      if ((i != 0) && (i%10 == 0))
      {
        OUTCHAR(' ');
      }

      acompare = tolower(firstlineorig[i]);
      bcompare = tolower(secondlineorig[i]);
      if (acompare == bcompare && acompare != (GtUchar) CONCRETEGAPSYMBOL &&
          gt_alphabet_encode(alphabet, acompare) != (GtUchar) WILDCARD)
      {
        OUTCHAR('|');
      } else
      {
        aorig = firstlineorig[i];
        borig = secondlineorig[i];
        CHECKORIGEQUAL(charequal,aorig,borig);
        if (charequal)
        {
          OUTCHAR('=');
        } else
        {
          OUTCHAR(' ');
        }
      }
    }
    OUTCHAR('\n');
  }
}

static void showeditoplineprotein(GtFile *outfp,
                                  GtUchar *genomicproteinline,
                                  GtUchar *referenceproteinline,
                                  unsigned long len,
                                  GtScoreMatrix *score_matrix,
                                  GtAlphabet *score_matrix_alphabet)
{
  unsigned long i;
  GtUchar genchar, refchar;
  int score;

  for (i = 0; i < len; i++)
  {
    genchar = genomicproteinline[i];
    refchar = referenceproteinline[i];

    /* an empty refchar implies an emtpy genchar */
    gt_assert(!(refchar == ' ' && genchar != ' '));

    /* output blank if necessary */
    if ((i != 0) && (i%10 == 0)) {
      OUTCHAR(' ');
    }

    /* handle the ``empty case'' */
    if (genchar == ' ') {
      OUTCHAR(' ');
    }
    else if (genchar == refchar) {
      OUTCHAR(EQUAL_AMINO_ACID_CHAR);
    }
    else {
      GtUchar code1, code2;
      int wcidx;

      /* output character depending on amino acid substitution score of the two
         characters
         XXX: handle output of stop codons appropriately */

      /* determine score */
      code1 = gt_alphabet_encode(score_matrix_alphabet, genchar);
      code2 = gt_alphabet_encode(score_matrix_alphabet, refchar);
      wcidx = gt_alphabet_size(score_matrix_alphabet) - 1;
      score = gt_score_matrix_get_score(score_matrix,
                                        code1 == WILDCARD ? wcidx : code1,
                                        code2 == WILDCARD ? wcidx : code2);

      /* output corresponding character depending on score */
      if (score > 0) {
        OUTCHAR(POSITIVE_SCORE_CHAR);
      }
      else if (score == 0) {
        OUTCHAR(ZERO_SCORE_CHAR);
      }
      else { /* score < 0 */
        OUTCHAR(NEGATIVE_SCORE_CHAR);
      }
    }
  }
}

/*
  The following function formats an alignment given by the
  \texttt{firstlinecompare} and the \texttt{secondlinecompare}.
  \texttt{numofcols} is the number of columns in the alignment.
  \texttt{linewidth} is the width according to which the alignment
  is formatted. \texttt{startfirst} and \texttt{startsecond} are
  the starting positions of the aligned sequences.
  The output goes to the file pointer \texttt{outfp}.
*/

static void formatalignment(GtFile *outfp, GtUchar *firstlineorig,
                            GtUchar *secondlineorig, unsigned long numofcols,
                            unsigned long linewidth, unsigned long startfirst,
                            unsigned long startsecond, unsigned long totalulen,
                            GtAlphabet *alphabet,
                            GtArrayShortIntronInfo *shortintroninfo,
                            bool reverse_subject_pos, bool wildcardimplosion)
{
  unsigned long len, i = 0;
  unsigned long tennerblocksadjustment = 0,
       firstinsertioncount    = 0,
       secondinsertioncount   = 0,
       currentshortintroninfo = 0,
       completeshortintronlen = 0;
  long j, numofblanks;
  unsigned long shortintronstart,
       shortintronend,
       shortintronlength;

  for (;;) {
    if ((numofcols - i) < linewidth)
    {
      len = numofcols -i;
      tennerblocksadjustment = (linewidth - len - 1)/10;
    }
    else
    {
      len = linewidth;
    }
    formatseqwithgaps(outfp,firstlineorig + i, len, &firstinsertioncount, false,
                      alphabet, wildcardimplosion);
    if (reverse_subject_pos)
    {
      OUTSTRINGNUM(NUMWIDTH+linewidth-len+tennerblocksadjustment,
                   totalulen + 1 - (i+startfirst+len-firstinsertioncount+
                                    completeshortintronlen));
    }
    else
    {
      OUTSTRINGNUM(NUMWIDTH+linewidth-len+tennerblocksadjustment,
                   i+startfirst+len-firstinsertioncount+completeshortintronlen);
    }
    OUTCHAR('\n');
    showeditopline(outfp, firstlineorig+i,secondlineorig+i,len, alphabet);
    formatseqwithgaps(outfp, secondlineorig + i, len, &secondinsertioncount,
                      false, alphabet, wildcardimplosion);

    OUTSTRINGNUM(NUMWIDTH+linewidth-len+tennerblocksadjustment,
                 i+startsecond+len-secondinsertioncount);

    OUTCHAR('\n');
    OUTCHAR('\n');

    i += len;
    if (i >= numofcols)
      break;
    OUTCHAR('\n');

    /* take care of short introns, if some remain */
    if (currentshortintroninfo < shortintroninfo->nextfreeShortIntronInfo) {
      /* current short intron starts after or at current position */
      gt_assert(shortintroninfo->spaceShortIntronInfo[currentshortintroninfo]
             .start >= i + completeshortintronlen);
      if (shortintroninfo->spaceShortIntronInfo[currentshortintroninfo].start
          == i + completeshortintronlen) {
        /* output short intron */
        shortintronstart = reverse_subject_pos
                           ? totalulen + 1 -
                             (i+startfirst-firstinsertioncount+
                              completeshortintronlen+1)
                           : i+startfirst-firstinsertioncount +
                             completeshortintronlen+1;
        shortintronend   = reverse_subject_pos
                           ? totalulen + 1 -
                             (i+startfirst-firstinsertioncount+
                              completeshortintronlen+
                              shortintroninfo->spaceShortIntronInfo
                              [currentshortintroninfo].length)
                           : i+startfirst-firstinsertioncount+
                             completeshortintronlen+
                             shortintroninfo->spaceShortIntronInfo
                             [currentshortintroninfo].length;
        shortintronlength = shortintroninfo->spaceShortIntronInfo
                            [currentshortintroninfo].length;

        numofblanks  = (long) linewidth - 33;
        numofblanks -= floor(log10((double) shortintronstart))+1;
        numofblanks -= floor(log10((double) shortintronend))+1;
        numofblanks -= floor(log10((double) shortintronlength))+1;
        numofblanks += (linewidth / 10 ) - 1;

        gt_file_xprintf(outfp, "// intron part %lu %lu (%lu n) not shown",
                        shortintronstart, shortintronend, shortintronlength);

        for (j = 0; j < numofblanks; j++) {
          OUTCHAR(' ');
        }

        gt_file_xprintf(outfp, "//\n\n");
        gt_file_xprintf(outfp, "//");
        for (j = 2; j < (long) linewidth - 2; j++) {
          if ((j != 0) && (j%10 == 0)) {
            OUTCHAR(' ');
          }
          OUTCHAR(CONCRETEINTRONSYMBOL);
        }
        gt_file_xprintf(outfp, "//\n\n\n");

        /* update short intron adjustment */
        completeshortintronlen += shortintroninfo->spaceShortIntronInfo
                                 [currentshortintroninfo].length;
        /* increase current intron */
        currentshortintroninfo++;
      }
    }
  }
  OUTCHAR('\n');
}

static void formatproteinalignment(GtFile *outfp,
                                   GtUchar *genomicdnaline,
                                   GtUchar *genomicproteinline,
                                   GtUchar *referenceproteinline,
                                   unsigned long numofcols,
                                   unsigned long linewidth,
                                   unsigned long startfirst,
                                   unsigned long startsecond,
                                   unsigned long totalulen,
                                   GtAlphabet *alphabet,
                                   GtScoreMatrix *score_matrix,
                                   GtAlphabet *score_matrix_alphabet,
                                   bool reverse_subject_pos,
                                   bool wildcardimplosion)
{
  unsigned long len, codon_remainder, i = 0,
                tennerblocksadjustment = 0,
                genomicdnainsertioncount = 0,
                referenceproteininsertioncount = 0,
                completeshortintronlen = 0;

  for (;;)
  {
    if ((numofcols - i) < linewidth)
    {
      len = numofcols - i;
      tennerblocksadjustment = (linewidth - len - 1)/10;
    }
    else
    {
      len = linewidth;
    }

    formatseqwithgaps(outfp, genomicdnaline + i, len, &genomicdnainsertioncount,
                      false, alphabet, wildcardimplosion);
    if (reverse_subject_pos)
    {
      OUTSTRINGNUM(NUMWIDTH+linewidth-len+tennerblocksadjustment,
                   totalulen + 1 - (i+startfirst+len-genomicdnainsertioncount+
                                    completeshortintronlen));
    }
    else
    {
      OUTSTRINGNUM(NUMWIDTH+linewidth-len+tennerblocksadjustment,
                   i+startfirst+len-genomicdnainsertioncount+
                   completeshortintronlen);
    }
    OUTCHAR('\n');
    formatseqwithgaps(outfp, genomicproteinline + i, len, NULL, false, alphabet,
                      wildcardimplosion);
    OUTCHAR('\n');
    showeditoplineprotein(outfp, genomicproteinline + i,
                          referenceproteinline + i, len, score_matrix,
                          score_matrix_alphabet);
    OUTCHAR('\n');
    formatseqwithgaps(outfp, referenceproteinline + i, len,
                      &referenceproteininsertioncount, true, alphabet,
                      wildcardimplosion);

    /* this is necessary for a correct output of protein positions */
    if ((i+len-referenceproteininsertioncount) % GT_CODON_LENGTH == 2)
      codon_remainder = 1;
    else
      codon_remainder = 0;
    OUTSTRINGNUM(NUMWIDTH+linewidth-len+tennerblocksadjustment,
                 startsecond+
                 ((i+len-referenceproteininsertioncount) / GT_CODON_LENGTH)+
                 codon_remainder);
    OUTCHAR('\n');
    OUTCHAR('\n');

    i += len;
    if (i >= numofcols)
      break;
    OUTCHAR('\n');
  }
  OUTCHAR('\n');
}

unsigned long gthfillthethreealignmentlines(GtUchar **firstline,
                                   GtUchar **secondline,
                                   GtUchar **thirdline,
                                   Editoperation *alignment,
                                   unsigned long lenalg,
                                   unsigned long indelcount,
                                   const GtUchar *genseqorig,
                                   unsigned long genseqlen,
                                   const GtUchar *refseqorig,
                                   unsigned long refseqlen,
                                   unsigned long translationschemenumber)
{
  unsigned long lengthofgenomicdnaline,
                lengthofgenomicproteinline,
                lengthofreferenceproteinline;

  /* set the lenght of the three output lines */
  lengthofgenomicdnaline       = genseqlen
                                 + MIN(indelcount, refseqlen * GT_CODON_LENGTH);
  lengthofgenomicproteinline   = genseqlen
                                 + MIN(indelcount, refseqlen * GT_CODON_LENGTH);
  lengthofreferenceproteinline = refseqlen * GT_CODON_LENGTH
                                 + MIN(indelcount, genseqlen);

  /* alloc space for the three output lines */
  *firstline = gt_malloc(lengthofgenomicdnaline * sizeof (GtUchar));
  *secondline =  gt_malloc(lengthofgenomicproteinline * sizeof (GtUchar));
  *thirdline = gt_malloc(lengthofreferenceproteinline * sizeof (GtUchar));

  /* fill the output lines */
  return filltheproteinlines(*firstline,
                             *secondline,
                             *thirdline,
                             lengthofgenomicdnaline,
                             lengthofgenomicproteinline,
                             lengthofreferenceproteinline,
                             genseqorig,
                             refseqorig,
                             alignment,
                             lenalg,
                             translationschemenumber);
}

void gthshowalignmentprotein(GtFile *outfp,
                             unsigned long linewidth,
                             Editoperation *alignment,
                             unsigned long lenalg,
                             unsigned long indelcount,
                             const GtUchar *genseqorig,
                             unsigned long genseqlen,
                             const GtUchar *refseqorig,
                             unsigned long refseqlen,
                             unsigned long startfirst,
                             unsigned long startsecond,
                             unsigned long totalulen,
                             GT_UNUSED unsigned long showintronmaxlen,
                             GtAlphabet *alphabet,
                             unsigned long translationschemenumber,
                             GtScoreMatrix *score_matrix,
                             GtAlphabet *score_matrix_alphabet,
                             bool reverse_subject_pos,
                             bool wildcardimplosion)
{
  unsigned long numofcols;
  GtUchar *genomicdnaline,
          *genomicproteinline,
          *referenceproteinline;
  GtArrayShortIntronInfo shortintroninfo;

  /* init */
  GT_INITARRAY(&shortintroninfo, ShortIntronInfo);

  numofcols = gthfillthethreealignmentlines(&genomicdnaline,
                                            &genomicproteinline,
                                            &referenceproteinline,
                                            alignment,
                                            lenalg,
                                            indelcount,
                                            genseqorig,
                                            genseqlen,
                                            refseqorig,
                                            refseqlen,
                                            translationschemenumber);

  /* output the three lines in a formated fashion */
  formatproteinalignment(outfp, genomicdnaline, genomicproteinline,
                         referenceproteinline, numofcols, linewidth, startfirst,
                         startsecond, totalulen, alphabet, score_matrix,
                         score_matrix_alphabet, reverse_subject_pos,
                         wildcardimplosion);

  /* free */
  GT_FREEARRAY(&shortintroninfo, ShortIntronInfo);
  gt_free(genomicdnaline);
  gt_free(genomicproteinline);
  gt_free(referenceproteinline);
}

unsigned long gthfillthetwoalignmentlines(GtUchar **firstline,
                                          GtUchar **secondline,
                                          const GtUchar *useq,
                                          unsigned long ulen,
                                          const GtUchar *vseq,
                                          unsigned long vlen,
                                          Editoperation *alignment,
                                          unsigned long lenalg,
                                          unsigned long linewidth,
                                          unsigned long showintronmaxlen,
                                          GtArrayShortIntronInfo
                                          *shortintroninfo,
                                          unsigned long indelcount)
{
  *firstline = gt_malloc(ulen + MIN(indelcount, vlen) * sizeof (GtUchar));
  *secondline = gt_malloc(vlen + MIN(indelcount, ulen) * sizeof (GtUchar));

  return fillthelines(*firstline, *secondline, useq, vseq, alignment, lenalg,
                      linewidth, showintronmaxlen, shortintroninfo);
}

/*EE
  The following function shows the \texttt{alignment} on
  the file pointer \texttt{outfp}. \texttt{lenalg} is the length
  of the alignment and \texttt{indelcount} is an upper bound on the number
  of insertions and deletions. If you do not exactly know this number then
  choose a value which is larger, e.g.\ MAX(ulen,vlen). The aligned
  sequences are \texttt{useqcompare} and \texttt{vseqcompare} of
  length \texttt{ulen} and \texttt{vlen}. The pointers
  \texttt{useqorig} and \texttt{vseqorig} refer to the original versions
  of the sequences.
  \texttt{startfirst} and \texttt{startsecond} are
  the starting positions of the aligned sequences.
  \texttt{linewidth} is the width the alignment is formatted for,
  If the argument \texttt{showintronmaxlen} is set to 0, introns are showed
  completely. Otherwise introns larger than \texttt{showintronmaxlen} are only
  shown partially.
*/

void gthshowalignmentdna(GtFile *outfp,
                         unsigned long linewidth,
                         Editoperation *alignment,
                         unsigned long lenalg,
                         unsigned long indelcount,
                         const GtUchar *useqorig,
                         unsigned long ulen,
                         const GtUchar *vseqorig,
                         unsigned long vlen,
                         unsigned long startfirst,
                         unsigned long startsecond,
                         unsigned long totalulen,
                         unsigned long showintronmaxlen,
                         GtAlphabet *alphabet,
                         bool reverse_subject_pos,
                         bool wildcardimplosion)
{
  GtUchar *firstlineorig, *secondlineorig;
  unsigned long numofcolsorig;
  GtArrayShortIntronInfo shortintroninfo;

  GT_INITARRAY(&shortintroninfo, ShortIntronInfo);

  numofcolsorig = gthfillthetwoalignmentlines(&firstlineorig, &secondlineorig,
                                              useqorig, ulen, vseqorig, vlen,
                                              alignment, lenalg, linewidth,
                                              showintronmaxlen,
                                              &shortintroninfo, indelcount);

  formatalignment(outfp, firstlineorig, secondlineorig, numofcolsorig,
                  linewidth, startfirst, startsecond, totalulen, alphabet,
                  &shortintroninfo, reverse_subject_pos, wildcardimplosion);
  gt_free(firstlineorig);
  gt_free(secondlineorig);
  GT_FREEARRAY(&shortintroninfo, ShortIntronInfo);
}
