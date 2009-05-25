/*
  Copyright (c) 2001-2003 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
  Copyright (c)      2009 Sascha Steinbiss <steinbiss@zbh.uni-hamburg.de>
  Copyright (c)      2009 Center for Bioinformatics, University of Hamburg

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

#include <string.h>
#include "core/assert_api.h"
#include "core/codon.h"
#include "core/ensure.h"
#include "core/ma.h"
#include "core/translator.h"
#include "core/unused_api.h"

#define GT_STANDARDTRANSLATIONTABLE  1
#define GT_AMINOACIDFAIL            -1

#define GT_NUMOFTRANSSCHEMES\
        ((unsigned int) (sizeof (schemetable)/sizeof (schemetable[0])))
#define GT_SIZEOFTRANSRANGE\
        ((unsigned int) (sizeof (transnum2index)/sizeof (transnum2index[0])))

#define GT_UNDEFTRANSNUM GT_NUMOFTRANSSCHEMES

#define GT_INCONSISTENT(BASE)\
        gt_error_set(err, "code=%lu with wildcard %c: inconsistent " \
                          "aminoacids %c and %c",\
                     (unsigned long) codeof2, wildcard, aa, newaa);\
        return GT_AMINOACIDFAIL

#define GT_ILLEGALCHAR(V)\
        gt_error_set(err, "illegal char %s='%c'(%lu)",#V,V,(unsigned long)(V));\
        return GT_AMINOACIDFAIL

#define GT_STARTAMINO 'M'

#define GT_T_BIT  ((unsigned char) 1)
#define GT_C_BIT (((unsigned char) 1) << 1)
#define GT_A_BIT (((unsigned char) 1) << 2)
#define GT_G_BIT (((unsigned char) 1) << 3)

#define GT_CASEWILDCARD\
          case 'n':\
          case 'N':\
          case 's':\
          case 'S':\
          case 'y':\
          case 'Y':\
          case 'w':\
          case 'W':\
          case 'r':\
          case 'R':\
          case 'k':\
          case 'K':\
          case 'v':\
          case 'V':\
          case 'b':\
          case 'B':\
          case 'd':\
          case 'D':\
          case 'h':\
          case 'H':\
          case 'm':\
          case 'M'

typedef struct
{
  const char *name;         /* the name of the translation */
  unsigned int identity;    /* the identity number */
  const char *aminos,       /* the amino acids in order */
             *startcodon;   /* the start codons */
} GtTranslationScheme;

/* according to http://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi */
static GtTranslationScheme schemetable[] = {
  {"Standard",
   (unsigned int) 1,
   "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
   "---M---------------M---------------M----------------------------"},
  {"Vertebrate Mitochondrial",
   (unsigned int) 2,
   "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSS**VVVVAAAADDEEGGGG",
   "--------------------------------MMMM---------------M------------"},
  {"Yeast Mitochondrial",
   (unsigned int) 3,
   "FFLLSSSSYY**CCWWTTTTPPPPHHQQRRRRIIMMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
   "-----------------------------------M----------------------------"},
  {"Mold Mitochondrial; Protozoan Mitochondrial; Coelenterate Mitochondrial; "
   "Mycoplasma; Spiroplasma",
   (unsigned int) 4,
   "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
   "--MM---------------M------------MMMM---------------M------------"},
  {"Invertebrate Mitochondrial",
   (unsigned int) 5,
   "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSSSSVVVVAAAADDEEGGGG",
   "---M----------------------------MMMM---------------M------------"},
  {"Ciliate Nuclear; Dasycladacean Nuclear; Hexamita Nuclear",
   (unsigned int) 6,
   "FFLLSSSSYYQQCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
   "-----------------------------------M----------------------------"},
  {"Echinoderm Mitochondrial",
   (unsigned int) 9,
   "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNNKSSSSVVVVAAAADDEEGGGG",
   "-----------------------------------M----------------------------"},
  {"Euplotid Nuclear",
   (unsigned int) 10,
   "FFLLSSSSYY**CCCWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
   "-----------------------------------M----------------------------"},
  {"Bacterial",
   (unsigned int) 11,
   "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
   "---M---------------M------------MMMM---------------M------------"},
  {"Alternative Yeast Nuclear",
   (unsigned int) 12,
   "FFLLSSSSYY**CC*WLLLSPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
   "-------------------M---------------M----------------------------"},
  {"Ascidian Mitochondrial",
   (unsigned int) 13,
   "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSSGGVVVVAAAADDEEGGGG",
   "-----------------------------------M----------------------------"},
  {"Flatworm Mitochondrial",
   (unsigned int) 14,
   "FFLLSSSSYYY*CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNNKSSSSVVVVAAAADDEEGGGG",
   "-----------------------------------M----------------------------"},
  {"Blepharisma Macronuclear",
   (unsigned int) 15,
   "FFLLSSSSYY*QCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
   "-----------------------------------M----------------------------"},
  {"Chlorophycean Mitochondrial",
   (unsigned int) 16,
   "FFLLSSSSYY*LCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
   "-----------------------------------M----------------------------"},
  {"Trematode Mitochondrial",
   (unsigned int) 21,
   "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNNKSSSSVVVVAAAADDEEGGGG",
   "-----------------------------------M----------------------------"},
  {"Scenedesmus Obliquus Mitochondrial",
   (unsigned int) 22,
   "FFLLSS*SYY*LCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
   "-----------------------------------M----------------------------"},
  {"Thraustochytrium Mitochondrial",
   (unsigned int) 23,
   "FF*LSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
   "--------------------------------M--M---------------M------------"}
};

static unsigned int transnum2index[] =
{
  GT_UNDEFTRANSNUM,
  0U,
  1U,
  2U,
  3U,
  4U,
  5U,
  GT_UNDEFTRANSNUM,
  GT_UNDEFTRANSNUM,
  6U,
  7U,
  8U,
  9U,
  10U,
  11U,
  12U,
  13U,
  GT_UNDEFTRANSNUM,
  GT_UNDEFTRANSNUM,
  GT_UNDEFTRANSNUM,
  GT_UNDEFTRANSNUM,
  14U,
  15U,
  16U
};

static GtTranslationScheme* getschemetable(unsigned int transnum, GtError *err)
{
  gt_assert(transnum < GT_SIZEOFTRANSRANGE);
  if (transnum2index[transnum] == GT_UNDEFTRANSNUM) {
    gt_error_set(err, "'%u' is not a valid translation table number!",
                 transnum);
    return NULL;
  }
  gt_assert(transnum != GT_UNDEFTRANSNUM);
  return schemetable + transnum2index[transnum];
}

/*
  The following table specifies a unsigned characters
  whose bits are set according to the characters encoded by a wildcard.
*/
static unsigned char wbitsvector[] =
{
  0 /* 0 */,
  0 /* 1 */,
  0 /* 2 */,
  0 /* 3 */,
  0 /* 4 */,
  0 /* 5 */,
  0 /* 6 */,
  0 /* 7 */,
  0 /* 8 */,
  0 /* 9 */,
  0 /* 10 */,
  0 /* 11 */,
  0 /* 12 */,
  0 /* 13 */,
  0 /* 14 */,
  0 /* 15 */,
  0 /* 16 */,
  0 /* 17 */,
  0 /* 18 */,
  0 /* 19 */,
  0 /* 20 */,
  0 /* 21 */,
  0 /* 22 */,
  0 /* 23 */,
  0 /* 24 */,
  0 /* 25 */,
  0 /* 26 */,
  0 /* 27 */,
  0 /* 28 */,
  0 /* 29 */,
  0 /* 30 */,
  0 /* 31 */,
  0 /* 32 */,
  0 /* 33 */,
  0 /* 34 */,
  0 /* 35 */,
  0 /* 36 */,
  0 /* 37 */,
  0 /* 38 */,
  0 /* 39 */,
  0 /* 40 */,
  0 /* 41 */,
  0 /* 42 */,
  0 /* 43 */,
  0 /* 44 */,
  0 /* 45 */,
  0 /* 46 */,
  0 /* 47 */,
  0 /* 48 */,
  0 /* 49 */,
  0 /* 50 */,
  0 /* 51 */,
  0 /* 52 */,
  0 /* 53 */,
  0 /* 54 */,
  0 /* 55 */,
  0 /* 56 */,
  0 /* 57 */,
  0 /* 58 */,
  0 /* 59 */,
  0 /* 60 */,
  0 /* 61 */,
  0 /* 62 */,
  0 /* 63 */,
  0 /* 64 */,
  0 /* 65 */,
  GT_C_BIT | GT_G_BIT | GT_T_BIT         /* [cgt] */ /* b */,
  0 /* 67 */,
  GT_A_BIT | GT_G_BIT | GT_T_BIT         /* [agt] */ /* d */,
  0 /* 69 */,
  0 /* 70 */,
  0 /* 71 */,
  GT_A_BIT | GT_C_BIT | GT_T_BIT         /* [act] */ /* h */,
  0 /* 73 */,
  0 /* 74 */,
  GT_G_BIT | GT_T_BIT                /* [gt] */ /* k */,
  0 /* 76 */,
  GT_A_BIT | GT_C_BIT                /* [ac] */ /* m */,
  GT_A_BIT | GT_C_BIT | GT_G_BIT | GT_T_BIT  /* [acgt] */ /* n */,
  0 /* 79 */,
  0 /* 80 */,
  0 /* 81 */,
  GT_A_BIT | GT_G_BIT                /* [ag] */ /* r */,
  GT_C_BIT | GT_G_BIT                /* [cg] */ /* s */,
  0 /* 84 */,
  0 /* 85 */,
  GT_A_BIT | GT_C_BIT | GT_G_BIT         /* [acg] */ /* v */,
  GT_A_BIT | GT_C_BIT                /* [at] */ /* w */,
  0 /* 88 */,
  GT_C_BIT | GT_T_BIT                /* [ct] */ /* y */,
  0 /* 90 */,
  0 /* 91 */,
  0 /* 92 */,
  0 /* 93 */,
  0 /* 94 */,
  0 /* 95 */,
  0 /* 96 */,
  0 /* 97 */,
  GT_C_BIT | GT_G_BIT | GT_T_BIT         /* [cgt] */ /* B */,
  0 /* 99 */,
  GT_A_BIT | GT_G_BIT | GT_T_BIT         /* [agt] */ /* D */,
  0 /* 101 */,
  0 /* 102 */,
  0 /* 103 */,
  GT_A_BIT | GT_C_BIT | GT_T_BIT         /* [act] */ /* H */,
  0 /* 105 */,
  0 /* 106 */,
  GT_G_BIT | GT_T_BIT                /* [gt] */ /* K */,
  0 /* 108 */,
  GT_A_BIT | GT_C_BIT                /* [ac] */ /* M */,
  GT_A_BIT | GT_C_BIT | GT_G_BIT | GT_T_BIT  /* [acgt] */ /* N */,
  0 /* 111 */,
  0 /* 112 */,
  0 /* 113 */,
  GT_A_BIT | GT_G_BIT                /* [ag] */ /* R */,
  GT_C_BIT | GT_G_BIT                /* [cg] */ /* S */,
  0 /* 116 */,
  0 /* 117 */,
  GT_A_BIT | GT_C_BIT | GT_G_BIT         /* [acg] */ /* V */,
  GT_A_BIT | GT_C_BIT                /* [at] */ /* W */,
  0 /* 120 */,
  GT_C_BIT | GT_T_BIT                /* [ct] */ /* Y */,
  0 /* 122 */,
  0 /* 123 */,
  0 /* 124 */,
  0 /* 125 */,
  0 /* 126 */,
  0 /* 127 */,
  0 /* 128 */,
  0 /* 129 */,
  0 /* 130 */,
  0 /* 131 */,
  0 /* 132 */,
  0 /* 133 */,
  0 /* 134 */,
  0 /* 135 */,
  0 /* 136 */,
  0 /* 137 */,
  0 /* 138 */,
  0 /* 139 */,
  0 /* 140 */,
  0 /* 141 */,
  0 /* 142 */,
  0 /* 143 */,
  0 /* 144 */,
  0 /* 145 */,
  0 /* 146 */,
  0 /* 147 */,
  0 /* 148 */,
  0 /* 149 */,
  0 /* 150 */,
  0 /* 151 */,
  0 /* 152 */,
  0 /* 153 */,
  0 /* 154 */,
  0 /* 155 */,
  0 /* 156 */,
  0 /* 157 */,
  0 /* 158 */,
  0 /* 159 */,
  0 /* 160 */,
  0 /* 161 */,
  0 /* 162 */,
  0 /* 163 */,
  0 /* 164 */,
  0 /* 165 */,
  0 /* 166 */,
  0 /* 167 */,
  0 /* 168 */,
  0 /* 169 */,
  0 /* 170 */,
  0 /* 171 */,
  0 /* 172 */,
  0 /* 173 */,
  0 /* 174 */,
  0 /* 175 */,
  0 /* 176 */,
  0 /* 177 */,
  0 /* 178 */,
  0 /* 179 */,
  0 /* 180 */,
  0 /* 181 */,
  0 /* 182 */,
  0 /* 183 */,
  0 /* 184 */,
  0 /* 185 */,
  0 /* 186 */,
  0 /* 187 */,
  0 /* 188 */,
  0 /* 189 */,
  0 /* 190 */,
  0 /* 191 */,
  0 /* 192 */,
  0 /* 193 */,
  0 /* 194 */,
  0 /* 195 */,
  0 /* 196 */,
  0 /* 197 */,
  0 /* 198 */,
  0 /* 199 */,
  0 /* 200 */,
  0 /* 201 */,
  0 /* 202 */,
  0 /* 203 */,
  0 /* 204 */,
  0 /* 205 */,
  0 /* 206 */,
  0 /* 207 */,
  0 /* 208 */,
  0 /* 209 */,
  0 /* 210 */,
  0 /* 211 */,
  0 /* 212 */,
  0 /* 213 */,
  0 /* 214 */,
  0 /* 215 */,
  0 /* 216 */,
  0 /* 217 */,
  0 /* 218 */,
  0 /* 219 */,
  0 /* 220 */,
  0 /* 221 */,
  0 /* 222 */,
  0 /* 223 */,
  0 /* 224 */,
  0 /* 225 */,
  0 /* 226 */,
  0 /* 227 */,
  0 /* 228 */,
  0 /* 229 */,
  0 /* 230 */,
  0 /* 231 */,
  0 /* 232 */,
  0 /* 233 */,
  0 /* 234 */,
  0 /* 235 */,
  0 /* 236 */,
  0 /* 237 */,
  0 /* 238 */,
  0 /* 239 */,
  0 /* 240 */,
  0 /* 241 */,
  0 /* 242 */,
  0 /* 243 */,
  0 /* 244 */,
  0 /* 245 */,
  0 /* 246 */,
  0 /* 247 */,
  0 /* 248 */,
  0 /* 249 */,
  0 /* 250 */,
  0 /* 251 */,
  0 /* 252 */,
  0 /* 253 */,
  0 /* 254 */,
  0 /* 255 */
};

/*
  The following function tries to solve the case, where
  the third base of a codon is a <wildcard>.
  It takes the code <codeof2> of the first two bases
  of a codon and checks if all characters encoded by <wildcard>
  lead to the same amino acid, when doing the translation.
  This means that the characters encoded by <wildcard> are
  equivalent.
*/
static char equivalentbits(const char *aminos,
                           unsigned int codeof2,
                           unsigned char wildcard,
                           GtError *err)
{
  unsigned char bits = wbitsvector[(int) wildcard];
  char aa = 0, newaa;
  bool aaundefined = true;

  if (bits & GT_T_BIT)
  {
    aa = aminos[codeof2 + GT_T_CODE];
    aaundefined = false;
  }
  if (bits & GT_C_BIT)
  {
    newaa = aminos[codeof2 + GT_C_CODE];
    if (aaundefined)
    {
      aa = newaa;
      aaundefined = false;
    } else
    {
      if (aa != newaa)
      {
        GT_INCONSISTENT('C');
      }
    }
  }
  if (bits & GT_A_BIT)
  {
    newaa = aminos[codeof2 + GT_A_CODE];
    if (aaundefined)
    {
      aa = newaa;
      aaundefined = false;
    } else
    {
      if (aa != newaa)
      {
        GT_INCONSISTENT('A');
      }
    }
  }
  if (bits & GT_G_BIT)
  {
    newaa = aminos[codeof2 + GT_G_CODE];
    if (aaundefined)
    {
      aa = newaa;
      aaundefined = false;
    } else
    {
      if (aa != newaa)
      {
        GT_INCONSISTENT('G');
      }
    }
  }
  if (aaundefined)
  {
    return GT_AMINOACIDFAIL;
  }
  return aa;
}

/*
  The following function finds the smallest character in the set
  of characters encoded by a wildcard. This set is given by the bit vector
  <bits>.
*/
static unsigned int smallestbase(unsigned char bits)
{
  if (bits & GT_T_BIT)
  {
    return GT_T_CODE;
  }
  if (bits & GT_C_BIT)
  {
    return GT_C_CODE;
  }
  if (bits & GT_A_BIT)
  {
    return GT_A_CODE;
  }
  if (bits & GT_G_BIT)
  {
    return GT_G_CODE;
  }
  fprintf(stderr, "%s: line %lu: This case is not supposed to occur\n",\
                  __FILE__, (unsigned long) __LINE__);\
  exit(EXIT_FAILURE);
  /*@ignore@*/
  return 0;
  /*@end@*/
}

static char codon2amino(const char *aminos, bool forward, unsigned char c0,
                        unsigned char c1,unsigned char c2,
                        unsigned int *coderet,
                        GtError *err)
{
  unsigned int code = 0;
  char aa;

  switch (c0)
  {
    case 'a':
    case 'A':
      if (forward)
      {
        code = (GT_A_CODE << 4);
      } else
      {
        code = 0;
      }
      break;
    case 'c':
    case 'C':
      if (forward)
      {
        code = (GT_C_CODE << 4);
      } else
      {
        code = (GT_G_CODE << 4);
      }
      break;
    case 'g':
    case 'G':
      if (forward)
      {
        code = (GT_G_CODE << 4);
      } else
      {
        code = (GT_C_CODE << 4);
      }
      break;
    case 't':
    case 'T':
    case 'u':
    case 'U':
      if (forward)
      {
        code = 0;
      } else
      {
        code = (GT_A_CODE << 4);
      }
      break;
    GT_CASEWILDCARD: /* delete this and the next line,
                        to inform about wildcards */
      code = (smallestbase(wbitsvector[(int) c0]) << 4);
      break;
    default:
      GT_ILLEGALCHAR(c0);
  }
  switch (c1)
  {
    case 'a':
    case 'A':
      if (forward)
      {
        code += (GT_A_CODE << 2);
      }
      break;
    case 'c':
    case 'C':
      if (forward)
      {
        code += (GT_C_CODE << 2);
      } else
      {
        code += (GT_G_CODE << 2);
      }
      break;
    case 'g':
    case 'G':
      if (forward)
      {
        code += (GT_G_CODE << 2);
      } else
      {
        code += (GT_C_CODE << 2);
      }
      break;
    case 't':
    case 'T':
    case 'u':
    case 'U':
      if (!forward)
      {
        code += (GT_A_CODE << 2);
      }
      break;
    GT_CASEWILDCARD: /* delete this and the next line,
                        to inform about wildcards */
      code += (smallestbase(wbitsvector[(int) c1]) << 2);
      break;
    default:
      GT_ILLEGALCHAR(c1);
  }
  switch (c2)
  {
    case 'a':
    case 'A':
      if (forward)
      {
        code += GT_A_CODE;
      }
      break;
    case 'c':
    case 'C':
      if (forward)
      {
        code += GT_C_CODE;
      } else
      {
        code += GT_G_CODE;
      }
      break;
    case 'g':
    case 'G':
      if (forward)
      {
        code += GT_G_CODE;
      } else
      {
        code += GT_C_CODE;
      }
      break;
    case 't':
    case 'T':
    case 'u':
    case 'U':
      if (!forward)
      {
        code += GT_A_CODE;
      }
      break;
    GT_CASEWILDCARD:
      aa = equivalentbits(aminos,code,c2, err);
      if (aa == GT_AMINOACIDFAIL)
      {
        /* no unique aminoacid => choose smallest base and compute aminos
           accordingly */
        code += smallestbase(wbitsvector[(int) c2]);
      } else
      {
        if (coderet != NULL)
        {
          *coderet = code;
        }
        return aa;
      }
      break;
    default:
      GT_ILLEGALCHAR(c2);
  }
  if (coderet != NULL)
    *coderet = code;
  return aminos[code];
}

struct GtTranslator {
  GtTranslationScheme *scheme;
  const char *dnaseq,
             *origseq;
  unsigned long dnalen,
                dnapos;
  unsigned long curframe;
};

GtTranslator* gt_translator_new()
{
  GtTranslator *tr;
  tr = gt_calloc(1, sizeof (GtTranslator));
  /* set standard translation table as default */
  tr->scheme = getschemetable(GT_STANDARDTRANSLATIONTABLE, NULL);
  gt_assert(tr->scheme); /* standard table must exist */
  return tr;
}

int gt_translator_set_translation_scheme(GtTranslator *tr,
                                         unsigned int transnum,
                                         GtError *err)
{
  gt_assert(tr && transnum < GT_NUMOFTRANSSCHEMES);
  tr->scheme = getschemetable(transnum, err);
  if (!tr->scheme)
    return -1;
  return 0;
}

int gt_translator_find_startcodon(GtTranslator *tr,
                                  const char *dnaseq,
                                  unsigned long dnalen,
                                  unsigned long *pos,
                                  GtError *err)
{
  unsigned int code;
  char *dnaptr, trans;
  gt_assert(tr && dnaseq && dnalen >= GT_CODON_LENGTH);

  for (dnaptr = (char*) dnaseq; dnaptr < dnaseq + dnalen - 3; dnaptr++) {
    trans = codon2amino(tr->scheme->aminos, true, *dnaptr, *(dnaptr+1),
                        *(dnaptr+2), &code, err);
    if (trans == GT_AMINOACIDFAIL)
      return -2;
    if (tr->scheme->startcodon[code] == GT_STARTAMINO) {
      *pos = (unsigned long) (dnaptr - dnaseq);
      return 0;
    }
  }
  gt_error_set(err, "cannot find start codon");
  return -1;
}

GtStrArray* gt_translator_get_translation_table_descriptions()
{
  unsigned long i;
  GtTranslationScheme *scheme;
  GtStr *str;
  GtStrArray *sa = gt_str_array_new();
  str = gt_str_new();
  for (i = 1; i < GT_SIZEOFTRANSRANGE; i++) {
    if (transnum2index[i] == GT_UNDEFTRANSNUM)
      continue;
    scheme = schemetable + transnum2index[i];
    gt_str_reset(str);
    gt_str_append_uint(str, scheme->identity);
    gt_str_append_cstr(str, ": ");
    gt_str_append_cstr(str, scheme->name);
    gt_str_array_add_cstr(sa, gt_str_get(str));
  }
  gt_str_delete(str);
  return sa;
}

int gt_translator_codon2amino(GtTranslator *tr,
                              char c1, char c2, char c3,
                              char *translated,
                              GtError *err)
{
  gt_assert(tr);
  gt_error_check(err);
  *translated = codon2amino(tr->scheme->aminos, true, c1, c2, c3, NULL, err);
  if (*translated == GT_AMINOACIDFAIL) {
    return -2;
  }
  return 0;
}

int gt_translator_start(GtTranslator *tr,
                        const char *dnaseq,
                        unsigned long dnalen,
                        char *translated,
                        unsigned int *frame,
                        GtError *err)
{
  gt_assert(tr && dnaseq && translated && frame && dnalen >= GT_CODON_LENGTH);
  gt_error_check(err);
  *translated = codon2amino(tr->scheme->aminos, true, *dnaseq, *(dnaseq+1),
                            *(dnaseq+2), NULL, err);
  *frame = tr->curframe % GT_CODON_LENGTH;
  gt_assert(*frame < GT_CODON_LENGTH);
  tr->curframe++;
  tr->curframe %= GT_CODON_LENGTH;
  if (*translated == GT_AMINOACIDFAIL)
    return -2;
  tr->dnaseq = tr->origseq = dnaseq;
  tr->dnaseq++;
  tr->dnalen = dnalen;
  return 0;
}

int gt_translator_next(GtTranslator *tr,
                       char *translated,
                       unsigned int *frame,
                       GtError *err)
{
  gt_assert(tr && translated && frame);
  gt_error_check(err);
  if (!tr->dnaseq || !tr->origseq) {
    gt_error_set(err, "No translation iteration was started yet");
    return -1;
  }
  if (tr->dnalen - (tr->dnaseq - tr->origseq) <= 2) {
    *translated = 0;
    return 0;
  }
  *translated = codon2amino(tr->scheme->aminos, true,
                            *tr->dnaseq,
                            *(tr->dnaseq+1),
                            *(tr->dnaseq+2), NULL, err);
  *frame = tr->curframe % GT_CODON_LENGTH;
  gt_assert(*frame < GT_CODON_LENGTH);
  tr->curframe++;
  tr->curframe %= GT_CODON_LENGTH;
  if (*translated == GT_AMINOACIDFAIL)
    return -2;
  tr->dnaseq++;
  return 0;
}

int gt_translator_translate_string(GtTranslator *tr,
                                   GtStr *protein,
                                   const char *dnaseq,
                                   unsigned long dnalen,
                                   unsigned int frame,
                                   GtError *err)
{
  char *dnaptr;
  char aa;
  gt_assert(protein && !gt_str_length(protein) && dnaseq
            && frame < GT_CODON_LENGTH && dnalen >= GT_CODON_LENGTH);
  gt_error_check(err);

  for (dnaptr = (char*) (dnaseq + frame);
      dnaptr < dnaseq + dnalen - 2;
      dnaptr += GT_CODON_LENGTH)
  {
    aa = codon2amino(tr->scheme->aminos, true, *dnaptr, *(dnaptr+1),
                     *(dnaptr+2), NULL, err);
    if (aa == GT_AMINOACIDFAIL)
      return (int) -2;
    gt_str_append_char(protein, aa);
  }
  return 0;
}

void gt_translator_delete(GtTranslator *tr)
{
  if (!tr) return;
  gt_free(tr);
}

int gt_translator_unit_test(GtError *err)
{
  int had_err = 0, test_errnum = 0, i = 0, j = 0, k = 0;
  GtTranslator *tr;
  GtError *test_err;
  const char *seq = "AGCTTTTCATTCTGACTGCAACGGGCAATATGTCTCTGTGT"
                    "GGATTAAAAAAAGAGTGTCTGATAGCAGCTTCTGAACTGGT"
                    "TACCTGCCGTGAGTAAATTAAAATTTTATTGACTTAGG";
  const char *no_startcodon = "AAAAAAAAAATCATCTCCCCATTTTTTT";
  const char *invalidseq  = "ZAGCTTTTCATTCTGACTGCAAATATGTCTCTGTGT";
  const char *invalidseq2 = "AGCTTTTCATTCTGACZTGCAAATATGTCTCTGTGT";
  char *bases = "AaCcGgTt";
  char translated;
  unsigned int frame;
  unsigned long pos = 0;
  GtStr *protein[3];
  GtStrArray *schemes;
  gt_error_check(err);

  tr = gt_translator_new();
  test_err = gt_error_new();
  protein[0] = gt_str_new();
  protein[1] = gt_str_new();
  protein[2] = gt_str_new();

  /* check uninitialised iterator use */
  test_errnum = gt_translator_next(tr, &translated, &frame, test_err);
  ensure(had_err, test_errnum && gt_error_is_set(test_err));

  /* do 3-frame translation */
  gt_error_unset(test_err);
  test_errnum = gt_translator_start(tr, seq, strlen(seq), &translated,
                                   &frame, test_err);
  while (!test_errnum && translated) {
    gt_str_append_char(protein[frame], translated);
    test_errnum = gt_translator_next(tr, &translated, &frame, test_err);
    ensure(had_err, !test_errnum && !gt_error_is_set(test_err));
  }

  /* check 3-frame translation */
  ensure(had_err, strcmp(gt_str_get(protein[0]),
                         "SFSF*LQRAICLCVD*KKSV**QLLNWLPAVSKLKFY*LR") == 0);
  ensure(had_err, strcmp(gt_str_get(protein[1]),
                         "AFHSDCNGQYVSVWIKKRVSDSSF*TGYLP*VN*NFIDL") == 0);
  ensure(had_err, strcmp(gt_str_get(protein[2]),
                         "LFILTATGNMSLCGLKKECLIAASELVTCRE*IKILLT*") == 0);

  /* check whole string translation */
  for (i=0;i<3;i++) {
    gt_str_reset(protein[i]);
    test_errnum = gt_translator_translate_string(tr, protein[i], seq,
                                                 strlen(seq), i, test_err);
    ensure(had_err, !test_errnum && !gt_error_is_set(test_err));
  }
  ensure(had_err, strcmp(gt_str_get(protein[0]),
                         "SFSF*LQRAICLCVD*KKSV**QLLNWLPAVSKLKFY*LR") == 0);
  ensure(had_err, strcmp(gt_str_get(protein[1]),
                         "AFHSDCNGQYVSVWIKKRVSDSSF*TGYLP*VN*NFIDL") == 0);
  ensure(had_err, strcmp(gt_str_get(protein[2]),
                         "LFILTATGNMSLCGLKKECLIAASELVTCRE*IKILLT*") == 0);

  /* check retrieval of table descriptions */
  schemes = gt_translator_get_translation_table_descriptions();
  ensure(had_err, gt_str_array_size(schemes) == GT_NUMOFTRANSSCHEMES);

  /* check switching translation scheme */
  test_errnum = gt_translator_set_translation_scheme(tr, 3, test_err);
  ensure(had_err, !test_errnum && !gt_error_is_set(test_err));

  /* check switching to invalid translation scheme */
  test_errnum = gt_translator_set_translation_scheme(tr, 7, test_err);
  ensure(had_err, test_errnum && gt_error_is_set(test_err));

  /* switch back to default translation scheme */
  gt_error_unset(test_err);
  test_errnum = gt_translator_set_translation_scheme(tr, 1, test_err);
  ensure(had_err, !test_errnum && !gt_error_is_set(test_err));

  /* check translation of sequence with invalid beginning */
  test_errnum = gt_translator_start(tr, invalidseq, strlen(invalidseq),
                                    &translated, &frame, test_err);
  ensure(had_err, test_errnum && gt_error_is_set(test_err));

  /* check translation of sequence with invalid character within */
  gt_error_unset(test_err);
  test_errnum = gt_translator_start(tr, invalidseq2, strlen(invalidseq2),
                                    &translated, &frame, test_err);
  while (!test_errnum && translated) {
    gt_str_append_char(protein[frame], translated);
    test_errnum = gt_translator_next(tr, &translated, &frame, test_err);
  }
  ensure(had_err, test_errnum && gt_error_is_set(test_err));

  /* find start codon -- positive */
  gt_error_unset(test_err);
  test_errnum = gt_translator_find_startcodon(tr, seq, strlen(seq), &pos,
                                              test_err);
  ensure(had_err, !test_errnum && !gt_error_is_set(test_err));
  ensure(had_err, pos == 11);
  /* find start codon -- fail */
  test_errnum = gt_translator_find_startcodon(tr, no_startcodon,
                                              strlen(no_startcodon),
                                              &pos, test_err);
  ensure(had_err, test_errnum && gt_error_is_set(test_err));

  /* check single codon translation */
  gt_error_unset(test_err);
  for (i=0; i<8; i++) {
    char c1 = bases[i];
    for (j=0; j<8; j++) {
      char c2 = bases[j];
      for (k=0; k<8; k++) {
        char c3 = bases[k], ret1, ret2;
        test_errnum = gt_translator_codon2amino(tr, c1, c2, c3, &ret1,
                                                test_err);
        ensure(had_err, !test_errnum && !gt_error_is_set(test_err));
        ret2 = codon2amino(tr->scheme->aminos, true, c1, c2, c3, NULL,
                           test_err);
        ensure(had_err, ret1 == ret2);
      }
    }
  }

  gt_translator_delete(tr);
  gt_str_delete(protein[0]);
  gt_str_delete(protein[1]);
  gt_str_delete(protein[2]);
  gt_str_array_delete(schemes);
  gt_error_delete(test_err);

  return had_err;
}
