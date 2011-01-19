/*
  Copyright (c) 2007-2009 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
  Copyright (c) 2007-2009 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2007-2009 Center for Bioinformatics, University of Hamburg

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
#include <ctype.h>
#include <stdlib.h>
#include <string.h>
#include "core/alphabet.h"
#include "core/chardef.h"
#include "core/cstr_api.h"
#include "core/fileutils_api.h"
#include "core/fa.h"
#include "core/gtdatapath.h"
#include "core/ma_api.h"
#include "core/mathsupport.h"
#include "core/str_array.h"
#include "core/xansi_api.h"

#define ALPHABET_GUESS_MAX_LENGTH       5000
#define ALPHABET_GUESS_PROTEIN_CHARS    "LIFEQPlifeqpXZ*-"

struct GtAlphabet {
  unsigned int domainsize,           /* size of domain of symbolmap */
               mapsize,              /* size of image of map, i.e. */
                                     /* mapping to [0..mapsize-1] */
               mappedwildcards,      /* number of mapped wildcards */
               bitspersymbol,        /* number of bits per symbol in
                                        bitspackedarray */
               reference_count;
  GtUchar wildcardshow,
          symbolmap[GT_MAXALPHABETCHARACTER+1], /* mapping of the symbols */
          *mapdomain,                        /* list of characters mapped */
          *characters;                       /* array of characters to show */
};

/*
  Some constants for the standard alphabet used. The name says it all.
*/

#define DNABASES                     "aAcCgGtTuU"
#define DNAWILDCARDS                 "nsywrkvbdhmNSYWRKVBDHM"
#define MAPSIZEDNA                   (GT_DNAALPHASIZE+1U)
#define DNAALPHABETDOMAIN            DNABASES DNAWILDCARDS
#define PROTEINUPPERAMINOACIDS       "LVIFKREDAGSTNQYWPHMC"
#define PROTEINALPHASIZE             20U
#define MAPSIZEPROTEIN               (PROTEINALPHASIZE+1U)
#define PROTEINWILDCARDS             "XUBZO*-"
#define PROTEINALPHABETDOMAIN        PROTEINUPPERAMINOACIDS PROTEINWILDCARDS

/*
  We use the following macro to access the \texttt{I}-th character of
  a line.
*/

#define LINE(I)          currentline[I]

/*
  We have developed a simple format to specify an alphabet
  and a corresponding alphabet transformation. This format specifies the
  characters of the alphabet (including wild card characters)
  and the pairs of symbols which are to be considered identical.
  The format is best explained by some examples. Consider a file
  containing the following lines:
  \begin{alltt}
  aA
  cC
  gG
  tTuU
  nsywrkvbdhmNSYWRKVBDHM
  \end{alltt}
  These line specify that the input sequence
  are allowed to contain the symbols \(a,c,g,t,u,n,s,y\)
  \(w,r,k,v,b,d,h,m\)
  in either lower or upper case. Moreover, the first four lines specify that
  \(a=A\), \(c=C\), \(g=G\), and \(t=T=u=U\). The last line specifies the
  wildcard symbols, which are replaced by unique symbols. Note that any
  wildcard symbol does not even match itself, if it occurs at different
  positions. Thus no false matches will be delivered. Note that however,
  a degenerate match may contain a wildcard, since this always leads
  to a mismatch.

  Consider a file containing the following lines:
  \begin{alltt}
  LVIF i
  KR +
  ED -
  AG s
  ST o
  NQ n
  YW a
  P p
  H h
  M m
  C c
  XBZ* x
  \end{alltt}
  This specifies the Protein alphabet
  \(L,V,I,F,K,R,E,D,A,G,S,T,N,Q,Y,W,P,H,M,C\) with some extra symbols
  \(X,U,B,Z,\ast\). All symbols occurring on the same line to the left of
  the first white space are considered to be pairwise equivalent. The symbol
  after the first white can be considered to be a comment. The symbols on
  the last line are considered to be wildcards. Again they are replaced
  by a unique character. This is the given parameter wildcard, if
  wildcard $>$ 0.
*/

/*
  The following function reads a file in the format as explained above.
  \texttt{mapfile} is the input filename. \texttt{fpin} is the corresponding
  file name. If the argument
  \texttt{wildcard} is larger than 0, then the characters  in the last
  line of the symbol mapping file are mapped to \texttt{wildcard}. Otherwise,
  they are mapped to the value one smaller than the line number they appear
  in (counting from 1). The result of the parsing is stored in
  \texttt{alpha}.
*/

static int read_symbolmap_from_lines(GtAlphabet *alpha,
                                     const GtStr *mapfile,
                                     const GtStrArray *lines,
                                     GtError *err)
{
  char cc;
  unsigned int cnum, allocateddomainsize = 0;
  unsigned long linecount, column;
  bool blankfound, ignore, preamble = true, haserr = false;
  const char *currentline;
  GtUchar chartoshow;

  gt_error_check(err);
  alpha->domainsize = alpha->mapsize = alpha->mappedwildcards = 0;
  for (cnum=0; cnum<=(unsigned int) GT_MAXALPHABETCHARACTER; cnum++)
  {
    alpha->symbolmap[cnum] = (GtUchar) UNDEFCHAR;
  }
  alpha->mapdomain = NULL;
  alpha->characters = gt_malloc(sizeof (GtUchar) *
                                (gt_str_array_size(lines)-1));
  for (linecount = 0; linecount < gt_str_array_size(lines); linecount++)
  {
    currentline = gt_str_array_get(lines,linecount);
    ignore = false;
    if (currentline != NULL && currentline[0] != '\0')
    {
      if (preamble)
      {
        if (LINE(0) == (GtUchar) '#')
        {
          ignore = true;
        } else
        {
          preamble = false;
        }
      }
      if (!ignore)
      {
        blankfound = false;
        for (column=0; LINE(column) != '\0'; column++)
        { /* for all chars in line */
          cc = LINE(column);
          if (ispunct((int) cc) || isalnum((int) cc))
          {
            if (alpha->symbolmap[(unsigned int) cc] != (GtUchar) UNDEFCHAR)
            {
              gt_error_set(err,"cannot map symbol '%c' to %u: "
                            "it is already mapped to %u",
                             cc,
                             alpha->mapsize,
                             (unsigned int) alpha->
                                            symbolmap[(unsigned int) cc]);
              haserr = true;
              break;
            }
            /* get same value */
            alpha->symbolmap[(unsigned int) cc] = (GtUchar) alpha->mapsize;
            if (alpha->domainsize >= allocateddomainsize)
            {
              allocateddomainsize += 8;
              alpha->mapdomain = gt_realloc(alpha->mapdomain,
                                        sizeof (GtUchar) * allocateddomainsize);
            }
            gt_assert(alpha->mapdomain != NULL);
            alpha->mapdomain[alpha->domainsize++] = (GtUchar) cc;
          } else
          {
            if (cc == (GtUchar) ' ')    /* first blank in line found */
            {
              blankfound = true;
              /*@innerbreak@*/ break;
            }
            gt_error_set(err,
                          "illegal character '%c' in line %lu of mapfile %s",
                          cc,linecount,gt_str_get(mapfile));
            haserr = true;
            break;
          }
        }
        if (haserr)
        {
          break;
        }
        if (blankfound)
        {
          if (isspace((int) LINE(column+1)))
          {
            gt_error_set(err,"illegal character '%c' at the end of "
                          "line %lu in mapfile %s",
                          LINE(column+1),linecount,gt_str_get(mapfile));
            haserr  = true;
            break;
          }
          /* use next character to display character */
          chartoshow = (GtUchar) LINE(column+1);
        } else
        {
          /* use first character of line to display character */
          chartoshow = (GtUchar) LINE(0);
        }
        if (linecount == gt_str_array_size(lines)-1)
        {
          alpha->wildcardshow = chartoshow;
        } else
        {
          alpha->characters[alpha->mapsize] = chartoshow;
        }
        alpha->mapsize++;
      }
    }
  }
  if (!haserr)
  {
    for (cnum=0;cnum<=(unsigned int) GT_MAXALPHABETCHARACTER; cnum++)
    {
      if (alpha->symbolmap[cnum] == (GtUchar) (alpha->mapsize - 1))
      {
        alpha->symbolmap[cnum] = (GtUchar) WILDCARD;
        alpha->mappedwildcards++;
      }
    }
  }
  /* there are mapsize-1 characters plus wildcard plus separator.
     hence there are mapsize+1 symbols in the range 0..mapsize.
     that is, mapsize is the largest symbol and we obtain */
  alpha->bitspersymbol = gt_determinebitspervalue((uint64_t) alpha->mapsize);
  return haserr ? -1 : 0;
}

/*
  The following function reads in a symbol map.
  \texttt{mapfile} is the input filename.
  If the argument
  \texttt{wildcard} is larger than 0, then the characters in the last
  line of the symbol mapping file are mapped to \texttt{wildcard}. Otherwise,
  they are mapped to \(i-1\) if they appear on line number \(i\)
  (counting from 1). The result of the parsing is stored in
  \texttt{alpha}.
*/

static int read_symbolmap(GtAlphabet *alpha,const GtStr *mapfile,GtError *err)
{
  bool haserr = false;
  GtStrArray *lines;

  gt_error_check(err);
  lines = gt_str_array_new_file(gt_str_get(mapfile));
  gt_assert(lines != NULL);
  if (read_symbolmap_from_lines(alpha,mapfile,lines,err) != 0)
  {
    haserr = true;
  }
  gt_str_array_delete(lines);
  return haserr ? -1 : 0;
}

static void assign_dna_symbolmap(GtUchar *symbolmap)
{
  unsigned int cnum;

  for (cnum=0; cnum<=(unsigned int) GT_MAXALPHABETCHARACTER; cnum++)
  {
    symbolmap[cnum] = (GtUchar) UNDEFCHAR;
  }
  symbolmap[(unsigned int) 'a'] = (GtUchar) 0;
  symbolmap[(unsigned int) 'A'] = (GtUchar) 0;
  symbolmap[(unsigned int) 'c'] = (GtUchar) 1;
  symbolmap[(unsigned int) 'C'] = (GtUchar) 1;
  symbolmap[(unsigned int) 'g'] = (GtUchar) 2;
  symbolmap[(unsigned int) 'G'] = (GtUchar) 2;
  symbolmap[(unsigned int) 't'] = (GtUchar) 3;
  symbolmap[(unsigned int) 'T'] = (GtUchar) 3;
  symbolmap[(unsigned int) 'u'] = (GtUchar) 3;
  symbolmap[(unsigned int) 'U'] = (GtUchar) 3;
  for (cnum=0; DNAWILDCARDS[cnum] != '\0'; cnum++)
  {
    symbolmap[(unsigned int) DNAWILDCARDS[cnum]] = (GtUchar) WILDCARD;
  }
}

GtAlphabet *gt_alphabet_clone(const GtAlphabet *alphabet)
{
  unsigned int i;
  GtAlphabet *newalpha;

  newalpha = gt_malloc(sizeof *newalpha);
  newalpha->domainsize = alphabet->domainsize;
  newalpha->mapsize = alphabet->mapsize;
  newalpha->mappedwildcards = alphabet->mappedwildcards;
  newalpha->wildcardshow = alphabet->wildcardshow;
  newalpha->reference_count = 0;
  for (i=0; i<=(unsigned int) GT_MAXALPHABETCHARACTER; i++)
  {
    newalpha->symbolmap[i] = alphabet->symbolmap[i];
  }
  for (i=0; i<newalpha->mapsize; i++)
  {
    newalpha->characters[i] = alphabet->characters[i];
  }
  for (i=0; i<newalpha->domainsize; i++)
  {
    newalpha->mapdomain[i] = alphabet->mapdomain[i];
  }
  return newalpha;
}

GtAlphabet* gt_alphabet_ref(GtAlphabet *alphabet)
{
  gt_assert(alphabet);
  alphabet->reference_count++;
  return alphabet;
}

void gt_alphabet_delete(GtAlphabet *alphabet)
{
  if (!alphabet) return;
  if (alphabet->reference_count) {
    alphabet->reference_count--;
    return;
  }
  gt_free(alphabet->mapdomain);
  gt_free(alphabet->characters);
  gt_free(alphabet);
}

void gt_alphabet_add_mapping(GtAlphabet *alphabet, const char *characters)
{
  size_t i, num_of_characters;
  gt_assert(alphabet && characters);
  num_of_characters = strlen(characters);
  gt_assert(num_of_characters);
  alphabet->mapdomain = gt_realloc(alphabet->mapdomain,
                                   (size_t) alphabet->domainsize
                                      + num_of_characters);
  memcpy(alphabet->mapdomain + alphabet->domainsize, characters,
         num_of_characters);
  alphabet->domainsize += num_of_characters;
  alphabet->symbolmap[(int) characters[0]] = (GtUchar) alphabet->mapsize;
  alphabet->characters = gt_realloc(alphabet->characters,
                                    (size_t) alphabet->domainsize);
  alphabet->characters[alphabet->mapsize] = (GtUchar) characters[0];
  for (i = 0; i < num_of_characters; i++)
    alphabet->symbolmap[(int) characters[i]] = (GtUchar) alphabet->mapsize;
  alphabet->mapsize++;
  alphabet->bitspersymbol = gt_determinebitspervalue((uint64_t)
                                                       alphabet->mapsize);
}

void gt_alphabet_add_wildcard(GtAlphabet *alphabet, char wildcard)
{
  gt_assert(alphabet);
  alphabet->mapdomain = gt_realloc(alphabet->mapdomain,
                                   (size_t) alphabet->domainsize + 1);
  alphabet->mapdomain[alphabet->domainsize] = (GtUchar) wildcard;
  alphabet->domainsize++;
  alphabet->symbolmap[(int) wildcard] = (GtUchar) WILDCARD;
  if (alphabet->wildcardshow == (GtUchar) UNDEFCHAR) {
    alphabet->wildcardshow = (GtUchar) wildcard;
    alphabet->mapsize++;
  }
  alphabet->mappedwildcards++;
}

/*
  The following function initializes the alphabet \texttt{alpha}
  in the same way as \texttt{read_symbolmap}, if it would be
  applied to a map file with the following content:
  \begin{alltt}
  aA
  cC
  gG
  tTuU
  nsywrkvbdhmNSYWRKVBDHM
  \end{alltt}
*/

static void assign_dna_alphabet(GtAlphabet *alpha)
{
  alpha->wildcardshow = (GtUchar) DNAWILDCARDS[0];
  alpha->mappedwildcards = (unsigned int) strlen(DNAWILDCARDS);
  alpha->domainsize = (unsigned int) strlen(DNAALPHABETDOMAIN);
  alpha->bitspersymbol = 3U; /* as we have to represent 4 + 2 characters */
  alpha->mapdomain = gt_malloc(sizeof (GtUchar) * alpha->domainsize);
  memcpy(alpha->mapdomain,(GtUchar *) DNAALPHABETDOMAIN,
         (size_t) alpha->domainsize);
  alpha->mapsize = MAPSIZEDNA;
  alpha->characters = gt_calloc((size_t) UCHAR_MAX+1, sizeof (GtUchar));
  memcpy(alpha->characters,"acgt",(size_t) (MAPSIZEDNA-1));
  alpha->characters[WILDCARD] = (GtUchar) DNAWILDCARDS[0];
  alpha->characters[MAPSIZEDNA-1] = (GtUchar) DNAWILDCARDS[0];
  assign_dna_symbolmap(alpha->symbolmap);
}

static void assignproteinsymbolmap(GtUchar *symbolmap)
{
  unsigned int cnum;

  for (cnum=0; cnum<=(unsigned int) GT_MAXALPHABETCHARACTER; cnum++)
  {
    symbolmap[cnum] = (GtUchar) UNDEFCHAR;
  }
  for (cnum=0; PROTEINUPPERAMINOACIDS[cnum] != '\0'; cnum++)
  {
    symbolmap[(unsigned int) PROTEINUPPERAMINOACIDS[cnum]] = (GtUchar) cnum;
  }
  for (cnum=0; PROTEINWILDCARDS[cnum] != '\0'; cnum++)
  {
    symbolmap[(unsigned int) PROTEINWILDCARDS[cnum]] = (GtUchar) WILDCARD;
  }
}

/*
  The following function initializes the alphabet \texttt{alpha}
  in the same way as \texttt{read_symbolmap}, if it would be
  applied to a map file with the following content:
  \begin{alltt}
  L
  V
  I
  F
  K
  R
  E
  D
  A
  G
  S
  T
  N
  Q
  Y
  W
  P
  H
  M
  C
  XUBZO*-
  \end{alltt}
  If the argument \texttt{wildcard} is 0, then the wildcard characters
  in the last line are mapped to 20. Otherwise they are mapped to
  the character \texttt{WILDCARD}, as defined in \texttt{chardef.h}
*/

static void assign_protein_alphabet(GtAlphabet *alpha)
{
  alpha->wildcardshow = (GtUchar) PROTEINWILDCARDS[0];
  alpha->domainsize = (unsigned int) strlen(PROTEINALPHABETDOMAIN);
  alpha->mappedwildcards = (unsigned int) strlen(PROTEINWILDCARDS);
  alpha->bitspersymbol = 5U; /* as we have to represent 20 + 2 characters */
  alpha->mapdomain = gt_malloc(sizeof (GtUchar) * alpha->domainsize);
  memcpy(alpha->mapdomain,
         (GtUchar *) PROTEINALPHABETDOMAIN,(size_t) alpha->domainsize);
  alpha->mapsize = MAPSIZEPROTEIN;
  alpha->characters = gt_calloc((size_t) UCHAR_MAX+1, sizeof (GtUchar));
  memcpy(alpha->characters,PROTEINUPPERAMINOACIDS,(size_t) (MAPSIZEPROTEIN-1));
  alpha->characters[WILDCARD] = (GtUchar) PROTEINWILDCARDS[0];
  alpha->characters[MAPSIZEPROTEIN-1] = (GtUchar) PROTEINWILDCARDS[0];
  assignproteinsymbolmap(alpha->symbolmap);
}

static int assign_protein_or_dna_alphabet(GtAlphabet *alpha,
                                          const GtStrArray *filenametab,
                                          GtError *err)
{
  int retval;

  gt_error_check(err);
  retval = gt_files_guess_if_protein_sequences(filenametab,err);
  if (retval < 0)
  {
    return -1;
  }
  if (retval == 1)
  {
    assign_protein_alphabet(alpha);
  } else
  {
    assign_dna_alphabet(alpha);
  }
  return 0;
}

GtAlphabet* gt_alphabet_new(bool isdna,
                            bool isprotein,
                            const GtStr *smapfile,
                            const GtStrArray *filenametab,
                            GtError *err)
{
  GtAlphabet *alpha;
  bool haserr = false;

  gt_error_check(err);

  alpha = gt_malloc(sizeof *alpha);
  alpha->reference_count = 0;
  alpha->characters = NULL;
  alpha->mapdomain = NULL;
  if (isdna)
  {
    assign_dna_alphabet(alpha);
  } else
  {
    if (isprotein)
    {
      assign_protein_alphabet(alpha);
    } else
    {
      gt_error_check(err);
      if (gt_str_length(smapfile) > 0)
      {
        GtStr *transpath = NULL;

        if (!gt_file_exists(gt_str_get(smapfile)))
        {
          GtStr *prog;
          const char *progname = gt_error_get_progname(err);

          gt_assert(progname != NULL);
          prog = gt_str_new();
          gt_str_append_cstr_nt(prog, progname,
                                gt_cstr_length_up_to_char(progname, ' '));
          transpath = gt_get_gtdata_path(gt_str_get(prog), err);
          gt_str_delete(prog);
          gt_str_append_cstr(transpath, "/trans/");
          gt_str_append_cstr(transpath, gt_str_get(smapfile));
        }
        if (read_symbolmap(alpha,
                           transpath == NULL ? smapfile : transpath,
                           err) != 0)
        {
          haserr = true;
        }
        gt_str_delete(transpath);
      } else
      {
        if (assign_protein_or_dna_alphabet(alpha,filenametab,err) != 0)
        {
          haserr = true;
        }
      }
    }
  }
  if (haserr)
  {
    gt_alphabet_delete(alpha);
    return NULL;
  }
  return alpha;
}

GtAlphabet* gt_alphabet_new_dna(void)
{
  GtAlphabet *a;
  a = gt_alphabet_new(true, false, NULL, NULL, NULL);
  gt_assert(a);
  return a;
}

GtAlphabet* gt_alphabet_new_protein(void)
{
  GtAlphabet *a;
  a = gt_alphabet_new(false, true, NULL, NULL, NULL);
  return a;
}

GtAlphabet* gt_alphabet_new_empty(void)
{
  GtAlphabet *a = gt_malloc(sizeof *a);
  a->domainsize = 0;
  a->mapsize = 0;
  a->mappedwildcards = 0;
  a->bitspersymbol = 0;
  a->reference_count = 0;
  a->wildcardshow = (GtUchar) UNDEFCHAR;
  memset(a->symbolmap, (int) UNDEFCHAR, (size_t) GT_MAXALPHABETCHARACTER+1);
  a->mapdomain = NULL;
  a->characters = NULL;
  return a;
}

GtAlphabet* gt_alphabet_guess(const char *sequence, unsigned long seqlen)
{
  unsigned long i;
  gt_assert(sequence && seqlen);
  for (i = 0;
       i < seqlen && i < (unsigned long) ALPHABET_GUESS_MAX_LENGTH;
        i++) {
    if (strchr(ALPHABET_GUESS_PROTEIN_CHARS, sequence[i]) != NULL)
      return gt_alphabet_new_protein();
  }
  return gt_alphabet_new_dna();
}

const GtUchar* gt_alphabet_symbolmap(const GtAlphabet *alphabet)
{
  gt_assert(alphabet);
  return alphabet->symbolmap;
}

unsigned int gt_alphabet_num_of_chars(const GtAlphabet *alphabet)
{
  gt_assert(alphabet);
  return alphabet->mapsize - 1;
}

unsigned int gt_alphabet_size(const GtAlphabet *alphabet)
{
  gt_assert(alphabet);
  return alphabet->mapsize;
}

const GtUchar* gt_alphabet_characters(const GtAlphabet *alphabet)
{
  gt_assert(alphabet);
  return alphabet->characters;
}

GtUchar gt_alphabet_wildcard_show(const GtAlphabet *alphabet)
{
  gt_assert(alphabet);
  return alphabet->wildcardshow;
}

unsigned int gt_alphabet_bits_per_symbol(const GtAlphabet *alphabet)
{
  gt_assert(alphabet);
  return alphabet->bitspersymbol;
}

void gt_alphabet_output(const GtAlphabet *alphabet, FILE *fpout)
{
  GtUchar chartoshow, currentcc, previouscc = 0, firstinline = 0;
  unsigned int cnum, linenum = 0;
  bool afternewline = true;
  gt_assert(alphabet && fpout);
  for (cnum=0; cnum < alphabet->domainsize; cnum++)
  {
    currentcc = alphabet->mapdomain[cnum];
    if (cnum > 0)
    {
      if (alphabet->symbolmap[currentcc] != alphabet->symbolmap[previouscc])
      {
        if (linenum < alphabet->mapsize-1)
        {
          chartoshow = alphabet->characters[linenum];
        } else
        {
          chartoshow = alphabet->wildcardshow;
        }
        if (firstinline != chartoshow)
        {
          fprintf(fpout," %c",(int) chartoshow);
        }
        gt_xfputc('\n',fpout);
        afternewline = true;
        linenum++;
      } else
      {
        afternewline = false;
      }
    }
    gt_xfputc((int) currentcc,fpout);
    if (afternewline)
    {
      firstinline = currentcc;
    }
    previouscc = currentcc;
  }
  if (linenum < alphabet->mapsize-1)
  {
    chartoshow = alphabet->characters[linenum];
  } else
  {
    chartoshow = alphabet->wildcardshow;
  }
  if (firstinline != chartoshow)
  {
    fprintf(fpout," %c",(int) chartoshow);
  }
  gt_xfputc((int) '\n',fpout);
}

void gt_alphabet_decode_seq_to_fp(const GtAlphabet *alphabet, FILE *fpout,
                                  const GtUchar *src, unsigned long len)
{
  unsigned long i;
  const GtUchar *characters;
  gt_assert(fpout);
  if (alphabet == NULL)
  {
    characters = (const GtUchar *) "acgt";
  } else
  {
    characters = alphabet->characters;
  }
  for (i = 0; i < len; i++)
  {
    gt_xfputc((int) characters[(int) src[i]],fpout);
  }
}

void gt_alphabet_printf_symbolstring(const GtAlphabet *alphabet,
                                     const GtUchar *w, unsigned long len)
{
  gt_alphabet_decode_seq_to_fp(alphabet, stdout, w, len);
}

static char converttoprettysymbol(const GtAlphabet *alphabet,
                                  GtUchar currentchar)
{
  char ret = '\0';
  if (alphabet == NULL)
  {
    ret = (char) currentchar;
  } else
  {

    if (currentchar == (GtUchar) WILDCARD)
    {
      ret = (char) alphabet->wildcardshow;
    } else
    {
      if (currentchar != (GtUchar) SEPARATOR)
      {
        gt_assert((unsigned int) currentchar < alphabet->mapsize-1);
        ret = (char) alphabet->characters[(int) currentchar];
      }
    }
  }
  gt_assert(ret != '\0');
  return ret;
}

void gt_alphabet_decode_seq_to_cstr(const GtAlphabet *alphabet, char *dest,
                                      const GtUchar *src, unsigned long len)
{
  unsigned long i;

  for (i = 0; i < len; i++)
  {
    dest[i] = converttoprettysymbol(alphabet, (GtUchar) src[i]);
  }
  dest[len] = '\0';
}

GtStr* gt_alphabet_decode_seq_to_str(const GtAlphabet *alphabet,
                                     const GtUchar *src,
                                     unsigned long len)
{
  char *buffer;
  GtStr *ret;
  gt_assert(alphabet && src);
  buffer = gt_malloc(sizeof (char) * len+1);
  gt_alphabet_decode_seq_to_cstr(alphabet, buffer, src, len);
  ret = gt_str_new_cstr(buffer);
  gt_free(buffer);
  return ret;
}

void gt_alphabet_echo_pretty_symbol(const GtAlphabet *alphabet, FILE *fpout,
                                    GtUchar currentchar)
{
  gt_xfputc((int) converttoprettysymbol(alphabet, currentchar), fpout);
}

GtUchar gt_alphabet_pretty_symbol(const GtAlphabet *alphabet,
                                  unsigned int currentchar)
{
  gt_assert(currentchar <= UCHAR_MAX);
  return (GtUchar) converttoprettysymbol(alphabet, (GtUchar) currentchar);
}

static unsigned int removelowercaseproteinchars(GtUchar *domainbuf,
                                                const GtAlphabet *alpha)
{
  unsigned int i, j = 0;

  for (i=0; i< alpha->domainsize - alpha->mappedwildcards; i++)
  {
    if (isalnum((int) alpha->mapdomain[i]) &&
        isupper((int) alpha->mapdomain[i]))
    {
      domainbuf[j++] = alpha->mapdomain[i];
    }
  }
  return j;
}

#define UNCAST(X) (*((const GtUchar *) (X)))

static int comparechar(const void *a,const void *b)
{
  if (UNCAST(a) < UNCAST(b))
  {
    return -1;
  }
  if (UNCAST(a) > UNCAST(b))
  {
    return 1;
  }
  return 0;
}

bool gt_alphabet_is_protein(const GtAlphabet *alphabet)
{
  GtAlphabet proteinalphabet;
  unsigned int i, reduceddomainsize1, reduceddomainsize2;
  bool isprot = false;
  GtUchar domainbuf1[GT_MAXALPHABETCHARACTER+1],
        domainbuf2[GT_MAXALPHABETCHARACTER+1];

  reduceddomainsize1 = removelowercaseproteinchars(&domainbuf1[0],alphabet);
  assign_protein_alphabet(&proteinalphabet);
  reduceddomainsize2 = removelowercaseproteinchars(&domainbuf2[0],
                                                   &proteinalphabet);
  if (reduceddomainsize1 == reduceddomainsize2)
  {
    qsort(&domainbuf1[0],(size_t) reduceddomainsize1,sizeof (char),comparechar);
    qsort(&domainbuf2[0],(size_t) reduceddomainsize2,sizeof (char),comparechar);
    for (i=0; i < reduceddomainsize2; i++)
    {
      if (domainbuf1[i] != domainbuf2[i])
      {
        isprot = false;
        break;
      }
    }
    isprot = true;
  } else
  {
    isprot = false;
  }
  gt_free(proteinalphabet.mapdomain);
  gt_free(proteinalphabet.characters);
  return isprot;
}

static bool check_symbolmap(const GtUchar *testsymbolmap,
                            const GtUchar *verifiedsymbolmap,
                            const char *testcharacters)
{
  unsigned int i;
  GtUchar cc1, cc2;

  for (i=0; testcharacters[i] != '\0'; i++)
  {
    cc1 = (GtUchar) testcharacters[i];
    if (isupper((int) cc1))
    {
      cc2 = (GtUchar) tolower((int) cc1);
    } else
    {
      gt_assert(islower((int) cc1));
      cc2 = (GtUchar) toupper((int) cc1);
    }
    if (testsymbolmap[cc1] != verifiedsymbolmap[cc1] &&
        testsymbolmap[cc2] != verifiedsymbolmap[cc2])
    {
      return false;
    }
  }
  return true;
}

bool gt_alphabet_is_dna(const GtAlphabet *alphabet)
{
  if (gt_alphabet_is_protein(alphabet))
  {
    return false;
  }
  if (alphabet->mapsize == MAPSIZEDNA)
  {
    GtUchar dnasymbolmap[GT_MAXALPHABETCHARACTER+1];

    assign_dna_symbolmap(&dnasymbolmap[0]);
    if (check_symbolmap(alphabet->symbolmap,&dnasymbolmap[0],"acgt"))
    {
      return true;
    }
  }
  return false;
}

bool gt_alphabet_valid_input(const GtAlphabet *alphabet, char c)
{
  if (alphabet->symbolmap[(int) c] != (GtUchar) UNDEFCHAR)
    return true;
  return false;
}

GtUchar gt_alphabet_encode(const GtAlphabet *alphabet, char c)
{
  gt_assert(alphabet);
  gt_assert(alphabet->symbolmap[(int) c] != (GtUchar) UNDEFCHAR);
  return alphabet->symbolmap[(int) c];
}

char gt_alphabet_decode(const GtAlphabet *alphabet, GtUchar c)
{
  gt_assert(alphabet);
  if (c == (GtUchar) alphabet->mapsize - 1)
    return (char) alphabet->wildcardshow;
  return converttoprettysymbol(alphabet, c);
}

void gt_alphabet_encode_seq(const GtAlphabet *alphabet, GtUchar *out,
                            const char *in, unsigned long length)
{
  unsigned long i;
  gt_assert(alphabet && out && in);
  for (i = 0; i < length; i++) {
    gt_assert(alphabet->symbolmap[(int) in[i]] != (GtUchar) UNDEFCHAR);
    out[i] = alphabet->symbolmap[(int) in[i]];
  }
}

GtAlphabet *gt_alphabet_new_from_file(const char *indexname, GtError *err)
{
  GtStr *tmpfilename;
  bool haserr = false;
  GtAlphabet *alpha;

  gt_error_check(err);
  tmpfilename = gt_str_new_cstr(indexname);
  gt_str_append_cstr(tmpfilename,GT_ALPHABETFILESUFFIX);
  alpha = gt_alphabet_new(false,false,tmpfilename,NULL,err);
  if (alpha == NULL)
  {
    haserr = true;
  }
  gt_str_delete(tmpfilename);
  if (haserr)
  {
    gt_alphabet_delete((GtAlphabet*) alpha);
    return NULL;
  }
  return alpha;
}

int gt_alphabet_to_file(const GtAlphabet *alpha, const char *indexname,
                        GtError *err)
{
  FILE *al1fp;
  bool haserr = false;

  gt_error_check(err);
  al1fp = gt_fa_fopen_with_suffix(indexname,GT_ALPHABETFILESUFFIX,"wb",err);
  if (al1fp == NULL)
  {
    haserr = true;
  }
  if (!haserr)
  {
    gt_alphabet_output(alpha,al1fp);
    gt_fa_xfclose(al1fp);
  }
  return haserr ? -1 : 0;
}
