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

#include <stdio.h>
#include <ctype.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdbool.h>
#include <limits.h>
#include <errno.h>
#include "libgtcore/chardef.h"
#include "libgtcore/cstr.h"
#include "libgtcore/error.h"
#include "libgtcore/fileutils.h"
#include "libgtcore/gtdatapath.h"
#include "libgtcore/str.h"
#include "libgtcore/strarray.h"
#include "libgtcore/symboldef.h"
#include "spacedef.h"
#include "qsorttype.h"
#include "alphadef.h"

#include "guessprot.pr"

 struct Alphabet                /* initial blank prevents select by skproto */
{
  unsigned int domainsize,           /* size of domain of symbolmap */
               mapsize,              /* size of image of map, i.e. */
                                     /* mapping to [0..mapsize-1] */
               mappedwildcards;      /* number of mapped wildcards */
  Uchar wildcardshow,
        symbolmap[UCHAR_MAX+1],      /* mapping of the symbols */
        *mapdomain,                  /* list of characters mapped */
        *characters;                 /* array of characters to show */
};

/*EE
  This file implements the datatype \texttt{alphabet}.
*/

/*
  Some constants for the standard alphabet used. The name says it all.
*/

#define DNABASES                     "aAcCgGtTuU"
#define DNAWILDCARDS                 "nsywrkvbdhmNSYWRKVBDHM"
#define MAPSIZEDNA                   5U
#define DNAALPHABETDOMAIN            DNABASES DNAWILDCARDS
#define PROTEINUPPERAMINOACIDS       "LVIFKREDAGSTNQYWPHMC"
#define MAPSIZEPROTEIN               21U
#define PROTEINWILDCARDS             "XUBZ*-"
#define PROTEINALPHABETDOMAIN        PROTEINUPPERAMINOACIDS PROTEINWILDCARDS

/*
  We use the following macro to access the \texttt{I}-th character of
  a line.
*/

#define LINE(I)          currentline[I]

/*EE
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

static int readsymbolmapfromlines(Alphabet *alpha,
                                  const Str *mapfile,
                                  const StrArray *lines,
                                  Error *err)
{
  char cc;
  unsigned int cnum, allocateddomainsize = 0;
  unsigned long linecount, column;
  bool blankfound, ignore, preamble = true, haserr = false;
  const char *currentline;
  Uchar chartoshow;

  error_check(err);
  alpha->domainsize = alpha->mapsize = alpha->mappedwildcards = 0;
  for (cnum=0; cnum<=UCHAR_MAX; cnum++)
  {
    alpha->symbolmap[cnum] = (Uchar) UNDEFCHAR;
  }
  alpha->mapdomain = NULL;
  ALLOCASSIGNSPACE(alpha->characters,NULL,Uchar,strarray_size(lines)-1);
  for (linecount = 0; linecount < strarray_size(lines); linecount++)
  {
    currentline = strarray_get(lines,linecount);
    ignore = false;
    if (currentline != NULL && currentline[0] != '\0')
    {
      if (preamble)
      {
        if (LINE(0) == (Uchar) '#')
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
            if (alpha->symbolmap[(unsigned int) cc] != (Uchar) UNDEFCHAR)
            {
              error_set(err,"cannot map symbol '%c' to %u: "
                            "it is already mapped to %u",
                             cc,
                             alpha->mapsize,
                             (unsigned int) alpha->
                                            symbolmap[(unsigned int) cc]);
              haserr = true;
              break;
            }
            /* get same value */
            alpha->symbolmap[(unsigned int) cc] = (Uchar) alpha->mapsize;
            if (alpha->domainsize >= allocateddomainsize)
            {
              allocateddomainsize += 8;
              ALLOCASSIGNSPACE(alpha->mapdomain,alpha->mapdomain,Uchar,
                               allocateddomainsize);
            }
            assert(alpha->mapdomain != NULL);
            alpha->mapdomain[alpha->domainsize++] = (Uchar) cc;
          } else
          {
            if (cc == (Uchar) ' ')    /* first blank in line found */
            {
              blankfound = true;
              /*@innerbreak@*/ break;
            }
            error_set(err,
                          "illegal character '%c' in line %lu of mapfile %s",
                          cc,linecount,str_get(mapfile));
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
            error_set(err,"illegal character '%c' at the end of "
                          "line %lu in mapfile %s",
                          LINE(column+1),linecount,str_get(mapfile));
            haserr  = true;
            break;
          }
          /* use next character to display character */
          chartoshow = (Uchar) LINE(column+1);
        } else
        {
          /* use first character of line to display character */
          chartoshow = (Uchar) LINE(0);
        }
        if (linecount == strarray_size(lines)-1)
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
    for (cnum=0;cnum<=UCHAR_MAX; cnum++)
    {
      if (alpha->symbolmap[cnum] == (Uchar) (alpha->mapsize - 1))
      {
        alpha->symbolmap[cnum] = (Uchar) WILDCARD;
        alpha->mappedwildcards++;
      }
    }
  }
  return haserr ? -1 : 0;
}

/*EE
  The following function reads in a symbol map.
  \texttt{mapfile} is the input filename.
  If the argument
  \texttt{wildcard} is larger than 0, then the characters in the last
  line of the symbol mapping file are mapped to \texttt{wildcard}. Otherwise,
  they are mapped to \(i-1\) if they appear on line number \(i\)
  (counting from 1). The result of the parsing is stored in
  \texttt{alpha}.
*/

static int readsymbolmap(Alphabet *alpha,const Str *mapfile,Error *err)
{
  bool haserr = false;
  StrArray *lines;

  error_check(err);
  lines = strarray_new_file(str_get(mapfile));
  assert(lines != NULL);
  if (readsymbolmapfromlines(alpha,mapfile,lines,err) != 0)
  {
    haserr = true;
  }
  strarray_delete(lines);
  return haserr ? -1 : 0;
}

static void assignDNAsymbolmap(Uchar *symbolmap)
{
  unsigned int cnum;

  for (cnum=0; cnum<=(unsigned int) UCHAR_MAX; cnum++)
  {
    symbolmap[cnum] = (Uchar) UNDEFCHAR;
  }
  symbolmap[(unsigned int) 'a'] = (Uchar) 0;
  symbolmap[(unsigned int) 'A'] = (Uchar) 0;
  symbolmap[(unsigned int) 'c'] = (Uchar) 1;
  symbolmap[(unsigned int) 'C'] = (Uchar) 1;
  symbolmap[(unsigned int) 'g'] = (Uchar) 2;
  symbolmap[(unsigned int) 'G'] = (Uchar) 2;
  symbolmap[(unsigned int) 't'] = (Uchar) 3;
  symbolmap[(unsigned int) 'T'] = (Uchar) 3;
  symbolmap[(unsigned int) 'u'] = (Uchar) 3;
  symbolmap[(unsigned int) 'U'] = (Uchar) 3;
  for (cnum=0; DNAWILDCARDS[cnum] != '\0'; cnum++)
  {
    symbolmap[(unsigned int) DNAWILDCARDS[cnum]] = (Uchar) WILDCARD;
  }
}

/*EE
  The following function initializes the alphabet \texttt{alpha}
  in the same way as \texttt{readsymbolmap}, if it would be
  applied to a map file with the following content:
  \begin{alltt}
  aA
  cC
  gG
  tTuU
  nsywrkvbdhmNSYWRKVBDHM
  \end{alltt}
*/

static void assignDNAalphabet(Alphabet *alpha)
{
  alpha->wildcardshow = (Uchar) DNAWILDCARDS[0];
  alpha->mappedwildcards = (unsigned int) strlen(DNAWILDCARDS);
  alpha->domainsize = (unsigned int) strlen(DNAALPHABETDOMAIN);
  ALLOCASSIGNSPACE(alpha->mapdomain,NULL,Uchar,alpha->domainsize);
  memcpy(alpha->mapdomain,(Uchar *) DNAALPHABETDOMAIN,
         (size_t) alpha->domainsize);
  alpha->mapsize = MAPSIZEDNA;
  ALLOCASSIGNSPACE(alpha->characters,NULL,char,MAPSIZEDNA-1);
  memcpy(alpha->characters,"acgt",(size_t) (MAPSIZEDNA-1));
  assignDNAsymbolmap(alpha->symbolmap);
}

static void assignproteinsymbolmap(Uchar *symbolmap)
{
  unsigned int cnum;

  for (cnum=0; cnum<=(unsigned int) UCHAR_MAX; cnum++)
  {
    symbolmap[cnum] = (Uchar) UNDEFCHAR;
  }
  for (cnum=0; PROTEINUPPERAMINOACIDS[cnum] != '\0'; cnum++)
  {
    symbolmap[(unsigned int) PROTEINUPPERAMINOACIDS[cnum]] = (Uchar) cnum;
  }
  for (cnum=0; PROTEINWILDCARDS[cnum] != '\0'; cnum++)
  {
    symbolmap[(unsigned int) PROTEINWILDCARDS[cnum]] = (Uchar) WILDCARD;
  }
}

/*EE
  The following function initializes the alphabet \texttt{alpha}
  in the same way as \texttt{readsymbolmap}, if it would be
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
  XUBZ*-
  \end{alltt}
  If the argument \texttt{wildcard} is 0, then the wildcard characters
  in the last line are mapped to 20. Otherwise they are mapped to
  the character \texttt{WILDCARD}, as defined in \texttt{chardef.h}
*/

static void assignProteinalphabet(Alphabet *alpha)
{
  alpha->wildcardshow = (Uchar) PROTEINWILDCARDS[0];
  alpha->domainsize = (unsigned int) strlen(PROTEINALPHABETDOMAIN);
  alpha->mappedwildcards = (unsigned int) strlen(PROTEINWILDCARDS);
  ALLOCASSIGNSPACE(alpha->mapdomain,NULL,Uchar,alpha->domainsize);
  memcpy(alpha->mapdomain,
         (Uchar *) PROTEINALPHABETDOMAIN,(size_t) alpha->domainsize);
  alpha->mapsize = MAPSIZEPROTEIN;
  ALLOCASSIGNSPACE(alpha->characters,NULL,char,MAPSIZEPROTEIN-1);
  memcpy(alpha->characters,PROTEINUPPERAMINOACIDS,(size_t) MAPSIZEPROTEIN-1);
  assignproteinsymbolmap(alpha->symbolmap);
}

static int assignProteinorDNAalphabet(Alphabet *alpha,
                                      const StrArray *filenametab,Error *err)
{
  int retval;

  error_check(err);
  retval = guessifproteinsequencestream(filenametab,err);
  if (retval < 0)
  {
    return -1;
  }
  if (retval == 1)
  {
    assignProteinalphabet(alpha);
  } else
  {
    assignDNAalphabet(alpha);
  }
  return 0;
}

void freeAlphabet(Alphabet **alpha)
{
  FREESPACE((*alpha)->mapdomain);
  FREESPACE((*alpha)->characters);
  FREESPACE(*alpha);
}

/*@null@*/ Alphabet *assigninputalphabet(bool isdna,
                                         bool isprotein,
                                         const Str *smapfile,
                                         const StrArray *filenametab,
                                         Error *err)
{
  Alphabet *alpha;
  bool haserr = false;

  error_check(err);
  ALLOCASSIGNSPACE(alpha,NULL,Alphabet,(size_t) 1);
  alpha->characters = NULL;
  alpha->mapdomain = NULL;
  if (isdna)
  {
    assignDNAalphabet(alpha);
  } else
  {
    if (isprotein)
    {
      assignProteinalphabet(alpha);
    } else
    {
      if (str_length(smapfile) > 0)
      {
        Str *transpath = NULL;

        if (!file_exists(str_get(smapfile)))
        {
          Str *prog;
          const char *progname = error_get_progname(err);

          assert(progname != NULL);
          prog = str_new();
          str_append_cstr_nt(prog, progname,
                             cstr_length_up_to_char(progname, ' '));
          transpath = gtdata_get_path(str_get(prog), err);
          str_delete(prog);
          str_append_cstr(transpath, "/trans/");
          str_append_cstr(transpath, str_get(smapfile));
        }
        if (readsymbolmap(alpha,
                          transpath == NULL ? smapfile : transpath,
                          err) != 0)
        {
          haserr = true;
        }
        str_delete(transpath);
      } else
      {
        if (assignProteinorDNAalphabet(alpha,filenametab,err) != 0)
        {
          haserr = true;
        }
      }
    }
  }
  if (haserr)
  {
    if (alpha != NULL)
    {
      freeAlphabet(&alpha);
    }
    return NULL;
  }
  return alpha;
}

const Uchar *getsymbolmapAlphabet(const Alphabet *alpha)
{
  return alpha->symbolmap;
}

unsigned int getnumofcharsAlphabet(const Alphabet *alpha)
{
  return alpha->mapsize-1;
}

unsigned int getmapsizeAlphabet(const Alphabet *alpha)
{
  return alpha->mapsize;
}

const Uchar *getcharactersAlphabet(const Alphabet *alpha)
{
  return alpha->characters;
}

void outputalphabet(FILE *fpout,const Alphabet *alpha)
{
  Uchar chartoshow, currentcc, previouscc = 0, firstinline = 0;
  unsigned int cnum, linenum = 0;
  bool afternewline = true;

  for (cnum=0; cnum < alpha->domainsize; cnum++)
  {
    currentcc = alpha->mapdomain[cnum];
    if (cnum > 0)
    {
      if (alpha->symbolmap[currentcc] != alpha->symbolmap[previouscc])
      {
        if (linenum < alpha->mapsize-1)
        {
          chartoshow = alpha->characters[linenum];
        } else
        {
          chartoshow = alpha->wildcardshow;
        }
        if (firstinline != chartoshow)
        {
          fprintf(fpout," %c",(int) chartoshow);
        }
        (void) putc('\n',fpout);
        afternewline = true;
        linenum++;
      } else
      {
        afternewline = false;
      }
    }
    (void) putc((int) currentcc,fpout);
    if (afternewline)
    {
      firstinline = currentcc;
    }
    previouscc = currentcc;
  }
  if (linenum < alpha->mapsize-1)
  {
    chartoshow = alpha->characters[linenum];
  } else
  {
    chartoshow = alpha->wildcardshow;
  }
  if (firstinline != chartoshow)
  {
    fprintf(fpout," %c",(int) chartoshow);
  }
  (void) putc((int) '\n',fpout);
}

/*
  Suppose the string \texttt{w} of length \texttt{wlen}
  was transformed according to the alphabet \texttt{alpha}.
  The following function shows each character in \texttt{w}
  as the characters specified in the transformation.
  The output goes to the given file pointer.
 */

void showsymbolstringgeneric(FILE *fpout,const Alphabet *alpha,
                             const Uchar *w,unsigned long wlen)
{
  unsigned long i;

  for (i = 0; i < wlen; i++)
  {
    (void) putc((int) alpha->characters[(int) w[i]],fpout);
  }
}

void echoprettysymbol(FILE *fpout,const Alphabet *alpha,Uchar currentchar)
{
  if (alpha == NULL)
  {
    (void) putc((int) currentchar,fpout);
  } else
  {
    if (currentchar == (Uchar) WILDCARD)
    {
      (void) putc((int) alpha->wildcardshow,fpout);
    } else
    {
      if (currentchar == (Uchar) SEPARATOR)
      {
        (void) fprintf(fpout,">\n");
      } else
      {
        assert((unsigned int) currentchar < alpha->mapsize-1);
        (void) putc((int) alpha->characters[(int) currentchar],fpout);
      }
    }
  }
}

Uchar getprettysymbol(const Alphabet *alpha,unsigned int currentchar)
{
   assert(currentchar < alpha->mapsize-1);
   return alpha->characters[currentchar];
}

/*
  The following function is a special case of the previous
  function showing the output on stdout.
*/

void showsymbolstring(const Alphabet *alpha,const Uchar *w,unsigned long wlen)
{
  showsymbolstringgeneric(stdout,alpha,w,wlen);
}

static unsigned int removelowercaseproteinchars(Uchar *domainbuf,
                                                const Alphabet *alpha)
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

#define UNCAST(X) (*((const Uchar *) (X)))

static Qsortcomparereturntype comparechar(const void *a,const void *b)
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

/*EE
  The following function checks if the given alphabet is the Protein
  alphabet with the aminoacids
  A, C, D, E, F, G, H, I, K, L, M, N, P, Q, R, S, T, V, W, Y written in
  lower or upper case.
*/

bool isproteinalphabet(const Alphabet *alpha)
{
  Alphabet proteinalphabet;
  unsigned int i, reduceddomainsize1, reduceddomainsize2;
  bool isprot = false;
  Uchar domainbuf1[UCHAR_MAX+1],
        domainbuf2[UCHAR_MAX+1];

  reduceddomainsize1 = removelowercaseproteinchars(&domainbuf1[0],alpha);
  assignProteinalphabet(&proteinalphabet);
  reduceddomainsize2 = removelowercaseproteinchars(&domainbuf2[0],
                                                   &proteinalphabet);
  if (reduceddomainsize1 == reduceddomainsize2)
  {
    qsort(&domainbuf1[0],(size_t) reduceddomainsize1,sizeof (char),
          (Qsortcomparefunction) comparechar);
    qsort(&domainbuf2[0],(size_t) reduceddomainsize2,sizeof (char),
          (Qsortcomparefunction) comparechar);
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
  FREESPACE(proteinalphabet.mapdomain);
  FREESPACE(proteinalphabet.characters);
  return isprot;
}

static bool checksymbolmap(const Uchar *testsymbolmap,
                           const Uchar *verifiedsymbolmap,
                           const char *testcharacters)
{
  unsigned int i;
  Uchar cc1, cc2;

  for (i=0; testcharacters[i] != '\0'; i++)
  {
    cc1 = (Uchar) testcharacters[i];
    if (isupper((int) cc1))
    {
      cc2 = (Uchar) tolower((int) cc1);
    } else
    {
      assert(islower((int) cc1));
      cc2 = (Uchar) toupper((int) cc1);
    }
    if (testsymbolmap[cc1] != verifiedsymbolmap[cc1] &&
        testsymbolmap[cc2] != verifiedsymbolmap[cc2])
    {
      return false;
    }
  }
  return true;
}

/*
  The following function checks if the given alphabet is the DNA
  alphabet with the bases A, C, G, T written in lower or upper case.
*/

bool isdnaalphabet(const Alphabet *alpha)
{
  if (isproteinalphabet(alpha))
  {
    return false;
  }
  if (alpha->mapsize == MAPSIZEDNA)
  {
    Uchar dnasymbolmap[UCHAR_MAX+1];

    assignDNAsymbolmap(&dnasymbolmap[0]);
    if (checksymbolmap(alpha->symbolmap,&dnasymbolmap[0],"acgt"))
    {
      return true;
    }
  }
  return false;
}
