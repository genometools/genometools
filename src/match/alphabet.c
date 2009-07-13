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
#include <stdbool.h>
#include <limits.h>
#include <errno.h>
#include "core/chardef.h"
#include "core/cstr.h"
#include "core/error.h"
#include "core/fileutils.h"
#include "core/gtdatapath.h"
#include "core/str.h"
#include "core/str_array.h"
#include "core/symboldef.h"
#include "core/mathsupport.h"
#include "spacedef.h"
#include "alphadef.h"

#include "guessprot.pr"

struct SfxAlphabet                /* initial blank prevents select by skproto */
{
  unsigned int domainsize,           /* size of domain of symbolmap */
               mapsize,              /* size of image of map, i.e. */
                                     /* mapping to [0..mapsize-1] */
               mappedwildcards,      /* number of mapped wildcards */
               bitspersymbol;        /* number of bits per symbol in
                                        bitspackedarray */
  GtUchar wildcardshow,
        symbolmap[MAXALPHABETCHARACTER+1], /* mapping of the symbols */
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
#define MAPSIZEDNA                   (DNAALPHASIZE+1U)
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

static int readsymbolmapfromlines(SfxAlphabet *alpha,
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
  for (cnum=0; cnum<=(unsigned int) MAXALPHABETCHARACTER; cnum++)
  {
    alpha->symbolmap[cnum] = (GtUchar) UNDEFCHAR;
  }
  alpha->mapdomain = NULL;
  ALLOCASSIGNSPACE(alpha->characters,NULL,GtUchar,gt_str_array_size(lines)-1);
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
              ALLOCASSIGNSPACE(alpha->mapdomain,alpha->mapdomain,GtUchar,
                               allocateddomainsize);
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
    for (cnum=0;cnum<=(unsigned int) MAXALPHABETCHARACTER; cnum++)
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

static int readsymbolmap(SfxAlphabet *alpha,const GtStr *mapfile,GtError *err)
{
  bool haserr = false;
  GtStrArray *lines;

  gt_error_check(err);
  lines = gt_str_array_new_file(gt_str_get(mapfile));
  gt_assert(lines != NULL);
  if (readsymbolmapfromlines(alpha,mapfile,lines,err) != 0)
  {
    haserr = true;
  }
  gt_str_array_delete(lines);
  return haserr ? -1 : 0;
}

static void assignDNAsymbolmap(GtUchar *symbolmap)
{
  unsigned int cnum;

  for (cnum=0; cnum<=(unsigned int) MAXALPHABETCHARACTER; cnum++)
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

/*EE
  The following function copies the alphabet \texttt{alpha2} into the alphabet
  \texttt{alpha1}.
*/

SfxAlphabet *gt_copyAlphabet(const SfxAlphabet *alpha2)
{
  unsigned int i;
  SfxAlphabet *alpha1;

  ALLOCASSIGNSPACE(alpha1,NULL,SfxAlphabet,(size_t) 1);
  alpha1->domainsize = alpha2->domainsize;
  alpha1->mapsize = alpha2->mapsize;
  alpha1->mappedwildcards = alpha2->mappedwildcards;
  alpha1->wildcardshow = alpha2->wildcardshow;
  for (i=0; i<=(unsigned int) MAXALPHABETCHARACTER; i++)
  {
    alpha1->symbolmap[i] = alpha2->symbolmap[i];
  }
  for (i=0; i<alpha1->mapsize; i++)
  {
    alpha1->characters[i] = alpha2->characters[i];
  }
  for (i=0; i<alpha1->domainsize; i++)
  {
    alpha1->mapdomain[i] = alpha2->mapdomain[i];
  }
  return alpha1;
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

static void assignDNAalphabet(SfxAlphabet *alpha)
{
  alpha->wildcardshow = (GtUchar) DNAWILDCARDS[0];
  alpha->mappedwildcards = (unsigned int) strlen(DNAWILDCARDS);
  alpha->domainsize = (unsigned int) strlen(DNAALPHABETDOMAIN);
  alpha->bitspersymbol = 3U; /* as we have to represent 4 + 2 characters */
  ALLOCASSIGNSPACE(alpha->mapdomain,NULL,GtUchar,alpha->domainsize);
  memcpy(alpha->mapdomain,(GtUchar *) DNAALPHABETDOMAIN,
         (size_t) alpha->domainsize);
  alpha->mapsize = MAPSIZEDNA;
  ALLOCASSIGNSPACE(alpha->characters,NULL,char,MAPSIZEDNA-1);
  memcpy(alpha->characters,"acgt",(size_t) (MAPSIZEDNA-1));
  assignDNAsymbolmap(alpha->symbolmap);
}

static void assignproteinsymbolmap(GtUchar *symbolmap)
{
  unsigned int cnum;

  for (cnum=0; cnum<=(unsigned int) MAXALPHABETCHARACTER; cnum++)
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

static void assignProteinalphabet(SfxAlphabet *alpha)
{
  alpha->wildcardshow = (GtUchar) PROTEINWILDCARDS[0];
  alpha->domainsize = (unsigned int) strlen(PROTEINALPHABETDOMAIN);
  alpha->mappedwildcards = (unsigned int) strlen(PROTEINWILDCARDS);
  alpha->bitspersymbol = 5U; /* as we have to represent 20 + 2 characters */
  ALLOCASSIGNSPACE(alpha->mapdomain,NULL,GtUchar,alpha->domainsize);
  memcpy(alpha->mapdomain,
         (GtUchar *) PROTEINALPHABETDOMAIN,(size_t) alpha->domainsize);
  alpha->mapsize = MAPSIZEPROTEIN;
  ALLOCASSIGNSPACE(alpha->characters,NULL,char,MAPSIZEPROTEIN-1);
  memcpy(alpha->characters,PROTEINUPPERAMINOACIDS,(size_t) (MAPSIZEPROTEIN-1));
  assignproteinsymbolmap(alpha->symbolmap);
}

static int assignProteinorDNAalphabet(SfxAlphabet *alpha,
                                      const GtStrArray *filenametab,
                                      GtError *err)
{
  int retval;

  gt_error_check(err);
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

void freeSfxAlphabet(SfxAlphabet **alpha)
{
  FREESPACE((*alpha)->mapdomain);
  FREESPACE((*alpha)->characters);
  FREESPACE(*alpha);
}

/*@null@*/ SfxAlphabet *assigninputalphabet(bool isdna,
                                         bool isprotein,
                                         const GtStr *smapfile,
                                         const GtStrArray *filenametab,
                                         GtError *err)
{
  SfxAlphabet *alpha;
  bool haserr = false;

  ALLOCASSIGNSPACE(alpha,NULL,SfxAlphabet,(size_t) 1);
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
        if (readsymbolmap(alpha,
                          transpath == NULL ? smapfile : transpath,
                          err) != 0)
        {
          haserr = true;
        }
        gt_str_delete(transpath);
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
      freeSfxAlphabet(&alpha);
    }
    return NULL;
  }
  return alpha;
}

const GtUchar *getsymbolmapAlphabet(const SfxAlphabet *alpha)
{
  return alpha->symbolmap;
}

unsigned int getnumofcharsAlphabet(const SfxAlphabet *alpha)
{
  return alpha->mapsize-1;
}

const GtUchar *getcharactersAlphabet(const SfxAlphabet *alpha)
{
  return alpha->characters;
}

GtUchar getwildcardshowAlphabet(const SfxAlphabet *alpha)
{
  return alpha->wildcardshow;
}

unsigned int getbitspersymbolAlphabet(const SfxAlphabet *alpha)
{
  return alpha->bitspersymbol;
}

void outputalphabet(FILE *fpout,const SfxAlphabet *alpha)
{
  GtUchar chartoshow, currentcc, previouscc = 0, firstinline = 0;
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

void fprintfsymbolstring(FILE *fpout,const SfxAlphabet *alpha,
                         const GtUchar *w,unsigned long wlen)
{
  unsigned long i;
  const GtUchar *characters;

  if (alpha == NULL)
  {
    characters = (const GtUchar *) "acgt";
  } else
  {
    characters = alpha->characters;
  }
  for (i = 0; i < wlen; i++)
  {
    (void) putc((int) characters[(int) w[i]],fpout);
  }
}

/*
  The following function is a special case of the previous
  function showing the output on stdout.
*/

void printfsymbolstring(const SfxAlphabet *alpha,
                        const GtUchar *w,unsigned long wlen)
{
  fprintfsymbolstring(stdout,alpha,w,wlen);
}

static char converttoprettysymbol(const SfxAlphabet *alpha, GtUchar currentchar)
{
  char ret = '\0';
  if (alpha == NULL)
  {
    ret = (char) currentchar;
  } else
  {

    if (currentchar == (GtUchar) WILDCARD)
    {
      ret = (char) alpha->wildcardshow;
    } else
    {
      if (currentchar != (GtUchar) SEPARATOR)
      {
        gt_assert((unsigned int) currentchar < alpha->mapsize-1);
        ret = (char) alpha->characters[(int) currentchar];
      }
    }
  }
  gt_assert(ret != '\0');
  return ret;
}

void sprintfsymbolstring(char *buffer,const SfxAlphabet *alpha,
                          const GtUchar *w,unsigned long wlen)
{
  unsigned long i;

  for (i = 0; i < wlen; i++)
  {
    buffer[i] = converttoprettysymbol(alpha, (GtUchar) w[i]);
  }
  buffer[wlen] = '\0';
}

void echoprettysymbol(FILE *fpout,const SfxAlphabet *alpha,GtUchar currentchar)
{
  (void) putc((int) converttoprettysymbol(alpha, currentchar), fpout);
}

GtUchar getprettysymbol(const SfxAlphabet *alpha,unsigned int currentchar)
{
  gt_assert(currentchar <= UCHAR_MAX);
  return (GtUchar) converttoprettysymbol(alpha, (GtUchar) currentchar);
}

static unsigned int removelowercaseproteinchars(GtUchar *domainbuf,
                                                const SfxAlphabet *alpha)
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

/*EE
  The following function checks if the given alphabet is the Protein
  alphabet with the aminoacids
  A, C, D, E, F, G, H, I, K, L, M, N, P, Q, R, S, T, V, W, Y written in
  lower or upper case.
*/

bool isproteinalphabet(const SfxAlphabet *alpha)
{
  SfxAlphabet proteinalphabet;
  unsigned int i, reduceddomainsize1, reduceddomainsize2;
  bool isprot = false;
  GtUchar domainbuf1[MAXALPHABETCHARACTER+1],
        domainbuf2[MAXALPHABETCHARACTER+1];

  reduceddomainsize1 = removelowercaseproteinchars(&domainbuf1[0],alpha);
  assignProteinalphabet(&proteinalphabet);
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
  FREESPACE(proteinalphabet.mapdomain);
  FREESPACE(proteinalphabet.characters);
  return isprot;
}

static bool checksymbolmap(const GtUchar *testsymbolmap,
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

/*
  The following function checks if the given alphabet is the DNA
  alphabet with the bases A, C, G, T written in lower or upper case.
*/

bool isdnaalphabet(const SfxAlphabet *alpha)
{
  if (isproteinalphabet(alpha))
  {
    return false;
  }
  if (alpha->mapsize == MAPSIZEDNA)
  {
    GtUchar dnasymbolmap[MAXALPHABETCHARACTER+1];

    assignDNAsymbolmap(&dnasymbolmap[0]);
    if (checksymbolmap(alpha->symbolmap,&dnasymbolmap[0],"acgt"))
    {
      return true;
    }
  }
  return false;
}
