/*
  Copyright (c) 2007 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
  Copyright (c) 2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#include <stdio.h>
#include <ctype.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdbool.h>
#include <limits.h>
#include "libgtcore/env.h"
#include "libgtcore/str.h"
#include "libgtcore/strarray.h"
#include "qsorttype.h"
#include "symboldef.h"
#include "arraydef.h"
#include "chardef.h"
#include "alphadef.h"

#include "guessprot.pr"
#include "scanpaths.pr"
#include "readnextline.pr"

 struct _Alphabet
{
  Uchar characters[UCHAR_MAX+1],     /* array of characters to show */
        mapdomain[UCHAR_MAX+1],      /* list of characters mapped */
        symbolmap[UCHAR_MAX+1];      /* mapping of the symbols */
  uint32_t domainsize,               /* size of domain of symbolmap */
           mapsize,                  /* size of image of map, i.e. */
                                     /* mapping to [0..mapsize-1] */
           mappedwildcards;          /* number of mapped wildcards */
};

/*EE
  This file implements the datatype \texttt{alphabet}.
*/

/*
  Some constants for the standard alphabet used. The name says it all.
*/

#define DNABASES                     "aAcCgGtTuU"
#define DNAWILDCARDS                 "nsywrkvbdhmNSYWRKVBDHM"
#define MAPSIZEDNA                   ((uint32_t) 5)
#define DNAALPHABETDOMAIN            DNABASES DNAWILDCARDS
#define PROTEINUPPERAMINOACIDS       "LVIFKREDAGSTNQYWPHMC"
#define PROTEINLOWERAMINOACIDS       "lvifkredagstnqywphmc"
#define MAPSIZEPROTEIN               ((uint32_t) 21)
#define PROTEINWILDCARDS             "XUBZ*-"
#define PROTEINALPHABETDOMAIN        PROTEINUPPERAMINOACIDS PROTEINWILDCARDS

/*
  We use the following macro to access the \texttt{I}-th character of
  a line.
*/

#define LINE(I)          line.spaceUchar[I]

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

static int readsymbolmapviafp(Alphabet *alpha,
                              Uchar wildcard,
                              const Str *mapfile,
                              FILE *fpin,
                              Env *env)
{
  Uchar cc;
  unsigned cnum, linecount = 0;
  ArrayUchar line;
  unsigned long column;
  bool blankfound, ignore, preamble = true;

  env_error_check(env);
  alpha->domainsize = alpha->mapsize = alpha->mappedwildcards = 0;
  for (cnum=0; cnum<=UCHAR_MAX; cnum++)
  {
    alpha->symbolmap[cnum] = (Uchar) UNDEFCHAR;
  }
  INITARRAY(&line,Uchar);
  while (true)
  {
    line.nextfreeUchar = 0;
    if (readnextline(fpin,&line,env) == EOF)
    {
      break;
    }
    linecount++;
    ignore = false;
    if (line.nextfreeUchar > 0)
    {
      assert(line.spaceUchar != NULL);
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
        for (column=0; column<line.nextfreeUchar; column++) 
        { /* for all chars in line */
          cc = LINE(column);
          if (ispunct((int) cc) || isalnum((int) cc))
          {
            if (alpha->symbolmap[(uint32_t) cc] != (Uchar) UNDEFCHAR)
            {
              env_error_set(env,"cannot map symbol '%c' to %u: "
                            "it is already mapped to %u",
                             cc,
                             alpha->mapsize,
                             (uint32_t) alpha->symbolmap[(uint32_t) cc]);
              return -1;
            }
            /* get same value */
            alpha->symbolmap[(uint32_t) cc] = (Uchar) alpha->mapsize;
            alpha->mapdomain[alpha->domainsize++] = cc;
          } else
          {
            if (cc == (Uchar) ' ')    /* first blank in line found */
            {
              blankfound = true;
              /*@innerbreak@*/ break;
            }
            env_error_set(env,
                          "illegal character '%c' in line %u of mapfile %s",
                          cc,linecount,str_get(mapfile));
            return -2;
          }
        }
        if (blankfound)
        {
          if (isspace((int) LINE(column+1)))
          {
            env_error_set(env,
                          "illegal character '%c' at the end of "
                          "line %u in mapfile %s",
                          LINE(column+1),linecount,str_get(mapfile));
            return -3;
          }
          /* use next character to display character */
          alpha->characters[alpha->mapsize++] = LINE(column+1);
        } else
        {
          /* use first character of line to display character */
          alpha->characters[alpha->mapsize++] = LINE(0);
        }
      }
    }
  }
  for (cnum=0;cnum<=UCHAR_MAX; cnum++)
  {
    if (alpha->symbolmap[cnum] == (Uchar) (alpha->mapsize - 1))
    {
      if (wildcard > 0)
      {
        alpha->symbolmap[cnum] = wildcard; /* modify mapping for wildcard */
      }
      alpha->mappedwildcards++;
    }
  }
  if (wildcard > 0)
  {
    alpha->characters[(uint32_t) wildcard] 
      = alpha->characters[alpha->mapsize-1];
  }
  FREEARRAY(&line,Uchar);
  return 0;
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

static int readsymbolmap(Alphabet *alpha,Uchar wildcard,
                         const Str *mapfile,Env *env)
{
  FILE *fpin;

  env_error_check(env);
  fpin = scanpathsforfile("MKVTREESMAPDIR",mapfile,env);
  if (fpin == NULL)
  {
    return -1;
  }
  if (readsymbolmapviafp(alpha,wildcard,mapfile,fpin,env) != 0)
  {
    return -2;
  }
  env_fa_xfclose(fpin,env);
  return 0;
}

static void assignDNAsymbolmap(Uchar *symbolmap)
{
  uint32_t cnum;

  for (cnum=0; cnum<=(uint32_t) UCHAR_MAX; cnum++)
  {
    symbolmap[cnum] = (Uchar) UNDEFCHAR;
  }
  symbolmap[(uint32_t) 'a'] = (Uchar) 0;
  symbolmap[(uint32_t) 'A'] = (Uchar) 0;
  symbolmap[(uint32_t) 'c'] = (Uchar) 1;
  symbolmap[(uint32_t) 'C'] = (Uchar) 1;
  symbolmap[(uint32_t) 'g'] = (Uchar) 2;
  symbolmap[(uint32_t) 'G'] = (Uchar) 2;
  symbolmap[(uint32_t) 't'] = (Uchar) 3;
  symbolmap[(uint32_t) 'T'] = (Uchar) 3;
  symbolmap[(uint32_t) 'u'] = (Uchar) 3;
  symbolmap[(uint32_t) 'U'] = (Uchar) 3;
  for (cnum=0; DNAWILDCARDS[cnum] != '\0'; cnum++)
  {
    symbolmap[(uint32_t) DNAWILDCARDS[cnum]] = (Uchar) WILDCARD;
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
  If the argument \texttt{wildcard} is 0, then the wildcard characters
  in the last line are mapped to 4. Otherwise they are mapped to
  the character \texttt{WILDCARD}, as defined in \texttt{chardef.h}
*/

static void assignDNAalphabet(Alphabet *alpha)
{
  alpha->domainsize = (uint32_t) strlen(DNAALPHABETDOMAIN);
  alpha->mappedwildcards = (uint32_t) strlen(DNAWILDCARDS);
  memcpy(alpha->mapdomain,
         (Uchar *) DNAALPHABETDOMAIN,
         (size_t) alpha->domainsize);
  alpha->mapsize = MAPSIZEDNA;
  memcpy(alpha->characters,"acgt",(size_t) (MAPSIZEDNA-1));
  alpha->characters[WILDCARD] = (Uchar) DNAWILDCARDS[0];
  alpha->characters[MAPSIZEDNA-1] = (Uchar) DNAWILDCARDS[0];
  assignDNAsymbolmap(alpha->symbolmap);
}

static void assignproteinsymbolmap(Uchar *symbolmap)
{
  uint32_t cnum;

  for (cnum=0; cnum<=(uint32_t) UCHAR_MAX; cnum++)
  {
    symbolmap[cnum] = (Uchar) UNDEFCHAR;
  }
  for (cnum=0; PROTEINUPPERAMINOACIDS[cnum] != '\0'; cnum++)
  {
    symbolmap[(uint32_t) PROTEINUPPERAMINOACIDS[cnum]] = (Uchar) cnum;
  }
  for (cnum=0; PROTEINWILDCARDS[cnum] != '\0'; cnum++)
  {
    symbolmap[(uint32_t) PROTEINWILDCARDS[cnum]] = (Uchar) WILDCARD;
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
  alpha->domainsize = (uint32_t) strlen(PROTEINALPHABETDOMAIN);
  alpha->mappedwildcards = (uint32_t) strlen(PROTEINWILDCARDS);
  memcpy(alpha->mapdomain,
         (Uchar *) PROTEINALPHABETDOMAIN,(size_t) alpha->domainsize);
  alpha->mapsize = MAPSIZEPROTEIN;
  memcpy(alpha->characters,PROTEINUPPERAMINOACIDS,(size_t) MAPSIZEPROTEIN-1);
  alpha->characters[WILDCARD] = (Uchar) PROTEINWILDCARDS[0];
  alpha->characters[MAPSIZEPROTEIN-1] = (Uchar) PROTEINWILDCARDS[0];
  assignproteinsymbolmap(alpha->symbolmap);
}

static int assignProteinorDNAalphabet(Alphabet *alpha,
                                      const StrArray *filenametab,Env *env)
{
  int retval = guessifproteinsequencestream(filenametab,env);
  if(retval < 0)
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

/*@null@*/ Alphabet *assigninputalphabet(bool isdna,
                                         bool isprotein,
                                         const Str *smapfile,
                                         const StrArray *filenametab,
                                         Env *env)
{
  Alphabet *alpha;

  env_error_check(env);
  ALLOCASSIGNSPACE(alpha,NULL,Alphabet,(size_t) 1);
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
        if (readsymbolmap(alpha,
                         (Uchar) WILDCARD,
                         smapfile,
                         env) != 0)
        {
          return NULL;
        }
      } else
      {
        if(assignProteinorDNAalphabet(alpha,filenametab,env) != 0)
        {
          return NULL;
        }
      }
    }
  }
  return alpha;
}

const Uchar *getsymbolmapAlphabet(const Alphabet *alpha)
{
  return alpha->symbolmap;
}

uint32_t getnumofcharsAlphabet(const Alphabet *alpha)
{
  return alpha->mapsize-1;
}

uint32_t getmapsizeAlphabet(const Alphabet *alpha)
{
  return alpha->mapsize;
}

const Uchar *getcharactersAlphabet(const Alphabet *alpha)
{
  return alpha->characters;
}

Uchar *copycharactersAlphabet(const Alphabet *alpha,Env *env)
{
  Uchar *characters;

  ALLOCASSIGNSPACE(characters,NULL,Uchar,alpha->domainsize);
  (void) memcpy(characters,alpha->characters,(size_t) alpha->domainsize);
  return characters;
}

void freeAlphabet(Alphabet **alpha,Env *env)
{
  FREESPACE(*alpha);
}

void outputalphabet(FILE *fpout,const Alphabet *alpha)
{
  Uchar currentcc, previouscc = 0, firstinline = 0;
  uint32_t cnum, linenum = 0;
  bool afternewline = true;

  for (cnum=0; cnum < alpha->domainsize; cnum++)
  {
    currentcc = alpha->mapdomain[cnum];
    if (cnum > 0)
    {
      if (alpha->symbolmap[currentcc] != alpha->symbolmap[previouscc])
      {
        if (firstinline != alpha->characters[linenum])
        {
          fprintf(fpout," %c",(int) alpha->characters[linenum]);
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
  if (firstinline != alpha->characters[linenum])
  {
    fprintf(fpout," %c",(int) alpha->characters[linenum]);
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
 
  for(i = 0; i < wlen; i++)
  {
    (void) putc((int) alpha->characters[(int) w[i]],fpout);
  }
}

/*
  The following function is a special case of the previous
  function showing the output on stdout.
*/

void showsymbolstring(const Alphabet *alpha,const Uchar *w,unsigned long wlen)
{
  showsymbolstringgeneric(stdout,alpha,w,wlen);
}

static uint32_t removelowercaseproteinchars(Uchar *domainbuf,
                                            const Alphabet *alpha)
{
  uint32_t i, j = 0;

  for(i=0; i< alpha->domainsize - alpha->mappedwildcards; i++)
  {
    if(isalnum((int) alpha->mapdomain[i]) && 
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
  if(UNCAST(a) < UNCAST(b))
  {
    return -1;
  }
  if(UNCAST(a) > UNCAST(b))
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
  uint32_t i, reduceddomainsize1, reduceddomainsize2;
  Uchar domainbuf1[UCHAR_MAX+1], 
        domainbuf2[UCHAR_MAX+1];

  reduceddomainsize1 = removelowercaseproteinchars(&domainbuf1[0],alpha);
  assignProteinalphabet(&proteinalphabet);
  reduceddomainsize2 = removelowercaseproteinchars(&domainbuf2[0],
                                                   &proteinalphabet);
  if(reduceddomainsize1 == reduceddomainsize2)
  {
    qsort(&domainbuf1[0],(size_t) reduceddomainsize1,sizeof(char),
          (Qsortcomparefunction) comparechar);
    qsort(&domainbuf2[0],(size_t) reduceddomainsize2,sizeof(char),
          (Qsortcomparefunction) comparechar);
    for(i=0; i < reduceddomainsize2; i++)
    {
      if(domainbuf1[i] != domainbuf2[i])
      {
        return false;
      }
    }
    return true;
  }
  return false;
}

static bool checksymbolmap(const Uchar *testsymbolmap,
                           const Uchar *verifiedsymbolmap,
                           const char *testcharacters)
{
  unsigned int i;
  Uchar cc1, cc2 = 0;

  for(i=0; testcharacters[i] != '\0'; i++)
  {
    cc1 = testcharacters[i];
    if(isupper((int) cc1))
    {
      cc2 = tolower((int) cc1);
    } else
    {
      if(islower((int) cc1))
      {
        cc2 = toupper((int) cc2);
      } else
      {
        fprintf(stderr,"checksymbolmap used for non-alphabet character %c\n",
                cc1);
        exit(EXIT_FAILURE);
      }
    }
    if(testsymbolmap[cc1] != verifiedsymbolmap[cc1] &&
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
  if(isproteinalphabet(alpha))
  {
    return false;
  }
  if(alpha->mapsize == MAPSIZEDNA)
  {
    Uchar dnasymbolmap[UCHAR_MAX+1];

    assignDNAsymbolmap(&dnasymbolmap[0]);
    if(checksymbolmap(alpha->symbolmap,&dnasymbolmap[0],"acgt"))
    {
      return true;
    }
  }
  return false;
}
