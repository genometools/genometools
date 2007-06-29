/*
  Copyright (c) 2007 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
  Copyright (c) 2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#ifndef TYPES_H
#define TYPES_H
#include <sys/types.h>
#include <limits.h>
#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <inttypes.h>

/*
  Show a boolean value as a string or as a character 0 or 1.
*/

#define SHOWBOOL(B) ((B) ? "true" : "false")
#define SHOWBIT(B)  ((B) ? '1' : '0')

/*
  This file contains some basic type definition.
*/

typedef uint8_t Uchar;         /* \Typedef{Uchar} */
typedef uint16_t Ushort;        /* \Typedef{Ushort} */

#define Uint32Const(N)   (N##U)  /* uint32_t constant */
#define Uint64Const(N)   (N##UL) /* uint64_t constant */

#ifdef S_SPLINT_S
#define Formatuint64_t "%lu"
#define Scanuint64_tcast(X) ((unsigned long *) (X))
#define PRINTSeqposcast(X)  ((unsigned long) (X))
#define PRINTuint64_tcast(X) ((unsigned long) (X))
#else
#define Formatuint64_t "%" PRIu64
#define Scanuint64_tcast(X) (X)
#define PRINTSeqposcast(X)  (X)
#define PRINTuint64_tcast(X) (X)
#endif

/*
  The following is the central case distinction to accomodate
  code for 32 bit integers and 64 bit integers.
*/

#ifdef _LP64

typedef uint64_t Seqpos;         /* \Typedef{Seqpos} */
#ifdef S_SPLINT_S
#define FormatSeqpos "%lu"
#else
#define FormatSeqpos   "%" PRIu64
#endif

#else

#ifdef S_SPLINT_S
#define FormatSeqpos "%lu"
#else
#define Formatuint64_t "%" PRIu64
#endif

#ifdef BIGNUM32
typedef uint64_t Seqpos;         /* \Typedef{Seqpos} */
#ifdef S_SPLINT_S
#define FormatSeqpos "%lu"
#else
#define FormatSeqpos "%" PRIu64
#endif
#else
typedef uint32_t Seqpos;         /* \Typedef{Seqpos} */
#ifdef S_SPLINT_S
#define FormatSeqpos "%lu"
#else
#define FormatSeqpos "%" PRIu32
#endif
#define Seqposequalsunsignedint
#endif /* BIGNUM32 */

#endif /* _LP64 */


/*
  Type of unsigned integer in \texttt{printf}.
*/

typedef unsigned long Showuint;     /* \Typedef{Showuint} */

/*
  Argument of a function from \texttt{ctype.h}.
*/

typedef int Ctypeargumenttype;      /* \Typedef{Ctypeargumenttype} */

/*
  type of second argument of fgets
*/

typedef int Fgetssizetype;      /* \Typedef{Fgetssizetype} */

/*
  Return type of \texttt{fgetc} and \texttt{getc}.
*/

typedef int Fgetcreturntype;        /* \Typedef{Fgetcreturntype} */

/*
  Type of first argument of \texttt{putc}.
*/

typedef int Fputcfirstargtype;      /* \Typedef{Fputcfirstargtype} */

/*
  Returntype of \texttt{putc}.
*/

typedef int Fputcreturntype;      /* \Typedef{Fputcreturntype} */

/*
  Return type of \texttt{strcmp}.
*/

typedef int Strcmpreturntype;       /* \Typedef{Strcmpreturntype} */

/*
  Type of a file descriptor.
*/

typedef int Filedesctype;           /* \Typedef{Filedesctype} */

/*
  Return type of \texttt{qsort} function.
*/

typedef int Qsortcomparereturntype; /* \Typedef{Qsortcomparefunction} */

/*
  Return type of \texttt{sprintf} function.
*/

typedef int Sprintfreturntype;     /* \Typedef{Sprintfreturntype}  */

/*
  Type of fieldwidth in \texttt{printf} format string.
*/

typedef int Fieldwidthtype;         /* \Typedef{Fieldwidthtype} */

/*
  Type of \texttt{argc}-parameter in main.
*/

typedef int Argctype;               /* \Typedef{Argctype} */

/*
  Return type of \texttt{getrlimit}
*/

typedef int Getrlimitreturntype;    /* \Typedef{Getrlimitreturntype} */

/*
  type of error flag for gzerror
*/

typedef int Gzerrorflagtype;   /* \Typedef{Gzerrorflagtype} */

/*
  This is the type of the third argument of the function gzread
*/

typedef unsigned int Gzreadthirdarg;  /* \Typedef{Gzreadthirdarg} */

/*
  This is the return type of fclose like functions
*/

typedef int Fclosereturntype;  /* \Typedef{Fclosereturntype} */

/*
  This is the return type of fprintf like functions
*/

typedef int Fprintfreturntype;   /* \Typedef{Fprintfreturntype} */

/*
  This is the return value of the function regcomp and regexec
*/

typedef int Regexreturntype;    /* \Typedef{Regexreturntype} */

/*
  This is the type of arguments for function sysconf.
*/

typedef int Sysconfargtype;         /* \Typedef{Sysconfargtype} */

/*
  This is the type for integer codes for strings of some fixed length
*/

/*
  The following structure stores information about special characters.
*/

typedef struct
{
  Seqpos specialcharacters,      /* total number of special syms */
         specialranges,          /* number of ranges with special syms */
         lengthofspecialprefix,  /* number of specials at start of sequence */
         lengthofspecialsuffix;  /* number of specials at end of sequence */
} Specialcharinfo;

/*
  The following is the type of the comparison function
  to be provided to the function \texttt{qsort}.
*/

typedef int (*Qsortcomparefunction)(const void *,const void *);

/*
  The following is the type of the function showing information in
  verbose mode.
*/

typedef void (*Showverbose)(char *);

typedef struct
{
  Seqpos uint0,
         uint1;
} PairSeqpos;                /* \Typedef{PairSeqpos} */

/*
  The following type stores an unsigned integer only of defined is True
*/

#endif
