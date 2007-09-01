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

#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include "libgtcore/env.h"

/*
  This file implements the function \texttt{getbasename} according to the
  specification in
  \begin{center}
  \texttt{http://www.unix-systems.org/onlinepubs/7908799/xsh/basename.html}
  \end{center}
  \texttt{getbasename} is equivalent to the function \texttt{basename} which
  is available on most unix systems, but in different libraries and
  slightly different functionality.

  \texttt{getbasename} takes the pathname pointed to by
  \texttt{path} and returns a pointer to the final
  component of the pathname, deleting any trailing
  \texttt{\CharQuote{/}} characters.

  If \texttt{path} consists entirely of the \texttt{\CharQuote{/}}
  character,  then \texttt{getbasename} returns a pointer to the string
  \StringQuote{/}.

  If \texttt{path} is a null pointer or points to an empty string,
  \texttt{getbasename} returns a pointer to the string
  \StringQuote{.}.

  \texttt{getbasename} function returns a pointer to static storage that
  may then be overwritten by a subsequent call to \texttt{getbasename}.

  \begin{center}
  \begin{tabular}{|l|l|}\hline
  \multicolumn{2}{|c|}{Example:}\\
  Input String&Output String\\\hline
  \StringQuote{/usr/lib}  &\StringQuote{lib}\\
  \StringQuote{/usr/}     &\StringQuote{usr}\\
  \StringQuote{/}         &\StringQuote{/}\\
  \StringQuote{}          &\StringQuote{.}\\
  \hline
  \end{tabular}
  \end{center}
*/

/*@null@*/ char *getbasename(const char *path,Env *env)
{
  char *sbuf, *c;
  unsigned char foundother = 0;
  size_t pathlen;

  env_error_check(env);
  pathlen = strlen(path);
  sbuf = env_ma_malloc(env,sizeof (char) * (pathlen+2));
  if (path == NULL || *path == '\0')
  {
    strcpy(sbuf,".");
    return sbuf;
  }
  strcpy(sbuf,path);
  for (c = sbuf + pathlen - 1; c >= sbuf; c--)
  {
    if (*c == '/')
    {
      if (foundother)
      {
        size_t i;
        c++;
        for (i=0; c[i] != '\0'; i++)
        {
          sbuf[i] = c[i];
        }
        sbuf[i] = '\0';
        break;
      }
      if (c > sbuf)
      {
        *c = '\0';
      }
    } else
    {
      foundother = 1U;
    }
  }
  return sbuf;
}
