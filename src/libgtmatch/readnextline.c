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

#include "arraydef.h"

/*EE
  The following functions read the next line from the given generic file pointer
  and stores it in the array \texttt{line}.
*/

int readnextline(FILE *fpin,ArrayUchar *line,Env *env)
{
  int cc;

  env_error_check(env);
  while (true)
  {
    cc = fgetc(fpin);
    if (cc == EOF)
    {
      return EOF;
    }
    if (cc == (int) '\n')
    {
      STOREINARRAY(line,Uchar,512,(Uchar) '\0');
      line->nextfreeUchar--;
      return 0;
    }
    STOREINARRAY(line,Uchar,512,(Uchar) cc);
  }
}
