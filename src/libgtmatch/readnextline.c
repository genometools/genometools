/*
  Copyright (c) 2007 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
  Copyright (c) 2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
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
