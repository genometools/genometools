#include "types.h"
#include "arraydef.h"

/*EE
  The following functions read the next line from the given generic file pointer
  and stores it in the array \texttt{line}.
*/

int readnextline(FILE *fpin,ArrayUchar *line,Env *env)
{
  Fgetcreturntype cc;

  while (true)
  {
    cc = fgetc(fpin);
    if (cc == EOF)
    {
      return EOF;
    }
    if (cc == (Fgetcreturntype) '\n')
    {
      STOREINARRAY(line,Uchar,512,(Uchar) '\0');
      line->nextfreeUchar--;
      return 0;
    }
    STOREINARRAY(line,Uchar,512,(Uchar) cc);
  }
}
