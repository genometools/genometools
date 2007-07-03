/*
  Copyright (c) 2007 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
  Copyright (c) 2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#include <ctype.h>
#include <stdbool.h>
#include "types.h"
#include "inputsymbol.h"
#include "genstream.h"

#include "genericstream.pr"

bool guessifproteinsequence(const Uchar *input,Seqpos inputlen)
{
  Uchar current;
  const Uchar *inputptr;
  Seqpos countnonbases = 0,
         countcharacters = 0,
         readnumoffirstcharacters,
         leastnumofnonbases;
  bool indesc = false;

  if (inputlen < (Seqpos) 1000)
  {
    readnumoffirstcharacters = inputlen;
  } else
  {
    readnumoffirstcharacters = (Seqpos) 1000;
  }
  leastnumofnonbases = readnumoffirstcharacters/10;
  for (inputptr = input; countnonbases < leastnumofnonbases &&
                        countcharacters < readnumoffirstcharacters &&
                        inputptr < input + inputlen; inputptr++)
  {
    current = *inputptr;
    if (indesc)
    {
      if (current == NEWLINESYMBOL)
      {
        indesc = false;
      }
    } else
    {
      if (current == FASTASEPARATOR)
      {
        indesc = true;
      } else
      {
        if (!isspace((Ctypeargumenttype) current))
        {
          countcharacters++;
          switch (current)
          {
            case 'L':
            case 'I':
            case 'F':
            case 'E':
            case 'Q':
            case 'P':
            case 'X':
            case 'Z':
              countnonbases++;
              break;
            default:
              break;
          }
        }
      }
    }
  }
  if (countnonbases >= leastnumofnonbases)
  {
    return true;
  }
  return false;
}

bool guessifproteinsequencestream(const char *inputfile)
{
  Fgetcreturntype currentchar;
  Seqpos countnonbases = 0,
         countcharacters = 0;
  bool indesc = false;
  Genericstream inputstream;

  opengenericstream(&inputstream,inputfile);
  for (;;)
  {
    if (inputstream.isgzippedstream)
    {
      currentchar = gzgetc(inputstream.stream.gzippedstream);
    } else
    {
      currentchar = fgetc(inputstream.stream.fopenstream);
    }
    if (indesc)
    {
      if (currentchar == NEWLINESYMBOL)
      {
        indesc = false;
      }
    } else
    {
      if (currentchar == FASTASEPARATOR)
      {
        indesc = true;
      } else
      {
        if (!isspace((Ctypeargumenttype) currentchar))
        {
          countcharacters++;
          switch (currentchar)
          {
            case 'L':
            case 'I':
            case 'F':
            case 'E':
            case 'Q':
            case 'P':
            case 'X':
            case 'Z':
              countnonbases++;
              break;
            default:
              break;
          }
        }
      }
    }
    if (countcharacters >= (Seqpos) 1000)
    {
      break;
    }
  }
  closegenericstream(&inputstream,inputfile);
  if (countnonbases >= countcharacters/10)
  {
    return true;
  }
  return false;
}
