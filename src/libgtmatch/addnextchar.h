/*
  Copyright (c) 2007 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
  Copyright (c) 2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#ifndef ADDNEXTCHAR_H
#define ADDNEXTCHAR_H
#include "divmodmul.h"

#define ADDNEXTCHAR(CODE,CC,NUMOFCHARS)\
        if ((NUMOFCHARS) == DNAALPHASIZE)\
        {\
          CODE = MULT4(CODE) | ((Uint) (CC));\
        } else\
        {\
          CODE = ((CODE) * (NUMOFCHARS)) + (Uint) (CC);\
        }

#endif
