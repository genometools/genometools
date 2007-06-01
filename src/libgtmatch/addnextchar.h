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
