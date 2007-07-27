/*
  Copyright (c) 2007 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
  Copyright (c) 2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#ifndef CODESPEC_H
#define CODESPEC_H

#define FROMCODE2SPECIALCODE(CODE,NUMOFCHARS)\
                            (((CODE) - ((NUMOFCHARS)-1)) / (NUMOFCHARS))

#endif
