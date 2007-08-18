/*
  Copyright (c) 2005-2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2005-2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#ifndef SFX_RI_DEF_H
#define SFX_RI_DEF_H

typedef struct Readintkeys Readintkeys;

#define SETREADINTKEYS(VALNAME,VAL,FORCEREAD)\
        setreadintkeys(riktab,VALNAME,VAL,sizeof (*(VAL)),FORCEREAD,env)

#endif
