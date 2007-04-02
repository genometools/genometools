/*
  Copyright (c) 2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#ifndef SAFEOP_H
#define SAFEOP_H

#define       safe_cast_to_long(ulong)\
              safe_cast_to_long_type(ulong, __FILE__, __LINE__)
long          safe_cast_to_long_type(unsigned long, const char*, unsigned int);

#define       safe_cast_to_ulong(slong)\
              safe_cast_to_ulong_type(slong, __FILE__, __LINE__)
unsigned long safe_cast_to_ulong_type(long, const char*, unsigned int);

#endif
