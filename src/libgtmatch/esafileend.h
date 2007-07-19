/*
  Copyright (c) 2007 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
  Copyright (c) 2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#ifndef ESAFILEEND_H
#define ESAFILEEND_H

#define ENCODEDSEQINDEX     ".esq"

/*
  The following defines the suffix of an indexname to be used for an
  reverse complemented index. The length is also defined
*/

#define RCMSUFFIX         ".rcm"

/*
  The following defines the suffix of an indexname to be used as
  a six frame index. The length is also defined
*/

#define SIXFRAMESUFFIX    ".6fr"

/*
  The following defines the suffix of a file to store the 
  table suftab.
*/

#define SUFTABSUFFIX  ".suf"

/*
  The following defines the suffix of a file to store the 
  table lcptab.
*/

#define LCPTABSUFFIX  ".lcp"

/*
  The following defines the suffix of a file to store the 
  large lcp table values.
*/

#define LARGELCPTABSUFFIX  ".llv"

/*
  The following defines the suffix of a file to store Burrows and
  Wheeler transform.
*/

#define BWTTABSUFFIX  ".bwt"

/*
  The following defines the suffix of a file to store bucket table.
*/

#define BCKTABSUFFIX  ".bck"

/*
  The following defines the suffix of a file to store alphabet.
*/

#define ALPHATABSUFFIX ".al1"

#endif
