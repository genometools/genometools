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

#ifndef ESAFILEEND_H
#define ESAFILEEND_H

/*
  The following defines the suffix of an indexname to be used for an
  reverse complemented index. The length is also defined.
*/

#define RCMSUFFIX         ".rcm"

/*
  The following defines the suffix of an indexname to be used as
  a six frame index. The length is also defined.
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

#define ALPHABETFILESUFFIX ".al1"

/*
  The following defines the suffix of a file to store sequence descriptions.
*/

#define DESTABSUFFIX ".des"

/*
  The following defines the suffix of the project file.
*/

#define PROJECTFILESUFFIX ".prj"

#endif
