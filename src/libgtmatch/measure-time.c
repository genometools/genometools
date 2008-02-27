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

#include <stdio.h>
#include <time.h>
#include "libgtcore/ma.h"
#include "measure-time-if.h"

 struct Measuretime
{
  clock_t startclock, overalltime;
  const char *eventdescription;
};

Measuretime *inittheclock(const char *event)
{
  Measuretime *mtime;

  mtime = ma_malloc(sizeof (Measuretime));
  mtime->startclock = clock();
  mtime->overalltime = 0;
  mtime->eventdescription = event;
  return mtime;
}

void deliverthetime(FILE *fp,Measuretime *mtime,const char *newevent)
{
  clock_t stopclock;

  stopclock = clock();
  fprintf(fp,"# TIME %s %.2f\n",mtime->eventdescription,
             (double) (stopclock-mtime->startclock)/(double) CLOCKS_PER_SEC);
  (void) fflush(fp);
  mtime->overalltime += (stopclock - mtime->startclock);
  if (newevent == NULL)
  {
    fprintf(fp,"# TIME overall %.2f\n",
                (double) mtime->overalltime/(double) CLOCKS_PER_SEC);
    (void) fflush(fp);
    ma_free(mtime);
  } else
  {
    mtime->startclock = stopclock;
    mtime->eventdescription = newevent;
  }
}
