/*
  Copyright (c) 2007 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
  Copyright (c) 2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#include "libgtcore/env.h"
#include "types.h"
#include "mapspec-def.h"

#define ASSIGNPTR2STARTPTR(TYPE)\
        if (mapspec->numofunits == 0)\
        {\
          *((TYPE **) mapspec->startptr) = NULL;\
        } else\
        {\
          voidptr = (((Uchar *) ptr) + byteoffset);\
          *((TYPE **) mapspec->startptr) = voidptr;\
        }

#define WRITEACTIONWITHTYPE(TYPE)\
        if (fwrite(*((TYPE **) mapspecptr->startptr),\
                   mapspecptr->sizeofunit,\
                   (size_t) mapspecptr->numofunits, fp) !=\
                    (size_t) mapspecptr->numofunits)\
        {\
          haserr = true;\
        }

#define STAMP\
        printf("STAMP(%lu,%s)\n",(Showuint) __LINE__,__FILE__);\
        (void) fflush(stdout)

static Uint64 expectedindexsize(const ArrayMapspecification *mapspectable)
{
  Uint64 sumup = 0;
  Mapspecification *mapspecptr;

  for (mapspecptr = mapspectable->spaceMapspecification;
       mapspecptr < mapspectable->spaceMapspecification +
                    mapspectable->nextfreeMapspecification; mapspecptr++)
  {
    sumup += (Uint64) mapspecptr->sizeofunit * (Uint64) mapspecptr->numofunits;
  }
  CHECKIFFITS32BITS(sumup);
  return sumup;
}

static void showmapspec(const Mapspecification *mapspec)
{
  printf("(%s,size=%lu,elems=%lu)",
           mapspec->name,
           (Showuint) mapspec->sizeofunit,
           (Showuint) mapspec->numofunits);
}

static int assigncorrecttype(Mapspecification *mapspec,
                             void *ptr,
                             Uint byteoffset,
                             Env *env)
{
  void *voidptr;
  bool haserr = false;

  env_error_check(env);
  switch (mapspec->typespec)
  {
    case UcharType:
      ASSIGNPTR2STARTPTR(Uchar);
      break;
    case UshortType:
      ASSIGNPTR2STARTPTR(Ushort);
      break;
    case UintType:
      ASSIGNPTR2STARTPTR(Uint);
      break;
    case PairUintType:
      ASSIGNPTR2STARTPTR(PairUint);
      break;
    case Uint64Type:
      ASSIGNPTR2STARTPTR(Uint64);
      break;
    case PairUint64Type:
      ASSIGNPTR2STARTPTR(PairUint64);
      break;
    default:
      env_error_set(env,"no assignment specification for size %lu",
               (Showuint) mapspec->sizeofunit);
      haserr = true;
  }
  if (!haserr)
  {
    printf("# assign pointer");
    showmapspec(mapspec);
    printf(" at byteoffset %lu\n",(Showuint) byteoffset);
  }
  return haserr ? -1 : 0;
}

int fillmapspecstartptr(Assignmapspec assignmapspec,
                        void **mappeduserptr,
                        void *assignmapinfo,
                        const char *tmpfilename,
                        Uint expectedsize,
                        Env *env)
{
  void *mapptr;
  Uint64 expected;
  Uint byteoffset;
  size_t numofbytes;
  ArrayMapspecification mapspectable;
  Mapspecification *mapspecptr;
  bool haserr = false;

  env_error_check(env);
  INITARRAY(&mapspectable,Mapspecification);
  assignmapspec(&mapspectable,assignmapinfo,env);
  mapptr = env_fa_mmap_read(env,tmpfilename, &numofbytes);
  if (mapptr == NULL)
  {
    env_error_set(env,"could not map datafile %s",tmpfilename);
    haserr = true;
  }
  *mappeduserptr = mapptr;
  if (!haserr)
  {
    if (assigncorrecttype(mapspectable.spaceMapspecification,
                          mapptr,0,env) != 0)
    {
      haserr = true;
    }
  }
  if (!haserr)
  {
    expected = expectedindexsize(&mapspectable);
    if ((Uint64) numofbytes != expected)
    {
      env_error_set(env,"%lu bytes read from %s, but %lu expected",
                         (Showuint) numofbytes,tmpfilename,(Showuint) expected);
      haserr = true;
    }
  }
  if (!haserr)
  {
    mapspecptr = mapspectable.spaceMapspecification;
    assert(mapspecptr != NULL);
    byteoffset = (Uint) mapspecptr->sizeofunit * mapspecptr->numofunits;
    for (mapspecptr++;
        mapspecptr < mapspectable.spaceMapspecification +
                     mapspectable.nextfreeMapspecification; mapspecptr++)
    {
      if (assigncorrecttype(mapspecptr,mapptr,byteoffset,env) != 0)
      {
        haserr = true;
        break;
      }
      byteoffset += mapspecptr->sizeofunit * mapspecptr->numofunits;
    }
  }
  if (!haserr)
  {
    if (expectedsize != byteoffset)
    {
      env_error_set(env,"mapping: expected size of index is %lu bytes, "
                        "but index has %lu bytes",
                        (Showuint) expectedsize,(Showuint) byteoffset);
      haserr = true;
    }
  }
  FREEARRAY(&mapspectable,Mapspecification);
  return haserr ? -1 : 0;
}

int flushtheindex2file(FILE *fp,
                       Assignmapspec assignmapspec,
                       void *assignmapinfo,
                       Uint expectedsize,
                       Env *env)
{
  ArrayMapspecification mapspectable;
  Mapspecification *mapspecptr;
  Uint byteoffset = 0;
  bool haserr = false;

  env_error_check(env);
  INITARRAY(&mapspectable,Mapspecification);
  assignmapspec(&mapspectable,assignmapinfo,env);
  assert(mapspectable.spaceMapspecification != NULL);
  for (mapspecptr = mapspectable.spaceMapspecification;
       mapspecptr < mapspectable.spaceMapspecification +
                    mapspectable.nextfreeMapspecification; 
       mapspecptr++)
  {
    printf("# flushtheindex2file");
    showmapspec(mapspecptr);
    printf(" at byteoffset %lu\n",(Showuint) byteoffset);
    if (mapspecptr->numofunits > 0)
    {
      switch (mapspecptr->typespec)
      {
        case UcharType:
          WRITEACTIONWITHTYPE(Uchar);
          break;
        case UshortType:
          WRITEACTIONWITHTYPE(Ushort);
          break;
        case UintType:
          WRITEACTIONWITHTYPE(Uint);
          break;
        case Uint64Type:
          WRITEACTIONWITHTYPE(Uint64);
          break;
        case PairUintType:
          WRITEACTIONWITHTYPE(PairUint);
          break;
        case PairUint64Type:
          WRITEACTIONWITHTYPE(PairUint64);
          break;
        default:
           env_error_set(env,"no map specification for size %lu",
                         (Showuint) mapspecptr->sizeofunit);
           haserr = true;
      }
    }
    if (haserr)
    {
      break;
    }
    byteoffset += mapspecptr->sizeofunit *
                  mapspecptr->numofunits;
  }
  if (!haserr)
  {
    if (expectedsize != byteoffset)
    {
      env_error_set(env,"flushindex: expected size of index is %lu bytes, "
                        "but index has %lu bytes",
                        (Showuint) expectedsize,(Showuint) byteoffset);
      haserr = true;
    }
  }
  FREEARRAY(&mapspectable,Mapspecification);
  return haserr ? -1 : 0;
}
