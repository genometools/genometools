#include "libgtcore/minmax.h"
#include "symboldef.h"
#include "chardef.h"
#include "arraydef.h"

#define COMPARESYMBOLS(A,B)\
        if((A) == (Uchar) SEPARATOR)\
        {\
          gl->ubound = uptr;\
        }\
        if((B) == (Uchar) SEPARATOR)\
        {\
          gl->vbound = vptr;\
        }\
        if((A) != (B) || ISSPECIAL(A))\
        {\
          break;\
        }

#define STOREFRONT(GL,F,V)            F = V
#define ROWVALUE(FVAL)                *(FVAL)

#define MINUSINFINITYFRONT(GL)        ((GL)->integermin)

/*
  A structure to store the global values.
*/

typedef long Frontvalue;

typedef struct
{
  long offset,   // absolute base address of the front
       width,    // width of the front
       left;     // left boundary (negative),
                 // -left is relative address of entry
} Frontspec;     // \Typedef{Frontspec}

DECLAREARRAYSTRUCT(Frontspec);

typedef struct
{
  const Uchar *useq, 
              *vseq,
              *ubound, 
              *vbound;
  long ulen, 
       vlen, 
       integermin;
  Frontvalue *frontspace;
} FrontResource;

#ifdef DEBUG
static void showfront(FrontResource *gl,Frontspec *fspec,long r)
{
  Frontvalue *fval;
  long k;

  for(fval = gl->frontspace + fspec->offset, k=fspec->left; 
      k < fspec->left + fspec->width; k++, fval++)
  {
    if(r <= 0 || k <= -r || k >= r)
    {
      printf("(k=%ld)%ld ",k,ROWVALUE(fval));
    } else
    {
      printf("(k=%ld)undef ",k);
    }
  }
  (void) putchar('\n');
}
#endif

static void frontspecparms(FrontResource *gl,Frontspec *fspec,long p,long r)
{
  if(r <= 0)
  {
    fspec->left = -p;
    fspec->width = p + p + 1;
  } else
  {
    fspec->left = MAX(-gl->ulen,-p);
    fspec->width = MIN(gl->vlen,p) - fspec->left + 1;
  }
  printf("p=%ld,offset=%ld,left=%ld,width=%ld\n",p,
                                                 fspec->offset,
                                                 fspec->left,
                                                 fspec->width);
}

static long accessfront(FrontResource *gl,Frontvalue *fval,
                        Frontspec *fspec,long k)
{
  if(fspec->left <= k && k < fspec->left + fspec->width)
  {
    return ROWVALUE(fval+k);
  }
  return MINUSINFINITYFRONT(gl);
}

/*
  The following function evaluates an entry \(front(k,p)\) in a 
  forward direction.
*/

static void evalentryforward(FrontResource *gl,
                             Frontvalue *fval,
                             Frontspec *fspec,
                             long k)
{
  long value, t;
  const Uchar *uptr, *vptr;
  Uchar a, b;
  Frontvalue *fptr;

  printf("evalentryforward(k=%ld)\n",k);
  fptr = gl->frontspace + fspec->offset - fspec->left;
  t = accessfront(gl,fptr,fspec,k) + 1;         // same diagonal
  printf("same: access(k=%ld)=%ld\n",k,t-1);

  value = accessfront(gl,fptr,fspec,k-1);       // diagonal below
  printf("below: access(k=%ld)=%ld\n",k-1,value);
  if(t < value)
  {
    t = value;
  }
  value = accessfront(gl,fptr,fspec,k+1) + 1;     // diagonal above
  printf("above: access(k=%ld)=%ld\n",k+1,value-1);
  if(t < value)
  {
    t = value;
  }
  printf("maximum: t=%ld\n",t);   // the maximum over three values
  if(t < 0 || t+k < 0)             // no negative value
  {
    STOREFRONT(gl,ROWVALUE(fval),MINUSINFINITYFRONT(gl));
  } else
  {
    uptr = gl->useq + t;
    vptr = gl->vseq + t + k;
    if(gl->ulen != 0 && gl->vlen != 0)  // only for nonempty strings
    {
      if(uptr == vptr)    // strings are equal
      {
        t = gl->ulen-1;
      } else
      {
        for(/* Nothing */; uptr < gl->ubound && vptr < gl->vbound;
            uptr++, vptr++)
        {
          a = *uptr;
          b = *vptr;
          COMPARESYMBOLS(a,b);
        }
        t = (long) (uptr - gl->useq);
      }
    }
    if(gl->useq + t > gl->ubound || gl->vseq + t + k > gl->vbound)
    {
      STOREFRONT(gl,ROWVALUE(fval),MINUSINFINITYFRONT(gl));
    } else
    {
      STOREFRONT(gl,ROWVALUE(fval),t);
    }
  }
}

/* 
  The following function evaluates a front in forward direction.
  It returns true if any of the returned values is at least 0.
*/

static bool evalfrontforward(FrontResource *gl,
                             Frontspec *prevfspec,
                             Frontspec *fspec,
                             long r)
{
  long k;
  bool defined = false; 
  Frontvalue *fval;

  for(fval = gl->frontspace + fspec->offset, k = fspec->left; 
      k < fspec->left + fspec->width; k++, fval++)
  {
    if(r <= 0 || k <= -r || k >= r)
    {
      evalentryforward(gl,fval,prevfspec,k);
      if(ROWVALUE(fval) >= 0)
      {
        defined = true;
      }
      printf("store front[k=%ld]=%ld ",k,ROWVALUE(fval));
      printf("at index %ld\n",(long) (fval-gl->frontspace));
    } else
    {
      printf("store front[k=%ld]=MINUSINFINITYFRONT ",k);
      printf("at index %ld\n",(long) (fval-gl->frontspace));
      STOREFRONT(gl,ROWVALUE(fval),MINUSINFINITYFRONT(gl));
    }
  }
  printf("frontvalues[r=%ld]=",r);
#ifdef DEBUG
  showfront(gl,fspec,r);
#endif
  return defined;
}

/*
  The following function evaluates the entry \(front(0,0)\) in a 
  forward direction.
*/

static void firstfrontforward(FrontResource *gl,Frontspec *fspec)
{
  Uchar a, b;
  const Uchar *uptr, *vptr;

  fspec->left = fspec->offset = 0;
  fspec->width = (long) 1;
  if(gl->ulen == 0 || gl->vlen == 0)
  {
    STOREFRONT(gl,ROWVALUE(&gl->frontspace[0]),0);
  } else
  {
    for(uptr = gl->useq, vptr = gl->vseq; 
        uptr < gl->ubound && 
        vptr < gl->vbound;
        uptr++, vptr++)
    {
      a = *uptr;
      b = *vptr;
      COMPARESYMBOLS(a,b);
    }
    STOREFRONT(gl,ROWVALUE(&gl->frontspace[0]),(long) (uptr - gl->useq));
  } 
  printf("forward front[0]=%ld\n",ROWVALUE(&gl->frontspace[0]));
}

long unitedistfrontSEPgeneric(bool withmaxdist,
                              unsigned long maxdist,
                              const Uchar *useq,
                              long ulen,
                              const Uchar *vseq,
                              long vlen,
                              Env *env)
{
  unsigned long currentallocated;
  FrontResource gl;
  Frontspec frontspecspace[2],
            *fspec, 
            *prevfspec;
  Frontvalue *fptr;
  long k, r, realdistance;

  printf("unitedistcheckSEPgeneric(ulen=%ld,vlen=%ld)\n",ulen,vlen);

  if(withmaxdist)
  {
    currentallocated = (maxdist+1) * (maxdist+1) + 1;
  } else
  {
    currentallocated = (unsigned long) 1;
  }
  ALLOCASSIGNSPACE(gl.frontspace,NULL,Frontvalue,currentallocated);
  gl.useq = useq;
  gl.vseq = vseq;
  gl.ubound = useq + ulen;
  gl.vbound = vseq + vlen;
  gl.ulen = ulen;
  gl.vlen = vlen;
  gl.integermin = -MAX(ulen,vlen);
  prevfspec = &frontspecspace[0];
  firstfrontforward(&gl,prevfspec);
  if(gl.ulen == gl.vlen && ROWVALUE(&gl.frontspace[0]) == gl.vlen)
  {
    realdistance = 0;
  } else
  {
    for(k=(long) 1, r=1-MIN(gl.ulen,gl.vlen); /* Nothing */ ; k++, r++)
    {
      if(withmaxdist && k > (long) maxdist)
      {
        realdistance = (long) (maxdist + 1);
        break;
      }
      if(prevfspec == &frontspecspace[0])
      {
        fspec = &frontspecspace[1];
      } else
      {
        fspec = &frontspecspace[0];
      }
      fspec->offset = prevfspec->offset + prevfspec->width;
      frontspecparms(&gl,fspec,k,r);
      while((unsigned long) (fspec->offset + fspec->width) >= currentallocated)
      {
        if(withmaxdist)
        {
          fprintf(stderr,"Not enough frontspace: "
                         "maxdist=%lu,p=%ld,offset=%ld,width=%ld\n",
                          maxdist,
                          k,
                          fspec->offset,
                          fspec->width);
          exit(EXIT_FAILURE);
        }
        currentallocated += (k+1);
        ALLOCASSIGNSPACE(gl.frontspace,gl.frontspace,
                         Frontvalue,currentallocated);
      }
      (void) evalfrontforward(&gl,prevfspec,fspec,r);
      fptr = gl.frontspace + fspec->offset - fspec->left;
      if(accessfront(&gl,fptr,fspec,(long) (vlen - ulen)) == ulen)
      {
        realdistance = k;
        break;
      }
      if(prevfspec == &frontspecspace[0])
      {
        prevfspec = &frontspecspace[1];
      } else
      {
        prevfspec = &frontspecspace[0];
      }
    }
    if(withmaxdist)
    {
      if(realdistance <= (long) maxdist)
      {
        printf("unitedistfrontSEP returns %ld\n",realdistance);
        FREESPACE(gl.frontspace);
        return realdistance;
      }
      printf("unitedistfrontSEP returns -1\n");
      FREESPACE(gl.frontspace);
      return (long) -1;
    }
  }
  printf("unitedistfrontSEP returns %ld\n",realdistance);
  FREESPACE(gl.frontspace);
  return realdistance;
}
