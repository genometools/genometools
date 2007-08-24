#include <limits.h>
#include "sarr-def.h"
#include "seqpos-def.h"
#include "symboldef.h"
#include "spacedef.h"

#define ABOVETOP  stackspace[nextfreeItvinfo]
#define TOP       stackspace[nextfreeItvinfo-1]
#define BELOWTOP  stackspace[nextfreeItvinfo-2]

#define INCSTACKSIZE  4

#define PUSHDFS(D,B)\
        stackspace[nextfreeItvinfo].depth = D;\
        stackspace[nextfreeItvinfo].lastisleafedge = B;\
        nextfreeItvinfo++

 DECLAREREADFUNCTION(Seqpos);

 DECLAREREADFUNCTION(Uchar);

 DECLAREREADFUNCTION(Largelcpvalue);

typedef struct Dfsinfo Dfsinfo;

typedef struct
{
  bool lastisleafedge;
  Seqpos depth;
  Dfsinfo *dfsinfo;
} Itvinfo;

static Itvinfo *allocItvinfo(Itvinfo *ptr,
                             unsigned long currentallocated,
                             unsigned long allocated,
                             Dfsinfo *(*allocateDfsinfo)(void *,Env *),
                             void *info,
                             Env *env)
{
  unsigned long i;

  ALLOCASSIGNSPACE(ptr,ptr,Itvinfo,allocated);
  if(allocateDfsinfo != NULL)
  {
    for(i=currentallocated; i<allocated; i++)
    {
      ptr[i].dfsinfo = allocateDfsinfo(info,env);
    }
  }
  return ptr;
}

static void freeItvinfo(Itvinfo *ptr,
                        unsigned long allocated,
                        void (*freeDfsinfo)(Dfsinfo *,void *,Env *),
                        void *info,
                        Env *env)
{
  unsigned long i;

  for(i=0; i<allocated; i++)
  {
    freeDfsinfo(ptr[i].dfsinfo,info,env);
  }
  FREESPACE(ptr);
}

int depthfirstesa(Suffixarray *suffixarray,
                  Uchar initialchar,
                  Dfsinfo *(*allocateDfsinfo)(void *,Env *),
                  void(*freeDfsinfo)(Dfsinfo *,void *,Env *),
                  int(*processleafedge)(bool,Seqpos,Dfsinfo *,
                                        Uchar,Seqpos,void *,
                                        Env *),
                  int(*processbranchedge)(bool,
                                          Seqpos,
                                          Dfsinfo *,
                                          void *,
                                          Env *),
                  int(*processcompletenode)(Dfsinfo *,void *,Env *),
                  int(*assignleftmostleaf)(Dfsinfo *,Seqpos,
                                           void *,Env *),
                  int(*assignrightmostleaf)(Dfsinfo *,Seqpos,Seqpos,
                                            Seqpos,void *,Env *),
                  void *info,
                  Env *env)
{
  int retval;
  Uchar tmpsmalllcpvalue;
  bool firstedge,
       firstrootedge;
  Seqpos previoussuffix, 
         previouslcp,
         currentindex,
         currentlcp = 0; /* May be necessary if the lcpvalue is used after the
                            outer while loop */
  unsigned long allocatedItvinfo,
                nextfreeItvinfo;
  Largelcpvalue tmpexception;
  Itvinfo *stackspace;
  Uchar leftchar;
  bool haserr = false;
  
  allocatedItvinfo = (unsigned long) INCSTACKSIZE;
  stackspace = allocItvinfo(NULL,
                            0,
                            allocatedItvinfo,
                            allocateDfsinfo,
                            info,
                            env);
  nextfreeItvinfo = 0;
  firstrootedge = true;
  if (!haserr)
  {
    PUSHDFS(0,true);
    if(assignleftmostleaf != NULL && 
       assignleftmostleaf(TOP.dfsinfo,0,info,env) != 0)
    {
      haserr = true;
    }
  }
  for(currentindex = 0; !haserr; currentindex++)
  {
    retval = readnextUcharfromstream(&tmpsmalllcpvalue,
                                     &suffixarray->lcptabstream,
                                     env);
    if (retval < 0)
    {
      haserr = true;
      break;
    }
    if (retval == 0)
    {
      env_error_set(env,"file %s: line %d: unexpected end of file when "
                        "reading lcptab",__FILE__,__LINE__);
      haserr = true;
      break;
    }
    if (tmpsmalllcpvalue == (Uchar) UCHAR_MAX)
    {
      retval = readnextLargelcpvaluefromstream(&tmpexception,
                                               &suffixarray->llvtabstream,
                                               env);
      if (retval < 0)
      {
        haserr = true;
        break;
      }
      if (retval == 0)
      {
        env_error_set(env,"file %s: line %d: unexpected end of file when "
                      "reading llvtab",__FILE__,__LINE__);
        haserr = true;
        break;
      }
      currentlcp = tmpexception.value;
    } else
    {
      currentlcp = (Seqpos) tmpsmalllcpvalue;
    }
    retval = readnextSeqposfromstream(&previoussuffix,
                                      &suffixarray->suftabstream,
                                      env);
    if (retval < 0)
    {
      haserr = true;
      break;
    }
    if (retval == 0)
    {
      env_error_set(env,"file %s: line %d: unexpected end of file when "
                    "reading suftab",__FILE__,__LINE__);
      haserr = true;
      break;
    }
    if (previoussuffix == 0)
    {
      leftchar = initialchar;
    } else
    {
      leftchar = getencodedchar(suffixarray->encseq,
                                previoussuffix-1,
                                suffixarray->readmode);
    }
    while(currentlcp < TOP.depth)    // splitting edge is reached 
    {
      if (TOP.lastisleafedge)       // last edge from top-node is leaf
      {                            // previoussuffix
        if (processleafedge != NULL &&
            processleafedge(false,TOP.depth,TOP.dfsinfo,leftchar,
                            previoussuffix,info,env) != 0)
        {
          haserr = true;
          break;
        }
      } else
      {
        if (processbranchedge != NULL &&
            processbranchedge(false,
                              ABOVETOP.depth,
                              ABOVETOP.dfsinfo,
                              info,
                              env) != 0)
        {
          haserr = true;
          break;
        }
      }
      if (assignrightmostleaf != NULL && 
          assignrightmostleaf(TOP.dfsinfo,
                              currentindex,
                              previoussuffix,
                              currentlcp,
                              info,env) != 0)
      {
        haserr = true;
        break; 
      }
      if (processcompletenode != NULL && 
          processcompletenode(TOP.dfsinfo,info,env) != 0)
      {
        haserr = true;
        break;
      }
      assert(nextfreeItvinfo > 0);
      nextfreeItvinfo--;
    }
    if (haserr)
    {
      break;
    }
    if (currentlcp == TOP.depth)
    {
      // add leaf edge to TOP-node
      if (firstrootedge && TOP.depth == 0)
      {
        firstedge = true;
        firstrootedge = false;
      } else
      {
        firstedge = false;
      }
      if (TOP.lastisleafedge)
      {
        if (processleafedge != NULL &&
            processleafedge(firstedge,TOP.depth,TOP.dfsinfo,
                            leftchar,previoussuffix,info,
                            env) != 0)
        {
          haserr = true;
          break;
        }
      } else
      {
        if (processbranchedge != NULL &&
            processbranchedge(firstedge,
                              ABOVETOP.depth,
                              ABOVETOP.dfsinfo,
                              info,
                              env) != 0)
        {
          haserr = true;
          break;
        }
        TOP.lastisleafedge = true;
      }
    } else
    {
      previouslcp = ABOVETOP.depth;
      if (nextfreeItvinfo >= allocatedItvinfo)
      {
        stackspace = allocItvinfo(stackspace,
                                  allocatedItvinfo,
                                  allocatedItvinfo+INCSTACKSIZE,
                                  allocateDfsinfo,
                                  info,
                                  env);
        allocatedItvinfo += INCSTACKSIZE;
      }
      PUSHDFS(currentlcp,true);
      if (BELOWTOP.lastisleafedge)
      {
        // replace leaf edge by internal Edge
        if(assignleftmostleaf != NULL &&
           assignleftmostleaf(TOP.dfsinfo,currentindex,info,env) != 0)
        {
          haserr = true;
          break;
        }
        if (processleafedge != NULL &&
            processleafedge(true,TOP.depth,TOP.dfsinfo,
                            leftchar,previoussuffix,
                            info,env) != 0)
        {
          haserr = true;
          break;
        }
        BELOWTOP.lastisleafedge = false;
      } else
      {
        if (processbranchedge != NULL && 
            processbranchedge(true,
                              previouslcp,
                              TOP.dfsinfo,
                              info,
                              env) != 0)
        {
          haserr = true;
          break;
        }
      }
    }
  }
  if (TOP.lastisleafedge)
  {
    if (previoussuffix == 0)
    {
      leftchar = initialchar;
    } else
    {
      leftchar = getencodedchar(suffixarray->encseq,
                                previoussuffix-1,
                                suffixarray->readmode);
    }
    retval = readnextSeqposfromstream(&previoussuffix,
                                      &suffixarray->suftabstream,
                                      env);
    if (retval < 0)
    {
      haserr = true;
    } else
    {
      if (retval == 0)
      {
        env_error_set(env,"file %s: line %d: unexpected end of file when "
                      "reading suftab",__FILE__,__LINE__);
        haserr = true;
      }
    }
    if (!haserr)
    {
      if (processleafedge != NULL &&
          processleafedge(false,TOP.depth,TOP.dfsinfo,
                          leftchar,previoussuffix,
                          info,env) != 0)
      {
        haserr = true;
      }
    }
    if (!haserr)
    {
      if (assignrightmostleaf != NULL && 
          assignrightmostleaf(TOP.dfsinfo,
                              currentindex,
                              previoussuffix,
                              currentlcp,
                              info,env) != 0)
      {
        haserr = true;
      }
    }
    if (!haserr)
    {
      if (processcompletenode != NULL && 
          processcompletenode(TOP.dfsinfo,info,env) != 0)
      {
        haserr = true;
      }
    }
  }
  freeItvinfo(stackspace,
              allocatedItvinfo,
              freeDfsinfo,
              info,
              env);
  return haserr ? -1 : 0;
}
