#include <limits.h>
#include "sarr-def.h"
#include "seqpos-def.h"
#include "symboldef.h"
#include "spacedef.h"

#define STATICSTACKSPACE 512

#define ABOVETOP  stackptr[nextfreeItvinfo]
#define TOP       stackptr[nextfreeItvinfo-1]
#define BELOWTOP  stackptr[nextfreeItvinfo-2]

#define PUSHDFS(D,B)\
        stackptr[nextfreeItvinfo].depth = D;\
        stackptr[nextfreeItvinfo].lastisleafedge = B;\
        nextfreeItvinfo++

 DECLAREREADFUNCTION(Seqpos);

 DECLAREREADFUNCTION(Uchar);

 DECLAREREADFUNCTION(Largelcpvalue);

typedef struct
{
  uint32_t start, 
           length;
} Listtype;

typedef struct Dfsinfo Dfsinfo;

typedef struct
{
  bool lastisleafedge;
  Seqpos depth;
  Dfsinfo *dfsinfo;
} Itvinfo;

int depthfirstvstreegeneric(
                  Suffixarray *suffixarray,
                  Readmode readmode,
                  Uchar initialchar,
                  int(*allocateextrastackelements)(Itvinfo *,unsigned long,
                                                   void *,Env *),
                  int(*reallocateextrastackelements)(Itvinfo *,
                                                     unsigned long,
                                                     unsigned long,
                                                     void *,Env *),
                  void(*freenodestackspace)(Itvinfo *,unsigned long,
                                            void *,Env *),
                  int(*processleafedge)(bool,Itvinfo *,Uchar,Seqpos,void *,
                                        Env *),
                  int(*processbranchedge)(bool,
                                          Seqpos,
                                          Itvinfo *,
                                          void *,
                                          Env *),
                  int(*processcompletenode)(Itvinfo *,void *,Env *),
                  int(*assignleftmostleaf)(Itvinfo *,Seqpos,
                                           void *,Env *),
                  int(*assignrightmostleaf)(Itvinfo *,Seqpos,Seqpos,
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
  Itvinfo *stackptr,
           stackspace[STATICSTACKSPACE];
  Uchar leftchar;
  bool haserr = false;

  allocatedItvinfo = (unsigned long) STATICSTACKSPACE;
  nextfreeItvinfo = 0;
  firstrootedge = true;
  stackptr = &stackspace[0];
  if (allocateextrastackelements != NULL &&
      allocateextrastackelements(stackptr,(unsigned long) STATICSTACKSPACE,
                                 info,env) != 0)
  {
    haserr = true;
  }
  if (!haserr)
  {
    PUSHDFS(0,true);
    if(assignleftmostleaf != NULL && 
       assignleftmostleaf(&TOP,0,info,env) != 0)
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
                                readmode);
    }
    while(currentlcp < TOP.depth)    // splitting edge is reached 
    {
      if (TOP.lastisleafedge)       // last edge from top-node is leaf
      {                            // previoussuffix
        if (processleafedge != NULL &&
            processleafedge(false,&TOP,leftchar,previoussuffix,info,env) != 0)
        {
          haserr = true;
          break;
        }
      } else
      {
        if (processbranchedge != NULL &&
            processbranchedge(false,
                              ABOVETOP.depth,
                              &ABOVETOP,
                              info,
                              env) != 0)
        {
          haserr = true;
          break;
        }
      }
      if (assignrightmostleaf != NULL && 
          assignrightmostleaf(&TOP,
                              currentindex,
                              previoussuffix,
                              currentlcp,
                              info,env) != 0)
      {
        haserr = true;
        break; 
      }
      if (processcompletenode != NULL && 
          processcompletenode(&TOP,info,env) != 0)
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
            processleafedge(firstedge,&TOP,leftchar,previoussuffix,info,
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
                              &ABOVETOP,
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
        if (allocatedItvinfo == (unsigned long) STATICSTACKSPACE)
        {
          ALLOCASSIGNSPACE(stackptr,NULL,Itvinfo,
                           nextfreeItvinfo + STATICSTACKSPACE);
          memcpy(stackptr,&stackspace[0],
                 (size_t) (sizeof(Itvinfo) * STATICSTACKSPACE));
          if (reallocateextrastackelements != NULL &&
              reallocateextrastackelements(stackptr,
                                           (unsigned long) STATICSTACKSPACE,
                                           nextfreeItvinfo,info,env) != 0)
          {
            haserr = true;
            break;
          }
        } else
        {
          ALLOCASSIGNSPACE(stackptr,stackptr,Itvinfo,
                           allocatedItvinfo+STATICSTACKSPACE);
          if (reallocateextrastackelements != NULL &&
              reallocateextrastackelements(stackptr,
                                           allocatedItvinfo,
                                           (unsigned long) STATICSTACKSPACE,
                                           info,
                                           env) != 0)
          {
            haserr = true;
            break;
          }
        }
        allocatedItvinfo += STATICSTACKSPACE;
      }
      PUSHDFS(currentlcp,true);
      if (BELOWTOP.lastisleafedge)
      {
        // replace leaf edge by internal Edge
        if(assignleftmostleaf != NULL &&
           assignleftmostleaf(&TOP,currentindex,info,env) != 0)
        {
          haserr = true;
          break;
        }
        if (processleafedge != NULL &&
            processleafedge(true,&TOP,leftchar,previoussuffix,info,env) != 0)
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
                              &TOP,
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
                                readmode);
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
          processleafedge(false,&TOP,leftchar,previoussuffix,info,env) != 0)
      {
        haserr = true;
      }
    }
    if (!haserr)
    {
      if (assignrightmostleaf != NULL && 
          assignrightmostleaf(&TOP,
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
          processcompletenode(&TOP,info,env) != 0)
      {
        haserr = true;
      }
    }
  }
  if (freenodestackspace != NULL)
  {
    freenodestackspace(stackptr,allocatedItvinfo,info,env);
  }
  if (allocatedItvinfo > (unsigned long) STATICSTACKSPACE)
  {
    FREESPACE(stackptr);
  }
  return haserr ? -1 : 0;
}
