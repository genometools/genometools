#include "seqpos-def.h"

#define STATICSTACKSPACE 512

#define ABOVETOP  stackptr[nextfreeNodeinfo]
#define TOP       stackptr[nextfreeNodeinfo-1]
#define BELOWTOP  stackptr[nextfreeNodeinfo-2]

#define PUSHDFS(D,B)\
        stackptr[nextfreeNodeinfo].depth = D;\
        stackptr[nextfreeNodeinfo].lastisleafedge = B;\
        nextfreeNodeinfo++

 DECLAREREADFUNCTION(Seqpos);

 DECLAREREADFUNCTION(Uchar);

 DECLAREREADFUNCTION(Largelcpvalue);

typedef struct
{
  uint32_t start, 
           length;
} Listtype;

typedef struct
{
  bool lastisleafedge;
  Uchar commonchar;
  Seqpos depth;
  uint32_t uniquecharposstart,
           uniquecharposlength; // uniquecharpos[start..start+len-1]
  Listtype *nodeposlist;
  Seqpos leftmostleaf,
         rightmostleaf;
} Nodeinfo;

static int depthfirstvstreegeneric(
                  Suffixarray *suffixarray,
                  Readmode readmode,
                  int(*allocateextrastackelements)(Nodeinfo *,uint32_t,
                                                   void *,Env *),
                  int(*reallocateextrastackelements)(Nodeinfo *,
                                                     uint32_t,uint32_t,
                                                     void *,Env *),
                  void(*freenodestackspace)(Nodeinfo *,uint32_t,void *,Env *),
                  int(*processleafedge)(bool,Nodeinfo *,Uchar,Seqpos,void *,
                                        Env *),
                  int(*processbranchedge)(bool,
                                          uint32_t,
                                          Seqpos,
                                          Seqpos,
                                          void *,
                                          Env *),
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
         currentlcp = 0, /* May be necessary if the lcpvalue is used after the
                            outer while loop */
         allocatedNodeinfo,
         nextfreeNodeinfo;
  Largelcpvalue tmpexception;
  Nodeinfo *stackptr,
           stackspace[STATICSTACKSPACE];
  Uchar leftchar;
  bool haserr = false;

  allocatedNodeinfo = (Uint) STATICSTACKSPACE;
  nextfreeNodeinfo = 0;
  firstrootedge = true;
  stackptr = &stackspace[0];
  if(allocateextrastackelements != NULL &&
     allocateextrastackelements(stackptr,STATICSTACKSPACE,info,env) != 0)
  {
    haserr = true;
  }
  if(!haserr)
  {
    PUSHDFS(FIRSTLCPVALUE,true);
    ASSIGNLEFTMOSTLEAF(TOP,firstsuftabindex);
  }
  for(currentindex = firstsuftabindex; !haserr; currentindex++)
  {
    retval = readnextUcharfromstream(&tmpsmalllcpvalue,
                                     suffixarray->lcptabstream,
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
                                               suffixarray->llvtabstream,
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
                                      suffixarray->suftabstream,
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
    if(previoussuffix == 0)
    {
      leftchar = (Uchar) INITIALCHAR;
    } else
    {
      leftchar = getencodedchar(suffixarray->encseqtable,
                                previoussuffix-1,
                                readmode);
    }
    while(currentlcp < TOP.depth)    // splitting edge is reached 
    {
      if(TOP.lastisleafedge)       // last edge from top-node is leaf
      {                            // previoussuffix
        if(processleaf != NULL &&
           processleafedge(false,&TOP,leftchar,previoussuffix) != 0)
        {
          haserr = true;
          break;
        }
      } else
      {
        if(processbranchedge != NULL &&
           processbranchedge(false,
                             ABOVETOP.depth,
                             ABOVETOP.leftmostleaf,
                             ABOVETOP.rightmostleaf,
                             info,
                             env) != 0)
        {
          haserr = true;
          break;
        }
      }
      ASSIGNRIGHTMOSTLEAF(TOP,
                          currentindex,
                          previoussuffix,
                          currentlcp);
      PROCESSCOMPLETENODE(&TOP);
      assert(nextfreeNodeinfo > 0);
      nextfreeNodeinfo--;
    }
    if(haserr)
    {
      break;
    }
    if(currentlcp == TOP.depth)
    {
      // add leaf edge to TOP-node
      if(firstrootedge && TOP.depth == 0)
      {
        firstedge = true;
        firstrootedge = false;
      } else
      {
        firstedge = false;
      }
      if(TOP.lastisleafedge)
      {
        if(processleafedge != NULL &&
           processleafedge(firstedge,&TOP,leftchar,previoussuffix,info,
                           env) != 0)
        {
          haserr = true;
          break;
        }
      } else
      {
        if(processbranchedge != NULL &&
          processbranchedge(firstedge,
                             ABOVETOP.depth,
                             ABOVETOP.leftmostleaf,
                             ABOVETOP.rightmostleaf,
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
      if(nextfreeNodeinfo >= allocatedNodeinfo)
      {
        if(allocatedNodeinfo == (Uint) STATICSTACKSPACE)
        {
          ALLOCASSIGNSPACE(stackptr,NULL,Nodeinfo,
                           nextfreeNodeinfo + STATICSTACKSPACE);
          memcpy(stackptr,&stackspace[0],
                 (size_t) (sizeof(Nodeinfo) * STATICSTACKSPACE));
          if(reallocateextrastackelements != NULL &&
             reallocateextrastackelements(stackptr,
                                          STATICSTACKSPACE,
                                          nextfreeNodeinfo,info,env) != 0)
          {
            haserr = true;
            break;
          }
        } else
        {
          ALLOCASSIGNSPACE(stackptr,stackptr,Nodeinfo,
                           allocatedNodeinfo+STATICSTACKSPACE);
          if(reallocateextrastackelements != NULL &&
             reallocateextrastackelements(stackptr,
                                          allocatedNodeinfo,
                                          STATICSTACKSPACE,
                                          info,
                                          env) != 0)
          {
            haserr = true;
            break;
          }
        }
        allocatedNodeinfo += STATICSTACKSPACE;
      }
      PUSHDFS(currentlcp,true);
      if(BELOWTOP.lastisleafedge)
      {
        // replace leaf edge by internal Edge
        ASSIGNLEFTMOSTLEAF(TOP,currentindex);
        PROCESSSPLITLEAFEDGE;
        if(processleafedge(true,&TOP,leftchar,previoussuffix) != 0)
        {
          haserr = true;
          break;
        }
        BELOWTOP.lastisleafedge = false;
      } else
      {
        if(processbranchedge(true,previouslcp,
                             TOP.leftmostleaf,
                             TOP.rightmostleaf,
                             info,
                             env) != 0)
        {
          haserr = true;
          break;
        }
      }
    }
  }
  if(TOP.lastisleafedge)
  {
    if(previoussuffix == 0)
    {
      leftchar = (Uchar) INITIALCHAR;
    } else
    {
      leftchar = getencodedchar(suffixarray->encseqtable,
                                previoussuffix-1,
                                readmode);
    }
    retval = readnextSeqposfromstream(&previoussuffix,
                                      suffixarray->suftabstream,
                                      env);
    if (retval < 0)
    {
      haserr = true;
    }
    if (retval == 0)
    {
      env_error_set(env,"file %s: line %d: unexpected end of file when "
                    "reading suftab",__FILE__,__LINE__);
      haserr = true;
    }
    if(processleafedge(false,&TOP,leftchar,previoussuffix) != 0)
    {
      haserr = true;
    }
    ASSIGNRIGHTMOSTLEAF(TOP,
                        currentindex,
                        previoussuffix,
                        currentlcp);
    PROCESSCOMPLETENODE(&TOP);
  }
  if(freenodestackspace != NULL)
  {
    freenodestackspace(stackptr,allocatedNodeinfo,info,env);
  }
  if(allocatedNodeinfo > (Uint) STATICSTACKSPACE)
  {
    FREESPACE(stackptr);
  }
  return 0;
}
