/* static int */
/* SfxIGetLCPCacheElem(sfxInterface *iface, Seqpos pos, Env *env) */
/* { */
/*   assert(iface && env); */
/*   assert(pos >= iface->LCPCacheStart); */
/*   if(pos < iface->LCPCacheStart + iface->LCPCacheLen) */
/*   { */
/*   } */
/*   else */
/*   { */
/*     /\* extend cache up to required position *\/ */
/*     Seqpos stillNeeded = findMinOpenRequest(iface, SFX_REQUEST_LCPTAB); */
/*     assert(stillNeeded >= iface->LCPCacheStart); */
    
/*   } */
  
/* } */

/* static int */
/* SfxICacheLCPValuesUpto(sfxInterface *iface, Seqpos pos, Env *env) */
/* { */
/*   size_t neededCacheSize = pos - iface->LCPCacheStart + 1, i; */
/*   Seqpos first, last; */
/*   assert(iface && env); */
/*   if(neededCacheSize > iface->LCPCacheSize) */
/*     iface->LCPCache = env_ma_realloc(env, iface->LCPCache, neededCacheSize */
/*                                      * sizeof(iface->LCPCache[0])); */
/*   first = iface->LCPCacheLen; */
/*   last = pos - iface->LCPCacheStart; */
/*   for(i = first; i <= last; ++i) */
/*   { */
/*     iface->LCPCache[i] */
/*       = nextLcpvalueiterator(iface->lvi, */
/*                              (iface->outfileinfo.pageoffset == 0?true:false), */
                             
/*   } */
/* } */
