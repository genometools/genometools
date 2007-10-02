#if 0
extern Seqpos
SRLSymbolCountInSeqRegion(struct seqRangeList *rangeList, Seqpos start,
                          Seqpos end, Symbol esym, seqRangeListSearchHint *hint)
{
  struct seqRange *p;
  if(rangeList->numRanges == 0)
    return 0;
  p = SRLFindPositionNext(rangeList, start, hint);
  if(p)
  {
    Seqpos symCount = 0;
    Seqpos s = MAX(start, p->startPos);
    Symbol sym = MRAEncMapSymbol(rangeList->alphabet, esym);
    struct seqRange *maxRange = rangeList->ranges + rangeList->numRanges - 1;
    while(s <= end)
    {
      if(p->sym == sym)
        symCount += MIN(p->startPos + p->len, end + 1) - s;
      if(p == maxRange)
        break;
      s = (++p)->startPos;
    }
    return symCount;
  }
  else
    return 0;
}

extern Seqpos
SRLAllSymbolsCountInSeqRegion(struct seqRangeList *rangeList, Seqpos start,
                              Seqpos end, seqRangeListSearchHint *hint)
{
  struct seqRange *p;
  if(rangeList->numRanges == 0)
    return 0;
  p = SRLFindPositionNext(rangeList, start, hint);
  if(p)
  {
    Seqpos symCount = 0;
    Seqpos s = MAX(start, p->startPos);
    struct seqRange *maxRange = rangeList->ranges + rangeList->numRanges - 1;
    while(s <= end)
    {
      symCount += MIN(p->startPos + p->len, end + 1) - s;
      if(p == maxRange)
        break;
      s = (++p)->startPos;
    }
    return symCount;
  }
  else
    return 0;
}
#else
