static GtUword ft_longest_common_twobit_twobit(
                                      GtFtSequenceObject *useq,
                                      GtUword ustart,
                                      GtFtSequenceObject *vseq,
                                      const GtUword vstart)
{
  if (ustart < useq->substringlength && vstart < vseq->substringlength)
  {
    GtUword uptr;int ustep;
    GtUword vptr;int vstep;
    GtUword minsubstringlength = useq->substringlength - ustart, matchlength;

    if (useq->read_seq_left2right)
    {
      uptr = useq->offset + ustart; ustep = 1;
    } else
    {
      uptr = useq->offset - ustart; ustep = -1;
    }
    if (vseq->read_seq_left2right)
    {
      vptr = vseq->offset + vstart; vstep = 1;
    } else
    {
      vptr = vseq->offset - vstart; vstep = -1;
    }
    if (vseq->substringlength - vstart < minsubstringlength)
    {
      minsubstringlength = vseq->substringlength - vstart;
    }
    if (vseq->dir_is_complement)
    {
      for (matchlength = 0; matchlength < minsubstringlength;
           uptr += ustep, vptr += vstep, matchlength++)
      {
        const GtUchar cu = gt_twobitencoding_char_at_pos(
                              useq->twobitencoding,
                              uptr);
        if (cu != GT_COMPLEMENTBASE(gt_twobitencoding_char_at_pos(
                              vseq->twobitencoding,
                              vptr)))
          break;
      }
    } else
    {
      for (matchlength = 0; matchlength < minsubstringlength;
           uptr += ustep, vptr += vstep, matchlength++)
      {
        const GtUchar cu = gt_twobitencoding_char_at_pos(
                              useq->twobitencoding,
                              uptr);
        if (cu != gt_twobitencoding_char_at_pos(
                              vseq->twobitencoding,
                              vptr))
          break;
      }
    }
    return matchlength;
  }
  return 0;
}

static GtUword ft_longest_common_twobit_encseq_reader(
                                      GtFtSequenceObject *useq,
                                      GtUword ustart,
                                      GtFtSequenceObject *vseq,
                                      const GtUword vstart)
{
  if (ustart < useq->substringlength && vstart < vseq->substringlength)
  {
    GtUword uptr;int ustep;
    GtUword vptr;
    GtUword minsubstringlength = useq->substringlength - ustart, matchlength;

    if (useq->read_seq_left2right)
    {
      uptr = useq->offset + ustart; ustep = 1;
    } else
    {
      uptr = useq->offset - ustart; ustep = -1;
    }
    vptr = vstart;
    if (vseq->substringlength - vstart < minsubstringlength)
    {
      minsubstringlength = vseq->substringlength - vstart;
    }
    if (vseq->dir_is_complement)
    {
      for (matchlength = 0; matchlength < minsubstringlength;
           uptr += ustep, vptr++, matchlength++)
      {
        const GtUchar cu = gt_twobitencoding_char_at_pos(
                              useq->twobitencoding,
                              uptr);
        if (cu != GT_COMPLEMENTBASE(gt_sequenceobject_esr_get(vseq,vptr)))
          break;
      }
    } else
    {
      for (matchlength = 0; matchlength < minsubstringlength;
           uptr += ustep, vptr++, matchlength++)
      {
        const GtUchar cu = gt_twobitencoding_char_at_pos(
                              useq->twobitencoding,
                              uptr);
        if (cu != gt_sequenceobject_esr_get(vseq,vptr))
          break;
      }
    }
    return matchlength;
  }
  return 0;
}

static GtUword ft_longest_common_twobit_encseq(
                                      GtFtSequenceObject *useq,
                                      GtUword ustart,
                                      GtFtSequenceObject *vseq,
                                      const GtUword vstart)
{
  if (ustart < useq->substringlength && vstart < vseq->substringlength)
  {
    GtUword uptr;int ustep;
    GtUword vptr;int vstep;
    GtUword minsubstringlength = useq->substringlength - ustart, matchlength;

    if (useq->read_seq_left2right)
    {
      uptr = useq->offset + ustart; ustep = 1;
    } else
    {
      uptr = useq->offset - ustart; ustep = -1;
    }
    if (vseq->read_seq_left2right)
    {
      vptr = vseq->offset + vstart; vstep = 1;
    } else
    {
      vptr = vseq->offset - vstart; vstep = -1;
    }
    if (vseq->substringlength - vstart < minsubstringlength)
    {
      minsubstringlength = vseq->substringlength - vstart;
    }
    if (vseq->dir_is_complement)
    {
      for (matchlength = 0; matchlength < minsubstringlength;
           uptr += ustep, vptr += vstep, matchlength++)
      {
        const GtUchar cu = gt_twobitencoding_char_at_pos(
                              useq->twobitencoding,
                              uptr);
        if (cu != GT_COMPLEMENTBASE(gt_encseq_get_encoded_char(vseq->encseq,
                    vptr,
                    GT_READMODE_FORWARD)))
          break;
      }
    } else
    {
      for (matchlength = 0; matchlength < minsubstringlength;
           uptr += ustep, vptr += vstep, matchlength++)
      {
        const GtUchar cu = gt_twobitencoding_char_at_pos(
                              useq->twobitencoding,
                              uptr);
        if (cu != gt_encseq_get_encoded_char(vseq->encseq,
                    vptr,
                    GT_READMODE_FORWARD))
          break;
      }
    }
    return matchlength;
  }
  return 0;
}

static GtUword ft_longest_common_twobit_bytes(
                                      GtFtSequenceObject *useq,
                                      GtUword ustart,
                                      GtFtSequenceObject *vseq,
                                      const GtUword vstart)
{
  if (ustart < useq->substringlength && vstart < vseq->substringlength)
  {
    GtUword uptr;int ustep;
    const GtUchar *vptr;int vstep;
    GtUword minsubstringlength = useq->substringlength - ustart, matchlength;

    if (useq->read_seq_left2right)
    {
      uptr = useq->offset + ustart; ustep = 1;
    } else
    {
      uptr = useq->offset - ustart; ustep = -1;
    }
    if (vseq->read_seq_left2right)
    {
      vptr = vseq->bytesequenceptr + vseq->offset + vstart; vstep = 1;
    } else
    {
      vptr = vseq->bytesequenceptr + vseq->offset - vstart; vstep = -1;
    }
    if (vseq->substringlength - vstart < minsubstringlength)
    {
      minsubstringlength = vseq->substringlength - vstart;
    }
    if (vseq->dir_is_complement)
    {
      for (matchlength = 0; matchlength < minsubstringlength;
           uptr += ustep, vptr += vstep, matchlength++)
      {
        const GtUchar cu = gt_twobitencoding_char_at_pos(
                              useq->twobitencoding,
                              uptr);
        if (cu != GT_COMPLEMENTBASE(*vptr))
          break;
      }
    } else
    {
      for (matchlength = 0; matchlength < minsubstringlength;
           uptr += ustep, vptr += vstep, matchlength++)
      {
        const GtUchar cu = gt_twobitencoding_char_at_pos(
                              useq->twobitencoding,
                              uptr);
        if (cu != *vptr)
          break;
      }
    }
    return matchlength;
  }
  return 0;
}

static GtUword ft_longest_common_encseq_reader_twobit(
                                      GtFtSequenceObject *useq,
                                      GtUword ustart,
                                      GtFtSequenceObject *vseq,
                                      const GtUword vstart)
{
  if (ustart < useq->substringlength && vstart < vseq->substringlength)
  {
    GtUword uptr;
    GtUword vptr;int vstep;
    GtUword minsubstringlength = useq->substringlength - ustart, matchlength;

    uptr = ustart;
    if (vseq->read_seq_left2right)
    {
      vptr = vseq->offset + vstart; vstep = 1;
    } else
    {
      vptr = vseq->offset - vstart; vstep = -1;
    }
    if (vseq->substringlength - vstart < minsubstringlength)
    {
      minsubstringlength = vseq->substringlength - vstart;
    }
    if (vseq->dir_is_complement)
    {
      for (matchlength = 0; matchlength < minsubstringlength;
           uptr++, vptr += vstep, matchlength++)
      {
        const GtUchar cu = gt_sequenceobject_esr_get(useq,uptr);
        if (cu != GT_COMPLEMENTBASE(gt_twobitencoding_char_at_pos(
                              vseq->twobitencoding,
                              vptr)))
          break;
      }
    } else
    {
      for (matchlength = 0; matchlength < minsubstringlength;
           uptr++, vptr += vstep, matchlength++)
      {
        const GtUchar cu = gt_sequenceobject_esr_get(useq,uptr);
        if (cu != gt_twobitencoding_char_at_pos(
                              vseq->twobitencoding,
                              vptr))
          break;
      }
    }
    return matchlength;
  }
  return 0;
}

static GtUword ft_longest_common_encseq_reader_encseq_reader(
                                      GtFtSequenceObject *useq,
                                      GtUword ustart,
                                      GtFtSequenceObject *vseq,
                                      const GtUword vstart)
{
  if (ustart < useq->substringlength && vstart < vseq->substringlength)
  {
    GtUword uptr;
    GtUword vptr;
    GtUword minsubstringlength = useq->substringlength - ustart, matchlength;

    uptr = ustart;
    vptr = vstart;
    if (vseq->substringlength - vstart < minsubstringlength)
    {
      minsubstringlength = vseq->substringlength - vstart;
    }
    if (vseq->dir_is_complement)
    {
      for (matchlength = 0; matchlength < minsubstringlength;
           uptr++, vptr++, matchlength++)
      {
        const GtUchar cu = gt_sequenceobject_esr_get(useq,uptr);
        if (cu != GT_COMPLEMENTBASE(gt_sequenceobject_esr_get(vseq,vptr)))
          break;
      }
    } else
    {
      for (matchlength = 0; matchlength < minsubstringlength;
           uptr++, vptr++, matchlength++)
      {
        const GtUchar cu = gt_sequenceobject_esr_get(useq,uptr);
        if (cu != gt_sequenceobject_esr_get(vseq,vptr))
          break;
      }
    }
    return matchlength;
  }
  return 0;
}

static GtUword ft_longest_common_encseq_reader_encseq(
                                      GtFtSequenceObject *useq,
                                      GtUword ustart,
                                      GtFtSequenceObject *vseq,
                                      const GtUword vstart)
{
  if (ustart < useq->substringlength && vstart < vseq->substringlength)
  {
    GtUword uptr;
    GtUword vptr;int vstep;
    GtUword minsubstringlength = useq->substringlength - ustart, matchlength;

    uptr = ustart;
    if (vseq->read_seq_left2right)
    {
      vptr = vseq->offset + vstart; vstep = 1;
    } else
    {
      vptr = vseq->offset - vstart; vstep = -1;
    }
    if (vseq->substringlength - vstart < minsubstringlength)
    {
      minsubstringlength = vseq->substringlength - vstart;
    }
    if (vseq->dir_is_complement)
    {
      for (matchlength = 0; matchlength < minsubstringlength;
           uptr++, vptr += vstep, matchlength++)
      {
        const GtUchar cu = gt_sequenceobject_esr_get(useq,uptr);
        if (cu != GT_COMPLEMENTBASE(gt_encseq_get_encoded_char(vseq->encseq,
                    vptr,
                    GT_READMODE_FORWARD)))
          break;
      }
    } else
    {
      for (matchlength = 0; matchlength < minsubstringlength;
           uptr++, vptr += vstep, matchlength++)
      {
        const GtUchar cu = gt_sequenceobject_esr_get(useq,uptr);
        if (cu != gt_encseq_get_encoded_char(vseq->encseq,
                    vptr,
                    GT_READMODE_FORWARD))
          break;
      }
    }
    return matchlength;
  }
  return 0;
}

static GtUword ft_longest_common_encseq_reader_bytes(
                                      GtFtSequenceObject *useq,
                                      GtUword ustart,
                                      GtFtSequenceObject *vseq,
                                      const GtUword vstart)
{
  if (ustart < useq->substringlength && vstart < vseq->substringlength)
  {
    GtUword uptr;
    const GtUchar *vptr;int vstep;
    GtUword minsubstringlength = useq->substringlength - ustart, matchlength;

    uptr = ustart;
    if (vseq->read_seq_left2right)
    {
      vptr = vseq->bytesequenceptr + vseq->offset + vstart; vstep = 1;
    } else
    {
      vptr = vseq->bytesequenceptr + vseq->offset - vstart; vstep = -1;
    }
    if (vseq->substringlength - vstart < minsubstringlength)
    {
      minsubstringlength = vseq->substringlength - vstart;
    }
    if (vseq->dir_is_complement)
    {
      for (matchlength = 0; matchlength < minsubstringlength;
           uptr++, vptr += vstep, matchlength++)
      {
        const GtUchar cu = gt_sequenceobject_esr_get(useq,uptr);
        if (cu != GT_COMPLEMENTBASE(*vptr))
          break;
      }
    } else
    {
      for (matchlength = 0; matchlength < minsubstringlength;
           uptr++, vptr += vstep, matchlength++)
      {
        const GtUchar cu = gt_sequenceobject_esr_get(useq,uptr);
        if (cu != *vptr)
          break;
      }
    }
    return matchlength;
  }
  return 0;
}

static GtUword ft_longest_common_encseq_twobit(
                                      GtFtSequenceObject *useq,
                                      GtUword ustart,
                                      GtFtSequenceObject *vseq,
                                      const GtUword vstart)
{
  if (ustart < useq->substringlength && vstart < vseq->substringlength)
  {
    GtUword uptr;int ustep;
    GtUword vptr;int vstep;
    GtUword minsubstringlength = useq->substringlength - ustart, matchlength;

    if (useq->read_seq_left2right)
    {
      uptr = useq->offset + ustart; ustep = 1;
    } else
    {
      uptr = useq->offset - ustart; ustep = -1;
    }
    if (vseq->read_seq_left2right)
    {
      vptr = vseq->offset + vstart; vstep = 1;
    } else
    {
      vptr = vseq->offset - vstart; vstep = -1;
    }
    if (vseq->substringlength - vstart < minsubstringlength)
    {
      minsubstringlength = vseq->substringlength - vstart;
    }
    if (vseq->dir_is_complement)
    {
      for (matchlength = 0; matchlength < minsubstringlength;
           uptr += ustep, vptr += vstep, matchlength++)
      {
        const GtUchar cu = gt_encseq_get_encoded_char(useq->encseq,
                    uptr,
                    GT_READMODE_FORWARD);
        if (cu != GT_COMPLEMENTBASE(gt_twobitencoding_char_at_pos(
                              vseq->twobitencoding,
                              vptr)))
          break;
      }
    } else
    {
      for (matchlength = 0; matchlength < minsubstringlength;
           uptr += ustep, vptr += vstep, matchlength++)
      {
        const GtUchar cu = gt_encseq_get_encoded_char(useq->encseq,
                    uptr,
                    GT_READMODE_FORWARD);
        if (cu != gt_twobitencoding_char_at_pos(
                              vseq->twobitencoding,
                              vptr))
          break;
      }
    }
    return matchlength;
  }
  return 0;
}

static GtUword ft_longest_common_encseq_encseq_reader(
                                      GtFtSequenceObject *useq,
                                      GtUword ustart,
                                      GtFtSequenceObject *vseq,
                                      const GtUword vstart)
{
  if (ustart < useq->substringlength && vstart < vseq->substringlength)
  {
    GtUword uptr;int ustep;
    GtUword vptr;
    GtUword minsubstringlength = useq->substringlength - ustart, matchlength;

    if (useq->read_seq_left2right)
    {
      uptr = useq->offset + ustart; ustep = 1;
    } else
    {
      uptr = useq->offset - ustart; ustep = -1;
    }
    vptr = vstart;
    if (vseq->substringlength - vstart < minsubstringlength)
    {
      minsubstringlength = vseq->substringlength - vstart;
    }
    if (vseq->dir_is_complement)
    {
      for (matchlength = 0; matchlength < minsubstringlength;
           uptr += ustep, vptr++, matchlength++)
      {
        const GtUchar cu = gt_encseq_get_encoded_char(useq->encseq,
                    uptr,
                    GT_READMODE_FORWARD);
        if (cu != GT_COMPLEMENTBASE(gt_sequenceobject_esr_get(vseq,vptr)))
          break;
      }
    } else
    {
      for (matchlength = 0; matchlength < minsubstringlength;
           uptr += ustep, vptr++, matchlength++)
      {
        const GtUchar cu = gt_encseq_get_encoded_char(useq->encseq,
                    uptr,
                    GT_READMODE_FORWARD);
        if (cu != gt_sequenceobject_esr_get(vseq,vptr))
          break;
      }
    }
    return matchlength;
  }
  return 0;
}

static GtUword ft_longest_common_encseq_encseq(
                                      GtFtSequenceObject *useq,
                                      GtUword ustart,
                                      GtFtSequenceObject *vseq,
                                      const GtUword vstart)
{
  if (ustart < useq->substringlength && vstart < vseq->substringlength)
  {
    GtUword uptr;int ustep;
    GtUword vptr;int vstep;
    GtUword minsubstringlength = useq->substringlength - ustart, matchlength;

    if (useq->read_seq_left2right)
    {
      uptr = useq->offset + ustart; ustep = 1;
    } else
    {
      uptr = useq->offset - ustart; ustep = -1;
    }
    if (vseq->read_seq_left2right)
    {
      vptr = vseq->offset + vstart; vstep = 1;
    } else
    {
      vptr = vseq->offset - vstart; vstep = -1;
    }
    if (vseq->substringlength - vstart < minsubstringlength)
    {
      minsubstringlength = vseq->substringlength - vstart;
    }
    if (vseq->dir_is_complement)
    {
      for (matchlength = 0; matchlength < minsubstringlength;
           uptr += ustep, vptr += vstep, matchlength++)
      {
        const GtUchar cu = gt_encseq_get_encoded_char(useq->encseq,
                    uptr,
                    GT_READMODE_FORWARD);
        if (cu != GT_COMPLEMENTBASE(gt_encseq_get_encoded_char(vseq->encseq,
                    vptr,
                    GT_READMODE_FORWARD)))
          break;
      }
    } else
    {
      for (matchlength = 0; matchlength < minsubstringlength;
           uptr += ustep, vptr += vstep, matchlength++)
      {
        const GtUchar cu = gt_encseq_get_encoded_char(useq->encseq,
                    uptr,
                    GT_READMODE_FORWARD);
        if (cu != gt_encseq_get_encoded_char(vseq->encseq,
                    vptr,
                    GT_READMODE_FORWARD))
          break;
      }
    }
    return matchlength;
  }
  return 0;
}

static GtUword ft_longest_common_encseq_bytes(
                                      GtFtSequenceObject *useq,
                                      GtUword ustart,
                                      GtFtSequenceObject *vseq,
                                      const GtUword vstart)
{
  if (ustart < useq->substringlength && vstart < vseq->substringlength)
  {
    GtUword uptr;int ustep;
    const GtUchar *vptr;int vstep;
    GtUword minsubstringlength = useq->substringlength - ustart, matchlength;

    if (useq->read_seq_left2right)
    {
      uptr = useq->offset + ustart; ustep = 1;
    } else
    {
      uptr = useq->offset - ustart; ustep = -1;
    }
    if (vseq->read_seq_left2right)
    {
      vptr = vseq->bytesequenceptr + vseq->offset + vstart; vstep = 1;
    } else
    {
      vptr = vseq->bytesequenceptr + vseq->offset - vstart; vstep = -1;
    }
    if (vseq->substringlength - vstart < minsubstringlength)
    {
      minsubstringlength = vseq->substringlength - vstart;
    }
    if (vseq->dir_is_complement)
    {
      for (matchlength = 0; matchlength < minsubstringlength;
           uptr += ustep, vptr += vstep, matchlength++)
      {
        const GtUchar cu = gt_encseq_get_encoded_char(useq->encseq,
                    uptr,
                    GT_READMODE_FORWARD);
        if (cu != GT_COMPLEMENTBASE(*vptr))
          break;
      }
    } else
    {
      for (matchlength = 0; matchlength < minsubstringlength;
           uptr += ustep, vptr += vstep, matchlength++)
      {
        const GtUchar cu = gt_encseq_get_encoded_char(useq->encseq,
                    uptr,
                    GT_READMODE_FORWARD);
        if (cu != *vptr)
          break;
      }
    }
    return matchlength;
  }
  return 0;
}

static GtUword ft_longest_common_bytes_twobit(
                                      GtFtSequenceObject *useq,
                                      GtUword ustart,
                                      GtFtSequenceObject *vseq,
                                      const GtUword vstart)
{
  if (ustart < useq->substringlength && vstart < vseq->substringlength)
  {
    const GtUchar *uptr;int ustep;
    GtUword vptr;int vstep;
    GtUword minsubstringlength = useq->substringlength - ustart, matchlength;

    if (useq->read_seq_left2right)
    {
      uptr = useq->bytesequenceptr + useq->offset + ustart; ustep = 1;
    } else
    {
      uptr = useq->bytesequenceptr + useq->offset - ustart; ustep = -1;
    }
    if (vseq->read_seq_left2right)
    {
      vptr = vseq->offset + vstart; vstep = 1;
    } else
    {
      vptr = vseq->offset - vstart; vstep = -1;
    }
    if (vseq->substringlength - vstart < minsubstringlength)
    {
      minsubstringlength = vseq->substringlength - vstart;
    }
    if (vseq->dir_is_complement)
    {
      for (matchlength = 0; matchlength < minsubstringlength;
           uptr += ustep, vptr += vstep, matchlength++)
      {
        const GtUchar cu = *uptr;
        if (cu != GT_COMPLEMENTBASE(gt_twobitencoding_char_at_pos(
                              vseq->twobitencoding,
                              vptr)))
          break;
      }
    } else
    {
      for (matchlength = 0; matchlength < minsubstringlength;
           uptr += ustep, vptr += vstep, matchlength++)
      {
        const GtUchar cu = *uptr;
        if (cu != gt_twobitencoding_char_at_pos(
                              vseq->twobitencoding,
                              vptr))
          break;
      }
    }
    return matchlength;
  }
  return 0;
}

static GtUword ft_longest_common_bytes_encseq_reader(
                                      GtFtSequenceObject *useq,
                                      GtUword ustart,
                                      GtFtSequenceObject *vseq,
                                      const GtUword vstart)
{
  if (ustart < useq->substringlength && vstart < vseq->substringlength)
  {
    const GtUchar *uptr;int ustep;
    GtUword vptr;
    GtUword minsubstringlength = useq->substringlength - ustart, matchlength;

    if (useq->read_seq_left2right)
    {
      uptr = useq->bytesequenceptr + useq->offset + ustart; ustep = 1;
    } else
    {
      uptr = useq->bytesequenceptr + useq->offset - ustart; ustep = -1;
    }
    vptr = vstart;
    if (vseq->substringlength - vstart < minsubstringlength)
    {
      minsubstringlength = vseq->substringlength - vstart;
    }
    if (vseq->dir_is_complement)
    {
      for (matchlength = 0; matchlength < minsubstringlength;
           uptr += ustep, vptr++, matchlength++)
      {
        const GtUchar cu = *uptr;
        if (cu != GT_COMPLEMENTBASE(gt_sequenceobject_esr_get(vseq,vptr)))
          break;
      }
    } else
    {
      for (matchlength = 0; matchlength < minsubstringlength;
           uptr += ustep, vptr++, matchlength++)
      {
        const GtUchar cu = *uptr;
        if (cu != gt_sequenceobject_esr_get(vseq,vptr))
          break;
      }
    }
    return matchlength;
  }
  return 0;
}

static GtUword ft_longest_common_bytes_encseq(
                                      GtFtSequenceObject *useq,
                                      GtUword ustart,
                                      GtFtSequenceObject *vseq,
                                      const GtUword vstart)
{
  if (ustart < useq->substringlength && vstart < vseq->substringlength)
  {
    const GtUchar *uptr;int ustep;
    GtUword vptr;int vstep;
    GtUword minsubstringlength = useq->substringlength - ustart, matchlength;

    if (useq->read_seq_left2right)
    {
      uptr = useq->bytesequenceptr + useq->offset + ustart; ustep = 1;
    } else
    {
      uptr = useq->bytesequenceptr + useq->offset - ustart; ustep = -1;
    }
    if (vseq->read_seq_left2right)
    {
      vptr = vseq->offset + vstart; vstep = 1;
    } else
    {
      vptr = vseq->offset - vstart; vstep = -1;
    }
    if (vseq->substringlength - vstart < minsubstringlength)
    {
      minsubstringlength = vseq->substringlength - vstart;
    }
    if (vseq->dir_is_complement)
    {
      for (matchlength = 0; matchlength < minsubstringlength;
           uptr += ustep, vptr += vstep, matchlength++)
      {
        const GtUchar cu = *uptr;
        if (cu != GT_COMPLEMENTBASE(gt_encseq_get_encoded_char(vseq->encseq,
                    vptr,
                    GT_READMODE_FORWARD)))
          break;
      }
    } else
    {
      for (matchlength = 0; matchlength < minsubstringlength;
           uptr += ustep, vptr += vstep, matchlength++)
      {
        const GtUchar cu = *uptr;
        if (cu != gt_encseq_get_encoded_char(vseq->encseq,
                    vptr,
                    GT_READMODE_FORWARD))
          break;
      }
    }
    return matchlength;
  }
  return 0;
}

static GtUword ft_longest_common_bytes_bytes(
                                      GtFtSequenceObject *useq,
                                      GtUword ustart,
                                      GtFtSequenceObject *vseq,
                                      const GtUword vstart)
{
  if (ustart < useq->substringlength && vstart < vseq->substringlength)
  {
    const GtUchar *uptr;int ustep;
    const GtUchar *vptr;int vstep;
    GtUword minsubstringlength = useq->substringlength - ustart, matchlength;

    if (useq->read_seq_left2right)
    {
      uptr = useq->bytesequenceptr + useq->offset + ustart; ustep = 1;
    } else
    {
      uptr = useq->bytesequenceptr + useq->offset - ustart; ustep = -1;
    }
    if (vseq->read_seq_left2right)
    {
      vptr = vseq->bytesequenceptr + vseq->offset + vstart; vstep = 1;
    } else
    {
      vptr = vseq->bytesequenceptr + vseq->offset - vstart; vstep = -1;
    }
    if (vseq->substringlength - vstart < minsubstringlength)
    {
      minsubstringlength = vseq->substringlength - vstart;
    }
    if (vseq->dir_is_complement)
    {
      for (matchlength = 0; matchlength < minsubstringlength;
           uptr += ustep, vptr += vstep, matchlength++)
      {
        const GtUchar cu = *uptr;
        if (cu != GT_COMPLEMENTBASE(*vptr))
          break;
      }
    } else
    {
      for (matchlength = 0; matchlength < minsubstringlength;
           uptr += ustep, vptr += vstep, matchlength++)
      {
        const GtUchar cu = *uptr;
        if (cu != *vptr)
          break;
      }
    }
    return matchlength;
  }
  return 0;
}

static GtUword ft_longest_common_twobit_twobit_wildcard(
                                      GtFtSequenceObject *useq,
                                      GtUword ustart,
                                      GtFtSequenceObject *vseq,
                                      const GtUword vstart)
{
  if (ustart < useq->substringlength && vstart < vseq->substringlength)
  {
    GtUword uptr;int ustep;
    GtUword vptr;int vstep;
    GtUword minsubstringlength = useq->substringlength - ustart, matchlength;

    if (useq->read_seq_left2right)
    {
      uptr = useq->offset + ustart; ustep = 1;
    } else
    {
      uptr = useq->offset - ustart; ustep = -1;
    }
    if (vseq->read_seq_left2right)
    {
      vptr = vseq->offset + vstart; vstep = 1;
    } else
    {
      vptr = vseq->offset - vstart; vstep = -1;
    }
    if (vseq->substringlength - vstart < minsubstringlength)
    {
      minsubstringlength = vseq->substringlength - vstart;
    }
    if (vseq->dir_is_complement)
    {
      for (matchlength = 0; matchlength < minsubstringlength;
           uptr += ustep, vptr += vstep, matchlength++)
      {
        const GtUchar cu = gt_twobitencoding_char_at_pos(
                              useq->twobitencoding,
                              uptr);
        if (cu == WILDCARD ||
            cu != GT_COMPLEMENTBASE(gt_twobitencoding_char_at_pos(
                              vseq->twobitencoding,
                              vptr)))
          break;
      }
    } else
    {
      for (matchlength = 0; matchlength < minsubstringlength;
           uptr += ustep, vptr += vstep, matchlength++)
      {
        const GtUchar cu = gt_twobitencoding_char_at_pos(
                              useq->twobitencoding,
                              uptr);
        if (cu == WILDCARD ||
            cu != gt_twobitencoding_char_at_pos(
                              vseq->twobitencoding,
                              vptr))
          break;
      }
    }
    return matchlength;
  }
  return 0;
}

static GtUword ft_longest_common_twobit_encseq_reader_wildcard(
                                      GtFtSequenceObject *useq,
                                      GtUword ustart,
                                      GtFtSequenceObject *vseq,
                                      const GtUword vstart)
{
  if (ustart < useq->substringlength && vstart < vseq->substringlength)
  {
    GtUword uptr;int ustep;
    GtUword vptr;
    GtUword minsubstringlength = useq->substringlength - ustart, matchlength;

    if (useq->read_seq_left2right)
    {
      uptr = useq->offset + ustart; ustep = 1;
    } else
    {
      uptr = useq->offset - ustart; ustep = -1;
    }
    vptr = vstart;
    if (vseq->substringlength - vstart < minsubstringlength)
    {
      minsubstringlength = vseq->substringlength - vstart;
    }
    if (vseq->dir_is_complement)
    {
      for (matchlength = 0; matchlength < minsubstringlength;
           uptr += ustep, vptr++, matchlength++)
      {
        const GtUchar cu = gt_twobitencoding_char_at_pos(
                              useq->twobitencoding,
                              uptr);
        if (cu == WILDCARD ||
            cu != GT_COMPLEMENTBASE(gt_sequenceobject_esr_get(vseq,vptr)))
          break;
      }
    } else
    {
      for (matchlength = 0; matchlength < minsubstringlength;
           uptr += ustep, vptr++, matchlength++)
      {
        const GtUchar cu = gt_twobitencoding_char_at_pos(
                              useq->twobitencoding,
                              uptr);
        if (cu == WILDCARD ||
            cu != gt_sequenceobject_esr_get(vseq,vptr))
          break;
      }
    }
    return matchlength;
  }
  return 0;
}

static GtUword ft_longest_common_twobit_encseq_wildcard(
                                      GtFtSequenceObject *useq,
                                      GtUword ustart,
                                      GtFtSequenceObject *vseq,
                                      const GtUword vstart)
{
  if (ustart < useq->substringlength && vstart < vseq->substringlength)
  {
    GtUword uptr;int ustep;
    GtUword vptr;int vstep;
    GtUword minsubstringlength = useq->substringlength - ustart, matchlength;

    if (useq->read_seq_left2right)
    {
      uptr = useq->offset + ustart; ustep = 1;
    } else
    {
      uptr = useq->offset - ustart; ustep = -1;
    }
    if (vseq->read_seq_left2right)
    {
      vptr = vseq->offset + vstart; vstep = 1;
    } else
    {
      vptr = vseq->offset - vstart; vstep = -1;
    }
    if (vseq->substringlength - vstart < minsubstringlength)
    {
      minsubstringlength = vseq->substringlength - vstart;
    }
    if (vseq->dir_is_complement)
    {
      for (matchlength = 0; matchlength < minsubstringlength;
           uptr += ustep, vptr += vstep, matchlength++)
      {
        const GtUchar cu = gt_twobitencoding_char_at_pos(
                              useq->twobitencoding,
                              uptr);
        if (cu == WILDCARD ||
            cu != GT_COMPLEMENTBASE(gt_encseq_get_encoded_char(vseq->encseq,
                    vptr,
                    GT_READMODE_FORWARD)))
          break;
      }
    } else
    {
      for (matchlength = 0; matchlength < minsubstringlength;
           uptr += ustep, vptr += vstep, matchlength++)
      {
        const GtUchar cu = gt_twobitencoding_char_at_pos(
                              useq->twobitencoding,
                              uptr);
        if (cu == WILDCARD ||
            cu != gt_encseq_get_encoded_char(vseq->encseq,
                    vptr,
                    GT_READMODE_FORWARD))
          break;
      }
    }
    return matchlength;
  }
  return 0;
}

static GtUword ft_longest_common_twobit_bytes_wildcard(
                                      GtFtSequenceObject *useq,
                                      GtUword ustart,
                                      GtFtSequenceObject *vseq,
                                      const GtUword vstart)
{
  if (ustart < useq->substringlength && vstart < vseq->substringlength)
  {
    GtUword uptr;int ustep;
    const GtUchar *vptr;int vstep;
    GtUword minsubstringlength = useq->substringlength - ustart, matchlength;

    if (useq->read_seq_left2right)
    {
      uptr = useq->offset + ustart; ustep = 1;
    } else
    {
      uptr = useq->offset - ustart; ustep = -1;
    }
    if (vseq->read_seq_left2right)
    {
      vptr = vseq->bytesequenceptr + vseq->offset + vstart; vstep = 1;
    } else
    {
      vptr = vseq->bytesequenceptr + vseq->offset - vstart; vstep = -1;
    }
    if (vseq->substringlength - vstart < minsubstringlength)
    {
      minsubstringlength = vseq->substringlength - vstart;
    }
    if (vseq->dir_is_complement)
    {
      for (matchlength = 0; matchlength < minsubstringlength;
           uptr += ustep, vptr += vstep, matchlength++)
      {
        const GtUchar cu = gt_twobitencoding_char_at_pos(
                              useq->twobitencoding,
                              uptr);
        if (cu == WILDCARD ||
            cu != GT_COMPLEMENTBASE(*vptr))
          break;
      }
    } else
    {
      for (matchlength = 0; matchlength < minsubstringlength;
           uptr += ustep, vptr += vstep, matchlength++)
      {
        const GtUchar cu = gt_twobitencoding_char_at_pos(
                              useq->twobitencoding,
                              uptr);
        if (cu == WILDCARD ||
            cu != *vptr)
          break;
      }
    }
    return matchlength;
  }
  return 0;
}

static GtUword ft_longest_common_encseq_reader_twobit_wildcard(
                                      GtFtSequenceObject *useq,
                                      GtUword ustart,
                                      GtFtSequenceObject *vseq,
                                      const GtUword vstart)
{
  if (ustart < useq->substringlength && vstart < vseq->substringlength)
  {
    GtUword uptr;
    GtUword vptr;int vstep;
    GtUword minsubstringlength = useq->substringlength - ustart, matchlength;

    uptr = ustart;
    if (vseq->read_seq_left2right)
    {
      vptr = vseq->offset + vstart; vstep = 1;
    } else
    {
      vptr = vseq->offset - vstart; vstep = -1;
    }
    if (vseq->substringlength - vstart < minsubstringlength)
    {
      minsubstringlength = vseq->substringlength - vstart;
    }
    if (vseq->dir_is_complement)
    {
      for (matchlength = 0; matchlength < minsubstringlength;
           uptr++, vptr += vstep, matchlength++)
      {
        const GtUchar cu = gt_sequenceobject_esr_get(useq,uptr);
        if (cu == WILDCARD ||
            cu != GT_COMPLEMENTBASE(gt_twobitencoding_char_at_pos(
                              vseq->twobitencoding,
                              vptr)))
          break;
      }
    } else
    {
      for (matchlength = 0; matchlength < minsubstringlength;
           uptr++, vptr += vstep, matchlength++)
      {
        const GtUchar cu = gt_sequenceobject_esr_get(useq,uptr);
        if (cu == WILDCARD ||
            cu != gt_twobitencoding_char_at_pos(
                              vseq->twobitencoding,
                              vptr))
          break;
      }
    }
    return matchlength;
  }
  return 0;
}

static GtUword ft_longest_common_encseq_reader_encseq_reader_wildcard(
                                      GtFtSequenceObject *useq,
                                      GtUword ustart,
                                      GtFtSequenceObject *vseq,
                                      const GtUword vstart)
{
  if (ustart < useq->substringlength && vstart < vseq->substringlength)
  {
    GtUword uptr;
    GtUword vptr;
    GtUword minsubstringlength = useq->substringlength - ustart, matchlength;

    uptr = ustart;
    vptr = vstart;
    if (vseq->substringlength - vstart < minsubstringlength)
    {
      minsubstringlength = vseq->substringlength - vstart;
    }
    if (vseq->dir_is_complement)
    {
      for (matchlength = 0; matchlength < minsubstringlength;
           uptr++, vptr++, matchlength++)
      {
        const GtUchar cu = gt_sequenceobject_esr_get(useq,uptr);
        if (cu == WILDCARD ||
            cu != GT_COMPLEMENTBASE(gt_sequenceobject_esr_get(vseq,vptr)))
          break;
      }
    } else
    {
      for (matchlength = 0; matchlength < minsubstringlength;
           uptr++, vptr++, matchlength++)
      {
        const GtUchar cu = gt_sequenceobject_esr_get(useq,uptr);
        if (cu == WILDCARD ||
            cu != gt_sequenceobject_esr_get(vseq,vptr))
          break;
      }
    }
    return matchlength;
  }
  return 0;
}

static GtUword ft_longest_common_encseq_reader_encseq_wildcard(
                                      GtFtSequenceObject *useq,
                                      GtUword ustart,
                                      GtFtSequenceObject *vseq,
                                      const GtUword vstart)
{
  if (ustart < useq->substringlength && vstart < vseq->substringlength)
  {
    GtUword uptr;
    GtUword vptr;int vstep;
    GtUword minsubstringlength = useq->substringlength - ustart, matchlength;

    uptr = ustart;
    if (vseq->read_seq_left2right)
    {
      vptr = vseq->offset + vstart; vstep = 1;
    } else
    {
      vptr = vseq->offset - vstart; vstep = -1;
    }
    if (vseq->substringlength - vstart < minsubstringlength)
    {
      minsubstringlength = vseq->substringlength - vstart;
    }
    if (vseq->dir_is_complement)
    {
      for (matchlength = 0; matchlength < minsubstringlength;
           uptr++, vptr += vstep, matchlength++)
      {
        const GtUchar cu = gt_sequenceobject_esr_get(useq,uptr);
        if (cu == WILDCARD ||
            cu != GT_COMPLEMENTBASE(gt_encseq_get_encoded_char(vseq->encseq,
                    vptr,
                    GT_READMODE_FORWARD)))
          break;
      }
    } else
    {
      for (matchlength = 0; matchlength < minsubstringlength;
           uptr++, vptr += vstep, matchlength++)
      {
        const GtUchar cu = gt_sequenceobject_esr_get(useq,uptr);
        if (cu == WILDCARD ||
            cu != gt_encseq_get_encoded_char(vseq->encseq,
                    vptr,
                    GT_READMODE_FORWARD))
          break;
      }
    }
    return matchlength;
  }
  return 0;
}

static GtUword ft_longest_common_encseq_reader_bytes_wildcard(
                                      GtFtSequenceObject *useq,
                                      GtUword ustart,
                                      GtFtSequenceObject *vseq,
                                      const GtUword vstart)
{
  if (ustart < useq->substringlength && vstart < vseq->substringlength)
  {
    GtUword uptr;
    const GtUchar *vptr;int vstep;
    GtUword minsubstringlength = useq->substringlength - ustart, matchlength;

    uptr = ustart;
    if (vseq->read_seq_left2right)
    {
      vptr = vseq->bytesequenceptr + vseq->offset + vstart; vstep = 1;
    } else
    {
      vptr = vseq->bytesequenceptr + vseq->offset - vstart; vstep = -1;
    }
    if (vseq->substringlength - vstart < minsubstringlength)
    {
      minsubstringlength = vseq->substringlength - vstart;
    }
    if (vseq->dir_is_complement)
    {
      for (matchlength = 0; matchlength < minsubstringlength;
           uptr++, vptr += vstep, matchlength++)
      {
        const GtUchar cu = gt_sequenceobject_esr_get(useq,uptr);
        if (cu == WILDCARD ||
            cu != GT_COMPLEMENTBASE(*vptr))
          break;
      }
    } else
    {
      for (matchlength = 0; matchlength < minsubstringlength;
           uptr++, vptr += vstep, matchlength++)
      {
        const GtUchar cu = gt_sequenceobject_esr_get(useq,uptr);
        if (cu == WILDCARD ||
            cu != *vptr)
          break;
      }
    }
    return matchlength;
  }
  return 0;
}

static GtUword ft_longest_common_encseq_twobit_wildcard(
                                      GtFtSequenceObject *useq,
                                      GtUword ustart,
                                      GtFtSequenceObject *vseq,
                                      const GtUword vstart)
{
  if (ustart < useq->substringlength && vstart < vseq->substringlength)
  {
    GtUword uptr;int ustep;
    GtUword vptr;int vstep;
    GtUword minsubstringlength = useq->substringlength - ustart, matchlength;

    if (useq->read_seq_left2right)
    {
      uptr = useq->offset + ustart; ustep = 1;
    } else
    {
      uptr = useq->offset - ustart; ustep = -1;
    }
    if (vseq->read_seq_left2right)
    {
      vptr = vseq->offset + vstart; vstep = 1;
    } else
    {
      vptr = vseq->offset - vstart; vstep = -1;
    }
    if (vseq->substringlength - vstart < minsubstringlength)
    {
      minsubstringlength = vseq->substringlength - vstart;
    }
    if (vseq->dir_is_complement)
    {
      for (matchlength = 0; matchlength < minsubstringlength;
           uptr += ustep, vptr += vstep, matchlength++)
      {
        const GtUchar cu = gt_encseq_get_encoded_char(useq->encseq,
                    uptr,
                    GT_READMODE_FORWARD);
        if (cu == WILDCARD ||
            cu != GT_COMPLEMENTBASE(gt_twobitencoding_char_at_pos(
                              vseq->twobitencoding,
                              vptr)))
          break;
      }
    } else
    {
      for (matchlength = 0; matchlength < minsubstringlength;
           uptr += ustep, vptr += vstep, matchlength++)
      {
        const GtUchar cu = gt_encseq_get_encoded_char(useq->encseq,
                    uptr,
                    GT_READMODE_FORWARD);
        if (cu == WILDCARD ||
            cu != gt_twobitencoding_char_at_pos(
                              vseq->twobitencoding,
                              vptr))
          break;
      }
    }
    return matchlength;
  }
  return 0;
}

static GtUword ft_longest_common_encseq_encseq_reader_wildcard(
                                      GtFtSequenceObject *useq,
                                      GtUword ustart,
                                      GtFtSequenceObject *vseq,
                                      const GtUword vstart)
{
  if (ustart < useq->substringlength && vstart < vseq->substringlength)
  {
    GtUword uptr;int ustep;
    GtUword vptr;
    GtUword minsubstringlength = useq->substringlength - ustart, matchlength;

    if (useq->read_seq_left2right)
    {
      uptr = useq->offset + ustart; ustep = 1;
    } else
    {
      uptr = useq->offset - ustart; ustep = -1;
    }
    vptr = vstart;
    if (vseq->substringlength - vstart < minsubstringlength)
    {
      minsubstringlength = vseq->substringlength - vstart;
    }
    if (vseq->dir_is_complement)
    {
      for (matchlength = 0; matchlength < minsubstringlength;
           uptr += ustep, vptr++, matchlength++)
      {
        const GtUchar cu = gt_encseq_get_encoded_char(useq->encseq,
                    uptr,
                    GT_READMODE_FORWARD);
        if (cu == WILDCARD ||
            cu != GT_COMPLEMENTBASE(gt_sequenceobject_esr_get(vseq,vptr)))
          break;
      }
    } else
    {
      for (matchlength = 0; matchlength < minsubstringlength;
           uptr += ustep, vptr++, matchlength++)
      {
        const GtUchar cu = gt_encseq_get_encoded_char(useq->encseq,
                    uptr,
                    GT_READMODE_FORWARD);
        if (cu == WILDCARD ||
            cu != gt_sequenceobject_esr_get(vseq,vptr))
          break;
      }
    }
    return matchlength;
  }
  return 0;
}

static GtUword ft_longest_common_encseq_encseq_wildcard(
                                      GtFtSequenceObject *useq,
                                      GtUword ustart,
                                      GtFtSequenceObject *vseq,
                                      const GtUword vstart)
{
  if (ustart < useq->substringlength && vstart < vseq->substringlength)
  {
    GtUword uptr;int ustep;
    GtUword vptr;int vstep;
    GtUword minsubstringlength = useq->substringlength - ustart, matchlength;

    if (useq->read_seq_left2right)
    {
      uptr = useq->offset + ustart; ustep = 1;
    } else
    {
      uptr = useq->offset - ustart; ustep = -1;
    }
    if (vseq->read_seq_left2right)
    {
      vptr = vseq->offset + vstart; vstep = 1;
    } else
    {
      vptr = vseq->offset - vstart; vstep = -1;
    }
    if (vseq->substringlength - vstart < minsubstringlength)
    {
      minsubstringlength = vseq->substringlength - vstart;
    }
    if (vseq->dir_is_complement)
    {
      for (matchlength = 0; matchlength < minsubstringlength;
           uptr += ustep, vptr += vstep, matchlength++)
      {
        const GtUchar cu = gt_encseq_get_encoded_char(useq->encseq,
                    uptr,
                    GT_READMODE_FORWARD);
        if (cu == WILDCARD ||
            cu != GT_COMPLEMENTBASE(gt_encseq_get_encoded_char(vseq->encseq,
                    vptr,
                    GT_READMODE_FORWARD)))
          break;
      }
    } else
    {
      for (matchlength = 0; matchlength < minsubstringlength;
           uptr += ustep, vptr += vstep, matchlength++)
      {
        const GtUchar cu = gt_encseq_get_encoded_char(useq->encseq,
                    uptr,
                    GT_READMODE_FORWARD);
        if (cu == WILDCARD ||
            cu != gt_encseq_get_encoded_char(vseq->encseq,
                    vptr,
                    GT_READMODE_FORWARD))
          break;
      }
    }
    return matchlength;
  }
  return 0;
}

static GtUword ft_longest_common_encseq_bytes_wildcard(
                                      GtFtSequenceObject *useq,
                                      GtUword ustart,
                                      GtFtSequenceObject *vseq,
                                      const GtUword vstart)
{
  if (ustart < useq->substringlength && vstart < vseq->substringlength)
  {
    GtUword uptr;int ustep;
    const GtUchar *vptr;int vstep;
    GtUword minsubstringlength = useq->substringlength - ustart, matchlength;

    if (useq->read_seq_left2right)
    {
      uptr = useq->offset + ustart; ustep = 1;
    } else
    {
      uptr = useq->offset - ustart; ustep = -1;
    }
    if (vseq->read_seq_left2right)
    {
      vptr = vseq->bytesequenceptr + vseq->offset + vstart; vstep = 1;
    } else
    {
      vptr = vseq->bytesequenceptr + vseq->offset - vstart; vstep = -1;
    }
    if (vseq->substringlength - vstart < minsubstringlength)
    {
      minsubstringlength = vseq->substringlength - vstart;
    }
    if (vseq->dir_is_complement)
    {
      for (matchlength = 0; matchlength < minsubstringlength;
           uptr += ustep, vptr += vstep, matchlength++)
      {
        const GtUchar cu = gt_encseq_get_encoded_char(useq->encseq,
                    uptr,
                    GT_READMODE_FORWARD);
        if (cu == WILDCARD ||
            cu != GT_COMPLEMENTBASE(*vptr))
          break;
      }
    } else
    {
      for (matchlength = 0; matchlength < minsubstringlength;
           uptr += ustep, vptr += vstep, matchlength++)
      {
        const GtUchar cu = gt_encseq_get_encoded_char(useq->encseq,
                    uptr,
                    GT_READMODE_FORWARD);
        if (cu == WILDCARD ||
            cu != *vptr)
          break;
      }
    }
    return matchlength;
  }
  return 0;
}

static GtUword ft_longest_common_bytes_twobit_wildcard(
                                      GtFtSequenceObject *useq,
                                      GtUword ustart,
                                      GtFtSequenceObject *vseq,
                                      const GtUword vstart)
{
  if (ustart < useq->substringlength && vstart < vseq->substringlength)
  {
    const GtUchar *uptr;int ustep;
    GtUword vptr;int vstep;
    GtUword minsubstringlength = useq->substringlength - ustart, matchlength;

    if (useq->read_seq_left2right)
    {
      uptr = useq->bytesequenceptr + useq->offset + ustart; ustep = 1;
    } else
    {
      uptr = useq->bytesequenceptr + useq->offset - ustart; ustep = -1;
    }
    if (vseq->read_seq_left2right)
    {
      vptr = vseq->offset + vstart; vstep = 1;
    } else
    {
      vptr = vseq->offset - vstart; vstep = -1;
    }
    if (vseq->substringlength - vstart < minsubstringlength)
    {
      minsubstringlength = vseq->substringlength - vstart;
    }
    if (vseq->dir_is_complement)
    {
      for (matchlength = 0; matchlength < minsubstringlength;
           uptr += ustep, vptr += vstep, matchlength++)
      {
        const GtUchar cu = *uptr;
        if (cu == WILDCARD ||
            cu != GT_COMPLEMENTBASE(gt_twobitencoding_char_at_pos(
                              vseq->twobitencoding,
                              vptr)))
          break;
      }
    } else
    {
      for (matchlength = 0; matchlength < minsubstringlength;
           uptr += ustep, vptr += vstep, matchlength++)
      {
        const GtUchar cu = *uptr;
        if (cu == WILDCARD ||
            cu != gt_twobitencoding_char_at_pos(
                              vseq->twobitencoding,
                              vptr))
          break;
      }
    }
    return matchlength;
  }
  return 0;
}

static GtUword ft_longest_common_bytes_encseq_reader_wildcard(
                                      GtFtSequenceObject *useq,
                                      GtUword ustart,
                                      GtFtSequenceObject *vseq,
                                      const GtUword vstart)
{
  if (ustart < useq->substringlength && vstart < vseq->substringlength)
  {
    const GtUchar *uptr;int ustep;
    GtUword vptr;
    GtUword minsubstringlength = useq->substringlength - ustart, matchlength;

    if (useq->read_seq_left2right)
    {
      uptr = useq->bytesequenceptr + useq->offset + ustart; ustep = 1;
    } else
    {
      uptr = useq->bytesequenceptr + useq->offset - ustart; ustep = -1;
    }
    vptr = vstart;
    if (vseq->substringlength - vstart < minsubstringlength)
    {
      minsubstringlength = vseq->substringlength - vstart;
    }
    if (vseq->dir_is_complement)
    {
      for (matchlength = 0; matchlength < minsubstringlength;
           uptr += ustep, vptr++, matchlength++)
      {
        const GtUchar cu = *uptr;
        if (cu == WILDCARD ||
            cu != GT_COMPLEMENTBASE(gt_sequenceobject_esr_get(vseq,vptr)))
          break;
      }
    } else
    {
      for (matchlength = 0; matchlength < minsubstringlength;
           uptr += ustep, vptr++, matchlength++)
      {
        const GtUchar cu = *uptr;
        if (cu == WILDCARD ||
            cu != gt_sequenceobject_esr_get(vseq,vptr))
          break;
      }
    }
    return matchlength;
  }
  return 0;
}

static GtUword ft_longest_common_bytes_encseq_wildcard(
                                      GtFtSequenceObject *useq,
                                      GtUword ustart,
                                      GtFtSequenceObject *vseq,
                                      const GtUword vstart)
{
  if (ustart < useq->substringlength && vstart < vseq->substringlength)
  {
    const GtUchar *uptr;int ustep;
    GtUword vptr;int vstep;
    GtUword minsubstringlength = useq->substringlength - ustart, matchlength;

    if (useq->read_seq_left2right)
    {
      uptr = useq->bytesequenceptr + useq->offset + ustart; ustep = 1;
    } else
    {
      uptr = useq->bytesequenceptr + useq->offset - ustart; ustep = -1;
    }
    if (vseq->read_seq_left2right)
    {
      vptr = vseq->offset + vstart; vstep = 1;
    } else
    {
      vptr = vseq->offset - vstart; vstep = -1;
    }
    if (vseq->substringlength - vstart < minsubstringlength)
    {
      minsubstringlength = vseq->substringlength - vstart;
    }
    if (vseq->dir_is_complement)
    {
      for (matchlength = 0; matchlength < minsubstringlength;
           uptr += ustep, vptr += vstep, matchlength++)
      {
        const GtUchar cu = *uptr;
        if (cu == WILDCARD ||
            cu != GT_COMPLEMENTBASE(gt_encseq_get_encoded_char(vseq->encseq,
                    vptr,
                    GT_READMODE_FORWARD)))
          break;
      }
    } else
    {
      for (matchlength = 0; matchlength < minsubstringlength;
           uptr += ustep, vptr += vstep, matchlength++)
      {
        const GtUchar cu = *uptr;
        if (cu == WILDCARD ||
            cu != gt_encseq_get_encoded_char(vseq->encseq,
                    vptr,
                    GT_READMODE_FORWARD))
          break;
      }
    }
    return matchlength;
  }
  return 0;
}

static GtUword ft_longest_common_bytes_bytes_wildcard(
                                      GtFtSequenceObject *useq,
                                      GtUword ustart,
                                      GtFtSequenceObject *vseq,
                                      const GtUword vstart)
{
  if (ustart < useq->substringlength && vstart < vseq->substringlength)
  {
    const GtUchar *uptr;int ustep;
    const GtUchar *vptr;int vstep;
    GtUword minsubstringlength = useq->substringlength - ustart, matchlength;

    if (useq->read_seq_left2right)
    {
      uptr = useq->bytesequenceptr + useq->offset + ustart; ustep = 1;
    } else
    {
      uptr = useq->bytesequenceptr + useq->offset - ustart; ustep = -1;
    }
    if (vseq->read_seq_left2right)
    {
      vptr = vseq->bytesequenceptr + vseq->offset + vstart; vstep = 1;
    } else
    {
      vptr = vseq->bytesequenceptr + vseq->offset - vstart; vstep = -1;
    }
    if (vseq->substringlength - vstart < minsubstringlength)
    {
      minsubstringlength = vseq->substringlength - vstart;
    }
    if (vseq->dir_is_complement)
    {
      for (matchlength = 0; matchlength < minsubstringlength;
           uptr += ustep, vptr += vstep, matchlength++)
      {
        const GtUchar cu = *uptr;
        if (cu == WILDCARD ||
            cu != GT_COMPLEMENTBASE(*vptr))
          break;
      }
    } else
    {
      for (matchlength = 0; matchlength < minsubstringlength;
           uptr += ustep, vptr += vstep, matchlength++)
      {
        const GtUchar cu = *uptr;
        if (cu == WILDCARD ||
            cu != *vptr)
          break;
      }
    }
    return matchlength;
  }
  return 0;
}

GtLongestCommonFunc ft_longest_common_func_tab[] =
{
  /* 0 */ ft_longest_common_twobit_twobit,
  /* 1 */ ft_longest_common_twobit_encseq_reader,
  /* 2 */ ft_longest_common_twobit_encseq,
  /* 3 */ ft_longest_common_twobit_bytes,
  /* 4 */ ft_longest_common_encseq_reader_twobit,
  /* 5 */ ft_longest_common_encseq_reader_encseq_reader,
  /* 6 */ ft_longest_common_encseq_reader_encseq,
  /* 7 */ ft_longest_common_encseq_reader_bytes,
  /* 8 */ ft_longest_common_encseq_twobit,
  /* 9 */ ft_longest_common_encseq_encseq_reader,
  /* 10 */ ft_longest_common_encseq_encseq,
  /* 11 */ ft_longest_common_encseq_bytes,
  /* 12 */ ft_longest_common_bytes_twobit,
  /* 13 */ ft_longest_common_bytes_encseq_reader,
  /* 14 */ ft_longest_common_bytes_encseq,
  /* 15 */ ft_longest_common_bytes_bytes,
  /* 16 */ ft_longest_common_twobit_twobit_wildcard,
  /* 17 */ ft_longest_common_twobit_encseq_reader_wildcard,
  /* 18 */ ft_longest_common_twobit_encseq_wildcard,
  /* 19 */ ft_longest_common_twobit_bytes_wildcard,
  /* 20 */ ft_longest_common_encseq_reader_twobit_wildcard,
  /* 21 */ ft_longest_common_encseq_reader_encseq_reader_wildcard,
  /* 22 */ ft_longest_common_encseq_reader_encseq_wildcard,
  /* 23 */ ft_longest_common_encseq_reader_bytes_wildcard,
  /* 24 */ ft_longest_common_encseq_twobit_wildcard,
  /* 25 */ ft_longest_common_encseq_encseq_reader_wildcard,
  /* 26 */ ft_longest_common_encseq_encseq_wildcard,
  /* 27 */ ft_longest_common_encseq_bytes_wildcard,
  /* 28 */ ft_longest_common_bytes_twobit_wildcard,
  /* 29 */ ft_longest_common_bytes_encseq_reader_wildcard,
  /* 30 */ ft_longest_common_bytes_encseq_wildcard,
  /* 31 */ ft_longest_common_bytes_bytes_wildcard
};
const int ft_longest_common_num_modes = 4;
const int ft_longest_common_func_first_wildcard = 16;
