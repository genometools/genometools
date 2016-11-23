static GtUword ft_longest_common_twobit_twobit(
                                      GtFtSequenceObject *useq,
                                      GtUword ustart,
                                      GtFtSequenceObject *vseq,
                                      const GtUword vstart)
{
  GtUword upos, vpos;

  for (upos = ustart, vpos = vstart;
       upos < useq->substringlength && vpos < vseq->substringlength;
       upos++, vpos++)
  {
    const GtUchar cu = gt_twobitencoding_char_at_pos(
                    useq->twobitencoding,
                    useq->read_seq_left2right
                                             ? useq->offset + upos
                                             : useq->offset - upos
                    );

    if (cu != (vseq->dir_is_complement
          ? GT_COMPLEMENTBASE(gt_twobitencoding_char_at_pos(
                    vseq->twobitencoding,
                    vseq->read_seq_left2right
                                             ? vseq->offset + vpos
                                             : vseq->offset - vpos
                    ))
          : gt_twobitencoding_char_at_pos(
                    vseq->twobitencoding,
                    vseq->read_seq_left2right
                                             ? vseq->offset + vpos
                                             : vseq->offset - vpos
                    )))
    {
      break;
    }
  }
  return upos - ustart;
}

static GtUword ft_longest_common_twobit_encseq(
                                      GtFtSequenceObject *useq,
                                      GtUword ustart,
                                      GtFtSequenceObject *vseq,
                                      const GtUword vstart)
{
  GtUword upos, vpos;

  for (upos = ustart, vpos = vstart;
       upos < useq->substringlength && vpos < vseq->substringlength;
       upos++, vpos++)
  {
    const GtUchar cu = gt_twobitencoding_char_at_pos(
                    useq->twobitencoding,
                    useq->read_seq_left2right
                                             ? useq->offset + upos
                                             : useq->offset - upos
                    );

    if (cu != (vseq->dir_is_complement
          ? GT_COMPLEMENTBASE(gt_encseq_get_encoded_char(vseq->encseq,
                    vseq->read_seq_left2right
                                             ? vseq->offset + vpos
                                             : vseq->offset - vpos,
                    GT_READMODE_FORWARD))
          : gt_encseq_get_encoded_char(vseq->encseq,
                    vseq->read_seq_left2right
                                             ? vseq->offset + vpos
                                             : vseq->offset - vpos,
                    GT_READMODE_FORWARD)))
    {
      break;
    }
  }
  return upos - ustart;
}

static GtUword ft_longest_common_twobit_bytes(
                                      GtFtSequenceObject *useq,
                                      GtUword ustart,
                                      GtFtSequenceObject *vseq,
                                      const GtUword vstart)
{
  GtUword upos, vpos;

  for (upos = ustart, vpos = vstart;
       upos < useq->substringlength && vpos < vseq->substringlength;
       upos++, vpos++)
  {
    const GtUchar cu = gt_twobitencoding_char_at_pos(
                    useq->twobitencoding,
                    useq->read_seq_left2right
                                             ? useq->offset + upos
                                             : useq->offset - upos
                    );

    if (cu != (vseq->dir_is_complement
          ? GT_COMPLEMENTBASE(vseq->bytesequenceptr[
                    vseq->read_seq_left2right
                                             ? vseq->offset + vpos
                                             : vseq->offset - vpos
                    ])
          : vseq->bytesequenceptr[
                    vseq->read_seq_left2right
                                             ? vseq->offset + vpos
                                             : vseq->offset - vpos
                    ]))
    {
      break;
    }
  }
  return upos - ustart;
}

static GtUword ft_longest_common_encseq_twobit(
                                      GtFtSequenceObject *useq,
                                      GtUword ustart,
                                      GtFtSequenceObject *vseq,
                                      const GtUword vstart)
{
  GtUword upos, vpos;

  for (upos = ustart, vpos = vstart;
       upos < useq->substringlength && vpos < vseq->substringlength;
       upos++, vpos++)
  {
    const GtUchar cu = gt_encseq_get_encoded_char(useq->encseq,
                    useq->read_seq_left2right
                                             ? useq->offset + upos
                                             : useq->offset - upos,
                    GT_READMODE_FORWARD);

    if (cu != (vseq->dir_is_complement
          ? GT_COMPLEMENTBASE(gt_twobitencoding_char_at_pos(
                    vseq->twobitencoding,
                    vseq->read_seq_left2right
                                             ? vseq->offset + vpos
                                             : vseq->offset - vpos
                    ))
          : gt_twobitencoding_char_at_pos(
                    vseq->twobitencoding,
                    vseq->read_seq_left2right
                                             ? vseq->offset + vpos
                                             : vseq->offset - vpos
                    )))
    {
      break;
    }
  }
  return upos - ustart;
}

static GtUword ft_longest_common_encseq_encseq(
                                      GtFtSequenceObject *useq,
                                      GtUword ustart,
                                      GtFtSequenceObject *vseq,
                                      const GtUword vstart)
{
  GtUword upos, vpos;

  for (upos = ustart, vpos = vstart;
       upos < useq->substringlength && vpos < vseq->substringlength;
       upos++, vpos++)
  {
    const GtUchar cu = gt_encseq_get_encoded_char(useq->encseq,
                    useq->read_seq_left2right
                                             ? useq->offset + upos
                                             : useq->offset - upos,
                    GT_READMODE_FORWARD);

    if (cu != (vseq->dir_is_complement
          ? GT_COMPLEMENTBASE(gt_encseq_get_encoded_char(vseq->encseq,
                    vseq->read_seq_left2right
                                             ? vseq->offset + vpos
                                             : vseq->offset - vpos,
                    GT_READMODE_FORWARD))
          : gt_encseq_get_encoded_char(vseq->encseq,
                    vseq->read_seq_left2right
                                             ? vseq->offset + vpos
                                             : vseq->offset - vpos,
                    GT_READMODE_FORWARD)))
    {
      break;
    }
  }
  return upos - ustart;
}

static GtUword ft_longest_common_encseq_bytes(
                                      GtFtSequenceObject *useq,
                                      GtUword ustart,
                                      GtFtSequenceObject *vseq,
                                      const GtUword vstart)
{
  GtUword upos, vpos;

  for (upos = ustart, vpos = vstart;
       upos < useq->substringlength && vpos < vseq->substringlength;
       upos++, vpos++)
  {
    const GtUchar cu = gt_encseq_get_encoded_char(useq->encseq,
                    useq->read_seq_left2right
                                             ? useq->offset + upos
                                             : useq->offset - upos,
                    GT_READMODE_FORWARD);

    if (cu != (vseq->dir_is_complement
          ? GT_COMPLEMENTBASE(vseq->bytesequenceptr[
                    vseq->read_seq_left2right
                                             ? vseq->offset + vpos
                                             : vseq->offset - vpos
                    ])
          : vseq->bytesequenceptr[
                    vseq->read_seq_left2right
                                             ? vseq->offset + vpos
                                             : vseq->offset - vpos
                    ]))
    {
      break;
    }
  }
  return upos - ustart;
}

static GtUword ft_longest_common_bytes_twobit(
                                      GtFtSequenceObject *useq,
                                      GtUword ustart,
                                      GtFtSequenceObject *vseq,
                                      const GtUword vstart)
{
  GtUword upos, vpos;

  for (upos = ustart, vpos = vstart;
       upos < useq->substringlength && vpos < vseq->substringlength;
       upos++, vpos++)
  {
    const GtUchar cu = useq->bytesequenceptr[
                    useq->read_seq_left2right
                                             ? useq->offset + upos
                                             : useq->offset - upos
                    ];

    if (cu != (vseq->dir_is_complement
          ? GT_COMPLEMENTBASE(gt_twobitencoding_char_at_pos(
                    vseq->twobitencoding,
                    vseq->read_seq_left2right
                                             ? vseq->offset + vpos
                                             : vseq->offset - vpos
                    ))
          : gt_twobitencoding_char_at_pos(
                    vseq->twobitencoding,
                    vseq->read_seq_left2right
                                             ? vseq->offset + vpos
                                             : vseq->offset - vpos
                    )))
    {
      break;
    }
  }
  return upos - ustart;
}

static GtUword ft_longest_common_bytes_encseq(
                                      GtFtSequenceObject *useq,
                                      GtUword ustart,
                                      GtFtSequenceObject *vseq,
                                      const GtUword vstart)
{
  GtUword upos, vpos;

  for (upos = ustart, vpos = vstart;
       upos < useq->substringlength && vpos < vseq->substringlength;
       upos++, vpos++)
  {
    const GtUchar cu = useq->bytesequenceptr[
                    useq->read_seq_left2right
                                             ? useq->offset + upos
                                             : useq->offset - upos
                    ];

    if (cu != (vseq->dir_is_complement
          ? GT_COMPLEMENTBASE(gt_encseq_get_encoded_char(vseq->encseq,
                    vseq->read_seq_left2right
                                             ? vseq->offset + vpos
                                             : vseq->offset - vpos,
                    GT_READMODE_FORWARD))
          : gt_encseq_get_encoded_char(vseq->encseq,
                    vseq->read_seq_left2right
                                             ? vseq->offset + vpos
                                             : vseq->offset - vpos,
                    GT_READMODE_FORWARD)))
    {
      break;
    }
  }
  return upos - ustart;
}

static GtUword ft_longest_common_bytes_bytes(
                                      GtFtSequenceObject *useq,
                                      GtUword ustart,
                                      GtFtSequenceObject *vseq,
                                      const GtUword vstart)
{
  GtUword upos, vpos;

  for (upos = ustart, vpos = vstart;
       upos < useq->substringlength && vpos < vseq->substringlength;
       upos++, vpos++)
  {
    const GtUchar cu = useq->bytesequenceptr[
                    useq->read_seq_left2right
                                             ? useq->offset + upos
                                             : useq->offset - upos
                    ];

    if (cu != (vseq->dir_is_complement
          ? GT_COMPLEMENTBASE(vseq->bytesequenceptr[
                    vseq->read_seq_left2right
                                             ? vseq->offset + vpos
                                             : vseq->offset - vpos
                    ])
          : vseq->bytesequenceptr[
                    vseq->read_seq_left2right
                                             ? vseq->offset + vpos
                                             : vseq->offset - vpos
                    ]))
    {
      break;
    }
  }
  return upos - ustart;
}

static GtUword ft_longest_common_twobit_twobit_wildcard(
                                      GtFtSequenceObject *useq,
                                      GtUword ustart,
                                      GtFtSequenceObject *vseq,
                                      const GtUword vstart)
{
  GtUword upos, vpos;

  for (upos = ustart, vpos = vstart;
       upos < useq->substringlength && vpos < vseq->substringlength;
       upos++, vpos++)
  {
    const GtUchar cu = gt_twobitencoding_char_at_pos(
                    useq->twobitencoding,
                    useq->read_seq_left2right
                                             ? useq->offset + upos
                                             : useq->offset - upos
                    );

    if ((cu == WILDCARD) || cu != (vseq->dir_is_complement
          ? GT_COMPLEMENTBASE(gt_twobitencoding_char_at_pos(
                    vseq->twobitencoding,
                    vseq->read_seq_left2right
                                             ? vseq->offset + vpos
                                             : vseq->offset - vpos
                    ))
          : gt_twobitencoding_char_at_pos(
                    vseq->twobitencoding,
                    vseq->read_seq_left2right
                                             ? vseq->offset + vpos
                                             : vseq->offset - vpos
                    )))
    {
      break;
    }
  }
  return upos - ustart;
}

static GtUword ft_longest_common_twobit_encseq_wildcard(
                                      GtFtSequenceObject *useq,
                                      GtUword ustart,
                                      GtFtSequenceObject *vseq,
                                      const GtUword vstart)
{
  GtUword upos, vpos;

  for (upos = ustart, vpos = vstart;
       upos < useq->substringlength && vpos < vseq->substringlength;
       upos++, vpos++)
  {
    const GtUchar cu = gt_twobitencoding_char_at_pos(
                    useq->twobitencoding,
                    useq->read_seq_left2right
                                             ? useq->offset + upos
                                             : useq->offset - upos
                    );

    if ((cu == WILDCARD) || cu != (vseq->dir_is_complement
          ? GT_COMPLEMENTBASE(gt_encseq_get_encoded_char(vseq->encseq,
                    vseq->read_seq_left2right
                                             ? vseq->offset + vpos
                                             : vseq->offset - vpos,
                    GT_READMODE_FORWARD))
          : gt_encseq_get_encoded_char(vseq->encseq,
                    vseq->read_seq_left2right
                                             ? vseq->offset + vpos
                                             : vseq->offset - vpos,
                    GT_READMODE_FORWARD)))
    {
      break;
    }
  }
  return upos - ustart;
}

static GtUword ft_longest_common_twobit_bytes_wildcard(
                                      GtFtSequenceObject *useq,
                                      GtUword ustart,
                                      GtFtSequenceObject *vseq,
                                      const GtUword vstart)
{
  GtUword upos, vpos;

  for (upos = ustart, vpos = vstart;
       upos < useq->substringlength && vpos < vseq->substringlength;
       upos++, vpos++)
  {
    const GtUchar cu = gt_twobitencoding_char_at_pos(
                    useq->twobitencoding,
                    useq->read_seq_left2right
                                             ? useq->offset + upos
                                             : useq->offset - upos
                    );

    if ((cu == WILDCARD) || cu != (vseq->dir_is_complement
          ? GT_COMPLEMENTBASE(vseq->bytesequenceptr[
                    vseq->read_seq_left2right
                                             ? vseq->offset + vpos
                                             : vseq->offset - vpos
                    ])
          : vseq->bytesequenceptr[
                    vseq->read_seq_left2right
                                             ? vseq->offset + vpos
                                             : vseq->offset - vpos
                    ]))
    {
      break;
    }
  }
  return upos - ustart;
}

static GtUword ft_longest_common_encseq_twobit_wildcard(
                                      GtFtSequenceObject *useq,
                                      GtUword ustart,
                                      GtFtSequenceObject *vseq,
                                      const GtUword vstart)
{
  GtUword upos, vpos;

  for (upos = ustart, vpos = vstart;
       upos < useq->substringlength && vpos < vseq->substringlength;
       upos++, vpos++)
  {
    const GtUchar cu = gt_encseq_get_encoded_char(useq->encseq,
                    useq->read_seq_left2right
                                             ? useq->offset + upos
                                             : useq->offset - upos,
                    GT_READMODE_FORWARD);

    if ((cu == WILDCARD) || cu != (vseq->dir_is_complement
          ? GT_COMPLEMENTBASE(gt_twobitencoding_char_at_pos(
                    vseq->twobitencoding,
                    vseq->read_seq_left2right
                                             ? vseq->offset + vpos
                                             : vseq->offset - vpos
                    ))
          : gt_twobitencoding_char_at_pos(
                    vseq->twobitencoding,
                    vseq->read_seq_left2right
                                             ? vseq->offset + vpos
                                             : vseq->offset - vpos
                    )))
    {
      break;
    }
  }
  return upos - ustart;
}

static GtUword ft_longest_common_encseq_encseq_wildcard(
                                      GtFtSequenceObject *useq,
                                      GtUword ustart,
                                      GtFtSequenceObject *vseq,
                                      const GtUword vstart)
{
  GtUword upos, vpos;

  for (upos = ustart, vpos = vstart;
       upos < useq->substringlength && vpos < vseq->substringlength;
       upos++, vpos++)
  {
    const GtUchar cu = gt_encseq_get_encoded_char(useq->encseq,
                    useq->read_seq_left2right
                                             ? useq->offset + upos
                                             : useq->offset - upos,
                    GT_READMODE_FORWARD);

    if ((cu == WILDCARD) || cu != (vseq->dir_is_complement
          ? GT_COMPLEMENTBASE(gt_encseq_get_encoded_char(vseq->encseq,
                    vseq->read_seq_left2right
                                             ? vseq->offset + vpos
                                             : vseq->offset - vpos,
                    GT_READMODE_FORWARD))
          : gt_encseq_get_encoded_char(vseq->encseq,
                    vseq->read_seq_left2right
                                             ? vseq->offset + vpos
                                             : vseq->offset - vpos,
                    GT_READMODE_FORWARD)))
    {
      break;
    }
  }
  return upos - ustart;
}

static GtUword ft_longest_common_encseq_bytes_wildcard(
                                      GtFtSequenceObject *useq,
                                      GtUword ustart,
                                      GtFtSequenceObject *vseq,
                                      const GtUword vstart)
{
  GtUword upos, vpos;

  for (upos = ustart, vpos = vstart;
       upos < useq->substringlength && vpos < vseq->substringlength;
       upos++, vpos++)
  {
    const GtUchar cu = gt_encseq_get_encoded_char(useq->encseq,
                    useq->read_seq_left2right
                                             ? useq->offset + upos
                                             : useq->offset - upos,
                    GT_READMODE_FORWARD);

    if ((cu == WILDCARD) || cu != (vseq->dir_is_complement
          ? GT_COMPLEMENTBASE(vseq->bytesequenceptr[
                    vseq->read_seq_left2right
                                             ? vseq->offset + vpos
                                             : vseq->offset - vpos
                    ])
          : vseq->bytesequenceptr[
                    vseq->read_seq_left2right
                                             ? vseq->offset + vpos
                                             : vseq->offset - vpos
                    ]))
    {
      break;
    }
  }
  return upos - ustart;
}

static GtUword ft_longest_common_bytes_twobit_wildcard(
                                      GtFtSequenceObject *useq,
                                      GtUword ustart,
                                      GtFtSequenceObject *vseq,
                                      const GtUword vstart)
{
  GtUword upos, vpos;

  for (upos = ustart, vpos = vstart;
       upos < useq->substringlength && vpos < vseq->substringlength;
       upos++, vpos++)
  {
    const GtUchar cu = useq->bytesequenceptr[
                    useq->read_seq_left2right
                                             ? useq->offset + upos
                                             : useq->offset - upos
                    ];

    if ((cu == WILDCARD) || cu != (vseq->dir_is_complement
          ? GT_COMPLEMENTBASE(gt_twobitencoding_char_at_pos(
                    vseq->twobitencoding,
                    vseq->read_seq_left2right
                                             ? vseq->offset + vpos
                                             : vseq->offset - vpos
                    ))
          : gt_twobitencoding_char_at_pos(
                    vseq->twobitencoding,
                    vseq->read_seq_left2right
                                             ? vseq->offset + vpos
                                             : vseq->offset - vpos
                    )))
    {
      break;
    }
  }
  return upos - ustart;
}

static GtUword ft_longest_common_bytes_encseq_wildcard(
                                      GtFtSequenceObject *useq,
                                      GtUword ustart,
                                      GtFtSequenceObject *vseq,
                                      const GtUword vstart)
{
  GtUword upos, vpos;

  for (upos = ustart, vpos = vstart;
       upos < useq->substringlength && vpos < vseq->substringlength;
       upos++, vpos++)
  {
    const GtUchar cu = useq->bytesequenceptr[
                    useq->read_seq_left2right
                                             ? useq->offset + upos
                                             : useq->offset - upos
                    ];

    if ((cu == WILDCARD) || cu != (vseq->dir_is_complement
          ? GT_COMPLEMENTBASE(gt_encseq_get_encoded_char(vseq->encseq,
                    vseq->read_seq_left2right
                                             ? vseq->offset + vpos
                                             : vseq->offset - vpos,
                    GT_READMODE_FORWARD))
          : gt_encseq_get_encoded_char(vseq->encseq,
                    vseq->read_seq_left2right
                                             ? vseq->offset + vpos
                                             : vseq->offset - vpos,
                    GT_READMODE_FORWARD)))
    {
      break;
    }
  }
  return upos - ustart;
}

static GtUword ft_longest_common_bytes_bytes_wildcard(
                                      GtFtSequenceObject *useq,
                                      GtUword ustart,
                                      GtFtSequenceObject *vseq,
                                      const GtUword vstart)
{
  GtUword upos, vpos;

  for (upos = ustart, vpos = vstart;
       upos < useq->substringlength && vpos < vseq->substringlength;
       upos++, vpos++)
  {
    const GtUchar cu = useq->bytesequenceptr[
                    useq->read_seq_left2right
                                             ? useq->offset + upos
                                             : useq->offset - upos
                    ];

    if ((cu == WILDCARD) || cu != (vseq->dir_is_complement
          ? GT_COMPLEMENTBASE(vseq->bytesequenceptr[
                    vseq->read_seq_left2right
                                             ? vseq->offset + vpos
                                             : vseq->offset - vpos
                    ])
          : vseq->bytesequenceptr[
                    vseq->read_seq_left2right
                                             ? vseq->offset + vpos
                                             : vseq->offset - vpos
                    ]))
    {
      break;
    }
  }
  return upos - ustart;
}

GtLongestCommonFunc ft_longest_common_func_tab[] =
{
  /* 0 */ ft_longest_common_all,
  /* 1 */ ft_longest_common_twobit_twobit,
  /* 2 */ ft_longest_common_twobit_encseq,
  /* 3 */ ft_longest_common_twobit_bytes,
  /* 4 */ ft_longest_common_encseq_twobit,
  /* 5 */ ft_longest_common_encseq_encseq,
  /* 6 */ ft_longest_common_encseq_bytes,
  /* 7 */ ft_longest_common_bytes_twobit,
  /* 8 */ ft_longest_common_bytes_encseq,
  /* 9 */ ft_longest_common_bytes_bytes,
  /* 10 */ ft_longest_common_twobit_twobit_wildcard,
  /* 11 */ ft_longest_common_twobit_encseq_wildcard,
  /* 12 */ ft_longest_common_twobit_bytes_wildcard,
  /* 13 */ ft_longest_common_encseq_twobit_wildcard,
  /* 14 */ ft_longest_common_encseq_encseq_wildcard,
  /* 15 */ ft_longest_common_encseq_bytes_wildcard,
  /* 16 */ ft_longest_common_bytes_twobit_wildcard,
  /* 17 */ ft_longest_common_bytes_encseq_wildcard,
  /* 18 */ ft_longest_common_bytes_bytes_wildcard
};
const int ft_longest_comon_func_first_wildcard = 9;
