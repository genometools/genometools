print ([[

File format for option -offsetfile:

The file supplied to option -offsetfile defines a mapping table named
``offsets''. It maps the sequence-regions given in the GFF3_file to offsets.
It can be defined as follows:

offsets = {
  chr1  = 1000,
  chr2  = 500
}

When this example is used, all features with seqid ``chr1'' will be offset by
1000 and all features with seqid ``chr2'' by 500.


File format for option -chseqids:

The file supplied to option -chseqids defines a mapping table named
``chseqids''. It maps the sequence-regions given in the GFF3_file to other names.
It can be defined as follows:

chseqids = {
  chr1  = "seq1",
  chr2  = "seq2"
}

When this example is used, all sequence ids ``chr1'' will be changed to ``seq1''
and all sequence ids ``chr2'' to ``seq2''.


If the options -offsetfile and -chseqids are used together, the offsets are
applied before the sequence ids are changed.]])
