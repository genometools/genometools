typedef struct
{
  const char *name;
  int rank;
  bool incolumn;
} GtSEdisplayStruct;

static GtSEdisplayStruct gt_display_arguments_table[] =
{
  {"alignment", 0, false},
  {"bitscore", 11, true},
  {"cigar", 1, true},
  {"evalue", 8, true},
  {"failed_seed", 5, false},
  {"fstperquery", 3, false},
  {"polinfo", 2, false},
  {"q.seqdesc", 10, true},
  {"s.seqdesc", 9, true},
  {"seed", 4, true},
  {"seed_in_algn", 6, true},
  {"seqlength", 7, true}
};

typedef enum
{
  Gt_Alignment_display,
  Gt_Cigar_display,
  Gt_Polinfo_display,
  Gt_Fstperquery_display,
  Gt_Seed_display,
  Gt_Failed_seed_display,
  Gt_Seed_in_algn_display,
  Gt_Seqlength_display,
  Gt_Evalue_display,
  Gt_S_seqdesc_display,
  Gt_Q_seqdesc_display,
  Gt_Bitscore_display
} GtSeedExtendDisplay_enum;

const char *gt_querymatch_display_help(void)
{
  return "specify what information about the matches to display\n"
         "alignment:    display alignment (possibly followed by \n"
         "              =<number> to specify width of alignment columns)\n"
         "cigar:        show cigar string representing alignment\n"
         "polinfo:      display polishing information for displayed \n"
         "              alignment\n"
         "fstperquery:  output only the first found match per query\n"
         "seed:         display the seed of the match, i.e. the length and \n"
         "              the start position of the seed in both instances\n"
         "failed_seed:  display the seed of the match that was extended, but \n"
         "              failed (after extension) the filter conditions\n"
         "seed_in_algn: display the seed in alignment\n"
         "seqlength:    display length of sequences in which the two \n"
         "              match-instances occur\n"
         "evalue:       display evalue\n"
         "s.seqdesc:    display sequence description of subject sequence\n"
         "q.seqdesc:    display sequence description of query sequence\n"
         "bitscore:     display bit score\n"
;}
