/*
  Copyright (c) 2005-2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2005-2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#include <gt.h>

#define DELIMITER         ','
#define FORWARDSTRANDCHAR '+'
#define REVERSESTRANDCHAR '-'

typedef struct {
  Str* id;
  bool forward;
  Array *exons; /* the exon ranges */
} SplicedAlignment;

static void initSplicedAlignment(SplicedAlignment *sa)
{
  assert(sa);
  sa->id = str_new();
  sa->forward = true;
  sa->exons = array_new(sizeof(Range));
}

static void parse_input_line(SplicedAlignment *alignment, const char *line,
                             unsigned long line_length)
{
  long leftpos, rightpos;
  unsigned long i = 0;
  Range exon;

#define CHECKLINELENGTH\
        if (i >= line_length) {                     \
          fprintf(stderr, "incomplete input line\n" \
                  "line=%s", line);                 \
          exit(EXIT_FAILURE);                       \
        }

  /* parsing id */
  for (;;) {
    CHECKLINELENGTH;
    if (line[i] == DELIMITER) {
      /* reference id has been saved, skip this character and break */
      i++;
      CHECKLINELENGTH;
      break;
    }
    else {
      /* save this character of the reference id */
      str_append_char(alignment->id, line[i]);
    }

    /* increase counter */
    i++;
  }

  /* parsing orientation */
  if (line[i] == FORWARDSTRANDCHAR)
    alignment->forward = true;
  else if (line[i] == REVERSESTRANDCHAR)
    alignment->forward = false;
  else {
    error("wrong formatted input line, orientation must be %c or %c\n"
          "line=%s", FORWARDSTRANDCHAR, REVERSESTRANDCHAR, line);
  }
  i++;
  CHECKLINELENGTH;

  if (line[i] != DELIMITER)
    error("incomplete input line\nline=%s", line);

  while (1) {
    if (line[i] == DELIMITER) {
      i++;
      CHECKLINELENGTH;
      if (sscanf(line+i, "%ld-%ld", &leftpos, &rightpos) != 2) {
        error("incomplete input line\nline=%s", line);
      }
      exon.start = leftpos;
      exon.end   = rightpos;

      /* save exon */
      array_add(alignment->exons, exon);
    }
    i++;
    if (i >= line_length)
      break;
  }

  /* alignment contains at least one exon */
  assert(array_size(alignment->exons));
}

static void parse_input_file(Array *spliced_alignments,
                             const char *file_name)
{
  FILE *input_file;
  SplicedAlignment sa;
  Str *line;

  line = str_new();
  input_file = xfopen(file_name, "r");

  while (str_read_next_line(line, input_file) != EOF) {
    /* init new spliced alignment */
    initSplicedAlignment(&sa);
    /* parse input line and save result in spliced alignment */
    parse_input_line(&sa, str_get(line), str_length(line));
    /* store spliced alignment */
    array_add(spliced_alignments, sa);
    /* reset array */
    str_reset(line);
  }

  xfclose(input_file);
  str_free(line);
}

static Range get_genomic_range(const void *sa)
{
  SplicedAlignment *alignment = (SplicedAlignment*) sa;
  Range range;
  assert(alignment);
  range.start = ((Range*) array_get_first(alignment->exons))->start;
  range.end   = ((Range*) array_get_last(alignment->exons))->end;
  return range;
}

static Strand get_strand(const void *sa)
{
  SplicedAlignment *alignment = (SplicedAlignment*) sa;
  if (alignment->forward)
    return STRAND_FORWARD;
  return STRAND_REVERSE;
}

static void get_exons(Array *exon_ranges, const void *sa)
{
  SplicedAlignment *alignment = (SplicedAlignment*) sa;
  assert(alignment);
  array_add_array(exon_ranges, alignment->exons);
}

static void process_splice_form(Array *spliced_alignments_in_form,
                                /*@unused@*/ const void *set_of_sas,
                                /*@unused@*/ unsigned long number_of_sas,
                                /*@unused@*/ size_t size_of_sa,
                                /*@unused@*/ void *userdata)
{
  unsigned long i;

  printf("contains [");
  for (i = 0; i < array_size(spliced_alignments_in_form); i++) {
    if (i)
      putchar(',');
    printf("%lu", *((unsigned long*) array_get(spliced_alignments_in_form, i)));
  }
  printf("]\n");
}

static int range_compare_long_first(Range range_a, Range range_b)
{
  assert(range_a.start <= range_a.end && range_b.start <= range_b.end);

  if ((range_a.start == range_b.start) && (range_a.end == range_b.end))
    return 0; /* range_a == range_b */

  if ((range_a.start < range_b.start) ||
      ((range_a.start == range_b.start) && (range_a.end > range_b.end)))
    return -1; /* range_a < range_b */

  return 1; /* range_a > range_b */
}

static int compare_spliced_alignment(const void *a, const void *b)
{
  SplicedAlignment *sa_a = (SplicedAlignment*) a,
                   *sa_b = (SplicedAlignment*) b;
  Range range_a, range_b;
  range_a.start = ((Range*) array_get_first(sa_a->exons))->start;
  range_a.end   = ((Range*) array_get_last(sa_a->exons))->end;
  range_b.start = ((Range*) array_get_first(sa_b->exons))->start;
  range_b.end   = ((Range*) array_get_last(sa_b->exons))->end;
  return range_compare_long_first(range_a, range_b);
}

int gt_consensus_sa(int argc, char *argv[], Error *err)
{
  Array *spliced_alignments;
  SplicedAlignment *sa;
  unsigned long i;
  error_check(err);

  /* check number of arguments */
  if (argc != 2) {
    fprintf(stderr, "Usage: %s spliced_alignment_file\n", argv[0]);
    fprintf(stderr, "Read file containing spliced alingments, compute "
                    "consensus spliced alignments, \n"
                    "and print them to stdout.\n");
    exit(EXIT_FAILURE); /* XXX */
  }

  /* parse input file and store resuilts in the spliced alignment array */
  spliced_alignments = array_new(sizeof(SplicedAlignment));
  parse_input_file(spliced_alignments, argv[1]);

  /* sort spliced alignments */
  qsort(array_get_space(spliced_alignments), array_size(spliced_alignments),
        sizeof(SplicedAlignment), compare_spliced_alignment);

  /* compute the consensus spliced alignments */
  consensus_sa(array_get_space(spliced_alignments),
               array_size(spliced_alignments), sizeof(SplicedAlignment),
               get_genomic_range, get_strand, get_exons, process_splice_form,
               NULL, NULL);

  /* free */
  for (i = 0; i < array_size(spliced_alignments); i++) {
    sa = array_get(spliced_alignments, i);
    str_free(sa->id);
    array_free(sa->exons);
  }
  array_free(spliced_alignments);

  return 0;
}
