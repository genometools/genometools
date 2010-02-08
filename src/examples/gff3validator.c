#include "genometools.h"

/* A simple GFF3 validator which validates all given GFF3 files (without
   checking types). */

int main(int argc, const char *argv[])
{
  GtNodeStream *gff3_in_stream;
  GtGenomeNode *gn;
  GtError *err;
  int had_err;

  if (gt_version_check(GT_MAJOR_VERSION, GT_MINOR_VERSION, GT_MICRO_VERSION)) {
    fprintf(stderr, "error: %s\n", gt_version_check(GT_MAJOR_VERSION,
                                                    GT_MINOR_VERSION,
                                                    GT_MICRO_VERSION));
    return EXIT_FAILURE;
  }

  /* initialize */
  gt_lib_init();

  /* create error object */
  err = gt_error_new();

  /* create GFF3 input stream (with ID attribute checking) */
  gff3_in_stream = gt_gff3_in_stream_new_unsorted(argc-1, argv+1);
  gt_gff3_in_stream_check_id_attributes((GtGFF3InStream*) gff3_in_stream);

  /* pull the features through the stream and free them afterwards */
  while (!(had_err = gt_node_stream_next(gff3_in_stream, &gn, err)) && gn)
    gt_genome_node_delete(gn);

  /* handle error */
  if (had_err)
    fprintf(stderr, "%s: error: %s\n", argv[0], gt_error_get(err));
  else
    printf("input is valid GFF3\n");

  /* free */
  gt_node_stream_delete(gff3_in_stream);
  gt_error_delete(err);

  if (had_err)
    return EXIT_FAILURE;
  return EXIT_SUCCESS;
}
