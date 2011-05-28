#include "genometools.h"

int main(int argc, const char *argv[])
{
  GtNodeStream *gff3_in_stream, *sort_stream, *gff3_out_stream;
  GtGenomeNode *gn;
  GtError *err;
  int had_err;

  gt_lib_init();
  err = gt_error_new();

  gff3_in_stream = gt_gff3_in_stream_new_unsorted(argc-1, argv+1);
  sort_stream = gt_sort_stream_new(gff3_in_stream);
  gff3_out_stream = gt_gff3_out_stream_new(sort_stream, NULL);

  while (!(had_err = gt_node_stream_next(gff3_out_stream, &gn, err)) && gn)
    gt_genome_node_delete(gn);

  if (had_err)
    fprintf(stderr, "%s: error: %s\n", argv[0], gt_error_get(err));

  gt_node_stream_delete(gff3_out_stream);
  gt_node_stream_delete(sort_stream);
  gt_node_stream_delete(gff3_in_stream);
  gt_error_delete(err);

  if (had_err)
    return EXIT_FAILURE;
  return EXIT_SUCCESS;
}
