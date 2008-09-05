#include "genometools.h"

static void handle_error(Error *err)
{
  fprintf(stderr, "error writing canvas %s\n", error_get(err));
  exit(EXIT_FAILURE);
}

int main(int argc, char *argv[])
{
  const char *style_file, *gff3_file, *output_file, *seqid;
  Style *style;
  FeatureIndex *feature_index;
  Range range;
  Diagram *diagram;
  Canvas *canvas;
  Error *err = error_new();

  if (argc != 4) {
    fprintf(stderr, "Usage: %s style_file GFF3_file output_file\n", argv[0]);
    return EXIT_FAILURE;
  }

  style_file = argv[1];
  gff3_file = argv[2];
  output_file = argv[3];

  /* create style */
  if (!(style = style_new(false, err)))
    handle_error(err);

  /* load style file */
  if (style_load_file(style, style_file, err))
    handle_error(err);

  /* create feature index */
  feature_index = feature_index_new();

  /* add GFF3 file to index */
  if (feature_index_add_gff3file(feature_index, gff3_file, err))
    handle_error(err);

  /* create diagram for first sequence ID in feature index */
  seqid = feature_index_get_first_seqid(feature_index);
  range = feature_index_get_range_for_seqid(feature_index, seqid);
  diagram = diagram_new(feature_index, seqid, &range, style);

  /* create canvas */
  canvas = canvas_new(style, GRAPHICS_PNG, 800 /* width */, NULL);

  /* sketch diagram on canvas */
  diagram_sketch(diagram, canvas);

  /* write canvas to file */
  if (canvas_to_file(canvas, output_file, err))
    handle_error(err);

  /* free */
  canvas_delete(canvas);
  diagram_delete(diagram);
  feature_index_delete(feature_index);
  style_delete(style);
  error_delete(err);

  return EXIT_SUCCESS;
}
