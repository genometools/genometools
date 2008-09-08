#include "genometools.h"

static void handle_error(GT_Error *err)
{
  fprintf(stderr, "error writing canvas %s\n", gt_error_get(err));
  exit(EXIT_FAILURE);
}

int main(int argc, char *argv[])
{
  const char *style_file, *gff3_file, *output_file, *seqid;
  GT_Style *style;
  GT_FeatureIndex *feature_index;
  GT_Range range;
  GT_Diagram *diagram;
  GT_Canvas *canvas;
  GT_Error *err = gt_error_new();

  if (argc != 4) {
    fprintf(stderr, "Usage: %s style_file GFF3_file output_file\n", argv[0]);
    return EXIT_FAILURE;
  }

  style_file = argv[1];
  gff3_file = argv[2];
  output_file = argv[3];

  /* create style */
  if (!(style = gt_style_new(false, err)))
    handle_error(err);

  /* load style file */
  if (gt_style_load_file(style, style_file, err))
    handle_error(err);

  /* create feature index */
  feature_index = gt_feature_index_new();

  /* add GFF3 file to index */
  if (gt_feature_index_add_gff3file(feature_index, gff3_file, err))
    handle_error(err);

  /* create diagram for first sequence ID in feature index */
  seqid = gt_feature_index_get_first_seqid(feature_index);
  gt_feature_index_get_range_for_seqid(feature_index, &range, seqid);
  diagram = gt_diagram_new(feature_index, seqid, &range, style);

  /* create canvas */
  canvas = gt_canvas_cairo_file_new(style, GRAPHICS_PNG, 800 /* width */, NULL);

  /* sketch diagram on canvas */
  gt_diagram_sketch(diagram, canvas);

  /* write canvas to file */
  if (gt_canvas_cairo_file_to_file((GT_CanvasCairoFile*) canvas, output_file,
                                   err)) {
    handle_error(err);
  }

  /* free */
  gt_canvas_delete(canvas);
  gt_diagram_delete(diagram);
  gt_feature_index_delete(feature_index);
  gt_style_delete(style);
  gt_error_delete(err);

  return EXIT_SUCCESS;
}
