#include "genometools.h"

static void handle_error(GtError *err)
{
  fprintf(stderr, "error: %s\n", gt_error_get(err));
  exit(EXIT_FAILURE);
}

int main(int argc, char *argv[])
{
  const char *style_file, *png_file, *gff3_file;
  char *seqid;
  GtStyle *style;
  GtFeatureIndex *feature_index;
  GtRange range;
  GtDiagram *diagram;
  GtLayout *layout;
  GtCanvas *canvas;
  unsigned long height;
  GtError *err;

  if (argc != 4) {
    fprintf(stderr, "Usage: %s style_file PNG_file GFF3_file\n", argv[0]);
    return EXIT_FAILURE;
  }

  style_file = argv[1];
  png_file = argv[2];
  gff3_file = argv[3];

  /* initialize */
  gt_lib_init();

  /* create error object */
  err = gt_error_new();

  /* create style */
  if (!(style = gt_style_new(err)))
    handle_error(err);

  /* load style file */
  if (gt_style_load_file(style, style_file, err))
    handle_error(err);

  /* create feature index */
  feature_index = gt_feature_index_memory_new();

  /* add GFF3 file to index */
  if (gt_feature_index_add_gff3file(feature_index, gff3_file, err))
    handle_error(err);

  /* create diagram for first sequence ID in feature index */
  if (!(seqid = gt_feature_index_get_first_seqid(feature_index, err))) {
    if (gt_error_is_set(err))
      handle_error(err);
  }
  if (gt_feature_index_get_range_for_seqid(feature_index, &range, seqid, err))
    handle_error(err);
  diagram = gt_diagram_new(feature_index, seqid, &range, style, err);
  gt_free(seqid);
  if (gt_error_is_set(err))
    handle_error(err);

  /* create layout with given width, determine resulting image height */
  layout = gt_layout_new(diagram, 600, style, err);
  if (!layout)
    handle_error(err);
  if (gt_layout_get_height(layout, &height, err))
    handle_error(err);

  /* create PNG canvas */
  canvas = gt_canvas_cairo_file_new(style, GT_GRAPHICS_PNG, 600, height,
                                    NULL, err);
  if (!canvas)
    handle_error(err);

  /* sketch layout on canvas */
  if (gt_layout_sketch(layout, canvas, err))
    handle_error(err);

  /* write canvas to file */
  if (gt_canvas_cairo_file_to_file((GtCanvasCairoFile*) canvas, png_file, err))
    handle_error(err);

  /* free */
  gt_canvas_delete(canvas);
  gt_layout_delete(layout);
  gt_diagram_delete(diagram);
  gt_feature_index_delete(feature_index);
  gt_style_delete(style);
  gt_error_delete(err);
  /* perform static data cleanup */
  gt_lib_clean();
  return EXIT_SUCCESS;
}
