#include "genometools.h"

static GtArray* create_example_features(void)
{
  GtArray *features;
  GtGenomeNode *gene, *exon, *intron; /* features */
  GtStr *seqid; /* holds the sequence id the features refer to */

  /* construct the example features */
  features = gt_array_new(sizeof (GtGenomeNode*));
  seqid = gt_str_new_cstr("chromosome_21");

  /* construct a gene on the forward strand with two exons */
  gene = gt_feature_node_new(seqid, "gene", 100, 900, GT_STRAND_FORWARD);
  exon = gt_feature_node_new(seqid, "exon", 100, 200, GT_STRAND_FORWARD);
  gt_feature_node_add_child((GtFeatureNode*) gene, (GtFeatureNode*) exon);
  intron = gt_feature_node_new(seqid, "intron", 201, 799, GT_STRAND_FORWARD);
  gt_feature_node_add_child((GtFeatureNode*) gene, (GtFeatureNode*) intron);
  exon = gt_feature_node_new(seqid, "exon", 800, 900, GT_STRAND_FORWARD);
  gt_feature_node_add_child((GtFeatureNode*) gene, (GtFeatureNode*) exon);

  /* store forward gene in feature array */
  gt_array_add(features, gene);

  /* construct a single-exon gene on the reverse strand
     (within the intron of the forward strand gene) */
  gene = gt_feature_node_new(seqid, "gene", 400, 600, GT_STRAND_REVERSE);
  exon = gt_feature_node_new(seqid, "exon", 400, 600, GT_STRAND_REVERSE);
  gt_feature_node_add_child((GtFeatureNode*) gene, (GtFeatureNode*) exon);

  /* store reverse gene in feature array */
  gt_array_add(features, gene);

  /* free */
  gt_str_delete(seqid);

  return features;
}

static void handle_error(GtError *err)
{
  fprintf(stderr, "error writing canvas %s\n", gt_error_get(err));
  exit(EXIT_FAILURE);
}

static void draw_example_features(GtArray *features, const char *style_file,
                                  const char *output_file)
{
  GtRange range = { 1, 1000 }; /* the genomic range to draw */
  GtStyle *style;
  GtDiagram *diagram;
  GtLayout *layout;
  GtCanvas *canvas;
  unsigned long height;
  GtError *err = gt_error_new();

  /* create style */
  if (!(style = gt_style_new(err)))
    handle_error(err);

  /* load style file */
  if (gt_style_load_file(style, style_file, err))
    handle_error(err);

  /* create diagram */
  diagram = gt_diagram_new_from_array(features, &range, style);

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
  if (gt_canvas_cairo_file_to_file((GtCanvasCairoFile*) canvas, output_file,
                                   err)) {
    handle_error(err);
  }

  /* free */
  gt_canvas_delete(canvas);
  gt_layout_delete(layout);
  gt_diagram_delete(diagram);
  gt_style_delete(style);
  gt_error_delete(err);
}

static void delete_example_features(GtArray *features)
{
  unsigned long i;
  for (i = 0; i < gt_array_size(features); i++)
    gt_genome_node_delete(*(GtGenomeNode**) gt_array_get(features, i));
  gt_array_delete(features);
}

int main(int argc, char *argv[])
{
  GtArray *features; /* stores the created example features */

  if (argc != 3) {
    fprintf(stderr, "Usage: %s style_file output_file\n", argv[0]);
    return EXIT_FAILURE;
  }

  gt_lib_init();

  features = create_example_features();

  draw_example_features(features, argv[1], argv[2]);

  delete_example_features(features);

  gt_lib_clean();
  return EXIT_SUCCESS;
}
