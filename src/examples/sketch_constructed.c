#include "genometools.h"

static GT_Array* create_example_features(void)
{
  GT_Array *features;
  GT_GenomeNode *forward_gene, *reverse_gene, *exon, *intron; /* features */
  GT_Str *seqid; /* holds the sequence id the features refer to */
  GT_FeatureTypeFactory *type_factory; /* used the create feature types */
  GT_GenomeFeatureType *type; /* hold a feature type */
  GT_Range range; /* used to define intervals on the genomic sequence */

  /* construct the example features */
  features = gt_array_new(sizeof (GT_GenomeNode*));
  type_factory = gt_feature_type_factory_any_new();
  seqid = gt_str_new_cstr("chromosome_21");

  /* construct a gene on the forward strand with two exons */
  type = gt_feature_type_factory_create_gft(type_factory, "gene");
  range.start = 100; range.end = 900;
  forward_gene = gt_genome_feature_new(seqid, type, range, GT_STRAND_FORWARD);

  type = gt_feature_type_factory_create_gft(type_factory, "exon");
  range.start = 100; range.end = 200;
  exon = gt_genome_feature_new(seqid, type, range, GT_STRAND_FORWARD);
  /* exon belongs to forward gene */
  gt_genome_node_is_part_of_genome_node(forward_gene, exon);

  type = gt_feature_type_factory_create_gft(type_factory, "intron");
  range.start = 201; range.end = 799;
  intron = gt_genome_feature_new(seqid, type, range, GT_STRAND_FORWARD);
  /* intron belongs to forward gene */
  gt_genome_node_is_part_of_genome_node(forward_gene, intron);

  type = gt_feature_type_factory_create_gft(type_factory, "exon");
  range.start = 800; range.end = 900;
  exon = gt_genome_feature_new(seqid, type, range, GT_STRAND_FORWARD);
  /* exon belongs to forward gene */
  gt_genome_node_is_part_of_genome_node(forward_gene, exon);

  /* store forward gene in feature array */
  gt_array_add(features, forward_gene);

  /* construt a single-exon gene on the reverse strand
     (within the intron of the forward strand gene) */
  type = gt_feature_type_factory_create_gft(type_factory, "gene");
  range.start = 400; range.end = 600;
  reverse_gene = gt_genome_feature_new(seqid, type, range, GT_STRAND_REVERSE);

  type = gt_feature_type_factory_create_gft(type_factory, "exon");
  range.start = 400; range.end = 600;
  exon = gt_genome_feature_new(seqid, type, range, GT_STRAND_REVERSE);
  /* exon belongs to reverse gene */
  gt_genome_node_is_part_of_genome_node(reverse_gene, exon);

  /* store reverse gene in feature array */
  gt_array_add(features, reverse_gene);

  /* free */
  gt_str_delete(seqid);
  gt_feature_type_factory_delete(type_factory);

  return features;
}

static void handle_error(GT_Error *err)
{
  fprintf(stderr, "error writing canvas %s\n", gt_error_get(err));
  exit(EXIT_FAILURE);
}

static void draw_example_features(GT_Array *features, const char *style_file,
                                  const char *output_file)
{
  GT_Range range = { 1, 1000 }; /* the genomic range to draw */
  GT_Style *style;
  GT_Diagram *diagram;
  GT_Canvas *canvas;
  GT_Error *err = gt_error_new();

  /* create style */
  if (!(style = gt_style_new(false, err)))
    handle_error(err);

  /* load style file */
  if (gt_style_load_file(style, style_file, err))
    handle_error(err);

  /* create diagram */
  diagram = gt_diagram_new_from_array(features, &range, style);

  /* create canvas */
  canvas = gt_canvas_new(style, GRAPHICS_PNG, 800 /* width */, NULL);

  /* sketch diagram on canvas */
  gt_diagram_sketch(diagram, canvas);

  /* write canvas to file */
  if (gt_canvas_to_file(canvas, output_file, err))
    handle_error(err);

  /* free */
  gt_canvas_delete(canvas);
  gt_diagram_delete(diagram);
  gt_style_delete(style);
  gt_error_delete(err);
}

static void delete_example_features(GT_Array *features)
{
  unsigned long i;
  for (i = 0; i < gt_array_size(features); i++)
    gt_genome_node_rec_delete(*(GT_GenomeNode**) gt_array_get(features, i));
  gt_array_delete(features);
}

int main(int argc, char *argv[])
{
  GT_Array *features; /* stores the created example features */

  if (argc != 3) {
    fprintf(stderr, "Usage: %s style_file output_file\n", argv[0]);
    return EXIT_FAILURE;
  }

  features = create_example_features();

  draw_example_features(features, argv[1], argv[2]);

  delete_example_features(features);

  return EXIT_SUCCESS;
}
