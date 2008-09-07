#include "genometools.h"

static Array* create_example_features(void)
{
  Array *features;
  GenomeNode *forward_gene, *reverse_gene, *exon, *intron; /* actual features */
  Str *seqid; /* holds the sequence id the features refer to */
  FeatureTypeFactory *type_factory; /* used the create feature types */
  GenomeFeatureType *type; /* hold a feature type */
  Range range; /* used to define intervals on the genomic sequence */

  /* construct the example features */
  features = array_new(sizeof (GenomeNode*));
  type_factory = feature_type_factory_any_new();
  seqid = str_new_cstr("chromosome_21");

  /* construct a gene on the forward strand with two exons */
  type = feature_type_factory_create_gft(type_factory, "gene");
  range.start = 100; range.end = 900;
  forward_gene = genome_feature_new(seqid, type, range, STRAND_FORWARD);

  type = feature_type_factory_create_gft(type_factory, "exon");
  range.start = 100; range.end = 200;
  exon = genome_feature_new(seqid, type, range, STRAND_FORWARD);
  /* exon belongs to forward gene */
  genome_node_is_part_of_genome_node(forward_gene, exon);

  type = feature_type_factory_create_gft(type_factory, "intron");
  range.start = 201; range.end = 799;
  intron = genome_feature_new(seqid, type, range, STRAND_FORWARD);
  /* intron belongs to forward gene */
  genome_node_is_part_of_genome_node(forward_gene, intron);

  type = feature_type_factory_create_gft(type_factory, "exon");
  range.start = 800; range.end = 900;
  exon = genome_feature_new(seqid, type, range, STRAND_FORWARD);
  /* exon belongs to forward gene */
  genome_node_is_part_of_genome_node(forward_gene, exon);

  /* store forward gene in feature array */
  array_add(features, forward_gene);

  /* construct a single-exon gene on the reverse strand
     (within the intron of the forward strand gene) */
  type = feature_type_factory_create_gft(type_factory, "gene");
  range.start = 400; range.end = 600;
  reverse_gene = genome_feature_new(seqid, type, range, STRAND_REVERSE);

  type = feature_type_factory_create_gft(type_factory, "exon");
  range.start = 400; range.end = 600;
  exon = genome_feature_new(seqid, type, range, STRAND_REVERSE);
  /* exon belongs to reverse gene */
  genome_node_is_part_of_genome_node(reverse_gene, exon);

  /* store reverse gene in feature array */
  array_add(features, reverse_gene);

  /* free */
  feature_type_factory_delete(type_factory);

  return features;
}

static void handle_error(Error *err)
{
  fprintf(stderr, "error writing canvas %s\n", error_get(err));
  exit(EXIT_FAILURE);
}

static void draw_example_features(Array *features, const char *style_file,
                                  const char *output_file)
{
  Range range = { 1, 1000 }; /* the genomic range to draw */
  Style *style;
  Diagram *diagram;
  Canvas *canvas;
  Error *err = error_new();

  /* create style */
  if (!(style = style_new(false, err)))
    handle_error(err);

  /* load style file */
  if (style_load_file(style, style_file, err))
    handle_error(err);

  /* create diagram */
  diagram = diagram_new_from_array(features, &range, style);

  /* create canvas */
  canvas = canvas_cairo_file_new(style, GRAPHICS_PNG, 800 /* width */, NULL);

  /* sketch diagram on canvas */
  diagram_sketch(diagram, canvas);

  /* write canvas to file */
  if (canvas_cairo_file_to_file((CanvasCairoFile*) canvas, output_file, err))
    handle_error(err);

  /* free */
  canvas_delete(canvas);
  diagram_delete(diagram);
  style_delete(style);
  error_delete(err);
}

static void delete_example_features(Array *features)
{
  unsigned long i;
  for (i = 0; i < array_size(features); i++)
    genome_node_rec_delete(*(GenomeNode**) array_get(features, i));
  array_delete(features);
}

int main(int argc, char *argv[])
{
  Array *features; /* stores the created example features */

  if (argc != 3) {
    fprintf(stderr, "Usage: %s style_file output_file\n", argv[0]);
    return EXIT_FAILURE;
  }

  features = create_example_features();

  draw_example_features(features, argv[1], argv[2]);

  delete_example_features(features);

  return EXIT_SUCCESS;
}
