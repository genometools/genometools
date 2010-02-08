#include "genometools.h"

/* This example shows a custom stream implementation (i.e., an implementation of
   the NodeStream interface). In this case, an input stream is implemented, but
   intermediate and output streams are implemented similarily. */

/* public custom stream definition */
typedef struct CustomStream CustomStream;

/* private custom stream definition */
struct CustomStream {
  const GtNodeStream parent_instance;
  bool emitted_standard_gene;
};

static const GtNodeStreamClass* custom_stream_class(void);

#define custom_stream_cast(NS)\
        gt_node_stream_cast(custom_stream_class(), NS)

static int custom_stream_next(GtNodeStream *ns, GtGenomeNode **gn,
                              GT_UNUSED GtError *err)
{
  CustomStream *cs = custom_stream_cast(ns);
  gt_error_check(err);
  if (!cs->emitted_standard_gene) {
    *gn = gt_feature_node_new_standard_gene();
    cs->emitted_standard_gene = true;
  }
  else
    *gn = NULL;
  return 0;
}

static const GtNodeStreamClass* custom_stream_class(void)
{
 static const GtNodeStreamClass *nsc = NULL;
 if (!nsc) {
   nsc = gt_node_stream_class_new(sizeof (CustomStream),
                                  NULL, /* we don't need a dedicated free
                                           method, because our CustomStream
                                           doesn't own object which have to be
                                           freed separately */
                                  custom_stream_next);
 }
 return nsc;
}

static GtNodeStream* custom_stream_new(void)
{
  GtNodeStream *ns = gt_node_stream_create(custom_stream_class(), false);
  CustomStream *cs = custom_stream_cast(ns);
  cs->emitted_standard_gene = false;
  return ns;
}

int main(GT_UNUSED int argc, const char *argv[])
{
  GtNodeStream *custom_stream, *gff3_out_stream;
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

  /* create custom input stream  */
  custom_stream = custom_stream_new();

  /* create GFF3 output stream */
  gff3_out_stream = gt_gff3_out_stream_new(custom_stream, NULL);

  /* pull the features through the stream */
  had_err = gt_node_stream_pull(gff3_out_stream, err);

  /* handle error */
  if (had_err)
    fprintf(stderr, "%s: error: %s\n", argv[0], gt_error_get(err));

  /* free */
  gt_node_stream_delete(gff3_out_stream);
  gt_node_stream_delete(custom_stream);
  gt_error_delete(err);

  if (had_err)
    return EXIT_FAILURE;
  return EXIT_SUCCESS;
}
