/* Public domain. */

#include "genometools.h" /* include the GenomeTools ``all-in-one'' header */

int main(GT_UNUSED int argc, GT_UNUSED char *argv[])
{
  if (gt_version_check(GT_MAJOR_VERSION, GT_MINOR_VERSION, GT_MICRO_VERSION)) {
    fprintf(stderr, "error: %s\n", gt_version_check(GT_MAJOR_VERSION,
                                                    GT_MINOR_VERSION,
                                                    GT_MICRO_VERSION));
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}
