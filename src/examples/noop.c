/* Public domain. */

#include "genometools.h" /* include the GenomeTools ``all-in-one'' header */

int main(GT_UNUSED int argc, GT_UNUSED char *argv[])
{
  gt_lib_init();
  if (gt_version_check(GT_MAJOR_VERSION, GT_MINOR_VERSION, GT_MICRO_VERSION)) {
    fprintf(stderr, "error: %s\n", gt_version_check(GT_MAJOR_VERSION,
                                                    GT_MINOR_VERSION,
                                                    GT_MICRO_VERSION));
    return EXIT_FAILURE;
  }

  if (gt_lib_clean())
    return GT_EXIT_PROGRAMMING_ERROR; /* programmer error */
  return EXIT_SUCCESS;
}
