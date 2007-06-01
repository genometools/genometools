#include <string.h>

#define GZIPSUFFIX        ".gz"
#define GZIPSUFFIXLENGTH  (sizeof (GZIPSUFFIX)-1)

unsigned char checkgzipsuffix(const char *filename)
{
  size_t filenamelength = strlen (filename);

  if (filenamelength < GZIPSUFFIXLENGTH ||
      strcmp (filename + filenamelength - GZIPSUFFIXLENGTH, GZIPSUFFIX) != 0)
  {
    return 0;
  }
  return (unsigned char) 1;
}
