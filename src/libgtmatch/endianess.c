#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>

bool islittleendian(void)
{
  int x = 1;

  if (*(char *) &x == 1)
  {
    return true;
  } else
  {
    return false;
  }
}
