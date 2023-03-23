#include "arraymalloc.h"
#include <stdlib.h>

void **arraymalloc2d(int nx, int ny, size_t typesize)
{
  int i;
  void **array2d;

  size_t mallocsize;

  // total memory requirements including pointers

  mallocsize = nx*sizeof(void *) + nx*ny*typesize;

  array2d = (void **) malloc(mallocsize);

  // set first pointer to first element of data

  array2d[0] = (void *) (array2d + nx);

  for(i=1; i < nx; i++)
    {
      // set other pointers to point at subsequent rows

      array2d[i] = (void *) (((char *) array2d[i-1]) + ny*typesize);
    }

  return array2d;
}
