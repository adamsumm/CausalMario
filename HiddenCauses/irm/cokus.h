#ifndef COKUS_H
#define COKUS_H
#include <stdlib.h>

typedef unsigned long uint32;
void seedMT(uint32);
uint32 randomMT(void);
#define myrand() (double) (((unsigned long) randomMT()) / 4294967296.)

#endif /* COKUS_H*/
