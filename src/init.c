
#include <stdio.h>
#include <stdlib.h>

#include "allvars.h"
#include "proto.h"

/*
 *
 */
void memory_alloc()
{
  Tcon_data = malloc(ndata_max*sizeof(double));
  Fcon_data = malloc(ndata_max*sizeof(double));
  Fcerrs_data = malloc(ndata_max*sizeof(double));
}
