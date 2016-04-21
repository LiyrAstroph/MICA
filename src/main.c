/*
 * PATCH
 * a non-Parameteric ApproaCH to constrain tthe transfer function in reverberation mapping
 * 
 * Yan-Rong Li, liyanrong@mail.ihep.ac.cn
 *
 * Reference: Li, Y.-R. et al. 2016, ApJ
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "allvars.h"
#include "proto.h"

int main(int argc, char** argv)
{
  printf("=========Start of RMPATCH==============\n");
  if(argc<2)
  {
    error_exit(1);
  }

  strcpy(fname_param, argv[1]);
  
  run();

  printf("========= End of RMPATCH ==============\n");
  return 0;
}