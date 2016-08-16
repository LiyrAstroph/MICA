/*
 * MICA (Multiple and Inhomogeneous Component Analysis)
 * A Non-parameteric ApproaCH to Constrain the Transfer Function in Reverberation Mapping
 * 
 * Yan-Rong Li, liyanrong@mail.ihep.ac.cn
 *
 * Reference: Li, Y.-R. et al. 2016, arXiv:1608.03741
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