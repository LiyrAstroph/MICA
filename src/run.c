#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "allvars.h"
#include "proto.h"

void run()
{
  read_param();
  init();
//  test();//
  if(flag_sim)
  {
    simulate_con(1.5, 80.0, 1.0, 10.0);
    simulate_line();
  }
  else
  {
    mcmc_con_run();
    mcmc_conline_run();
  }
}