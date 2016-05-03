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
    simulate_con(0.6, 100.0, 1.0, 10.0);
    simulate_line();
  }
  else
  {
    mcmc_con_run();
    mcmc_conline_run();
  }
}