#include <stdio.h>
#include <stdlib.h>

#include "allvars.h"
#include "proto.h"

void error_exit(int code)
{
  switch(code)
  {
    case 1:
      printf("no parameter file.\n");
      exit(1);

    case 2:
      printf("cannot open file %s.\n", str_error_exit);
      exit(2);

    case 3:
      printf("cannot recognize the parameter %s.\n", str_error_exit);
      exit(3);

    case 4:
      printf("file %s does not exist.\n", str_error_exit);
      exit(4);

    case 5:
      printf("end of parameter file %s.\n", fname_param);
      exit(5);

    case 6:
      printf("error in reading file %s.\n", str_error_exit);
      exit(6);

    default:
      printf("cannot recognize the error code.\n");
      exit(-1);
  }
}