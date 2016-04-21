
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "allvars.h"
#include "proto.h"

void init()
{
  nall = ncon + nline;
  memory_alloc();
  read_data();
}

/*
 *
 */
void memory_alloc()
{
  workspace_ipiv = malloc(ndata_max*sizeof(double));
  workspace = malloc(10*nall_data*sizeof(double) + 10);

  Tcon_data = malloc(ndata_max*sizeof(double));
  if(Tcon_data==NULL)
  {
    strcpy(str_error_exit, "Tcon_data");
    error_exit(7);
  }

  Fcon_data = malloc(ndata_max*sizeof(double));
  if(Fcon_data==NULL)
  {
    strcpy(str_error_exit, "Fcon_data");
    error_exit(7);
  }

  Fcerrs_data = malloc(ndata_max*sizeof(double));
  if(Fcerrs_data==NULL)
  {
    strcpy(str_error_exit, "Fcerrs_data");
    error_exit(7);
  }

  Tline_data = malloc(ndata_max*sizeof(double));
  if(Tline_data==NULL)
  {
    strcpy(str_error_exit, "Tline_data");
    error_exit(7);
  }

  Fline_data = malloc(ndata_max*sizeof(double));
  if(Fline_data==NULL)
  {
    strcpy(str_error_exit, "Fline_data");
    error_exit(7);
  }

  Flerrs_data = malloc(ndata_max*sizeof(double));
  if(Flerrs_data==NULL)
  {
    strcpy(str_error_exit, "Flerrs_data");
    error_exit(7);
  }

  Smat = array_malloc(nall_data*nall_data);
  Nmat = array_malloc(nall_data*nall_data);
  Cmat = array_malloc(nall_data*nall_data);
  ICmat = array_malloc(nall_data*nall_data);
  ICvmat = array_malloc(nall_data*nall_data);

  Tmat1 = array_malloc(nall_data*nall_data);
  Tmat2 = array_malloc(nall_data*nall_data);
  Tmat3 = array_malloc(nall_data*nall_data);
  Tmat4 = array_malloc(nall_data*nall_data);

  Tcon = array_malloc(ncon);
  Fcon = array_malloc(ncon);
  Fcerrs = array_malloc(ncon);
  Tline = array_malloc(nline);
  Fline = array_malloc(nline);
  Flerrs = array_malloc(nline);

}
