#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#include "allvars.h"
#include "proto.h"

/*
 * read parameter file
 */

void read_param()
{
  FILE *fp;
  char buf[200], buf1[100], buf2[100], buf3[100], buf4[100];

  fp = fopen(fname_param, "r");

  if(fp==NULL)
  {
    strcpy(str_error_exit, fname_param);
    error_exit(2);
  }

//*******************************************
// read file name for continuum

  buf1[0]='%';
  while(buf1[0]=='%')
  {
    if(feof(fp))
    {
      error_exit(5);
    }

    fgets(buf, 200, fp);
    if(sscanf(buf, "%s%s%s", buf1, buf2, buf3)<1)
    {
      buf1[0]='%';
    }
  }
  strcpy(fname_con, buf1);
  printf("%s\n", fname_con);
  if((access(fname_con, 0))!=0)
  {
    strcpy(str_error_exit, fname_con);
    error_exit(4);
  }
  

//*******************************************
// read file name for line

  buf1[0]='%';
  while(buf1[0]=='%')
  {
    if(feof(fp))
    {
      error_exit(5);
    }

    fgets(buf, 200, fp);
    if(sscanf(buf, "%s%s%s", buf1, buf2, buf3)<1)
    {
      buf1[0]='%';
    }
  }
  strcpy(fname_line, buf1);
  printf("%s\n", fname_line);
  if((access(fname_line, 0))!=0)
  {
    strcpy(str_error_exit, fname_line);
    error_exit(4);
  }

//*******************************************
// read range for time lag to be solved.
  buf1[0]='%';
  while(buf1[0]=='%')
  {
    if(feof(fp))
    {
      error_exit(5);
    }

    fgets(buf, 200, fp);
    if(sscanf(buf, "%s%s%s%s", buf1, buf2, buf3, buf4)<2)
    {
      buf1[0]='%';
    }
  }
  if(strcmp(buf1, "taulim")!=0)
  {
    strcat(buf1, ". expecting taulim");
    strcpy(str_error_exit, buf1);
    error_exit(3);
  }
  tau_lim_low=atof(buf2);
  tau_lim_up=atof(buf3);
  printf("%f %f\n", tau_lim_low, tau_lim_up);

//*******************************************
// read number nc.
  buf1[0]='%';
  while(buf1[0]=='%')
  {
    if(feof(fp))
    {
      error_exit(5);
    }

    fgets(buf, 200, fp);
    if(sscanf(buf, "%s%s%s", buf1, buf2, buf3)<1)
    {
      buf1[0]='%';
    }
  }
  if(strcmp(buf1, "nc")!=0)
  {
    strcat(buf1, ". expecting nc");
    strcpy(str_error_exit, buf1);
    error_exit(3);
  }
  nc=atoi(buf2);
  printf("%d\n", nc);

//*******************************************
// read number nmcmc.
  buf1[0]='%';
  while(buf1[0]=='%')
  {
    if(feof(fp))
    {
      error_exit(5);
    }

    fgets(buf, 200, fp);
    if(sscanf(buf, "%s%s%s", buf1, buf2, buf3)<1)
    {
      buf1[0]='%';
    }
  }
  if(strcmp(buf1, "nmcmc")!=0)
  {
    strcat(buf1, ". expecting nmcmc");
    strcpy(str_error_exit, buf1);
    error_exit(3);
  }
  nmcmc=atoi(buf2);
  printf("%d\n", nmcmc);

//*******************************************
// read number nbuilt.
  buf1[0]='%';
  while(buf1[0]=='%')
  {
    if(feof(fp))
    {
      error_exit(5);
    }

    fgets(buf, 200, fp);
    if(sscanf(buf, "%s%s%s", buf1, buf2, buf3)<1)
    {
      buf1[0]='%';
    }
  }
  if(strcmp(buf1, "nbuilt")!=0)
  {
    strcat(buf1, ". expecting nbuilt");
    strcpy(str_error_exit, buf1);
    error_exit(3);
  }
  nbuilt=atoi(buf2);
  printf("%d\n", nbuilt);

//*******************************************
// read flag for detrending.
  buf1[0]='%';
  while(buf1[0]=='%')
  {
    if(feof(fp))
    {
      error_exit(5);
    }

    fgets(buf, 200, fp);
    if(sscanf(buf, "%s%s%s", buf1, buf2, buf3)<1)
    {
      buf1[0]='%';
    }
  }
  if(strcmp(buf1, "flag_detrend")!=0)
  {
    strcat(buf1, ". expecting flag_detrend");
    strcpy(str_error_exit, buf1);
    error_exit(3);
  }
  flag_detrend=atoi(buf2);
  printf("%d\n", flag_detrend);  

//*******************************************
// read flag for sim.
  buf1[0]='%';
  while(buf1[0]=='%')
  {
    if(feof(fp))
    {
      error_exit(5);
    }

    fgets(buf, 200, fp);
    if(sscanf(buf, "%s%s%s", buf1, buf2, buf3)<1)
    {
      buf1[0]='%';
    }
  }
  if(strcmp(buf1, "flag_sim")!=0)
  {
    strcat(buf1, ". expecting flag_sim");
    strcpy(str_error_exit, buf1);
    error_exit(3);
  }
  flag_sim=atoi(buf2);
  printf("%d\n", flag_sim);    
}

/*
 * read data files
 */

void read_data()
{
  FILE *fp;
  int i;
  char buf[200];

//*************************************
// read continuum data

  fp = fopen(fname_con, "r");
  if(fp==NULL)
  {
    strcpy(str_error_exit, fname_con);
    error_exit(2);
  }

  i = 0;
  while(!feof(fp))
  {
    if(fgets(buf, 200, fp)==NULL)
      break;
    if(sscanf(buf, "%lf %lf %lf", &Tcon_data[i], &Fcon_data[i], &Fcerrs_data[i]) < 3)
    {
      strcpy(str_error_exit, fname_con);
      error_exit(6);
    }
    i++;

    if(i > ndata_max)
    {
      strcpy(str_error_exit, fname_con);
      error_exit(8);
    }
  }
  ncon_data = i;
  printf("ncon_data: %d\n", ncon_data);
  fclose(fp);

//*************************************
// read line data

  fp = fopen(fname_line, "r");
  if(fp==NULL)
  {
    strcpy(str_error_exit, fname_line);
    error_exit(2);
  }

  i = 0;
  while(!feof(fp))
  {
    if(fgets(buf, 200, fp)==NULL)
      break;
    if(sscanf(buf, "%lf %lf %lf", &Tline_data[i], &Fline_data[i], &Flerrs_data[i]) < 3)
    {
      strcpy(str_error_exit, fname_line);
      error_exit(6);
    }
    i++;

    if(i > ndata_max)
    {
      strcpy(str_error_exit, fname_line);
      error_exit(8);
    }
  }
  nline_data = i;
  printf("nline_data: %d\n", nline_data); 
  fclose(fp); 
}
