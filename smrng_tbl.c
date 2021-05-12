/*
 *  This program tabulates the upper quantiles
 *    of the Studentised maximum range distribution.
 *
 *  command format: smrng_tbl k_end alpha [index [nrng]]
 *
 *  Arguments
 *    k_end:   k = 2, ..., k_end.
 *               If k_end > 100,
 *               k = 2, ..., 20, 50, 100, 200, 500, 1000.
 *    alpha:   upper probability
 *    [index]: If index==2, df runs from 1 to 40.
 *    [nrng]:  number of independent ranges
 *
 *  Required functions:
 *    extern double smrng_lq()
 *      extern double smrng_lp()
 *        extern double rng_lp()
 *          extern double nrml_p()
 *    static void line(int i)
 *
 *  Include files:
 *    <stdio.h>
 *    <stdlib.h>
 *    <math.h>
 *
 *  Note
 *    The table can be stored in a file by piping such as
 *      ./smrng_tbl 20 0.05 2 10 > smrng05.txt
 *
 *  Stored in:
 *    smrng_tbl.c
 *
 *  History
 *    2002-10-04: Last modified.
 *    2018-11-10: Created for the new version.
 *    2019-04-26: k_end > 100
 *    2021-05-12: Studentised maximum range
 *
 *  Coded by Tetsuhisa Miwa.
 */


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define EPS (1.0e-8)

extern double smrng_lq(double p, int k, int df, int nrng,
                       double xeps, double peps, int *itr);

static void line(int i)
{
  for( ; i > 0; i--)
    printf("-");
  printf("\n");
}

int main(int argc, char **argv)
{
  double  alpha, xeps, peps, q;
  int     kupper[5]={50, 100, 200, 500, 1000}, k[99], ke, j;
  int     index=1, nrng=1, df[106], i, itr, itrmax=0;
  FILE    *fout;

  if(argc < 3) {
    printf("command format: smrng_tbl k_end alpha [index [nrng]]\n");
    exit(1);
  }

  ke = atoi(argv[1]) - 2; // end value of k
  if(ke <= 98) {
    for(j=0; j <= ke; j++)
      k[j] = j + 2;
  } else {
    ke = 23;
    for(j=0; j <= 18; j++)
      k[j] = j + 2;
    for( ; j <= ke; j++)
      k[j] = kupper[j - 19];
  }

  alpha = atof(argv[2]);
  xeps = EPS;
  peps = alpha*EPS;

  if(argc >= 4) {
    index = atoi(argv[3]);
    if(index != 1)  // index value should be 1 or 2
      index = 2;
  }
  for(i=0; i < 20*index; i++)
    df[i] = i + 1;
  for(i=0; i < 5; i++)
    df[i + 20*index] = 120*index/(5 - i);
  df[5 + 20*index] = 0;

  if(argc >= 5)
    nrng = atoi(argv[4]);
  
  printf("The Studentised maximum range upper quantiles\n"
         "q(k, df, no.ranges=%4i; alpha=%5.2lf)\n", nrng, alpha);
  line(7*ke + 12);
  printf(" df  k->%3i", k[0]);
  for(j=1; j <= ke; j++)
    printf("%7i", k[j]);
  printf("\n");
  line(7*ke + 12);

  for(i=0; i < 6+20*index; i++){
    if(df[i] == 0)
      printf("Inf  ");
    else
      printf("%3i  ", df[i]);

    for(j=0; j <= ke; j++){
      q = smrng_lq(1.0-alpha, k[j], df[i], nrng, xeps, peps, &itr);
      if(q < 100.0)
        printf("%7.3lf", q);
      else
        printf("%7.2lf", q);
      if(itr > itrmax)
        itrmax = itr;
    }
    printf("\n");

    if((i+1)%10==0)
      line(7*ke+12);
    if((i+1)==20 && index==2){
      printf(" df  k->%3i", k[0]);
      for(j=1; j <= ke; j++)
        printf("%7i", k[j]);
      printf("\n");
      line(7*ke+12);
    }
  }
  line(7*ke+12);

  printf("max.iterations = %5i\n", itrmax);
}
