/*
 *  Test program for smrng_lq().
 *    Command format: ./smrng_lq_tst k df alpha [nrng [xeps]]
 *
 *  Required functions:
 *    extern double smrng_lq()
 *      extern double smrng_lp()
 *        extern double rng_lp()
 *          extern double nrml_p()
 *
 */


#include <stdio.h>
#include <stdlib.h>
#include <math.h>

extern double smrng_lq(double p, int k, int df, int nrng,
                       double xeps, double peps, int *itr);

int main(int argc, char **argv)
{
  int k, df, itr, nrng=1;
  double x, x0, x1, alpha, xeps=1.0e-8, peps;

  if(argc < 4) {
    printf("Command format: smrng_lq_tst k df alpha [nrng [xeps]]\n");
    exit (1);
  }
  k = atoi(argv[1]);
  df = atoi(argv[2]);
  alpha = atof(argv[3]);
  if(argc >= 5)
    nrng = atoi(argv[4]);
  if(argc >= 6)
    xeps = atof(argv[5]);
  peps = alpha*xeps;

  x = smrng_lq(1.0 - alpha, k, df, nrng, xeps, peps, &itr);
  printf("itr = %4d, quantile = %20.16g\n", itr, x);

  // Interpolation between df=240 and infinity.
  if(df > 240)
    {
      x0 = smrng_lq(1.0-alpha, k, 0, nrng, xeps, peps, &itr);
      x1 = smrng_lq(1.0-alpha, k, 240, nrng, xeps, peps, &itr);
      x = (x1-x0) * (240.0/df) + x0;
      printf("Interpolation in 1/df\n"
             "itr = %4d, quantile = %20.16g\n", itr, x);
    }
  exit (0);
}
