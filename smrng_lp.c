/*
 *  double smrng_lp(double q, int k, int df, int nrng)
 *    returns lower probability of
 *    the Studentised maximum range distribution.
 *
 *  Arguments
 *    q:    Studentised maximum range value
 *    k:    number of treatments for each range
 *    df:   error degrees of freedom (df<=0 means df=infinity)
 *    nrng: number of independent ranges
 *
 *  Required functions
 *    extern double rng_lp()
 *    static double rupper()
 *    static double rlower()
 *    static double chi2u()
 *    static double chi2l()
 *    static double coef()
 *    static double f()
 *
 *  Include files
 *    <math.h>
 *
 *  References
 *    Copenhaver, M. D. and B. Holland (1988).
 *      Computation of the distribution of the maximum Studentized
 *      range statistic with application to multiple significance
 *      testing of simple effects,
 *      J. Statist. Comput. Siml., vol. 30, 1 -- 15.
 *
 *  Note
 *    1) The 40-node Gauss-Legendre quadrature is used.
 *    2) The accuracy is of order e-11 or more (I hope).
 *    3) This accuracy is not guaranteed for k > 1000 or nrng > 100.
 *    4) Integrates twice if ru/q < su (ru: upper limit of max range).
 *
 *  Stored in
 *   smrng_lp.c
 *
 *  History
 *    c. 1994:    First written in Fortran for Studentised range.
 *    2018-11-02: Created for the new version.
 *    2021-05-10: Consider maximum of several ranges.
 *
 *  License
 *    GPLv3 (Free and No Warranty)
 *    https://www.gnu.org/licenses/
 *
 *  Coded by Tetsuhisa Miwa.
 */


#include <math.h>
#define LOGSQRTPI 0.572364942924700087071713675676529356  // log(sqrt(pi))

extern double rng_lp(double r, int k);

/* Upper limit of max range with approx upper prob=0.5e-13.
 */
static double rupper(int k, int nrng)
{
  double  rn1, rn100, y;

  rn1 = 0.42*pow(log(k - 0.5), 0.9) + 10.465;
  if(nrng <= 1)
    return(rn1);
  rn100 = 0.2866*pow(log(k - 0.9), 1.05) + 11.451;
  y = 0.2273*(rn100 - rn1)*pow(log((double)nrng), 0.97) + rn1;
  return(y);
}

/* Lower limit of max range with approx lower prob=0.5e-13.
 */
static double rlower(int k, int nrng)
{
  double  z1, z100, dk, bk, x1, x100, x, z;

  if(k <= 40) {
    z1 = -27.12/pow(log(k + 0.5), 2.1) + 1.8800;
    if(nrng <= 1)
      return(exp(z1));
    z100 = -5.749/pow(log((double)k), 0.12) + 6.4651;
    dk = (k < 8) ? 2.934/(k + 1.0) + 0.522 : 0.86 - 0.0015*k;
    bk = (k < 8) ? 7.88/(k + 2.0) + 0.112 : 16.875/(k + 10.0) - 0.0375;
    x1 = 1.0/pow(log(1.0 + dk), bk);
    x100 = 1.0/pow(log(100.0 + dk), bk);
    x = 1.0/pow(log(nrng + dk), bk);
    z = exp((z100 - z1)/(x100 - x1)*(x - x1) + z1);
  }
  else {
    z1 = 449.4*pow(log(k + 10.0), 0.012) - 455.6678;
    if(nrng <= 1)
      return(z1);
    z100 = 3.149*pow(log(k + 1.0), 0.48) - 1.2017;
    bk = (k <= 55) ? -0.08478*log((double)k) + 0.5738 
      : 0.03220*log((double)k) + 0.1050;
    x1 = pow(log(2.0), bk);
    x100 = pow(log(101.0), bk);
    x = pow(log(nrng + 1.0), bk);
    z = (z100 - z1)/(x100 - x1)*(x - x1) + z1;
  }
  return(z);
}

/* Upper limit for chi^2(df) with approx upper prob=0.5e-13.
 */
static double chi2u(int df)
{
  double  first[5]={56.73, 61.26, 65.01, 68.38, 71.50};
  double  w, z, ddf=2.0/9.0/df;

  if(df <= 5)
    return(first[df-1]);
  if(df <= 20)
    w = 7.391 - 3.050/df + 5.208/(df*df);
  else
    w = 7.441 - 5.209/df + 29.27/(df*df);
  // Wilson-Hilferty approximation.
  z = df*pow(w*sqrt(ddf) + (1.0 - ddf), 3.0);
  return(z);
}

/* Lower limit for chi^2(df) with approx lower prob=0.5e-13.
 */
static double chi2l(int df)
{
  double  first[5]={3.926e-27, 1.0e-13, 3.281e-09, 6.324e-07, 1.546e-05};
  double  w, z, ddf=2.0/9.0/df;

  if(df <= 5)
    return(first[df-1]);
  if(df <= 20) {
    // Log approximation.
    w = -8.645 - 70.72/df + 77.47/(df*df);
    z = df*exp(w/sqrt(0.5*df)-1.0/df);
  }
  else {
    // Wilson-Hilferty approximation.
    w = -7.451 + 10.07/df + 82.83/(df*df);
    z = df*pow(w*sqrt(ddf) + (1.0 - ddf), 3.0);
  }    
  return(z);
}

/* Coefficient of chi distribution (Note: not chi^2 distribution).
 */
static double coef(int df)
{
  int n;
  double g = (df%2 == 1) ? LOGSQRTPI : 0.0;

  for(n=df-2; n > 0; n -= 2)
    g += log(0.5*n);
  return (2.0 * exp(0.5*df*(log(0.5*df)-1.0) - g));
}

/* Integrand function
 */
static double f(double s, int k, int df, int nrng, double q, int isw)
{
  double y=exp((df - 1.0)*log(s) + 0.5*df*(1.0 - s*s));
  if(isw == 0)
    return (y);
  else
    return (y*pow(rng_lp(s*q, k), (double)nrng));
}


double smrng_lp(double q, int k, int df, int nrng)
{
  // 40 nodes and weights for Gauss-Legendre quadrature.
  const double nd[20]={
    0.998237709710559200349622702420586492,
    0.990726238699457006453054352221372155,
    0.977259949983774262663370283712903807,
    0.957916819213791655804540999452759285,
    0.932812808278676533360852166845205716,
    0.902098806968874296728253330868493104,
    0.865959503212259503820781808354619964,
    0.824612230833311663196320230666098774,
    0.778305651426519387694971545506494848,
    0.727318255189927103280996451754930549,
    0.671956684614179548379354514961494110,
    0.612553889667980237952612450230694877,
    0.549467125095128202075931305529517970,
    0.483075801686178712908566574244823005,
    0.413779204371605001524879745803713683,
    0.341994090825758473007492481179194310,
    0.268152185007253681141184344808596183,
    0.192697580701371099715516852065149895,
    0.116084070675255208483451284408024114,
    0.0387724175060508219331934440246232947
  };
  const double wt[20]={
    0.00452127709853319125847173287818533273,
    0.0104982845311528136147421710672796524,
    0.0164210583819078887128634848823639273,
    0.0222458491941669572615043241842085732,
    0.0279370069800234010984891575077210773,
    0.0334601952825478473926781830864108490,
    0.0387821679744720176399720312904461623,
    0.0438709081856732719916746860417154958,
    0.0486958076350722320614341604481463881,
    0.0532278469839368243549964797722605046,
    0.0574397690993915513666177309104259856,
    0.0613062424929289391665379964083985959,
    0.0648040134566010380745545295667527300,
    0.0679120458152339038256901082319239860,
    0.0706116473912867796954836308552868324,
    0.0728865823958040590605106834425178359,
    0.0747231690579682642001893362613246732,
    0.0761103619006262423715580759224948230,
    0.0770398181642479655883075342838102485,
    0.0775059479784248112637239629583263270
  };
  double  sl, su, cnst, rlq, ruq, sll, x;
  double  p=0.0, p1, cntr, wdth;
  int     isw=0, i;

  if(q <= 0.0)
    return(0.0);
  // df = infinity
  if(df <= 0)
    return(pow(rng_lp(q, k), (double)nrng));

  // Upper and lower integral limits
  sl = sqrt(chi2l(df)/df);
  su = sqrt(chi2u(df)/df);
  cnst = coef(df);

  // Lower limit of max range.
  rlq = rlower(k, nrng)/q;
  if(rlq >= su)
    return(0.0);
  if(rlq > sl)
    sl = rlq;

  // Upper limit of max range.
  ruq = rupper(k, nrng)/q;
  if(ruq <= sl)
    return(1.0);

  // If ru/q < su, then integrate twice:
  //   1) \int_{sl}^{ru/q}
  //   2) \int_{ru/q}^{su}, where rng_p(s*q)=1.0
  // First integrate the latter (with isw=0).
  if(ruq < su) {
    sll = sl;
    sl = ruq;
  }
  else
    isw = 1;

  for( ; isw < 2; isw++) {
    p1 = 0.0;
    cntr = 0.5*(sl+su);
    wdth = 0.5*(su-sl);
    for(i=0; i < 20; i++) {
      x = wdth*nd[i];
      p1 += wt[i] * (f(cntr-x, k, df, nrng, q, isw)
                     + f(cntr+x, k, df, nrng, q, isw));
    }
    p += wdth*p1;

    if(isw == 0) {
      su = ruq;
      sl = sll;
    }
  }

  return (cnst*p);
}
