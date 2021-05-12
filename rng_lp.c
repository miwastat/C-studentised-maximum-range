/*
 *  double rng_lp(double r, int k)
 *    returns lower probability of the range distribution.
 *
 *  Arguments
 *    r: range value
 *    k: number of treatments
 *
 *  Required functions
 *    extern double nrml_p()
 *    static double nrml_ip()
 *    static double ulim()
 *    static double f()
 *
 *  Include files
 *    <math.h>
 *
 *  Note
 *    1) The 20-node Gauss-Legendre quadrature is used.
 *    2) The accuracy is of order e-12 (I hope).
 *    3) This accuracy is not guaranteed for k > 1000.
 *
 *  References
 *    H. O. Hartley (1942). Biometrika, 32, 309-310.
 *    
 *  Stored in 
 *    rng_lp.c
 *
 *  History
 *    c. 1994:    First written in Fortran.
 *    2018-11-01: Created with ulim() function.
 *    2019-04-23: Modified for new version.
 *    2021-05-08: Last modified.
 *
 *  License
 *    GPLv3 (Free and No Warranty)
 *    https://www.gnu.org/licenses/
 *
 *  Coded by Tetsuhisa Miwa.
 */


#include <math.h>
#define BORDER  3.7
#define CNST0   0.398942280401432677939946059934381868  // 1/sqrt(2*pi)
#define MAX(X, Y)  ((X < Y) ? Y : X)
#define MIN(X, Y)  ((X < Y) ? X : Y)

extern double nrml_p(double u, int upper);

/* Normal probability in interval (a, b).
 */
static double nrml_ip(double a, double b)
{
  double  border=(BORDER);

  if(a >= b)
    return(0.0);

  if(a > border)
    return(nrml_p(a, 1) - nrml_p(b, 1));
  else if(b < -border)
    return(nrml_p(b, 0) - nrml_p(a, 0));
  else
    return(nrml_p(b, 2) - nrml_p(a, 2));
}

/* Upper integral limit for Hartley's formula.
 * The limit depends on both r and k.
 */
static double ulim(double r, int k)
{
  double ulim13, rmin, rmin10, a1, a2, a3, d1, d2, w, z;

  // If k > 1000, use the value for k=1000.
  if(k > 1000) 
    k = 1000;

  // Approximate upper limit at r=13.
  ulim13 = 1.403*sqrt(log((double)k) + 28.127);

  // Calculate approximate rmin(k).
  // Return 0.0 if r <= rmin(k).
  w = log((double)k);
  rmin = exp(2.3641 - 4.669/w - 9.499/(w*w) - 13.293/(w*w*w));
  if(r <= rmin)
    return(0.0);

  // Upper integral limit depending on whether k <= 10 or k > 10.
  if(k <= 10) {
    d1 = 0.02173*log(8.7/(k - 1.3));
    d2 = 8.4 + 0.2*k;
    z = MIN(1.0, MAX(0.0, d1*(d2 - r)) + 0.199 + 0.134*r - 0.00500*r*r);
  } else {
    rmin10 = 0.07856;
    a1 = (k < 30) ? 8.889*log(k - 3.0) + 24.70 : 54.0;
    a2 = (k < 30) ? 0.06873*log(k - 7.0) + 0.9245 : 1.14;
    if(k < 22)
      a3 = -0.6031*log(k + 6.0) + 1.6877;
    else
      a3 = (k <= 35) ? -0.31 : 0.308*log(k - 5.0) - 1.3576;
    w = a1*pow((r - rmin + rmin10)/(42.0 - rmin + rmin10), a2) + a3;
    z = (w > 9.0) ? 1.0 : 0.199 + 0.134*w - 0.00500*w*w;
  }
  return(ulim13*z);
}

/* Integrand function
 */
static double f(double x, double r, int k)
{
  double y = exp(-0.5*x*x) * pow(nrml_ip(x - r, x), k - 1);
  return(y);
}

double rng_lp(double r, int k)
{
  // 20 nodes and weights for Gauss-Legendre quadrature.
  const double nd[10]={
    0.993128599185094924786122388471320278,
    0.963971927277913791267666131197277222,
    0.912234428251325905867752441203298113,
    0.839116971822218823394529061701520685,
    0.746331906460150792614305070355641590,
    0.636053680726515025452836696226285937,
    0.510867001950827098004364050955250998,
    0.373706088715419560672548177024927237,
    0.227785851141645078080496195368574625,
    0.0765265211334973337546404093988382110
  };
  const double wt[10]={
    0.0176140071391521183118619623518528164,
    0.0406014298003869413310399522749321099,
    0.0626720483341090635695065351870416064,
    0.0832767415767047487247581432220462061,
    0.101930119817240435036750135480349876,
    0.118194531961518417312377377711382287,
    0.131688638449176626898494499748163135,
    0.142096109318382051329298325067164933,
    0.149172986472603746787828737001969437,
    0.152753387130725850698084331955097593
  };

  double  xu, p=0.0, cntr, wdth, x;
  int     ix;

  if(r <= 0.0)
    return(0.0);

  // Normal probability.
  if(k == 2)
    return(2.0*nrml_p(r/sqrt(2.0), 2));
  
  // Upper integral limit.
  xu = ulim(r, k);

  // 2nd term of Hartley's formula.
  if(xu > 0.5*r) {
    wdth = 0.5*(xu - 0.5*r);
    cntr = 0.5*(xu + 0.5*r);
    for(ix=0; ix < 10; ix++) {
      x = wdth*nd[ix];
      p += wt[ix] * (f(cntr - x, r, k) + f(cntr + x, r, k));
    }
    p *= 2.0*k*wdth*(CNST0);
  }

  // Add 1st term.
  p += pow(2.0*nrml_p(0.5*r, 2), (double)k);
  return(p);
}
