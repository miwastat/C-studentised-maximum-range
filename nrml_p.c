/*
 *  double nrml_p(double u, int upper)
 *    returns lower, upper or central probability
 *    of standard normal distribution.
 *
 *  Arguments
 *    u:     normal deviate
 *    upper: upper==0 -> lower probability
 *           upper==1 -> upper probability
 *           upper==2 -> central probability
 *                       probability from 0 to u (negative for u < 0.0)
 *
 *  Required functions
 *    None
 *
 *  Include files
 *    <math.h>
 *
 *  References
 *    Yamauti, Ziro (ed).
 *      Statistical Tables and Formulas with Computer Applications,
 *      Japanese Standards Association, 1972 (in Japanese).
 *
 *  Note
 *    (1) |u| >  border: Laplace's approximation
 *        |u| <= border: Shenton's approximation
 *        Accuracy
 *          term  border  absolute error
 *            28     3.7    10^(-16)
 *            24     3.3    10^(-14)
 *            20     3.0    10^(-12)
 *            16     2.7    10^(-10)
 *            12     2.4    10^( -8)
 *            10     2.2    10^( -6)
 *
 *    (2) We can use system's erf() if available (and faster).
 *        erf(x)  = 2/sqrt(pi) \int_0^x exp(-t^2) dt
 *                = Pr{ |X(mu=0, sigma=sqrt(1/2))| < x }
 *        erfc(x) = 1.0 - erf(x)
 *        lower_prob = u > 0.0 ? 0.5 + 0.5 * erf(|u|/sqrt(2))
 *                             : 0.5 - 0.5 * erf(|u|/sqrt(2))
 *        or simply  = 0.5 + 0.5 * erf(u/sqrt(2)) if erf(-x) = -erf(x)
 *
 *    (3) If we want upper probability for large u, avoid
 *          upper_prob = 1.0 - lower_prob
 *        which can cause loss of significant digits. Use upper=1.
 *        For small values of u
 *          Pr{0.0 < U < u} = lower_prob - 0.5
 *        also causes loss of significant digits. Use upper=2.
 *
 *  Stored in
 *    nrml_p.c
 *
 *  History:
 *    c. 1988:    First written in Fortran.
 *    2013-05-29: Last modified in the old version.
 *    2017-02-08: TERM and BORDER are fixed.
 *                lower, upper or central probability is specified.
 *    2021-05-07: Last modified.
 *
 *  License
 *    GPLv3 (Free and No Warranty)
 *    https://www.gnu.org/licenses/
 *
 *  Coded by Tetsuhisa Miwa.
 */

#include <math.h>
#define TERM    28
#define BORDER  3.7
#define CNST0   0.398942280401432677939946059934381868  // 1/sqrt(2*pi)

double nrml_p(double u, int upper)
{
  int     term=(TERM), sw=-1;
  double  border=(BORDER), w=fabs(u), uu=u*u, p=0.0;
  double  dnrml=(CNST0) * exp(-0.5*u*u);

  if(w > border) {
    // Laplace's approximation for large |u|.
    for( ; term > 0; term--)
      p = term/(w + p);
    p = dnrml/(w + p);

    if(upper >= 2)
      p = (u > 0.0) ?  0.5 - p : -0.5 + p;
    else if(u > 0.0 && upper != 1 || u < 0.0 && upper == 1)
        p = 1.0 - p;
  }
  else {
    // Shenton's approximation for small |u|.
    for( ; term > 0; term--, sw = -sw)
      p = term*uu / (2.0*term + 1.0 + sw*p);
    p = dnrml*w / (1.0 - p);

    if(upper >= 2)
      p = (u < 0.0) ? -p : p;
    else if(u > 0.0 && upper != 1 || u < 0.0 && upper == 1)
      p = 0.5 + p;
    else
      p = 0.5 - p;
  }
  return(p);
}
